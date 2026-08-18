"""
GATE 3.5 - Task 2 (RS241 extension): B3 (free per-host lookup embedding,
diagnostic) evaluated on the three RS241-derived hosts (S. enterica,
V. natriegens, C. glutamicum) -- Gate 3 only tested the 3 primary hosts.
These are the hosts where the free-embedding fallback question is most
acute, since they are genuinely never seen at all (not even held out from a
shared fold structure -- RS241 has no fold structure, it is a wholly
separate reserved evaluation set, mechanically excluded from training,
audit_leakage.py check 2).

DATA GAP FOUND AND DOCUMENTED (not silently worked around): RS241 lists 241
oligo IDs with expression values (data/rs241.parquet), but the released
supplementary tables only contain REGULATORY SEQUENCE TEXT for 207 of those
241 IDs -- the remaining 34 are not reconstructable from any released table
(checked against both the processed data/three_host_library.parquet and the
raw Supplementary Data Table 1 directly). This was not previously
documented because Gate 3 never needed RS241 sequence text (B3 was
primary-hosts-only). Concretely this reduces every RS241 usable-N figure by
~14-15% for the sequence-based baselines specifically (85-90% of each
config's usable N has a recoverable sequence) -- reported per
host/readout/config below, not silently substituted.

Uses the two frozen configs from data/splits/rs241_configs.json:
  PRIMARY:   train on EC+BS+PA (matches the main benchmark's host set)
  SECONDARY: train on EC+PA only (drops weak-coverage B. subtilis, higher N)
No fold rotation applies to the RS241 side (it has no folds) -- primary-host
training uses ALL of that host's frozen-split data (all 5 folds pooled),
since none of it is being tested against here (the test set is RS241
sequences, structurally disjoint, confirmed by audit_leakage.py check 2).
"""
import numpy as np
import pandas as pd
import json
import sys
from pathlib import Path
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import matthews_corrcoef, roc_auc_score, mean_squared_error
import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
b2mod = import_module("31_baseline2_calibration")

DATA = Path(__file__).resolve().parent.parent / "data"
CACHE = DATA / "baseline_cache"
OUT = Path(__file__).resolve().parent.parent / "out" / "baselines"
OUT.mkdir(parents=True, exist_ok=True)

PRIMARY_HOSTS = ["EC", "BS", "PA"]
RS241_HOSTS = ["SE", "VN", "CG"]
LOGREG_C = 1.0
RIDGE_ALPHA = 10.0
SEQ_LEN = 165
BASE_TO_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}


def one_hot(seq):
    arr = np.zeros((4, SEQ_LEN), dtype=np.float32)
    for i, b in enumerate(seq):
        idx = BASE_TO_IDX.get(b)
        if idx is not None:
            arr[idx, i] = 1.0
    return arr


def load_primary_full(host):
    """Full pooled data for a primary host, ALL folds (no test-fold holdout
    -- see docstring: nothing from the primary hosts' own folds is being
    tested against in this evaluation)."""
    d = dict(np.load(CACHE / f"{host}_baseline_data.npz", allow_pickle=True))
    d["kmer"] = b2mod.build_kmer_features(d["onehot"])
    return d


def load_rs241_with_sequence():
    rs241 = pd.read_parquet(DATA / "rs241.parquet")
    rs241["id"] = rs241["id"].astype(str)
    lib = pd.read_parquet(DATA / "three_host_library.parquet")
    lib["OLIGO ID"] = lib["OLIGO ID"].astype(str)
    seq_map = dict(zip(lib["OLIGO ID"], lib["Regulatory Sequence"]))
    rs241["has_sequence"] = rs241["id"].isin(seq_map)
    rs241["sequence"] = rs241["id"].map(seq_map)
    n_total = len(rs241)
    n_avail = int(rs241["has_sequence"].sum())
    print(f"RS241 sequence availability: {n_avail}/{n_total} ids have recoverable regulatory sequence text")
    rs241_seq = rs241[rs241["has_sequence"]].copy()
    onehot = np.stack([one_hot(s) for s in rs241_seq["sequence"]]).astype(np.float32)
    kmer = b2mod.build_kmer_features(onehot)
    rs241_seq = rs241_seq.reset_index(drop=True)
    return rs241_seq, kmer, n_total, n_avail


def run_config(config_name, train_hosts, held_out, readout, primary_data, rs241_seq, rs241_kmer):
    X_parts, y_active_parts, y_strength_parts, active_mask_parts, host_ind_parts = [], [], [], [], []
    for i, h in enumerate(train_hosts):
        d = primary_data[h]
        if readout == "transcription":
            usable, active, strength = d["tx_usable"], d["tx_active"], d["tx_norm"]
            usable_reg = usable
        else:
            usable, active, strength = d["tl_usable"], d["tl_active"], d["protein_log10"]
            usable_reg = d["tl_usable_for_regression"]
        idx = np.where(usable)[0]
        X_parts.append(d["kmer"][idx])
        y_active_parts.append(active[idx].astype(int))
        y_strength_parts.append(strength[idx])
        active_mask_parts.append(active[idx] & usable_reg[idx])
        host_ind = np.zeros((len(idx), len(train_hosts)), dtype=np.float32)
        host_ind[:, i] = 1.0
        host_ind_parts.append(host_ind)

    X_seq = np.concatenate(X_parts)
    host_ind = np.concatenate(host_ind_parts)
    X_train = np.concatenate([X_seq, host_ind], axis=1)
    y_active_train = np.concatenate(y_active_parts)
    y_strength_train = np.concatenate(y_strength_parts)
    active_mask_train = np.concatenate(active_mask_parts)

    clf = LogisticRegression(C=LOGREG_C, max_iter=300).fit(X_train, y_active_train) \
        if len(set(y_active_train)) > 1 else None
    reg = None
    if active_mask_train.sum() >= 3:
        Xa, ya = X_train[active_mask_train], y_strength_train[active_mask_train]
        if readout == "transcription":
            ya = np.log1p(np.clip(ya, 0, None))
        if len(set(np.round(ya, 6))) > 1:
            reg = Ridge(alpha=RIDGE_ALPHA).fit(Xa, ya)

    # RS241 held-out host test set: usable per the frozen config, AND sequence-available
    tag = "EC_BS_PA" if set(train_hosts) == {"EC", "BS", "PA"} else "EC_PA"
    ro_tag = "tx" if readout == "transcription" else "tl"
    usable_col = f"usable_heldout_{held_out}_train_{tag}_{ro_tag}"
    usable_mask = rs241_seq[usable_col].fillna(False).values.astype(bool)
    test_idx = np.where(usable_mask)[0]
    if len(test_idx) == 0:
        return {"held_out_host": held_out, "config": config_name, "n_test": 0,
                "note": "no usable+sequence-available rows for this held-out host/readout/config"}

    X_seq_test = rs241_kmer[test_idx]
    fallback_ind = np.full((len(test_idx), len(train_hosts)), 1.0 / len(train_hosts), dtype=np.float32)
    X_test = np.concatenate([X_seq_test, fallback_ind], axis=1)

    val_col = f"tx_log2_{held_out}" if readout == "transcription" else f"tl_log10_{held_out}"
    y_strength_test = rs241_seq.iloc[test_idx][val_col].values.astype(np.float32)

    if readout == "transcription":
        y_active_test = (y_strength_test > 0).astype(int)  # log2 fold-change > 0 == active/induced
    else:
        # translation: no floor-value correction attempted for RS241 (out of scope,
        # small N) -- active defined the SAME way Gate 1/2 did for RS241 usability,
        # i.e. any usable measurement counts; this is a simplification vs the primary
        # hosts' floor-corrected definition, disclosed here explicitly.
        y_active_test = np.ones(len(test_idx), dtype=int)

    pred_prob = clf.predict_proba(X_test)[:, 1] if clf is not None else np.full(len(X_test), np.nan)
    pred_class = (pred_prob > 0.5).astype(int) if clf is not None else np.zeros(len(X_test), dtype=int)
    mcc = matthews_corrcoef(y_active_test, pred_class) if clf is not None and len(set(y_active_test)) > 1 else None
    try:
        auc = roc_auc_score(y_active_test, pred_prob) if clf is not None and len(set(y_active_test)) > 1 else None
    except Exception:
        auc = None

    rho = None
    if reg is not None and len(test_idx) > 1:
        y_reg_test = y_strength_test
        if readout == "transcription":
            y_reg_test = np.log1p(np.clip(y_reg_test, 0, None)) if (y_reg_test >= 0).all() else y_reg_test
        y_pred = reg.predict(X_test)
        if len(set(np.round(y_reg_test, 6))) > 1:
            rho = spearmanr(y_reg_test, y_pred).correlation

    return {"held_out_host": held_out, "config": config_name, "train_hosts": train_hosts,
            "n_train": int(len(X_train)), "n_test": int(len(test_idx)),
            "mcc": float(mcc) if mcc is not None and mcc == mcc else None,
            "auc": float(auc) if auc is not None and auc == auc else None,
            "spearman_rho": float(rho) if rho is not None and rho == rho else None}


def main():
    primary_data = {h: load_primary_full(h) for h in PRIMARY_HOSTS}
    rs241_seq, rs241_kmer, n_total, n_avail = load_rs241_with_sequence()

    configs = {"PRIMARY": ["EC", "BS", "PA"], "SECONDARY": ["EC", "PA"]}
    results = {}
    for config_name, train_hosts in configs.items():
        print(f"\n=== Config {config_name}: train on {train_hosts} ===")
        results[config_name] = {}
        for held_out in RS241_HOSTS:
            results[config_name][held_out] = {}
            for readout in ["transcription", "translation"]:
                r = run_config(config_name, train_hosts, held_out, readout, primary_data, rs241_seq, rs241_kmer)
                print(f"  held_out={held_out} {readout}: n_test={r.get('n_test')}, "
                      f"mcc={r.get('mcc')}, auc={r.get('auc')}, rho={r.get('spearman_rho')}")
                results[config_name][held_out][readout] = r

    output = {"sequence_availability": {"n_total_rs241": n_total, "n_sequence_available": n_avail,
                                         "pct_available": round(100 * n_avail / n_total, 1)},
              "results": results}
    with open(OUT / "baseline3_rs241_extension.json", "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nWrote {OUT / 'baseline3_rs241_extension.json'}")


if __name__ == "__main__":
    main()
