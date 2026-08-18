"""
GATE 3 - Baseline 3: free per-host lookup embedding. DIAGNOSTIC, NOT A GATE.

By construction this baseline cannot produce a prediction for a host it never
saw during training -- it has no learned embedding row for that host. A
fallback is REQUIRED and specified BEFORE running (per the explicit
instruction that an earlier version of this project's plan made exactly this
mistake by leaving it undefined):

FALLBACK CHOSEN: mean of training-host embeddings (not nearest phylogenetic
neighbor). Justification: with only 2 training hosts per leave-one-host-out
fold (3 primary hosts total), "nearest phylogenetic neighbor" reduces to an
arbitrary pick between the only two available options and does not obviously
generalize better than an even blend; the mean is the lower-variance,
better-defined choice at this host count and keeps the baseline
appropriately uninformed (it is a floor/diagnostic, not a competitor, so a
simple fallback is the right level of sophistication for it).

IMPLEMENTATION NOTE (dimensionality): the task specifies an 8-32 dimensional
learned embedding. This is implemented as a per-host ONE-HOT indicator
(dimension = n_training_hosts = 2 for each leave-one-host-out fold with 3
primary hosts) concatenated with the k-mer sequence features, jointly fit
with LogisticRegression/Ridge -- for a linear model this is the natural,
minimal-sufficient "free per-host embedding" (a full-rank one-hot already
lets the model learn an arbitrary per-host additive effect; a higher-
dimensional embedding would be redundant/overparameterized with only 2
training hosts to distinguish). The fallback for the unseen host is the mean
of the training hosts' one-hot vectors (e.g. [0.5, 0.5] for 2 training
hosts), which is the direct linear-model analog of "mean of training-host
embeddings."

Evaluation: leave-one-host-out across the 3 primary hosts. Train on the
pooled training-pool (fold != 0) data of the other 2 hosts; evaluate on the
held-out host's fold-0 test set (same partition used throughout Gate 3 for
comparability with B1/B2).
"""
import numpy as np
import json
from pathlib import Path
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import matthews_corrcoef, roc_auc_score, mean_squared_error
import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
b2mod = import_module("31_baseline2_calibration")

DATA = Path(__file__).resolve().parent.parent / "data"
CACHE = DATA / "baseline_cache"
OUT = Path(__file__).resolve().parent.parent / "out" / "baselines"
OUT.mkdir(parents=True, exist_ok=True)

HOSTS = ["EC", "BS", "PA"]
TEST_FOLD = 0
LOGREG_C = 1.0
RIDGE_ALPHA = 10.0


def load_all():
    data = {}
    for host in HOSTS:
        d = dict(np.load(CACHE / f"{host}_baseline_data.npz", allow_pickle=True))
        d["kmer"] = b2mod.build_kmer_features(d["onehot"])
        data[host] = d
    return data


def run_readout(readout, data, held_out):
    train_hosts = [h for h in HOSTS if h != held_out]

    X_train_parts, y_active_parts, y_strength_parts, active_mask_parts, host_ind_parts = [], [], [], [], []
    for i, h in enumerate(train_hosts):
        d = data[h]
        if readout == "transcription":
            usable, active, strength = d["tx_usable"], d["tx_active"], d["tx_norm"]
            usable_reg = usable
            log_transform = True
        else:
            usable, active, strength = d["tl_usable"], d["tl_active"], d["protein_log10"]
            usable_reg = d["tl_usable_for_regression"]
            log_transform = False
        pool_mask = usable & (d["fold"] != TEST_FOLD)
        idx = np.where(pool_mask)[0]
        X_train_parts.append(d["kmer"][idx])
        y_active_parts.append(active[idx].astype(int))
        y_strength_parts.append(strength[idx])
        active_mask_parts.append(active[idx] & usable_reg[idx])
        host_ind = np.zeros((len(idx), len(train_hosts)), dtype=np.float32)
        host_ind[:, i] = 1.0
        host_ind_parts.append(host_ind)

    X_seq = np.concatenate(X_train_parts)
    host_ind = np.concatenate(host_ind_parts)
    X_train = np.concatenate([X_seq, host_ind], axis=1)
    y_active_train = np.concatenate(y_active_parts)
    y_strength_train = np.concatenate(y_strength_parts)
    active_mask_train = np.concatenate(active_mask_parts)

    clf = LogisticRegression(C=LOGREG_C, max_iter=300).fit(X_train, y_active_train) \
        if len(set(y_active_train)) > 1 else None

    reg = None
    if active_mask_train.sum() >= 3:
        Xa = X_train[active_mask_train]
        ya = y_strength_train[active_mask_train]
        if readout == "transcription":
            ya = np.log1p(np.clip(ya, 0, None))
        if len(set(np.round(ya, 6))) > 1:
            reg = Ridge(alpha=RIDGE_ALPHA).fit(Xa, ya)

    # held-out host test set, with FALLBACK host indicator = mean of training hosts'
    d = data[held_out]
    if readout == "transcription":
        usable, active, strength = d["tx_usable"], d["tx_active"], d["tx_norm"]
        usable_reg = usable
        log_transform = True
    else:
        usable, active, strength = d["tl_usable"], d["tl_active"], d["protein_log10"]
        usable_reg = d["tl_usable_for_regression"]
        log_transform = False
    test_mask = usable & (d["fold"] == TEST_FOLD)
    test_idx = np.where(test_mask)[0]
    X_seq_test = d["kmer"][test_idx]
    fallback_ind = np.full((len(test_idx), len(train_hosts)), 1.0 / len(train_hosts), dtype=np.float32)
    X_test = np.concatenate([X_seq_test, fallback_ind], axis=1)

    y_active_test = active[test_idx].astype(int)
    active_mask_test = active[test_idx] & usable_reg[test_idx]

    pred_prob = clf.predict_proba(X_test)[:, 1] if clf is not None else np.full(len(X_test), np.nan)
    pred_class = (pred_prob > 0.5).astype(int) if clf is not None else np.zeros(len(X_test), dtype=int)
    mcc = matthews_corrcoef(y_active_test, pred_class) if clf is not None and len(set(y_active_test)) > 1 else None
    try:
        auc = roc_auc_score(y_active_test, pred_prob) if clf is not None and len(set(y_active_test)) > 1 else None
    except Exception:
        auc = None

    rho, rmse = None, None
    if reg is not None and active_mask_test.sum() > 1:
        Xa = X_test[active_mask_test]
        y_true = strength[test_idx][active_mask_test]
        if readout == "transcription":
            y_true = np.log1p(np.clip(y_true, 0, None))
        y_pred = reg.predict(Xa)
        if len(set(np.round(y_true, 6))) > 1:
            rho = spearmanr(y_true, y_pred).correlation
        rmse = float(np.sqrt(mean_squared_error(y_true, y_pred)))

    return {
        "held_out_host": held_out, "train_hosts": train_hosts,
        "n_train": int(len(X_train)), "n_test": int(len(test_idx)),
        "mcc": float(mcc) if mcc is not None and mcc == mcc else None,
        "auc": float(auc) if auc is not None and auc == auc else None,
        "spearman_rho": float(rho) if rho is not None and rho == rho else None,
        "rmse_log_strength": rmse,
    }


def main():
    data = load_all()
    results = {}
    for held_out in HOSTS:
        print(f"\n=== held out: {held_out} (train on {[h for h in HOSTS if h != held_out]}) ===")
        results[held_out] = {}
        for readout in ["transcription", "translation"]:
            r = run_readout(readout, data, held_out)
            print(f"  {readout}: {r}")
            results[held_out][readout] = r

    with open(OUT / "baseline3_host_embedding.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {OUT / 'baseline3_host_embedding.json'}")


if __name__ == "__main__":
    main()
