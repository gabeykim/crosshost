"""
GATE 5 - Task 4: train models for RS241 zero-shot evaluation (S. enterica,
V. natriegens, C. glutamicum -- the only hosts with NO training data at all,
the purest cross-host test in the project).

Gate 4's LOHO models don't apply here: they were all trained on exactly 2 of
the 3 primary hosts (holding the 3rd out). For RS241 hosts, none of the 3
primary hosts is held out -- both frozen configs use different pools:
  PRIMARY:   train on ALL of EC+BS+PA (all 5 folds each, no exclusion needed
             -- audit_leakage.py check 2 already confirms RS241 sequences are
             structurally disjoint from the primary-host fold assignment)
  SECONDARY: train on ALL of EC+PA only (drops weak-coverage B. subtilis,
             the deeper-phylogenetic-distance variant per the frozen config)

5 SEEDS per (config, variant) instead of Gate 4's "5 folds" -- there is no
natural data partition here (every seed sees identical training data), so
this is draw-level/training-stochasticity uncertainty, not partition-level
fold uncertainty. Disclosed as such, not conflated with Gate 4's fold
resolution. 2 configs x 2 variants x 5 seeds = 20 fits.

RS241 sequence availability: 207/241 recoverable (Gate 3.5 finding) --
evaluation is restricted to those 207 throughout.
"""
import sys
import time
import json
import numpy as np
import pandas as pd
import torch
from pathlib import Path
from scipy.stats import spearmanr
from sklearn.metrics import matthews_corrcoef, roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
m = import_module("40_film_cnn_model")

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"
MODELS_DIR = OUT / "models"
MODELS_DIR.mkdir(parents=True, exist_ok=True)
DEVICE = "mps" if torch.backends.mps.is_available() else "cpu"

CONFIGS = {"PRIMARY": ["EC", "BS", "PA"], "SECONDARY": ["EC", "PA"]}
RS241_HOSTS = ["SE", "VN", "CG"]
N_SEEDS = 5
EPOCHS = 35
PATIENCE = 8
BATCH_SIZE = 1024
SEQ_LEN = 165
BASE_TO_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}


def one_hot(seq):
    arr = np.zeros((4, SEQ_LEN), dtype=np.float32)
    for i, b in enumerate(seq):
        idx = BASE_TO_IDX.get(b)
        if idx is not None:
            arr[idx, i] = 1.0
    return arr


def load_rs241_eval_arrays():
    """Returns rs241_seq_df (207 rows, sequence-available only), onehot, kmer-free
    (not needed here -- CNN uses raw onehot)."""
    rs241 = pd.read_parquet(DATA / "rs241.parquet")
    rs241["id"] = rs241["id"].astype(str)
    lib = pd.read_parquet(DATA / "three_host_library.parquet")
    lib["OLIGO ID"] = lib["OLIGO ID"].astype(str)
    seq_map = dict(zip(lib["OLIGO ID"], lib["Regulatory Sequence"]))
    rs241["has_sequence"] = rs241["id"].isin(seq_map)
    rs241["sequence"] = rs241["id"].map(seq_map)
    rs241_seq = rs241[rs241["has_sequence"]].copy().reset_index(drop=True)
    onehot = np.stack([one_hot(s) for s in rs241_seq["sequence"]]).astype(np.float32)
    return rs241_seq, onehot


def evaluate_rs241(model, host_vec_row, rs241_seq, onehot, held_out, config_tag):
    results = {}
    hv = np.tile(host_vec_row, (len(onehot), 1)).astype(np.float32)
    preds = m.predict_film(model, onehot, hv)
    for readout, active_col_prefix, strength_col, usable_suffix in [
        ("transcription", "tx", f"tx_log2_{held_out}", "tx"),
        ("translation", "tl", f"tl_log10_{held_out}", "tl"),
    ]:
        usable_col = f"usable_heldout_{held_out}_train_{config_tag}_{usable_suffix}"
        usable_mask = rs241_seq[usable_col].fillna(False).values.astype(bool)
        n_test = int(usable_mask.sum())
        if n_test < 3:
            results[readout] = {"n_test": n_test, "mcc": None, "auc": None, "spearman_rho": None}
            continue
        y_strength = rs241_seq.loc[usable_mask, strength_col].values.astype(np.float64)
        pred_strength = preds[f"{active_col_prefix}_strength"][usable_mask]
        rho = spearmanr(y_strength, pred_strength).correlation if len(set(np.round(y_strength, 6))) > 1 else None

        if readout == "transcription":
            y_active = (y_strength > 0).astype(int)
            prob = 1 / (1 + np.exp(-preds["tx_active_logit"][usable_mask]))
            pred_class = (prob > 0.5).astype(int)
            mcc = matthews_corrcoef(y_active, pred_class) if len(set(y_active)) > 1 else None
            try:
                auc = roc_auc_score(y_active, prob) if len(set(y_active)) > 1 else None
            except Exception:
                auc = None
        else:
            # translation active/inactive not computable for RS241 (no floor-correction
            # available for this reserved set) -- disclosed simplification, matches
            # Gate 3.5's scripts/37 B3-RS241 extension precedent
            mcc, auc = None, None

        results[readout] = {
            "n_test": n_test,
            "mcc": float(mcc) if mcc is not None and mcc == mcc else None,
            "auc": float(auc) if auc is not None and auc == auc else None,
            "spearman_rho": float(rho) if rho is not None and rho == rho else None,
        }
    return results


def main():
    t_start = time.time()
    g_z, gcols, p_z, pcols, imputed = m.load_host_features()
    assert imputed == []
    variants = {"genomic": (g_z, len(gcols)), "physiology": (p_z, len(pcols))}

    rs241_seq, onehot = load_rs241_eval_arrays()
    print(f"RS241 sequence-available rows: {len(rs241_seq)} / 241")

    all_results = {}
    n_fits_done = 0
    n_fits_total = len(CONFIGS) * len(variants) * N_SEEDS
    for config_name, train_hosts in CONFIGS.items():
        config_tag = "_".join(train_hosts)  # matches usable_heldout_*_train_EC_BS_PA_* / _EC_PA_* column naming
        all_results[config_name] = {}
        for variant_name, (feat_df, host_dim) in variants.items():
            all_results[config_name][variant_name] = {"per_seed": {}}
            seq_tr, hv_tr, targets_tr, masks_tr, _ = m.build_pooled_arrays(
                train_hosts, feat_df, fold_filter_fn=None)  # ALL folds -- no primary host held out here
            print(f"\n=== config={config_name} variant={variant_name}: pool N={len(seq_tr)} ===")

            for seed in range(N_SEEDS):
                t0 = time.time()
                model, epochs_run = m.train_film_model(
                    seq_tr, hv_tr, targets_tr, masks_tr, host_dim=host_dim,
                    epochs=EPOCHS, batch_size=BATCH_SIZE, seed=seed, device=DEVICE,
                    use_film=True, patience=PATIENCE, verbose=False)
                elapsed = time.time() - t0
                n_fits_done += 1

                seed_result = {}
                for held_out in RS241_HOSTS:
                    host_vec_row = feat_df.loc[held_out].values
                    seed_result[held_out] = evaluate_rs241(model, host_vec_row, rs241_seq, onehot, held_out, config_tag)

                all_results[config_name][variant_name]["per_seed"][seed] = {
                    "eval": seed_result, "epochs_run": epochs_run, "wall_clock_sec": elapsed, "n_train": int(len(seq_tr))
                }
                print(f"[{n_fits_done}/{n_fits_total}] config={config_name} variant={variant_name} seed={seed}: "
                      f"{elapsed:.1f}s, epochs_run={epochs_run}")
                for h in RS241_HOSTS:
                    print(f"    {h}: {seed_result[h]}")

                model_path = MODELS_DIR / f"rs241_{config_name}_{variant_name}_seed{seed}.pt"
                torch.save(model.state_dict(), model_path)

                with open(OUT / "gate5_rs241_results.json", "w") as f:
                    json.dump(all_results, f, indent=2, default=str)

    print(f"\nTotal wall time: {time.time()-t_start:.1f}s")
    with open(OUT / "gate5_rs241_results.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"Wrote {OUT / 'gate5_rs241_results.json'}")


if __name__ == "__main__":
    main()
