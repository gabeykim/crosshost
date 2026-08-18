"""
GATE 3.5 - Task 2: full 5-fold rotation for B1-B4 (Gate 3 used a single
fixed test fold throughout). Re-runs each baseline with every fold in turn
as the held-out test set, using the SAME logic as the Gate 3 scripts
(imported and called directly, with each module's TEST_FOLD global
monkey-patched per fold rather than re-implemented, to guarantee identical
methodology to Gate 3 -- only the fold varies).

Reports, per baseline x readout x host: mean/std/min/max ACROSS folds.
For B2 specifically, also decomposes fold-to-fold variance vs draw-to-draw
variance (the two are conflated in Gate 3's single-fold report) -- this is
the number that determines whether Gate 3's single-fold results are
trustworthy on their own.

Uses the k-mer/linear substitute throughout (not the CNN), exactly as Gate 3
did, specifically because it is ~4 orders of magnitude faster and makes a
full 5x rotation (2,100+ B2 fits alone) affordable in seconds rather than
requiring the CNN's multi-hour-per-fit cost (see scripts/38 for that cost,
measured separately for Gate 4/7 planning).
"""
import numpy as np
import pandas as pd
import json
import time
import sys
from pathlib import Path
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
b1mod = import_module("30_baseline1_mean_majority")
b2mod = import_module("31_baseline2_calibration")
b3mod = import_module("32_baseline3_host_embedding")

DATA = Path(__file__).resolve().parent.parent / "data"
CACHE = DATA / "baseline_cache"
OUT = Path(__file__).resolve().parent.parent / "out" / "baselines"
OUT.mkdir(parents=True, exist_ok=True)

HOSTS = ["EC", "BS", "PA"]
N_FOLDS = 5
N_VALUES = [0, 10, 30, 100, 300, 1000, 3000]


def load_primary_data():
    data = {}
    for host in HOSTS:
        d = dict(np.load(CACHE / f"{host}_baseline_data.npz", allow_pickle=True))
        d["kmer"] = b2mod.build_kmer_features(d["onehot"])
        data[host] = d
    return data


def run_b1_all_folds(data):
    results = {host: {"transcription": [], "translation": []} for host in HOSTS}
    for test_fold in range(N_FOLDS):
        b1mod.TEST_FOLD = test_fold
        for host in HOSTS:
            d = data[host]
            tx = b1mod.evaluate_readout(d, "tx_usable", "tx_active", "tx_norm")
            tl = b1mod.evaluate_readout(d, "tl_usable", "tl_active", "protein_log10",
                                         usable_reg_key="tl_usable_for_regression")
            results[host]["transcription"].append({"fold": test_fold, **tx})
            results[host]["translation"].append({"fold": test_fold, **tl})
    b1mod.TEST_FOLD = 0  # restore default
    return results


def run_b2_all_folds(data):
    results = {host: {"transcription": {}, "translation": {}} for host in HOSTS}
    for test_fold in range(N_FOLDS):
        b2mod.TEST_FOLD = test_fold
        for host in HOSTS:
            d = data[host]
            for readout in ["transcription", "translation"]:
                res, meta = b2mod.run_host_readout(host, readout, d, d["kmer"])
                for N in N_VALUES:
                    results[host][readout].setdefault(N, {})[test_fold] = {"draws": res[N], "meta": meta}
    b2mod.TEST_FOLD = 0
    return results


def run_b3_all_folds(data):
    results = {host: {"transcription": [], "translation": []} for host in HOSTS}
    for test_fold in range(N_FOLDS):
        b3mod.TEST_FOLD = test_fold
        for held_out in HOSTS:
            for readout in ["transcription", "translation"]:
                r = b3mod.run_readout(readout, data, held_out)
                results[held_out][readout].append({"fold": test_fold, **r})
    b3mod.TEST_FOLD = 0
    return results


def run_b4_all_folds():
    df = pd.read_parquet(DATA / "three_host_library.parquet")
    df["OLIGO ID"] = df["OLIGO ID"].astype(str)
    fold_df = pd.read_parquet(DATA / "splits" / "fold_assignment_FINAL.parquet")
    id_to_fold = dict(zip(fold_df["OLIGO ID"], fold_df["fold"]))
    df["fold"] = df["OLIGO ID"].map(id_to_fold)
    df = df[df["fold"].notna()].copy()
    df["fold"] = df["fold"].astype(int)

    three_host = pd.read_parquet(DATA / "three_host.parquet")
    three_host["OLIGO ID"] = three_host["OLIGO ID"].astype(str)
    df = df.merge(three_host[["OLIGO ID"] + [f"{h}_tx_usable" for h in HOSTS] +
                              [f"{h}_tx_active" for h in HOSTS] + [f"{h}_tx_norm" for h in HOSTS]],
                   on="OLIGO ID", suffixes=("", "_dup"))

    results = {host: {"transcription": [], "translation": []} for host in HOSTS}
    for test_fold in range(N_FOLDS):
        for host in HOSTS:
            test_mask = (df["fold"] == test_fold) & df[f"{host}_tx_usable"]
            active_mask = test_mask & df[f"{host}_tx_active"]
            sigma_score = df.loc[active_mask, f"{host}_best_sigma70_match_score"]
            tx_val = df.loc[active_mask, f"{host}_tx_norm"]
            valid = sigma_score.notna() & tx_val.notna()
            log_tx = np.log1p(np.clip(tx_val[valid], 0, None))
            rho_tx = spearmanr(sigma_score[valid], log_tx).correlation if valid.sum() > 1 else None

            protein_col = f"{host}_protein (log10)"
            dg = df.loc[test_mask, f"{host}_delta_G"]
            protein = df.loc[test_mask, protein_col]
            valid_tl = dg.notna() & protein.notna()
            rho_tl = spearmanr(dg[valid_tl], protein[valid_tl]).correlation if valid_tl.sum() > 1 else None

            results[host]["transcription"].append({"fold": test_fold, "n": int(valid.sum()),
                "spearman_rho": float(rho_tx) if rho_tx is not None and rho_tx == rho_tx else None})
            results[host]["translation"].append({"fold": test_fold, "n": int(valid_tl.sum()),
                "spearman_rho": float(rho_tl) if rho_tl is not None and rho_tl == rho_tl else None})
    return results


def fold_summary(per_fold_values):
    vals = [v for v in per_fold_values if v is not None]
    if not vals:
        return {"mean": None, "std": None, "min": None, "max": None, "n_folds": 0}
    return {"mean": float(np.mean(vals)), "std": float(np.std(vals)),
            "min": float(np.min(vals)), "max": float(np.max(vals)), "n_folds": len(vals)}


def compile_b1_summary(b1_results):
    out = {}
    for host in HOSTS:
        out[host] = {}
        for readout in ["transcription", "translation"]:
            rows = b1_results[host][readout]
            out[host][readout] = {
                "mcc": fold_summary([r["mcc"] for r in rows]),
                "auc": fold_summary([r["auc"] for r in rows]),
                "spearman_rho": fold_summary([r["spearman_rho"] for r in rows]),
                "per_fold": rows,
            }
    return out


def compile_b3_summary(b3_results):
    out = {}
    for host in HOSTS:
        out[host] = {}
        for readout in ["transcription", "translation"]:
            rows = b3_results[host][readout]
            out[host][readout] = {
                "mcc": fold_summary([r["mcc"] for r in rows]),
                "auc": fold_summary([r["auc"] for r in rows]),
                "spearman_rho": fold_summary([r["spearman_rho"] for r in rows]),
                "per_fold": [{"fold": r["fold"], "mcc": r["mcc"], "auc": r["auc"],
                              "spearman_rho": r["spearman_rho"], "n_test": r["n_test"]} for r in rows],
            }
    return out


def compile_b4_summary(b4_results):
    out = {}
    for host in HOSTS:
        out[host] = {}
        for readout in ["transcription", "translation"]:
            rows = b4_results[host][readout]
            out[host][readout] = {
                "spearman_rho": fold_summary([r["spearman_rho"] for r in rows]),
                "per_fold": rows,
            }
    return out


def compile_b2_summary_and_variance(b2_results):
    """For each host/readout/N: fold-mean-of-draws per fold, then
    fold-to-fold std of those fold-means (FOLD variance) vs the average
    within-fold draw-to-draw std (DRAW variance). Ratio > 1 means fold
    identity matters more than which random draw you happened to get --
    the methodologically important number this task asked for."""
    summary = {}
    variance_comparison = {}
    for host in HOSTS:
        summary[host] = {}
        variance_comparison[host] = {}
        for readout in ["transcription", "translation"]:
            summary[host][readout] = {}
            variance_comparison[host][readout] = {}
            for N in N_VALUES:
                per_fold_entry = b2_results[host][readout][N]
                fold_means_rho, fold_stds_rho = [], []
                fold_means_mcc = []
                for test_fold in range(N_FOLDS):
                    draws = per_fold_entry[test_fold]["draws"]
                    rhos = [dr["spearman_rho"] for dr in draws if dr["spearman_rho"] is not None]
                    mccs = [dr["mcc"] for dr in draws if dr["mcc"] is not None]
                    if rhos:
                        fold_means_rho.append(float(np.mean(rhos)))
                        fold_stds_rho.append(float(np.std(rhos)))
                    if mccs:
                        fold_means_mcc.append(float(np.mean(mccs)))
                summary[host][readout][N] = {
                    "rho_across_folds": fold_summary(fold_means_rho),
                    "mcc_across_folds": fold_summary(fold_means_mcc),
                }
                if fold_means_rho:
                    fold_to_fold_std = float(np.std(fold_means_rho))
                    draw_to_draw_std_avg = float(np.mean(fold_stds_rho)) if fold_stds_rho else None
                    ratio = (fold_to_fold_std / draw_to_draw_std_avg) if draw_to_draw_std_avg and draw_to_draw_std_avg > 1e-9 else None
                    variance_comparison[host][readout][N] = {
                        "fold_to_fold_std_of_rho": fold_to_fold_std,
                        "mean_within_fold_draw_to_draw_std_of_rho": draw_to_draw_std_avg,
                        "ratio_fold_std_over_draw_std": ratio,
                        "fold_variance_dominates": ratio is not None and ratio > 1.0,
                    }
    return summary, variance_comparison


def main():
    t0 = time.time()
    print("Loading primary-host baseline cache + k-mer features...")
    data = load_primary_data()

    print("\n=== Running B1 (mean/majority) across 5 folds ===")
    b1_raw = run_b1_all_folds(data)
    b1_summary = compile_b1_summary(b1_raw)
    for host in HOSTS:
        for readout in ["transcription", "translation"]:
            s = b1_summary[host][readout]
            print(f"  {host} {readout}: mcc={s['mcc']['mean']:.3f}±{s['mcc']['std']:.3f}  "
                  f"auc={s['auc']['mean']:.3f}±{s['auc']['std']:.3f}")

    print(f"\n=== Running B2 (calibration curve) across 5 folds ({time.time()-t0:.1f}s elapsed) ===")
    b2_raw = run_b2_all_folds(data)
    b2_summary, b2_variance = compile_b2_summary_and_variance(b2_raw)
    for host in HOSTS:
        for readout in ["transcription", "translation"]:
            for N in [100, 3000]:
                vc = b2_variance[host][readout].get(N)
                if vc:
                    print(f"  {host} {readout} N={N}: fold_std={vc['fold_to_fold_std_of_rho']:.4f}, "
                          f"draw_std={vc['mean_within_fold_draw_to_draw_std_of_rho']}, "
                          f"ratio={vc['ratio_fold_std_over_draw_std']}")

    print(f"\n=== Running B3 (free host embedding, LOHO) across 5 folds ({time.time()-t0:.1f}s elapsed) ===")
    b3_raw = run_b3_all_folds(data)
    b3_summary = compile_b3_summary(b3_raw)
    for host in HOSTS:
        for readout in ["transcription", "translation"]:
            s = b3_summary[host][readout]
            print(f"  held_out={host} {readout}: rho={s['spearman_rho']['mean']:.3f}±{s['spearman_rho']['std']:.3f}")

    print(f"\n=== Running B4 (biophysical) across 5 folds ({time.time()-t0:.1f}s elapsed) ===")
    b4_raw = run_b4_all_folds()
    b4_summary = compile_b4_summary(b4_raw)
    for host in HOSTS:
        for readout in ["transcription", "translation"]:
            s = b4_summary[host][readout]
            print(f"  {host} {readout}: rho={s['spearman_rho']['mean']:.3f}±{s['spearman_rho']['std']:.3f}")

    print(f"\nTotal wall time: {time.time()-t0:.1f}s")

    output = {
        "n_folds": N_FOLDS,
        "b1_mean_majority": {"summary": b1_summary},
        "b2_calibration": {"summary": b2_summary, "variance_comparison": b2_variance},
        "b3_free_host_embedding": {"summary": b3_summary},
        "b4_biophysical": {"summary": b4_summary},
    }
    with open(OUT / "fold_rotation_all_baselines.json", "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"Wrote {OUT / 'fold_rotation_all_baselines.json'}")

    # also dump raw B2 per-draw data separately (large) for full reproducibility
    with open(OUT / "fold_rotation_b2_raw.json", "w") as f:
        json.dump(b2_raw, f, indent=2, default=str)
    print(f"Wrote {OUT / 'fold_rotation_b2_raw.json'}")


if __name__ == "__main__":
    main()
