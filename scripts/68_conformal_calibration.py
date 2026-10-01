"""
GATE 7 - Task 2: split-conformal prediction intervals (80%/90% target
coverage) and classifier calibration (reliability diagrams, ECE) for the
sequence-only model and the B2 per-host-N baseline. Scoped down per
CHARTER_AMENDMENTS.md's post-Gate-5 revision -- no deep-ensemble program,
just honest calibration figures on the baseline suite a benchmark should
ship regardless of the demonstration model's own outcome.

METHOD -- split conformal (Lei et al. 2018 / Vovk et al. 2005), the standard,
distribution-free approach: for a calibration set of size n with residuals
|y_i - pred_i|, the (1-alpha) prediction interval half-width is the
ceil((n+1)(1-alpha))/n empirical quantile of those residuals (the standard
finite-sample correction, not the naive quantile). Applied symmetrically
around each test point's prediction. Coverage is then checked empirically
on a DISJOINT test set never used for calibration.

SEQUENCE-ONLY: uses the already-trained, already-saved LOHO fold checkpoints
(out/models/seqonly_loho_{host}_fold{f}.pt) -- no retraining needed. For
held-out host X, fold f: the model was trained ONLY on the other two hosts
(any fold), so ALL of X's data is legitimately unseen. TEST = X's fold==f
rows (matching every other reported LOHO evaluation in this project, for
direct comparability). CALIBRATION = X's fold!=f rows (never trained on by
this fold's model either, and disjoint from TEST by construction).

B2 BASELINE: within-host by design (no cross-host structure to exploit).
TRAIN = N=3000 draw from host X's train pool (fold!=0, matching Gate 3's
TEST_FOLD=0 convention). CALIBRATION = the remainder of that train pool
(fold!=0, rows NOT drawn into the N=3000 training sample). TEST = fold==0.
N=3000 chosen as the representative "practitioner has measured a lot of
parts" scenario (matches H-MAIN's baseline arm).
"""
import sys
import json
import numpy as np
import pandas as pd
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.linear_model import LogisticRegression, Ridge
import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
m40 = import_module("40_film_cnn_model")
so = import_module("58_sequence_only_model")
b2 = import_module("31_baseline2_calibration")

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"
FIGS = OUT / "figures"
MODELS_DIR = OUT / "models"
RESULTS.mkdir(parents=True, exist_ok=True)
FIGS.mkdir(parents=True, exist_ok=True)

HOSTS = ["EC", "BS", "PA"]
READOUTS = ["transcription", "translation"]
N_FOLDS = 5
TEST_FOLD = 0
N_TRAIN_B2 = 3000
ALPHAS = {0.80: 0.20, 0.90: 0.10}  # target coverage -> alpha
N_ECE_BINS = 10


def conformal_halfwidth(calib_resid, target_coverage):
    alpha = ALPHAS[target_coverage]
    n = len(calib_resid)
    if n == 0:
        return None
    q_level = min(1.0, np.ceil((n + 1) * (1 - alpha)) / n)
    return float(np.quantile(calib_resid, q_level))


def ece(probs, labels, n_bins=N_ECE_BINS):
    bins = np.linspace(0, 1, n_bins + 1)
    bin_ids = np.digitize(probs, bins[1:-1])
    total_n = len(probs)
    ece_val = 0.0
    bin_stats = []
    for b in range(n_bins):
        mask = bin_ids == b
        n_b = mask.sum()
        if n_b == 0:
            bin_stats.append({"bin_lo": bins[b], "bin_hi": bins[b + 1], "n": 0, "mean_pred": None, "empirical_rate": None})
            continue
        mean_pred = float(probs[mask].mean())
        emp_rate = float(labels[mask].mean())
        ece_val += (n_b / total_n) * abs(mean_pred - emp_rate)
        bin_stats.append({"bin_lo": float(bins[b]), "bin_hi": float(bins[b + 1]), "n": int(n_b),
                           "mean_pred": mean_pred, "empirical_rate": emp_rate})
    return float(ece_val), bin_stats


READOUT_SPEC = [("transcription", "tx_active", "tx_strength", "tx"), ("translation", "tl_active", "tl_strength", "tl")]


def seqonly_conformal():
    """One forward pass per fold (both readouts share the same model output
    dict, scripts/58's SequenceOnlyCNN.forward returns all 4 heads at once)
    -- the first version of this script wastefully reloaded and re-ran
    inference per readout (2x the necessary model loads/forward passes),
    found and fixed after the process ran far longer than expected on CPU.

    DEVICE: MPS, now that predict_seqonly chunks internally at 1024 rows
    (scripts/58, Gate 7 fix). The earlier kIOGPUCommandBufferCallbackErrorOutOfMemory
    was caused by the UNBATCHED path trying to allocate the entire ~23,232-row
    calibration set on MPS in one shot; chunked at 1024 rows, benchmarked at
    ~1.3s regardless of row count in this range (scripts/audit_batching_fix.py)
    -- ~15,000x faster than the original unbatched-on-CPU estimate. Explicit
    del + torch.mps.empty_cache() per fold as a defensive margin on this 8GB
    unified-memory machine, since 15 models are loaded across the loop."""
    device = "mps" if torch.backends.mps.is_available() else "cpu"
    results = {}
    reliability_pooled = {h: {r: {"probs": [], "labels": []} for r in READOUTS} for h in HOSTS}
    for host in HOSTS:
        results[host] = {r: {} for r in READOUTS}
        per_fold_cov = {r: {0.80: [], 0.90: []} for r in READOUTS}
        d = dict(np.load(DATA / "baseline_cache" / f"{host}_baseline_data.npz", allow_pickle=True))
        targets, masks = m40.build_target_arrays(d)
        fold = d["fold"]
        for f in range(N_FOLDS):
            model = so.SequenceOnlyCNN().to(device)
            model.load_state_dict(torch.load(MODELS_DIR / f"seqonly_loho_{host}_fold{f}.pt", map_location=device))
            calib_idx = np.where(fold != f)[0]
            test_idx = np.where(fold == f)[0]
            preds_calib = so.predict_seqonly(model, d["onehot"][calib_idx])
            preds_test = so.predict_seqonly(model, d["onehot"][test_idx])
            print(f"  seqonly {host} fold {f}: n_calib={len(calib_idx)}, n_test={len(test_idx)}, device={device}")

            for readout, active_key, strength_key, active_prefix in READOUT_SPEC:
                calib_strength_mask = masks[strength_key][calib_idx]
                test_strength_mask = masks[strength_key][test_idx]
                y_calib = targets[strength_key][calib_idx][calib_strength_mask]
                pred_calib = preds_calib[f"{active_prefix}_strength"][calib_strength_mask]
                y_test = targets[strength_key][test_idx][test_strength_mask]
                pred_test = preds_test[f"{active_prefix}_strength"][test_strength_mask]
                resid_calib = np.abs(y_calib - pred_calib)

                for target_cov in [0.80, 0.90]:
                    hw = conformal_halfwidth(resid_calib, target_cov)
                    if hw is None or len(y_test) == 0:
                        continue
                    covered = np.abs(y_test - pred_test) <= hw
                    per_fold_cov[readout][target_cov].append({"fold": f, "empirical_coverage": float(covered.mean()),
                                                       "n_test": int(len(y_test)), "n_calib": int(len(y_calib)), "halfwidth": hw})

                test_active_mask = masks[active_key][test_idx]
                probs = 1 / (1 + np.exp(-preds_test[f"{active_prefix}_active_logit"][test_active_mask]))
                labels = targets[active_key][test_idx][test_active_mask]
                reliability_pooled[host][readout]["probs"].append(probs)
                reliability_pooled[host][readout]["labels"].append(labels)

            del model
            if device == "mps":
                torch.mps.empty_cache()

        for readout, *_ in READOUT_SPEC:
            results[host][readout] = {
                str(tc): {"per_fold": per_fold_cov[readout][tc],
                          "pooled_empirical_coverage": float(np.mean([x["empirical_coverage"] for x in per_fold_cov[readout][tc]])) if per_fold_cov[readout][tc] else None}
                for tc in [0.80, 0.90]
            }
    return results, reliability_pooled


def b2_conformal():
    results = {}
    reliability_pooled = {h: {r: {"probs": [], "labels": []} for r in READOUTS} for h in HOSTS}
    for host in HOSTS:
        results[host] = {}
        d = dict(np.load(DATA / "baseline_cache" / f"{host}_baseline_data.npz", allow_pickle=True))
        kmer = b2.build_kmer_features(d["onehot"])
        fold = d["fold"]
        targets, masks = m40.build_target_arrays(d)
        train_pool_idx = np.where(fold != TEST_FOLD)[0]
        test_idx = np.where(fold == TEST_FOLD)[0]

        rng = np.random.default_rng(42)
        n_sample = min(N_TRAIN_B2, len(train_pool_idx))
        train_sub = rng.choice(train_pool_idx, size=n_sample, replace=False)
        calib_idx = np.setdiff1d(train_pool_idx, train_sub)

        for readout, active_key, strength_key in [
            ("transcription", "tx_active", "tx_strength"),
            ("translation", "tl_active", "tl_strength"),
        ]:
            y_active_train = targets[active_key][train_sub]
            clf = LogisticRegression(C=b2.LOGREG_C, max_iter=300).fit(kmer[train_sub], y_active_train)

            reg_train_mask = masks[strength_key][train_sub]
            reg = None
            if reg_train_mask.sum() > 5 and len(set(np.round(targets[strength_key][train_sub][reg_train_mask], 6))) > 1:
                reg = Ridge(alpha=b2.RIDGE_ALPHA).fit(kmer[train_sub][reg_train_mask], targets[strength_key][train_sub][reg_train_mask])

            per_fold_cov = {0.80: [], 0.90: []}
            if reg is not None:
                calib_strength_mask = masks[strength_key][calib_idx]
                test_strength_mask = masks[strength_key][test_idx]
                y_calib = targets[strength_key][calib_idx][calib_strength_mask]
                pred_calib = reg.predict(kmer[calib_idx][calib_strength_mask])
                y_test = targets[strength_key][test_idx][test_strength_mask]
                pred_test = reg.predict(kmer[test_idx][test_strength_mask])
                resid_calib = np.abs(y_calib - pred_calib)
                for target_cov in [0.80, 0.90]:
                    hw = conformal_halfwidth(resid_calib, target_cov)
                    if hw is None or len(y_test) == 0:
                        continue
                    covered = np.abs(y_test - pred_test) <= hw
                    per_fold_cov[target_cov].append({"fold": TEST_FOLD, "empirical_coverage": float(covered.mean()),
                                                       "n_test": int(len(y_test)), "n_calib": int(len(y_calib)), "halfwidth": hw})

            test_active_mask = masks[active_key][test_idx]
            probs = clf.predict_proba(kmer[test_idx][test_active_mask])[:, 1]
            labels = targets[active_key][test_idx][test_active_mask]
            reliability_pooled[host][readout]["probs"].append(probs)
            reliability_pooled[host][readout]["labels"].append(labels)

            results[host][readout] = {
                str(tc): {"per_fold": per_fold_cov[tc],
                          "pooled_empirical_coverage": float(np.mean([x["empirical_coverage"] for x in per_fold_cov[tc]])) if per_fold_cov[tc] else None}
                for tc in [0.80, 0.90]
            }
    return results, reliability_pooled


def compute_reliability_and_ece(reliability_pooled, system_name):
    out = {}
    for host in HOSTS:
        out[host] = {}
        for readout in READOUTS:
            probs = np.concatenate(reliability_pooled[host][readout]["probs"]) if reliability_pooled[host][readout]["probs"] else np.array([])
            labels = np.concatenate(reliability_pooled[host][readout]["labels"]) if reliability_pooled[host][readout]["labels"] else np.array([])
            if len(probs) == 0:
                out[host][readout] = {"ece": None, "n": 0, "bins": []}
                continue
            ece_val, bins = ece(probs, labels)
            out[host][readout] = {"ece": ece_val, "n": int(len(probs)), "bins": bins}
    return out


def fig_reliability_diagrams(seqonly_rel, b2_rel):
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    for col, host in enumerate(HOSTS):
        for row, readout in enumerate(READOUTS):
            ax = axes[row, col]
            for rel, label, color in [(seqonly_rel, "sequence-only", "#4c72b0"), (b2_rel, "B2 (N=3000)", "#dd8452")]:
                bins = rel[host][readout]["bins"]
                xs = [b["mean_pred"] for b in bins if b["n"] > 0]
                ys = [b["empirical_rate"] for b in bins if b["n"] > 0]
                ax.plot(xs, ys, "o-", color=color, label=f"{label} (ECE={rel[host][readout]['ece']:.3f})" if rel[host][readout]["ece"] is not None else label)
            ax.plot([0, 1], [0, 1], "k--", linewidth=0.8, label="perfect calibration")
            ax.set_xlim(0, 1); ax.set_ylim(0, 1)
            ax.set_title(f"{host} {readout}", fontsize=10)
            if row == 1:
                ax.set_xlabel("mean predicted P(active)")
            if col == 0:
                ax.set_ylabel("empirical active rate")
            ax.legend(fontsize=6, loc="upper left")
    fig.suptitle("Reliability diagrams, zero-shot active/inactive classifier", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIGS / "gate7_reliability_diagrams.png", dpi=150)
    print(f"Wrote {FIGS / 'gate7_reliability_diagrams.png'}")


def fig_coverage(seqonly_cov, b2_cov):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for ax, target_cov in zip(axes, [0.80, 0.90]):
        x = np.arange(len(HOSTS) * len(READOUTS))
        labels = [f"{h}\n{'tx' if r == 'transcription' else 'tl'}" for h in HOSTS for r in READOUTS]
        seq_vals = [seqonly_cov[h][r][str(target_cov)]["pooled_empirical_coverage"] for h in HOSTS for r in READOUTS]
        b2_vals = [b2_cov[h][r][str(target_cov)]["pooled_empirical_coverage"] for h in HOSTS for r in READOUTS]
        width = 0.35
        ax.bar(x - width / 2, seq_vals, width, label="sequence-only", color="#4c72b0", edgecolor="black", linewidth=0.5)
        ax.bar(x + width / 2, b2_vals, width, label="B2 (N=3000)", color="#dd8452", edgecolor="black", linewidth=0.5)
        ax.axhline(target_cov, color="red", linestyle="--", linewidth=1, label=f"target ({int(target_cov*100)}%)")
        ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
        ax.set_title(f"Target coverage {int(target_cov*100)}%")
        ax.set_ylim(0, 1.05)
    axes[0].set_ylabel("empirical coverage")
    axes[1].legend(fontsize=8)
    fig.suptitle("Split-conformal empirical coverage vs target", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(FIGS / "gate7_conformal_coverage.png", dpi=150)
    print(f"Wrote {FIGS / 'gate7_conformal_coverage.png'}")


def main():
    print("Running sequence-only conformal calibration...")
    seqonly_cov, seqonly_rel_raw = seqonly_conformal()
    print("Running B2 (N=3000) conformal calibration...")
    b2_cov, b2_rel_raw = b2_conformal()

    seqonly_rel = compute_reliability_and_ece(seqonly_rel_raw, "sequence_only")
    b2_rel = compute_reliability_and_ece(b2_rel_raw, "b2")

    print("\n=== COVERAGE SUMMARY ===")
    rows = []
    for system_name, cov in [("sequence_only", seqonly_cov), ("b2_n3000", b2_cov)]:
        for host in HOSTS:
            for readout in READOUTS:
                for tc in [0.80, 0.90]:
                    v = cov[host][readout][str(tc)]
                    print(f"  {system_name} {host} {readout} target={tc}: empirical={v['pooled_empirical_coverage']}")
                    rows.append({"system": system_name, "host": host, "readout": readout, "target_coverage": tc,
                                 "empirical_coverage": v["pooled_empirical_coverage"]})

    print("\n=== ECE SUMMARY ===")
    for system_name, rel in [("sequence_only", seqonly_rel), ("b2_n3000", b2_rel)]:
        for host in HOSTS:
            for readout in READOUTS:
                print(f"  {system_name} {host} {readout}: ECE={rel[host][readout]['ece']}, n={rel[host][readout]['n']}")

    output = {"seqonly_conformal_coverage": seqonly_cov, "b2_conformal_coverage": b2_cov,
              "seqonly_reliability": seqonly_rel, "b2_reliability": b2_rel}
    with open(RESULTS / "gate7_conformal_calibration.json", "w") as f:
        json.dump(output, f, indent=2, default=str)
    pd.DataFrame(rows).to_csv(RESULTS / "gate7_conformal_calibration.csv", index=False)
    print(f"\nWrote gate7_conformal_calibration.{{json,csv}}")

    fig_reliability_diagrams(seqonly_rel, b2_rel)
    fig_coverage(seqonly_cov, b2_cov)


if __name__ == "__main__":
    main()
