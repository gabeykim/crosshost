"""
GATE 4 - Task 4: CNN fold-variance spot check. Gate 3.5 measured fold-to-fold
vs draw-to-draw variance using the k-mer/linear substitute only -- the CNN's
own fold variance was never measured. This script measures it directly for
one representative configuration: B. subtilis held out, transcription,
genomic features, N=100 (the exact H-MAIN-relevant cell), reusing the
fine-tuning results already produced by scripts/43 (N=100, head_only
mechanism, 10 draws, all 5 folds) rather than re-running anything --
scripts/43's calibration grid already contains everything this check needs.

Compares CNN fold-to-fold std (across the 5 folds' per-fold mean-of-10-draws)
against the k-mer model's fold-to-fold std for the identical configuration
(BS, transcription, N=100), read from Gate 3.5's
out/baselines/fold_rotation_all_baselines.json.
"""
import json
import numpy as np
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "out"


def main():
    with open(OUT / "gate4_calibration_curves.json") as f:
        cal = json.load(f)

    cnn_fold_means = []
    cnn_draw_stds = []
    for test_fold_str, cell_data in cal["BS"]["genomic"].items():
        draws = cell_data["N100_head_only"]
        rhos = [d["eval"]["transcription"]["spearman_rho"] for d in draws
                if d["eval"]["transcription"]["spearman_rho"] is not None]
        if rhos:
            cnn_fold_means.append(float(np.mean(rhos)))
            cnn_draw_stds.append(float(np.std(rhos)))

    cnn_fold_to_fold_std = float(np.std(cnn_fold_means)) if len(cnn_fold_means) > 1 else None
    cnn_draw_to_draw_std_avg = float(np.mean(cnn_draw_stds)) if cnn_draw_stds else None

    with open(OUT / "baselines" / "fold_rotation_all_baselines.json") as f:
        kmer = json.load(f)
    kmer_vc = kmer["b2_calibration"]["variance_comparison"]["BS"]["transcription"]["100"]

    result = {
        "configuration": "BS held out, transcription, genomic features, N=100",
        "cnn": {
            "fold_means": cnn_fold_means,
            "fold_to_fold_std": cnn_fold_to_fold_std,
            "mean_within_fold_draw_to_draw_std": cnn_draw_to_draw_std_avg,
            "ratio_fold_std_over_draw_std": (cnn_fold_to_fold_std / cnn_draw_to_draw_std_avg
                                              if cnn_draw_to_draw_std_avg and cnn_draw_to_draw_std_avg > 1e-9 else None),
        },
        "kmer_linear_baseline_gate3_5": kmer_vc,
    }

    cnn_ratio = result["cnn"]["ratio_fold_std_over_draw_std"]
    kmer_ratio = kmer_vc["ratio_fold_std_over_draw_std"]
    if cnn_ratio is not None and kmer_ratio is not None:
        if cnn_fold_to_fold_std > 2 * kmer_vc["fold_to_fold_std_of_rho"]:
            result["decision_needed"] = True
            result["decision_needed_reason"] = (
                f"CNN fold-to-fold std ({cnn_fold_to_fold_std:.4f}) is more than 2x the k-mer model's "
                f"({kmer_vc['fold_to_fold_std_of_rho']:.4f}) for the identical configuration -- the H-MAIN "
                f"decision rule (5 folds, t-interval) may be under-powered for the CNN specifically.")
        else:
            result["decision_needed"] = False
            result["decision_needed_reason"] = "CNN fold variance is comparable to or smaller than the k-mer model's for this configuration."
    else:
        result["decision_needed"] = None
        result["decision_needed_reason"] = "Could not compute both ratios (insufficient valid draws)."

    print(json.dumps(result, indent=2, default=str))
    with open(OUT / "gate4_cnn_fold_variance_spotcheck.json", "w") as f:
        json.dump(result, f, indent=2, default=str)
    print(f"\nWrote {OUT / 'gate4_cnn_fold_variance_spotcheck.json'}")


if __name__ == "__main__":
    main()
