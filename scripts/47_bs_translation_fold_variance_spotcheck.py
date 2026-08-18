"""
GATE 4.5 - Task 3.2: CNN fold-variance spot-check for B. subtilis TRANSLATION
(the other H-MAIN arm -- Gate 4 only checked transcription). Same protocol as
Gate 4's Task 4 (scripts/44): B. subtilis held out, genomic features, N=100,
frozen-trunk (head_only) mechanism -- the pre-registered primary. Uses the
existing Gate 4 calibration grid (genomic variant, unaffected by Gate 4.5's
physiology retrain), no retraining needed.
"""
import json
import numpy as np
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "out"


def main():
    with open(OUT / "gate4_calibration_curves.json") as f:
        cal = json.load(f)

    fold_means, draw_stds = [], []
    for fold in range(5):
        draws = cal["BS"]["genomic"][str(fold)]["N100_head_only"]
        rhos = [d["eval"]["translation"]["spearman_rho"] for d in draws if d["eval"]["translation"]["spearman_rho"] is not None]
        if rhos:
            fold_means.append(float(np.mean(rhos)))
            draw_stds.append(float(np.std(rhos)))

    cnn_fold_std = float(np.std(fold_means))
    cnn_draw_std = float(np.mean(draw_stds))
    cnn_ratio = cnn_fold_std / cnn_draw_std if cnn_draw_std > 1e-9 else None

    with open(OUT / "baselines" / "fold_rotation_all_baselines.json") as f:
        kmer = json.load(f)
    kmer_vc = kmer["b2_calibration"]["variance_comparison"]["BS"]["translation"]["100"]

    with open(OUT / "gate4_cnn_fold_variance_spotcheck.json") as f:
        transcription_spotcheck = json.load(f)

    result = {
        "configuration": "BS held out, TRANSLATION, genomic features, N=100, head_only mechanism",
        "cnn": {
            "fold_means": fold_means,
            "fold_to_fold_std": cnn_fold_std,
            "mean_within_fold_draw_to_draw_std": cnn_draw_std,
            "ratio_fold_std_over_draw_std": cnn_ratio,
        },
        "kmer_linear_baseline_gate3_5": kmer_vc,
        "comparison_to_transcription_arm": {
            "transcription_cnn_ratio": transcription_spotcheck["cnn"]["ratio_fold_std_over_draw_std"],
            "translation_cnn_ratio": cnn_ratio,
            "translation_worse": cnn_ratio > transcription_spotcheck["cnn"]["ratio_fold_std_over_draw_std"],
        },
    }

    ratio_vs_kmer = cnn_ratio / kmer_vc["ratio_fold_std_over_draw_std"] if kmer_vc["ratio_fold_std_over_draw_std"] else None
    result["ratio_vs_kmer_baseline"] = ratio_vs_kmer
    result["decision_needed"] = ratio_vs_kmer is not None and ratio_vs_kmer > 2.0
    result["decision_needed_reason"] = (
        f"CNN's fold/draw ratio for BS translation ({cnn_ratio:.3f}) is {ratio_vs_kmer:.2f}x the k-mer "
        f"baseline's ({kmer_vc['ratio_fold_std_over_draw_std']:.3f}), and is itself HIGHER than the "
        f"transcription arm's ratio ({transcription_spotcheck['cnn']['ratio_fold_std_over_draw_std']:.3f}) "
        f"checked in Gate 4 -- the translation arm of H-MAIN, not just transcription, shows CNN fold "
        f"variance substantially exceeding baseline-level variance. Both H-MAIN arms now show elevated "
        f"fold variance; the pre-registered decision rule's robustness should be treated with this in mind."
        if ratio_vs_kmer and ratio_vs_kmer > 2.0 else
        "Ratio versus k-mer baseline within the 2x threshold."
    )

    print(json.dumps(result, indent=2))
    with open(OUT / "gate4_5_bs_translation_fold_variance_spotcheck.json", "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nWrote {OUT / 'gate4_5_bs_translation_fold_variance_spotcheck.json'}")


if __name__ == "__main__":
    main()
