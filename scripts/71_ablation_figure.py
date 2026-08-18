"""GATE 7 - Task D deliverable figure: feature-group ablation deltas."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"
FIGS = OUT / "figures"

GROUP_LABELS = {"sigma_factor": "sigma-factor", "anti_SD": "anti-SD", "tAI_codon": "tAI/codon",
                "RNAP": "RNAP", "chaperone_heme": "chaperone/heme"}


def main():
    d = json.load(open(RESULTS / "gate7_host_feature_ablation.json"))
    summary = d["summary"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
    for ax, readout in zip(axes, ["transcription", "translation"]):
        base_std = summary["baseline"][readout]["std"]
        groups = list(GROUP_LABELS.keys())
        deltas = [summary["ablation_delta"][readout][g]["delta"] for g in groups]
        colors = ["#c44e52" if abs(v) >= base_std else "#4c72b0" for v in deltas]
        ax.bar(range(len(groups)), deltas, color=colors, edgecolor="black", linewidth=0.5)
        ax.axhline(0, color="black", linewidth=0.8)
        ax.axhspan(-base_std, base_std, color="gray", alpha=0.2, label=f"±1 baseline fold-std ({base_std:.3f})")
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels([GROUP_LABELS[g] for g in groups], rotation=30, ha="right", fontsize=9)
        ax.set_title(f"{readout}\nbaseline rho={summary['baseline'][readout]['mean']:.3f}")
        ax.legend(fontsize=8)
    axes[0].set_ylabel("delta rho (ablated - baseline)")
    fig.suptitle("Gate 7 Task 3: genomic feature-group ablation, B. subtilis held out, 5 folds\n"
                  "(mean-ablation at inference time; shaded band = noise floor)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig(FIGS / "gate7_feature_ablation.png", dpi=150)
    print(f"Wrote {FIGS / 'gate7_feature_ablation.png'}")


if __name__ == "__main__":
    main()
