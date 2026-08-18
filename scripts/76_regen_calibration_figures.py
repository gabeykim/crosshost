"""
GATE 8 - Task 2e: producing scripts for the calibration-curve figures Gate 7's
provenance triage flagged as FIGURE-REGENERABLE (underlying data reproducible,
plotting code not previously committed). Covers 10 of the 14: Gate 3's
single-fold curves (2), Gate 3.5's fold-resolved curves (2), Gate 4's
per-host genomic/physiology curves (6).
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "out"
N_VALUES = [0, 10, 30, 100, 300, 1000, 3000]
HOSTS = ["EC", "BS", "PA"]


def curve_from_single_fold(readout):
    """Gate 3's original single-fixed-test-fold B2 curve."""
    d = json.load(open(OUT / "baselines" / "baseline2_calibration_raw.json"))
    fig, ax = plt.subplots(figsize=(7, 5))
    for host in HOSTS:
        means, los, his = [], [], []
        for n in N_VALUES:
            draws = [x["spearman_rho"] for x in d[host]["results"][readout][str(n)] if x["spearman_rho"] is not None]
            if draws:
                means.append(np.mean(draws)); los.append(np.mean(draws) - np.std(draws)); his.append(np.mean(draws) + np.std(draws))
            else:
                means.append(np.nan); los.append(np.nan); his.append(np.nan)
        ax.errorbar(N_VALUES, means, yerr=[np.array(means) - np.array(los), np.array(his) - np.array(means)],
                    marker="o", label=host, capsize=3)
    ax.set_xscale("symlog")
    ax.set_xlabel("N calibration examples")
    ax.set_ylabel("Spearman rho")
    ax.set_title(f"Gate 3: B2 calibration curve, {readout} (single fixed test fold)")
    ax.legend()
    fig.tight_layout()
    out_path = OUT / "baselines" / f"calibration_curve_{readout}.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Wrote {out_path}")


def curve_fold_resolved(readout):
    """Gate 3.5's fold-rotated B2 curve, mean +/- fold-to-fold std."""
    d = json.load(open(OUT / "baselines" / "fold_rotation_b2_raw.json"))
    fig, ax = plt.subplots(figsize=(7, 5))
    for host in HOSTS:
        means, stds = [], []
        for n in N_VALUES:
            fold_means = []
            for f in range(5):
                draws = [x["spearman_rho"] for x in d[host][readout][str(n)][str(f)]["draws"] if x["spearman_rho"] is not None]
                if draws:
                    fold_means.append(np.mean(draws))
            means.append(np.mean(fold_means) if fold_means else np.nan)
            stds.append(np.std(fold_means) if fold_means else np.nan)
        ax.errorbar(N_VALUES, means, yerr=stds, marker="o", label=host, capsize=3)
    ax.set_xscale("symlog")
    ax.set_xlabel("N calibration examples")
    ax.set_ylabel("Spearman rho (fold-to-fold std shown)")
    ax.set_title(f"Gate 3.5: B2 calibration curve, {readout}, fold-resolved (5-fold rotation)")
    ax.legend()
    fig.tight_layout()
    out_path = OUT / "baselines" / f"calibration_curve_{readout}_fold_resolved.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Wrote {out_path}")


def gate4_curve(host, readout):
    """Gate 4's per-host genomic-vs-physiology calibration curve."""
    d = json.load(open(OUT / "gate4_calibration_curves.json"))
    d_loho = json.load(open(OUT / "gate4_loho_results.json"))
    cells = [(10, "head_only"), (30, "head_only"), (100, "head_only"),
             (300, "head_only"), (300, "top_conv"), (1000, "top_conv"), (3000, "top_conv")]
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for variant, color in [("genomic", "#4c72b0"), ("physiology", "#c44e52")]:
        xs, means = [0], []
        zeroshot = [d_loho[host][variant]["per_fold"][str(f)]["eval"][readout]["spearman_rho"] for f in range(5)]
        zeroshot = [v for v in zeroshot if v is not None]
        means.append(np.mean(zeroshot) if zeroshot else np.nan)
        for n, mech in cells:
            vals = []
            for f in range(5):
                draws = d[host][variant][str(f)][f"N{n}_{mech}"]
                vals.extend([x["eval"][readout]["spearman_rho"] for x in draws if x["eval"][readout]["spearman_rho"] is not None])
            if n not in xs:
                xs.append(n)
                means.append(np.mean(vals) if vals else np.nan)
        ax.plot(xs, means, marker="o", label=variant, color=color)
    ax.set_xscale("symlog")
    ax.set_xlabel("N calibration examples")
    ax.set_ylabel("Spearman rho")
    ax.set_title(f"Gate 4: {host} {readout}, genomic vs physiology calibration curve")
    ax.legend()
    fig.tight_layout()
    out_path = OUT / "figures" / f"calibration_curve_{host}_{readout}.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    for readout in ["transcription", "translation"]:
        curve_from_single_fold(readout)
        curve_fold_resolved(readout)
    for host in HOSTS:
        for readout in ["transcription", "translation"]:
            gate4_curve(host, readout)
