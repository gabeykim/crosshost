"""
GATE 5 - figures: H-MAIN comparison (model vs baseline, 90% CI error bars),
genomic-vs-physiology (H-SCIENCE), performance vs extrapolation distance.
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)


def err(ci):
    return [[ci["mean"] - ci["lower"]], [ci["upper"] - ci["mean"]]]


def plot_hmain():
    with open(RESULTS / "gate5_hmain_results.json") as f:
        d = json.load(f)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, arm, baseline_key, baseline_label in [
        (axes[0], "H_MAIN_TX", "N3000", "baseline @ N=3000"),
        (axes[1], "H_MAIN_TL", None, "baseline saturation ceiling (N≈300)"),
    ]:
        data = d[arm]
        if arm == "H_MAIN_TX":
            baseline_ci = data["baseline"]["N3000"]
        else:
            baseline_ci = data["baseline_ceiling"]

        labels, means, lowers, uppers, colors = [], [], [], [], []
        labels.append(baseline_label)
        means.append(baseline_ci["mean"]); lowers.append(baseline_ci["mean"]-baseline_ci["lower"]); uppers.append(baseline_ci["upper"]-baseline_ci["mean"])
        colors.append("black")

        for variant, color in [("genomic", "tab:blue"), ("physiology", "tab:orange")]:
            for mech, marker_suffix in [("head_only", " (frozen, primary)"), ("top_conv", " (top-conv, suppl.)")]:
                cell = data["model"][variant][mech]
                ci = cell["model_ci"]
                labels.append(f"{variant}{marker_suffix}")
                means.append(ci["mean"]); lowers.append(ci["mean"]-ci["lower"]); uppers.append(ci["upper"]-ci["mean"])
                colors.append(color)

        y = np.arange(len(labels))
        ax.barh(y, means, xerr=[lowers, uppers], color=colors, alpha=0.75, capsize=4)
        ax.set_yticks(y)
        ax.set_yticklabels(labels)
        ax.axvline(baseline_ci["mean"], color="black", linestyle="--", alpha=0.5)
        ax.set_xlabel("Spearman ρ (90% bootstrap CI)")
        ax.set_title(f"{arm.replace('_','-')}: model @ N=100 vs {baseline_label}\n(B. subtilis held out)")
        ax.grid(alpha=0.3, axis="x")
    plt.tight_layout()
    out_path = FIG / "gate5_hmain_comparison.png"
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Wrote {out_path}")


def plot_hscience():
    with open(RESULTS / "gate5_hscience_results.json") as f:
        d = json.load(f)
    hosts = ["EC", "BS", "PA"]
    host_names = {"EC": "E. coli", "BS": "B. subtilis", "PA": "P. aeruginosa"}
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    x = np.arange(len(hosts))
    width = 0.35
    for ax, readout in zip(axes, ["transcription", "translation"]):
        g_means, g_err_lo, g_err_hi = [], [], []
        p_means, p_err_lo, p_err_hi = [], [], []
        for host in hosts:
            g = d["full_grid"][host][readout]["head_only"]["genomic"]["spearman_rho"]
            p = d["full_grid"][host][readout]["head_only"]["physiology"]["spearman_rho"]
            g_means.append(g["mean"]); g_err_lo.append(g["mean"]-g["lower"]); g_err_hi.append(g["upper"]-g["mean"])
            p_means.append(p["mean"]); p_err_lo.append(p["mean"]-p["lower"]); p_err_hi.append(p["upper"]-p["mean"])
        ax.bar(x - width/2, g_means, width, yerr=[g_err_lo, g_err_hi], capsize=4, label="genomic (37-D)", color="tab:blue")
        ax.bar(x + width/2, p_means, width, yerr=[p_err_lo, p_err_hi], capsize=4, label="physiology (6-D)", color="tab:orange")
        ax.set_xticks(x); ax.set_xticklabels([host_names[h] for h in hosts])
        ax.axhline(0, color="black", linewidth=0.8)
        ax.set_ylabel("Spearman ρ, N=100, frozen-trunk (90% CI)")
        ax.set_title(f"{readout.capitalize()}: genomic vs physiology (H-SCIENCE)")
        ax.legend(); ax.grid(alpha=0.3, axis="y")
    plt.tight_layout()
    out_path = FIG / "gate5_hscience_comparison.png"
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Wrote {out_path}")


def plot_extrapolation():
    with open(RESULTS / "gate5_hscience_results.json") as f:
        d = json.load(f)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for ax, variant in zip(axes, ["genomic", "physiology"]):
        for readout, marker in [("transcription", "o"), ("translation", "s")]:
            rows = d["performance_vs_extrapolation_distance"][variant][readout]["rows"]
            xs = [r["extrapolation_distance"] for r in rows]
            ys = [r["performance_rho"] for r in rows]
            labels = [r["host"] for r in rows]
            ax.scatter(xs, ys, marker=marker, s=90, label=readout)
            for x_, y_, lab in zip(xs, ys, labels):
                ax.annotate(lab, (x_, y_), textcoords="offset points", xytext=(6, 4), fontsize=9)
        ax.set_xlabel("Extrapolation distance (z-score units beyond training-pair range)")
        ax.set_ylabel("Spearman ρ, N=100, frozen-trunk")
        ax.set_title(f"{variant} features: performance vs extrapolation distance\n(n=3 hosts -- no statistical claim, descriptive only)")
        ax.legend(); ax.grid(alpha=0.3)
    plt.tight_layout()
    out_path = FIG / "gate5_extrapolation_distance.png"
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    plot_hmain()
    plot_hscience()
    plot_extrapolation()
