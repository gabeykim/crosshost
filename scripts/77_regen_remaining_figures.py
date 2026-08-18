"""
GATE 8 - Task 2e: producing scripts for the remaining 3 of Gate 7's 14
FIGURE-REGENERABLE orphans. `gate5_5_performance_vs_ceiling.png` (the 14th)
is deliberately NOT regenerated -- it visualized the retired percent-of-
ceiling metric (out/GATE7_MEMO.md Task 1) and is superseded by
out/figures/gate8_ceiling_anomaly_resolution.png (scripts/73), which has
carried a committed script since Gate 8.
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"
FIGS = OUT / "figures"
HOSTS = ["EC", "BS", "PA"]
RS241_HOSTS = ["SE", "VN", "CG"]


def gate5_5_fourway_figure():
    df = pd.read_csv(RESULTS / "gate5_5_fourway_comparison.csv")
    systems = ["sequence_only", "free_embedding_B3", "genomic", "physiology"]
    colors = {"sequence_only": "#4c72b0", "free_embedding_B3": "#dd8452", "genomic": "#55a868", "physiology": "#c44e52"}
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for col, host in enumerate(HOSTS):
        for row, readout in enumerate(["transcription", "translation"]):
            ax = axes[row, col]
            sub = df[(df.host == host) & (df.readout == readout) & (df.eval_point == "N0_zeroshot")]
            xs = [s for s in systems if s in sub.system.values]
            means = [sub[sub.system == s].rho_mean.values[0] for s in xs]
            los = [sub[sub.system == s].rho_lower90.values[0] for s in xs]
            his = [sub[sub.system == s].rho_upper90.values[0] for s in xs]
            errs = [[m - l for m, l in zip(means, los)], [h - m for m, h in zip(means, his)]]
            ax.bar(range(len(xs)), means, yerr=errs, color=[colors[s] for s in xs], capsize=3, edgecolor="black")
            ax.set_xticks(range(len(xs)))
            ax.set_xticklabels([s.replace("_", "\n") for s in xs], fontsize=7)
            ax.set_title(f"{host} {readout}", fontsize=10)
            ax.axhline(0, color="gray", linewidth=0.5)
    fig.suptitle("Gate 5.5: four-way comparison, zero-shot LOHO, 90% bootstrap CI")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out_path = FIGS / "gate5_5_fourway_comparison.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Wrote {out_path}")


def gate5_5_rs241_fourway_figure():
    df = pd.read_csv(RESULTS / "gate5_5_rs241_fourway.csv")
    systems = ["seqonly", "genomic", "physiology"]
    colors = {"seqonly": "#4c72b0", "genomic": "#55a868", "physiology": "#c44e52"}
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for i, config in enumerate(["PRIMARY", "SECONDARY"]):
        for j, readout in enumerate(["transcription", "translation"]):
            ax = axes[j, i]
            sub = df[(df.config == config) & (df.readout == readout)]
            x = np.arange(len(RS241_HOSTS))
            width = 0.25
            for k, s in enumerate(systems):
                means = [sub[sub.host == h][f"{s}_mean"].values[0] for h in RS241_HOSTS]
                los = [sub[sub.host == h][f"{s}_lower90"].values[0] for h in RS241_HOSTS]
                his = [sub[sub.host == h][f"{s}_upper90"].values[0] for h in RS241_HOSTS]
                errs = [[m - l for m, l in zip(means, los)], [h - m for m, h in zip(means, his)]]
                ax.bar(x + (k - 1) * width, means, width=width, yerr=errs, label=s, color=colors[s], capsize=2, edgecolor="black")
            ax.set_xticks(x)
            ax.set_xticklabels(RS241_HOSTS)
            ax.set_title(f"{config} {readout}", fontsize=10)
            if i == 0 and j == 0:
                ax.legend(fontsize=8)
    fig.suptitle("Gate 5.5: RS241 four-way comparison (sequence-only/genomic/physiology)")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out_path = FIGS / "gate5_5_rs241_fourway.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Wrote {out_path}")


def gate5_rs241_zeroshot_figure():
    d = json.load(open(RESULTS / "gate5_rs241_results_summary.json"))
    perf = d["performance"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), sharey=True)
    for ax, config in zip(axes, ["PRIMARY", "SECONDARY"]):
        x = np.arange(len(RS241_HOSTS))
        width = 0.35
        for k, variant in enumerate(["genomic", "physiology"]):
            means = [perf[config][variant][h]["transcription"]["rho_ci"]["mean"] for h in RS241_HOSTS]
            los = [perf[config][variant][h]["transcription"]["rho_ci"]["lower"] for h in RS241_HOSTS]
            his = [perf[config][variant][h]["transcription"]["rho_ci"]["upper"] for h in RS241_HOSTS]
            errs = [[m - l for m, l in zip(means, los)], [h - m for m, h in zip(means, his)]]
            ax.bar(x + (k - 0.5) * width, means, width=width, yerr=errs, label=variant, capsize=3, edgecolor="black")
        ax.set_xticks(x)
        ax.set_xticklabels(RS241_HOSTS)
        ax.set_title(config)
        ax.legend(fontsize=9)
    axes[0].set_ylabel("Spearman rho (transcription, zero-shot)")
    fig.suptitle("Gate 5: RS241 zero-shot performance, genomic vs physiology, seed-resolved 90% CI")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    out_path = FIGS / "gate5_rs241_zeroshot.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    gate5_5_fourway_figure()
    gate5_5_rs241_fourway_figure()
    gate5_rs241_zeroshot_figure()
