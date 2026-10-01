"""
GATE 6 - Task 3 deliverable figures: full model comparison (zero-shot LOHO,
all systems), and percentage-of-ceiling across all models (corrected, see
scripts/64 module docstring for why the previous gate5_5_per_host_ceiling.json
was wrong).
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
FIGS.mkdir(parents=True, exist_ok=True)

HOSTS = ["EC", "BS", "PA"]
READOUTS = ["transcription", "translation"]
# The three conditioning mechanisms Section 3.5 cites (concat, per-host-heads
# averaged and nearest) live in gate8_5_conditioning_mechanisms.csv, not in
# gate6_full_comparison.csv, so Figure 5 did not plot the values the text cited
# (Gate 26 T4). The two tables are on identical footing -- sequence_only 0.367130
# and film_genomic/genomic 0.371166 agree exactly at EC transcription, same 5
# folds, same 90% bootstrap scheme -- so they can share a panel.
SYSTEMS = ["sequence_only", "free_embedding_B3", "genomic", "physiology",
           "concat", "perhost_heads_avg", "perhost_heads_nearest",
           "dnabert2", "promogen2"]
MECHANISM_SYSTEMS = {"concat", "perhost_heads_avg", "perhost_heads_nearest"}
SYSTEM_LABELS = {"sequence_only": "Sequence-only\nCNN (214K)", "free_embedding_B3": "Free host\nembedding",
                  "genomic": "Genomic\n+ FiLM", "physiology": "Physiology\n+ FiLM",
                  "dnabert2": "DNABERT-2\n(117M)", "promogen2": "PromoGen2\n(148M)",
                  "concat": "Concat-\nenation", "perhost_heads_avg": "Per-host\nheads (avg)",
                  "perhost_heads_nearest": "Per-host\nheads (near)",
                  "evo2": "Evo 2\n(40B)"}
COLORS = {"sequence_only": "#4c72b0", "free_embedding_B3": "#dd8452", "genomic": "#55a868",
          "physiology": "#c44e52", "dnabert2": "#8172b3", "promogen2": "#937860", "evo2": "#ccb974",
          "concat": "#64b5cd", "perhost_heads_avg": "#8c8c8c", "perhost_heads_nearest": "#b5b5b5"}


def fig_full_comparison():
    df = pd.read_csv(RESULTS / "gate6_full_comparison.csv")
    mech = pd.read_csv(RESULTS / "gate8_5_conditioning_mechanisms.csv")
    mech = mech[mech.system.isin(MECHANISM_SYSTEMS)].copy()
    mech["eval_point"] = "N0_zeroshot"
    df = pd.concat([df, mech], ignore_index=True)
    systems = [s for s in SYSTEMS if s in df.system.unique()] + (["evo2"] if "evo2" in df.system.unique() else [])
    # One y-scale across all six panels (Gate 25 A3): per-panel scaling made weak
    # B. subtilis performance look comparable to E. coli, the inverse of the point.
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharey=True)
    for col, host in enumerate(HOSTS):
        for row, readout in enumerate(READOUTS):
            ax = axes[row, col]
            sub0 = df[(df.host == host) & (df.readout == readout) & (df.eval_point == "N0_zeroshot")]
            x = np.arange(len(systems))
            means = [sub0[sub0.system == s].rho_mean.values[0] if len(sub0[sub0.system == s]) else np.nan for s in systems]
            lowers = [sub0[sub0.system == s].rho_lower90.values[0] if len(sub0[sub0.system == s]) else np.nan for s in systems]
            uppers = [sub0[sub0.system == s].rho_upper90.values[0] if len(sub0[sub0.system == s]) else np.nan for s in systems]
            errs = [[m - l for m, l in zip(means, lowers)], [u - m for m, u in zip(means, uppers)]]
            colors = [COLORS[s] for s in systems]
            ax.bar(x, means, yerr=errs, color=colors, capsize=3, edgecolor="black", linewidth=0.5)
            ax.axhline(0, color="gray", linewidth=0.5)
            ax.set_xticks(x)
            ax.set_xticklabels([SYSTEM_LABELS[s] for s in systems], fontsize=6, rotation=30, ha="right")
            ax.set_title(f"{host} {readout}", fontsize=10)
            if col == 0:
                ax.set_ylabel("Spearman rho (N=0 zero-shot)")
    fig.suptitle("Full model comparison, zero-shot leave-one-host-out, 90% bootstrap CI\n"
                 "(one y-scale across all panels)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIGS / "gate6_full_comparison.png", dpi=150)
    print(f"Wrote {FIGS / 'gate6_full_comparison.png'}")


def fig_pct_ceiling():
    ceiling = json.load(open(RESULTS / "gate5_5_per_host_ceiling_CORRECTED.json"))
    systems = ["sequence_only", "dnabert2", "promogen2"]
    if any("evo2" in ceiling[h][r]["systems"] for h in HOSTS for r in READOUTS):
        systems.append("evo2")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for ax, readout in zip(axes, READOUTS):
        x = np.arange(len(HOSTS))
        width = 0.8 / len(systems)
        for i, s in enumerate(systems):
            vals = [ceiling[h][readout]["systems"].get(s, {}).get("pct_of_max_ceiling", np.nan) for h in HOSTS]
            ax.bar(x + i * width - 0.4 + width / 2, vals, width=width, label=SYSTEM_LABELS[s].replace("\n", " "),
                   color=COLORS[s], edgecolor="black", linewidth=0.5)
        ax.axhline(100, color="red", linestyle="--", linewidth=1, label="ceiling (100%)")
        ax.set_xticks(x)
        ax.set_xticklabels(HOSTS)
        ax.set_title(readout)
        ax.set_ylabel("% of cross-host measurement-correlation ceiling")
    axes[1].legend(fontsize=8, loc="upper right")
    fig.suptitle("Gate 6 (corrected): percent-of-ceiling by system -- NOT the previously-reported\n"
                  "'remarkably consistent ~50%' (see out/GATE6_MEMO.md for the correction)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig(FIGS / "gate6_pct_of_ceiling_corrected.png", dpi=150)
    print(f"Wrote {FIGS / 'gate6_pct_of_ceiling_corrected.png'}")


if __name__ == "__main__":
    fig_full_comparison()
    fig_pct_ceiling()
