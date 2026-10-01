"""
GATE 8.6 - figures: (1) shift model vs reference-only vs shuffled-sequence,
per winning cell -- the regression-to-the-mean control's headline result.
(2) mechanism fold-to-fold variance, FiLM vs concat vs per-host-heads vs
sequence-only, EC transcription (the cell FiLM's instability was measured
on) plus all-cell summary.
"""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"
FIGS = OUT / "figures"
FIGS.mkdir(parents=True, exist_ok=True)


def fig_shift_controls():
    df = pd.read_csv(RESULTS / "gate8_6_shift_controls_decisive_comparison.csv")
    df = df.sort_values("effect_size_rho_orig_minus_refonly")
    # readout[:2] rendered BOTH "transcription" and "translation" as "tr", making two
    # distinct cells look like duplicates and the 11 groups look like 12 (Gate 25 A2).
    RD = {"transcription": "tx", "translation": "tl"}
    labels = [f"{r.pair}\n{RD[r.readout]} {r.variant.replace('_condition','')}" for r in df.itertuples()]
    x = np.arange(len(df))
    w = 0.27

    fig, ax = plt.subplots(figsize=(13, 6))
    ax.bar(x - w, df["refonly_rho_mean"], width=w, label="reference-value only (OLS)", color="#888888")
    ax.bar(x, df["orig_rho_mean"], width=w, label="sequence + reference (original)", color="#2b6cb0")
    ax.bar(x + w, df["shuffled_rho_mean"], width=w, label="shuffled-sequence + reference", color="#c05621")
    for i, row in enumerate(df.itertuples()):
        if row.beats_refonly_distinguishably:
            ax.annotate("survives", (i, max(row.orig_rho_mean, row.refonly_rho_mean) + 0.03),
                        ha="center", fontsize=9, fontweight="bold", color="#2f855a")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Spearman rho (test fold, mean of 5)")
    ax.set_title("Regression-to-the-mean control on the 11 originally-reported winning cells\n"
                 "Reference-value alone matches or beats the original model in 10 of 11; "
                 "1 of 11 survives")
    ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIGS / "gate8_6_shift_regression_to_mean_control.png", dpi=150)
    plt.close(fig)
    print(f"Wrote {FIGS / 'gate8_6_shift_regression_to_mean_control.png'}")


def fig_mechanism_variance():
    df = pd.read_csv(RESULTS / "gate8_6_mechanism_fold_variance.csv")
    mechs = ["sequence_only", "film_genomic", "concat", "perhost_heads_avg"]
    colors = {"sequence_only": "#718096", "film_genomic": "#c53030", "concat": "#2b6cb0", "perhost_heads_avg": "#2f855a"}
    labels = {"sequence_only": "sequence-only", "film_genomic": "FiLM", "concat": "concatenation", "perhost_heads_avg": "per-host heads (avg)"}

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    cells = df[["host", "readout"]].drop_duplicates().reset_index(drop=True)
    cell_labels = [f"{r.host}\n{ {'transcription': 'tx', 'translation': 'tl'}[r.readout] }"
                   for r in cells.itertuples()]
    x = np.arange(len(cells))
    w = 0.2
    for i, mech in enumerate(mechs):
        sub = df[df.mechanism == mech].set_index(["host", "readout"])
        vals = [sub.loc[(r.host, r.readout), "fold_to_fold_std"] if (r.host, r.readout) in sub.index else 0
                for r in cells.itertuples()]
        axes[0].bar(x + (i - 1.5) * w, vals, width=w, label=labels[mech], color=colors[mech])
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(cell_labels, fontsize=8)
    axes[0].set_ylabel("fold-to-fold std of zero-shot Spearman rho")
    axes[0].set_title("Fold-to-fold instability, all 4 systems, all 6 (host, readout) cells")
    axes[0].legend(fontsize=8)

    ec_tx = df[(df.host == "EC") & (df.readout == "transcription")].set_index("mechanism")
    for mech in mechs:
        axes[1].plot(range(5), [None] * 5, alpha=0)  # placeholder for spacing
    import json
    film_vals = json.load(open(OUT / "gate4_loho_results.json"))["EC"]["genomic"]["per_fold"]
    concat_vals = json.load(open(OUT / "gate8_5_concat_loho_results.json"))["EC"]["per_fold"]
    perhost_vals = json.load(open(OUT / "gate8_5_perhost_heads_loho_results.json"))["EC"]["per_fold"]
    seqonly_vals = json.load(open(OUT / "gate5_5_seqonly_loho_results.json"))["EC"]["per_fold"]
    folds = list(range(5))
    axes[1].plot(folds, [seqonly_vals[str(f)]["eval"]["transcription"]["spearman_rho"] for f in folds],
                 "o-", color=colors["sequence_only"], label=labels["sequence_only"])
    axes[1].plot(folds, [film_vals[str(f)]["eval"]["transcription"]["spearman_rho"] for f in folds],
                 "o-", color=colors["film_genomic"], label=labels["film_genomic"])
    axes[1].plot(folds, [concat_vals[str(f)]["eval"]["transcription"]["spearman_rho"] for f in folds],
                 "o-", color=colors["concat"], label=labels["concat"])
    axes[1].plot(folds, [perhost_vals[str(f)]["eval_avg"]["transcription"]["spearman_rho"] for f in folds],
                 "o-", color=colors["perhost_heads_avg"], label=labels["perhost_heads_avg"])
    axes[1].set_xlabel("fold")
    axes[1].set_ylabel("zero-shot Spearman rho")
    axes[1].set_title("EC transcription: per-fold trajectory\n(FiLM's swing is visible directly, not just its std)")
    axes[1].legend(fontsize=8)
    axes[1].set_xticks(folds)

    fig.tight_layout()
    fig.savefig(FIGS / "gate8_6_mechanism_fold_variance.png", dpi=150)
    plt.close(fig)
    print(f"Wrote {FIGS / 'gate8_6_mechanism_fold_variance.png'}")


if __name__ == "__main__":
    fig_shift_controls()
    fig_mechanism_variance()
