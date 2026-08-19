"""GATE 10 - required figures: 10x10 correlation matrix, modality comparison, GC-vs-activity."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"
FIGS = OUT / "figures"
FIGS.mkdir(parents=True, exist_ok=True)

PHYLO_ORDER = ["Ec", "Ef", "Se", "Ko", "Pa", "Pp", "Vn", "Bs", "Ll", "Cg"]
SPECIES_FULL = {
    "Ec": "E. coli", "Ef": "E. fergusonii", "Se": "S. enterica", "Ko": "K. oxytoca",
    "Pa": "P. agglomerans", "Pp": "P. putida", "Vn": "V. natriegens", "Bs": "B. subtilis",
    "Cg": "C. glutamicum", "Ll": "L. lactis",
}


def fig_correlation_matrix():
    m = pd.read_csv(RESULTS / "gate10_crosshost_correlation_matrix_phylo_ordered.csv", index_col=0)
    labels = [SPECIES_FULL[c] for c in m.columns]
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(m.values, vmin=0.5, vmax=1.0, cmap="viridis")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=9)
    for i in range(len(labels)):
        for j in range(len(labels)):
            v = m.values[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    color="white" if v < 0.8 else "black", fontsize=7)
    ax.axvline(6.5, color="red", linewidth=1.5)
    ax.axhline(6.5, color="red", linewidth=1.5)
    ax.axvline(8.5, color="orange", linewidth=1)
    ax.axhline(8.5, color="orange", linewidth=1)
    ax.set_title("DRAFTS cell-free cross-host Spearman correlation (all 45 pairs)\n"
                 "red line: Proteobacteria | Firmicutes+Actinobacteria boundary\n"
                 "Contrast with Johns et al. in-vivo: EC-PA rho~0.75, any BS pair rho~0.16-0.26")
    fig.colorbar(im, ax=ax, label="Spearman rho")
    fig.tight_layout()
    fig.savefig(FIGS / "gate10_crosshost_correlation_matrix.png", dpi=150)
    plt.close(fig)
    print(f"Wrote {FIGS / 'gate10_crosshost_correlation_matrix.png'}")


def fig_modality_comparison():
    d = json.load(open(RESULTS / "gate10_modality_comparison.json"))
    a = d["comparison_a_own_computation"]
    b = d["comparison_b_published_and_recomputed_rs234"]["b2_recomputed_from_drafts_own_released_data_NOT_a_published_quote"]["per_species"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel 1: EC-BS cross-host correlation, cell-free vs in-vivo
    labels = ["Cell-free\n(DRAFTS)\nEC-BS", "In-vivo (Johns)\nEC-BS\n(112-seq overlap,\nn=15, unreliable)",
              "In-vivo (Johns)\nEC-BS\n(full library,\nn=3668)"]
    vals = [a["cell_free_DRAFTS_EC_BS"]["spearman_rho"], a["in_vivo_Johns_EC_BS_on_same_112seq_overlap"]["spearman_rho"],
            a["in_vivo_Johns_EC_BS_full_library_reference"]["spearman_rho"]]
    colors = ["#2b6cb0", "#cbd5e0", "#c53030"]
    axes[0].bar(labels, vals, color=colors)
    axes[0].set_ylabel("Spearman rho, EC-BS cross-host correlation")
    axes[0].set_title("Comparison A: cell-free vs in-vivo,\nsame host pair (EC-BS)")
    axes[0].set_ylim(0, 0.7)
    for i, v in enumerate(vals):
        axes[0].text(i, v + 0.01, f"{v:.3f}", ha="center", fontsize=9)

    # Panel 2: RS234 per-species in-vitro vs in-vivo, DRAFTS's own paired design
    species = list(b.keys())
    rhos = [b[s]["spearman_rho_invitro_vs_invivo"] for s in species]
    axes[1].bar([SPECIES_FULL[s] for s in species], rhos, color="#2f855a")
    axes[1].set_ylabel("Spearman rho, in-vitro vs in-vivo (same species)")
    axes[1].set_title("Comparison B: within-species modality agreement\n(RS234, DRAFTS's own paired design)")
    axes[1].set_ylim(0, 1.0)
    axes[1].tick_params(axis="x", rotation=45)
    for i, v in enumerate(rhos):
        axes[1].text(i, v + 0.02, f"{v:.2f}", ha="center", fontsize=8)

    fig.suptitle("Is the cell-free host effect the same object as the in-vivo host effect?", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIGS / "gate10_modality_comparison.png", dpi=150)
    plt.close(fig)
    print(f"Wrote {FIGS / 'gate10_modality_comparison.png'}")


def fig_gc_vs_activity():
    df = pd.read_parquet(DATA / "drafts.parquet")
    fig, axes = plt.subplots(2, 5, figsize=(18, 7), sharex=True)
    species = PHYLO_ORDER
    for ax, sp in zip(axes.flat, species):
        usable = df[f"{sp}_usable"]
        gc = df.loc[usable, "gc_pct"]
        log_tx = np.log10(df.loc[usable, f"{sp}_tx"].clip(lower=1e-6))
        ax.scatter(gc, log_tx, s=3, alpha=0.3, color="#2b6cb0")
        ax.set_title(SPECIES_FULL[sp], fontsize=10)
        ax.set_xlabel("source GC%")
        ax.set_ylabel("log10(tx)")
    fig.suptitle("Source-genome GC content vs. recipient-host transcription activity, all 10 species\n"
                 "(all show a strong negative relationship -- the GC/phylum confound, Task 4.4)", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIGS / "gate10_gc_vs_activity.png", dpi=150)
    plt.close(fig)
    print(f"Wrote {FIGS / 'gate10_gc_vs_activity.png'}")


if __name__ == "__main__":
    fig_correlation_matrix()
    fig_modality_comparison()
    fig_gc_vs_activity()
