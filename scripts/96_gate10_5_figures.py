"""GATE 10.5 - figures: GC vs activity per host (in vivo), raw vs GC-controlled correlation matrices."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
g95 = import_module("95_gc_confound_invivo")

OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"
FIGS = OUT / "figures"
FIGS.mkdir(parents=True, exist_ok=True)

HOSTS = ["EC", "BS", "PA"]
HOST_FULL = {"EC": "E. coli", "BS": "B. subtilis", "PA": "P. aeruginosa"}


def fig_gc_vs_activity():
    df = g95.load_data()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    for ax, host in zip(axes, HOSTS):
        mask = df[f"{host}_tx_usable"] & df[f"{host}_tx_active"]
        ax.scatter(df.loc[mask, "gc_pct"], df.loc[mask, f"{host}_tx_norm"], s=2, alpha=0.15, color="#2b6cb0")
        ax.set_title(HOST_FULL[host], fontsize=11)
        ax.set_xlabel("source GC%")
        ax.set_ylabel("tx_norm" if host == "EC" else "")
    fig.suptitle("In-vivo (Johns et al.): source-genome GC% vs. transcription activity\n"
                 "(compare to DRAFTS cell-free, rho range -0.74 to -0.49 in all 10 hosts)", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIGS / "gate10_5_gc_vs_activity_invivo.png", dpi=150)
    plt.close(fig)
    print(f"Wrote {FIGS / 'gate10_5_gc_vs_activity_invivo.png'}")


def fig_correlation_comparison():
    d = json.load(open(RESULTS / "gate10_5_gc_confound.json"))
    t2 = d["task2_raw_vs_gc_controlled_correlations"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, readout in zip(axes, ["transcription", "translation"]):
        pairs = ["EC_BS", "EC_PA", "BS_PA"]
        raw = [t2[readout][p]["rho_raw"] for p in pairs]
        gc = [t2[readout][p]["rho_gc_controlled_residual"] for p in pairs]
        x = np.arange(len(pairs))
        w = 0.35
        ax.bar(x - w / 2, raw, width=w, label="raw", color="#c53030")
        ax.bar(x + w / 2, gc, width=w, label="GC-controlled", color="#2b6cb0")
        ax.set_xticks(x)
        ax.set_xticklabels(["EC-BS", "EC-PA", "BS-PA"])
        ax.set_ylabel("Spearman rho")
        ax.set_title(readout, fontsize=11)
        ax.set_ylim(0, 0.85)
        ax.legend(fontsize=9)
        for i, (rv, gv) in enumerate(zip(raw, gc)):
            ax.text(i - w / 2, rv + 0.02, f"{rv:.2f}", ha="center", fontsize=8)
            ax.text(i + w / 2, gv + 0.02, f"{gv:.2f}", ha="center", fontsize=8)
    fig.suptitle("In-vivo cross-host correlation: raw vs. GC-controlled\n"
                 "EC-PA vs. BS-pairs contrast survives GC control in both readouts", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIGS / "gate10_5_raw_vs_gc_controlled.png", dpi=150)
    plt.close(fig)
    print(f"Wrote {FIGS / 'gate10_5_raw_vs_gc_controlled.png'}")


def fig_phylum_stratified():
    d = json.load(open(RESULTS / "gate10_5_gc_confound.json"))
    t4 = d["task4_phylum_stratified"]["transcription"]
    phyla = ["Proteobacteria", "Firmicutes", "Actinobacteria", "Bacteroidetes"]
    pairs = ["EC_BS", "EC_PA", "BS_PA"]

    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.arange(len(phyla))
    w = 0.25
    colors = {"EC_BS": "#c53030", "EC_PA": "#2f855a", "BS_PA": "#805ad5"}
    for i, pair in enumerate(pairs):
        vals = [t4[pair].get(ph, {}).get("spearman_rho") for ph in phyla]
        vals = [v if v is not None else 0 for v in vals]
        ax.bar(x + (i - 1) * w, vals, width=w, label=pair.replace("_", "-"), color=colors[pair])
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(phyla, rotation=15)
    ax.set_ylabel("Spearman rho (transcription)")
    ax.set_title("Cross-host correlation WITHIN each source-genome phylum stratum\n"
                 "(controls for phylum composition without any functional-form assumption)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGS / "gate10_5_phylum_stratified.png", dpi=150)
    plt.close(fig)
    print(f"Wrote {FIGS / 'gate10_5_phylum_stratified.png'}")


if __name__ == "__main__":
    fig_gc_vs_activity()
    fig_correlation_comparison()
    fig_phylum_stratified()
