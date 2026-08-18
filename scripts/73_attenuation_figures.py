"""GATE 8 - Task 1 deliverable figures: raw vs. disattenuated cross-host
correlation, and the >100%-ceiling anomaly resolution."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"
FIGS = OUT / "figures"

RELIABILITIES = [0.5, 0.7, 0.9]
PAIRS = ["EC_PA", "EC_BS", "BS_PA"]
PAIR_LABELS = {"EC_PA": "EC-PA", "EC_BS": "EC-BS", "BS_PA": "BS-PA"}


def fig_disattenuation():
    d = json.load(open(RESULTS / "gate8_attenuation_analysis.json"))
    scenario2 = d["scenario2_symmetric_sensitivity_grid"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), sharey=True)
    for ax, readout in zip(axes, ["transcription", "translation"]):
        x = np.arange(len(PAIRS))
        obs = [scenario2[readout][p]["rho_observed"] for p in PAIRS]
        ax.bar(x - 0.3, obs, width=0.15, label="observed", color="#555555", edgecolor="black")
        colors = ["#4c72b0", "#55a868", "#c44e52"]
        for i, rel in enumerate(RELIABILITIES):
            vals = [scenario2[readout][p][f"reliability_{rel}"]["rho_corrected"] for p in PAIRS]
            ax.bar(x - 0.3 + (i + 1) * 0.15, vals, width=0.15, label=f"corrected @ reliability={rel}", color=colors[i], edgecolor="black")
        ax.set_xticks(x)
        ax.set_xticklabels([PAIR_LABELS[p] for p in PAIRS])
        ax.set_title(readout)
        ax.set_ylim(0, 1.05)
    axes[0].set_ylabel("Spearman rho")
    axes[1].legend(fontsize=8, loc="upper right")
    fig.suptitle("Gate 8 Task 1: raw vs. disattenuation-corrected cross-host correlation\n"
                  "(EC-PA gap survives BS-pairs at every tested reliability)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig(FIGS / "gate8_disattenuation.png", dpi=150)
    print(f"Wrote {FIGS / 'gate8_disattenuation.png'}")


def fig_ceiling_resolution():
    d = json.load(open(RESULTS / "gate8_attenuation_analysis.json"))
    scenario2 = d["scenario2_symmetric_sensitivity_grid"]
    # BS ceiling = max(EC_BS, BS_PA); tx ceiling pair = EC_BS, tl ceiling pair = BS_PA (established Gate 7)
    seqonly_bs = {"transcription": 0.263, "translation": 0.291}  # Gate 5.5 sequence-only, verified Gate 7 script 66
    ceiling_pair = {"transcription": "EC_BS", "translation": "BS_PA"}

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    for ax, readout in zip(axes, ["transcription", "translation"]):
        pair = ceiling_pair[readout]
        obs_ceiling = scenario2[readout][pair]["rho_observed"]
        corrected_ceilings = [scenario2[readout][pair][f"reliability_{r}"]["rho_corrected"] for r in RELIABILITIES]
        x = list(range(len(RELIABILITIES) + 1))
        labels = ["observed\n(raw)"] + [f"corrected\n@ rel={r}" for r in RELIABILITIES]
        vals = [obs_ceiling] + corrected_ceilings
        colors = ["#c44e52"] + ["#4c72b0" if v >= seqonly_bs[readout] else "#c44e52" for v in corrected_ceilings]
        ax.bar(x, vals, color=colors, edgecolor="black")
        ax.axhline(seqonly_bs[readout], color="black", linestyle="--", linewidth=1.5,
                   label=f"sequence-only model rho ({seqonly_bs[readout]})")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_title(f"{readout}\nceiling pair: {PAIR_LABELS[pair]}")
        ax.set_ylabel("Spearman rho")
        ax.legend(fontsize=8)
    fig.suptitle("Gate 8 Task 1: the >100%-of-ceiling anomaly, resolved\n"
                  "(blue = corrected ceiling now at/above model performance; red = model still exceeds)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    fig.savefig(FIGS / "gate8_ceiling_anomaly_resolution.png", dpi=150)
    print(f"Wrote {FIGS / 'gate8_ceiling_anomaly_resolution.png'}")


if __name__ == "__main__":
    fig_disattenuation()
    fig_ceiling_resolution()
