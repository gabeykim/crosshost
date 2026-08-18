"""
GATE 7 - Task 1 consequence: regenerate the cross-host measurement
correlation figure via a committed script (the original
out/figures/gate5_5_crosshost_measurement_correlation.png had no producing
script -- flagged by scripts/audit_provenance.py). The underlying numbers
(out/results/gate5_5_crosshost_measurement_correlation.json) were
independently re-verified reproducible bit-for-bit via scripts/57 during
Gate 6 and again here; only the plotting code was missing, not the data.

This figure now carries the FULL interpretive weight Task 1 assigns it,
since the percent-of-ceiling metric is RETIRED (see scripts/66 and
out/GATE7_MEMO.md) -- raw cross-host measurement correlation stands alone
as the descriptive finding.
"""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "out"
FIGS = OUT / "figures"
FIGS.mkdir(parents=True, exist_ok=True)

PAIRS = ["EC_PA", "EC_BS", "BS_PA"]
PAIR_LABELS = {"EC_PA": "E. coli –\nP. aeruginosa", "EC_BS": "E. coli –\nB. subtilis", "BS_PA": "B. subtilis –\nP. aeruginosa"}


def main():
    corr = json.load(open(OUT / "results" / "gate5_5_crosshost_measurement_correlation.json"))
    fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
    for ax, readout in zip(axes, ["transcription", "translation"]):
        rhos = [corr[readout][p]["rho"] for p in PAIRS]
        ns = [corr[readout][p]["n"] for p in PAIRS]
        colors = ["#4c72b0" if p == "EC_PA" else "#c44e52" for p in PAIRS]
        bars = ax.bar(range(len(PAIRS)), rhos, color=colors, edgecolor="black", linewidth=0.5)
        for i, (rho, n) in enumerate(zip(rhos, ns)):
            ax.text(i, rho + 0.02, f"n={n}", ha="center", fontsize=8)
        ax.set_xticks(range(len(PAIRS)))
        ax.set_xticklabels([PAIR_LABELS[p] for p in PAIRS], fontsize=9)
        ax.set_ylim(0, 0.9)
        ax.set_title(readout)
        ax.axhline(0, color="gray", linewidth=0.5)
    axes[0].set_ylabel("Spearman rho, raw measurements\n(sequences active in both hosts)")
    fig.suptitle("Cross-host measurement correlation, computed directly from raw data\n"
                  "(no model involved) -- the standalone finding, percent-of-ceiling framing retired (Gate 7)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig(FIGS / "gate7_crosshost_measurement_correlation.png", dpi=150)
    print(f"Wrote {FIGS / 'gate7_crosshost_measurement_correlation.png'}")


if __name__ == "__main__":
    main()
