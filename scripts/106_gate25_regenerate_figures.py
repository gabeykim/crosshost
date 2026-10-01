#!/usr/bin/env python3
"""Gate 25 Block A: regenerate all eight main-text figures into out/PREPRINT/figures/.

Two things this fixes beyond the individual figure defects:

1. Figure 3 is rebuilt from the CURRENT analysis.  The committed file showed the
   Gate 10 value 0.597 -- DRAFTS EC-BS on the 82-sequence Johns overlap, computed
   before Gate 15 defined the restriction regimes -- which appears nowhere in the
   manuscript.  Panel A now shows the four values the text actually reports.

2. Until now NO script produced out/PREPRINT/figures/*.png.  The producing scripts
   write out/figures/gate*.png, and Gate 12 renamed copies into the preprint
   directory by hand.  audit_provenance.py passed only because the Gate 20 embed
   script mentions every preprint filename -- a false pass.  This script closes
   that: it runs each producer and copies its output to the preprint name, so the
   preprint figures have a real generating script.

Run: python3 scripts/106_gate25_regenerate_figures.py
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "out" / "results"
FIGS = ROOT / "out" / "figures"
PREPRINT = ROOT / "out" / "PREPRINT" / "figures"

# producer script -> {source filename in out/figures: preprint filename}
PRODUCERS = {
    "73_attenuation_figures": {"gate8_disattenuation.png": "Figure1_disattenuation.png"},
    "96_gate10_5_figures": {"gate10_5_raw_vs_gc_controlled.png": "Figure2a_gc_control.png",
                            "gate10_5_phylum_stratified.png": "Figure2b_phylum_stratified.png"},
    "55_gate5_figures": {"gate5_hmain_comparison.png": "Figure4_hmain_kill_gate.png"},
    "65_gate6_figures": {"gate6_full_comparison.png": "Figure5_conditioning_mechanisms.png"},
    "68_conformal_calibration": {"gate7_reliability_diagrams.png": "Figure6_calibration_collapse.png"},
    "87_gate8_6_figures": {"gate8_6_mechanism_fold_variance.png": "Figure7_film_instability.png",
                           "gate8_6_shift_regression_to_mean_control.png":
                               "Figure8_shift_prediction_retraction.png"},
}

RS234_FULL = {"Ec": "E. coli", "Se": "S. enterica", "Pp": "P. putida", "Vn": "V. natriegens",
              "Ko": "K. oxytoca", "Bs": "B. subtilis", "Cg": "C. glutamicum"}


def figure3():
    """Panel A: the four values Section 3.2 reports, under matched restriction.
    Panel B: DRAFTS's own within-species in-vitro/in-vivo agreement (RS234)."""
    gc = pd.read_csv(RESULTS / "gate25_modality_gc_control.csv")
    g = {(r.modality, r.regime): r for r in gc.itertuples()}
    order = [("cell_free", "pooled"), ("in_vivo", "pooled"),
             ("cell_free", "co_active"), ("in_vivo", "co_active")]
    vals = [g[k].rho_raw for k in order]
    ns = [g[k].n for k in order]
    labels = ["Cell-free\n(DRAFTS)\npooled", "In vivo\n(Johns)\npooled",
              "Cell-free\n(DRAFTS)\nco-active", "In vivo\n(Johns)\nco-active"]
    colors = ["#2b6cb0", "#c53030", "#2b6cb0", "#c53030"]

    b = json.load(open(RESULTS / "gate10_modality_comparison.json"))
    per_sp = (b["comparison_b_published_and_recomputed_rs234"]
               ["b2_recomputed_from_drafts_own_released_data_NOT_a_published_quote"]["per_species"])

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    bars = ax.bar(np.arange(4), vals, color=colors, edgecolor="black", width=0.62)
    for rect, v, n in zip(bars, vals, ns):
        ax.annotate(f"{v:.3f}\nn={n:,}", (rect.get_x() + rect.get_width() / 2, v),
                    ha="center", va="bottom", fontsize=9)
    ax.set_xticks(np.arange(4))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylim(0, 0.85)
    ax.set_ylabel("Spearman rho, E. coli-B. subtilis transcription")
    ax.set_title("A. Same pair, same restriction, two modalities\n"
                 "(pooled: no difference; co-active: 2.6-fold gap)", fontsize=11)

    ax = axes[1]
    sp = list(RS234_FULL)
    rv = [per_sp[s]["spearman_rho_invitro_vs_invivo"] for s in sp]
    rn = [per_sp[s]["n"] for s in sp]
    bars = ax.bar(np.arange(len(sp)), rv, color="#2f855a", edgecolor="black")
    for rect, v, n in zip(bars, rv, rn):
        ax.annotate(f"{v:.2f}\nn={n}", (rect.get_x() + rect.get_width() / 2, v),
                    ha="center", va="bottom", fontsize=7)
    ax.set_xticks(np.arange(len(sp)))
    ax.set_xticklabels([RS234_FULL[s] for s in sp], rotation=30, ha="right", fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Spearman rho, in vitro vs in vivo (same species)")
    ax.set_title("B. Within-species agreement across modalities\n"
                 "(RS234, DRAFTS's own paired design, 7 species)", fontsize=11)

    fig.suptitle("Cell-free reproduces within-host behaviour but not host-to-host difference",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = PREPRINT / "Figure3_drafts_modality_comparison.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Wrote {out}")
    return {"panel_a": dict(zip([f"{m}_{r}" for m, r in order], zip(vals, ns))),
            "panel_b_n_species": len(sp)}


def run_producer(mod_name: str) -> bool:
    """Run a producer as a subprocess so its __main__ block executes.

    Importing the module and calling main() does NOT work here: only one of the
    six producers defines main(), so an import-and-call runner silently does
    nothing for the other five and then copies their stale output -- which is
    exactly the false pass this script exists to remove.
    """
    path = ROOT / "scripts" / f"{mod_name}.py"
    r = subprocess.run([sys.executable, str(path)], cwd=ROOT,
                       capture_output=True, text=True)
    if r.returncode != 0:
        print(f"FAIL: {mod_name} exited {r.returncode}\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
        return False
    return True


def main() -> int:
    if not (RESULTS / "gate25_modality_gc_control.csv").is_file():
        print("FAIL: run scripts/105_gate25_modality_gc_control.py first")
        return 1

    started = time.time()
    copied = []
    for mod_name, mapping in PRODUCERS.items():
        print(f"--- {mod_name} ---")
        if not run_producer(mod_name):
            return 1
        for src, dst in mapping.items():
            s, d = FIGS / src, PREPRINT / dst
            if not s.is_file():
                print(f"FAIL: producer did not write {s}")
                return 1
            # Freshness assertion: a producer that silently did nothing leaves an
            # old file behind, and copying it would look like a successful rebuild.
            if s.stat().st_mtime < started:
                print(f"FAIL: {src} was not rewritten by {mod_name} (stale output)")
                return 1
            shutil.copyfile(s, d)
            copied.append(dst)
            print(f"  {src} -> {dst}")

    figure3()
    copied.append("Figure3_drafts_modality_comparison.png")

    expected = {f"Figure{i}" for i in range(1, 9)}
    got = {c.split("_")[0] for c in copied} | {"Figure2"}
    missing = expected - got
    if missing:
        print(f"FAIL: no output for {sorted(missing)}")
        return 1
    print(f"\nregenerated {len(copied)} files covering all 8 main figures")
    return 0


if __name__ == "__main__":
    sys.exit(main())
