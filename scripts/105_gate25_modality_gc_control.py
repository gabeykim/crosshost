#!/usr/bin/env python3
"""Gate 25 Block C: GC partial-correlation control on the modality contrast.

Gate 10.5 applied a source-genome GC control to the in-vivo pairs of Section 3.1
but never to the 0.677-vs-0.258 cell-free/in-vivo comparison the title rests on.
This script applies it, using the same two methods Gate 10.5 used -- the closed-
form partial-Spearman formula and an independent rank-residual regression -- so
the numbers are directly comparable to the ones already in the manuscript.

Masks replicate scripts/98 (cell-free) and scripts/80 (in vivo) exactly; the
script asserts it reproduces the four published raw values before controlling for
anything, and refuses to continue if it does not.

No training, no manuscript edits.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "out" / "results"

INVIVO_COACTIVE_RETENTION = 3668 / 14088


def rho(x, y):
    return float(spearmanr(x, y).correlation)


def rank_residual(a, g):
    ra, rg = rankdata(a), rankdata(g)
    slope, intercept = np.polyfit(rg, ra, 1)
    return ra - (slope * rg + intercept)


def partial_formula(r_xy, r_xz, r_yz):
    denom = np.sqrt((1 - r_xz ** 2) * (1 - r_yz ** 2))
    return None if denom < 1e-9 else (r_xy - r_xz * r_yz) / denom


def controlled(x, y, g):
    """GC-controlled Spearman, both ways. Returns (formula, residual, r_xg, r_yg)."""
    r_xy, r_xg, r_yg = rho(x, y), rho(x, g), rho(y, g)
    f = partial_formula(r_xy, r_xg, r_yg)
    r = rho(rank_residual(x, g), rank_residual(y, g))
    return f, r, r_xg, r_yg


def cellfree_sets():
    d = pd.read_parquet(DATA / "drafts.parquet")
    gc = d["gc_pct"]
    out = {}
    for sp in ("Ec", "Bs"):
        reason = d.get(f"{sp}_unusable_reason", pd.Series([""] * len(d))).fillna("")
        usable = d[f"{sp}_usable"].astype(bool)          # RNA detected
        no_rna = reason == "no_RNA_counts"
        activity = d[f"{sp}_tx"].copy()
        activity[no_rna] = 0.0                            # pooled imputation
        out[sp] = {"usable": usable, "dna_adequate": usable | no_rna,
                   "tx": d[f"{sp}_tx"], "activity": activity}
    A, B = out["Ec"], out["Bs"]
    pooled = A["dna_adequate"] & B["dna_adequate"] & A["activity"].notna() & B["activity"].notna()
    coact = A["usable"] & B["usable"] & A["tx"].notna() & B["tx"].notna()
    return {
        "pooled": (A["activity"][pooled], B["activity"][pooled], gc[pooled]),
        "co_active": (A["tx"][coact], B["tx"][coact], gc[coact]),
    }


def invivo_sets():
    t = pd.read_parquet(DATA / "three_host.parquet")
    gc = t["gc_content"]
    pooled = t["EC_tx_usable"] & t["BS_tx_usable"]
    coact = pooled & t["EC_tx_active"] & t["BS_tx_active"]
    return {
        "pooled": (t["EC_tx_norm"][pooled], t["BS_tx_norm"][pooled], gc[pooled]),
        "co_active": (t["EC_tx_norm"][coact], t["BS_tx_norm"][coact], gc[coact]),
    }


EXPECTED = {("cell_free", "pooled"): (862, 0.616), ("cell_free", "co_active"): (807, 0.677),
            ("in_vivo", "pooled"): (14088, 0.655), ("in_vivo", "co_active"): (3668, 0.258)}


def main() -> int:
    sets = {"cell_free": cellfree_sets(), "in_vivo": invivo_sets()}
    rows, res = [], {}

    for modality, regimes in sets.items():
        for regime, (x, y, g) in regimes.items():
            n = len(x)
            raw = rho(x, y)
            exp_n, exp_rho = EXPECTED[(modality, regime)]
            if n != exp_n or abs(raw - exp_rho) > 0.0005:
                print(f"FAIL: {modality}/{regime} reproduces n={n}, rho={raw:.4f}; "
                      f"published n={exp_n}, rho={exp_rho}")
                return 1
            f, r, r_xg, r_yg = controlled(x, y, g)
            rows.append({"modality": modality, "regime": regime, "n": n,
                         "rho_raw": round(raw, 4),
                         "rho_gc_controlled_formula": round(f, 4),
                         "rho_gc_controlled_residual": round(r, 4),
                         "gc_vs_EC_activity": round(r_xg, 4),
                         "gc_vs_BS_activity": round(r_yg, 4),
                         "pct_change": round(100 * (f - raw) / raw, 1)})
            res[f"{modality}_{regime}"] = rows[-1]

    cf, iv = res["cell_free_co_active"], res["in_vivo_co_active"]
    cfp, ivp = res["cell_free_pooled"], res["in_vivo_pooled"]
    ratio_raw = cf["rho_raw"] / iv["rho_raw"]
    ratio_ctl = cf["rho_gc_controlled_formula"] / iv["rho_gc_controlled_formula"]
    ratio_pooled_raw = cfp["rho_raw"] / ivp["rho_raw"]
    ratio_pooled_ctl = cfp["rho_gc_controlled_formula"] / ivp["rho_gc_controlled_formula"]

    summary = {
        "question": "Does the cell-free/in-vivo modality contrast survive a source-genome GC "
                    "partial-correlation control?",
        "method": "Spearman partial correlation controlling for source-genome GC%, by the "
                  "closed-form formula and independently by rank-residual regression, the two "
                  "methods used in Gate 10.5.",
        "note_separate_populations": "The two modalities are measured on different sequence "
                                     "populations (DRAFTS 1,047-sequence library vs the Johns "
                                     "three-host library), so GC is controlled within each "
                                     "modality separately; this is not a paired control.",
        "co_active": {"cell_free": cf, "in_vivo": iv,
                      "ratio_raw": round(ratio_raw, 2),
                      "ratio_gc_controlled": round(ratio_ctl, 2)},
        "pooled": {"cell_free": cfp, "in_vivo": ivp,
                   "ratio_raw": round(ratio_pooled_raw, 2),
                   "ratio_gc_controlled": round(ratio_pooled_ctl, 2)},
        "verdict": ("SURVIVES" if ratio_ctl >= 2.0 else
                    "ATTENUATED" if ratio_ctl >= 1.5 else "DOES NOT SURVIVE"),
    }

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "gate25_modality_gc_control.csv", index=False)
    with open(RESULTS / "gate25_modality_gc_control.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    print(df.to_string(index=False))
    print()
    print(f"co-active ratio   raw {ratio_raw:.2f}x  ->  GC-controlled {ratio_ctl:.2f}x")
    print(f"pooled    ratio   raw {ratio_pooled_raw:.2f}x  ->  GC-controlled {ratio_pooled_ctl:.2f}x")
    print(f"VERDICT: {summary['verdict']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
