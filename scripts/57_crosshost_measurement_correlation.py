"""
GATE 5.5 - Task 2.1 (PRIORITY): cross-host measurement correlation, computed
directly from the raw data -- no model involved. For sequences with usable,
active measurements in BOTH hosts of a pair, how correlated are the
measurements themselves? This is the ceiling any cross-host model (sequence-
only or host-conditioned) can plausibly reach without host-specific
information, and it should have been computed at Gate 2.

Uses the SAME usability/activity definitions as every model evaluation in
this project (Gate 3's floor-correction for translation, "Spearman on
actives" convention) so the ceiling is directly comparable to model rho
values reported elsewhere.
"""
import pandas as pd
import numpy as np
import json
from pathlib import Path
from scipy.stats import spearmanr

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

HOSTS = ["EC", "BS", "PA"]
PAIRS = [("EC", "BS"), ("EC", "PA"), ("BS", "PA")]


def main():
    df = pd.read_parquet(DATA / "three_host.parquet")

    # need floor-corrected translation active/usable-for-regression per host,
    # matching scripts/28's definitions exactly -- recompute here directly from
    # three_host.parquet using the same floor-detection logic (fold-local not
    # needed here since this is a DATA property, not a model-training statistic --
    # we are correlating measurements, not fitting anything)
    tl_active_col = {}
    tl_usable_reg_col = {}
    for host in HOSTS:
        usable = df[f"{host}_tl_usable"].values.astype(bool)
        protein = df[f"{host}_protein_log10"].values.astype(np.float32)
        vals = protein[usable]
        rounded = np.round(vals, 2)
        u, c = np.unique(rounded, return_counts=True)
        mode_val, mode_frac = u[np.argmax(c)], c[np.argmax(c)] / len(vals)
        is_floor = mode_frac > 0.05
        if is_floor:
            at_floor = usable & (np.round(protein, 2) == mode_val)
            tl_usable_reg_col[host] = usable & ~at_floor
            tl_active_col[host] = usable & (np.round(protein, 2) > mode_val)
        else:
            tl_usable_reg_col[host] = usable.copy()
            tl_active_col[host] = usable & (protein > np.median(vals))
        print(f"{host}: tl floor={mode_val:.2f} (frac={mode_frac:.3f}), tl_active_n={tl_active_col[host].sum()}")

    results = {"transcription": {}, "translation": {}}

    print("\n=== TRANSCRIPTION: cross-host measurement correlation (on actives in BOTH hosts) ===")
    for h1, h2 in PAIRS:
        both_active = (df[f"{h1}_tx_usable"] & df[f"{h1}_tx_active"] &
                       df[f"{h2}_tx_usable"] & df[f"{h2}_tx_active"]).values
        n = int(both_active.sum())
        v1 = df.loc[both_active, f"{h1}_tx_norm"].values
        v2 = df.loc[both_active, f"{h2}_tx_norm"].values
        rho = spearmanr(v1, v2).correlation if n > 1 else None
        results["transcription"][f"{h1}_{h2}"] = {"n": n, "rho": float(rho) if rho == rho else None}
        print(f"  {h1} vs {h2}: n={n}, rho={rho:.4f}" if rho == rho else f"  {h1} vs {h2}: n={n}, rho=undefined")

    print("\n=== TRANSLATION: cross-host measurement correlation (on floor-corrected actives in BOTH hosts) ===")
    for h1, h2 in PAIRS:
        both_active = tl_active_col[h1] & tl_active_col[h2] & tl_usable_reg_col[h1] & tl_usable_reg_col[h2]
        n = int(both_active.sum())
        v1 = df.loc[both_active, f"{h1}_protein_log10"].values
        v2 = df.loc[both_active, f"{h2}_protein_log10"].values
        rho = spearmanr(v1, v2).correlation if n > 1 else None
        results["translation"][f"{h1}_{h2}"] = {"n": n, "rho": float(rho) if rho == rho else None}
        print(f"  {h1} vs {h2}: n={n}, rho={rho:.4f}" if rho == rho else f"  {h1} vs {h2}: n={n}, rho=undefined")

    with open(RESULTS / "gate5_5_crosshost_measurement_correlation.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {RESULTS / 'gate5_5_crosshost_measurement_correlation.json'}")


if __name__ == "__main__":
    main()
