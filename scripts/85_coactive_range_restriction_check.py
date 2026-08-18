"""
GATE 8.6 - Task 2b: is the co-active-restriction "halving" of cross-host
correlation (Gate 8.5 Task 1C) an N or range-restriction artifact?

Range restriction is a classical, mechanical way to shrink a correlation
coefficient: restricting to a narrower slice of the X (or Y) distribution
mechanically lowers the observed correlation even with NO change in the
underlying relationship, because Spearman/Pearson both depend on the
spread of ranks/values present in the sample. Co-active restriction could
plausibly do exactly this -- "active" sequences by construction exclude the
zero/near-zero tail, narrowing the range on both axes.

This script reports, per pair/readout, the RANGE (IQR and full range) of
each host's own values in both the pooled and co-active-restricted subsets,
so a reader can see directly whether the restricted subset's distribution
is meaningfully narrower -- not just smaller in N.
"""
import pandas as pd
import numpy as np
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"

HOSTS = ["EC", "BS", "PA"]
PAIRS = [("EC", "BS"), ("EC", "PA"), ("BS", "PA")]


def spread_stats(vals):
    vals = np.asarray(vals, dtype=np.float64)
    if len(vals) < 2:
        return {"n": len(vals), "iqr": None, "std": None, "range": None}
    q25, q75 = np.percentile(vals, [25, 75])
    return {"n": int(len(vals)), "iqr": float(q75 - q25), "std": float(vals.std()),
            "range": float(vals.max() - vals.min())}


def main():
    df = pd.read_parquet(DATA / "three_host.parquet")

    tl_active_col, tl_usable_reg_col = {}, {}
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

    rows = []
    for h1, h2 in PAIRS:
        # TRANSCRIPTION
        usable_both = (df[f"{h1}_tx_usable"] & df[f"{h2}_tx_usable"]).values
        active_both = (df[f"{h1}_tx_usable"] & df[f"{h1}_tx_active"] &
                        df[f"{h2}_tx_usable"] & df[f"{h2}_tx_active"]).values
        for h in (h1, h2):
            pooled_stats = spread_stats(df.loc[usable_both, f"{h}_tx_norm"])
            active_stats = spread_stats(df.loc[active_both, f"{h}_tx_norm"])
            rows.append({"pair": f"{h1}_{h2}", "readout": "transcription", "host_measured": h,
                         "pooled_n": pooled_stats["n"], "pooled_iqr": pooled_stats["iqr"],
                         "pooled_std": pooled_stats["std"], "pooled_range": pooled_stats["range"],
                         "coactive_n": active_stats["n"], "coactive_iqr": active_stats["iqr"],
                         "coactive_std": active_stats["std"], "coactive_range": active_stats["range"],
                         "iqr_ratio_coactive_over_pooled": (active_stats["iqr"] / pooled_stats["iqr"])
                             if pooled_stats["iqr"] and pooled_stats["iqr"] > 0 else None})

        # TRANSLATION
        usable_both_tl = (tl_usable_reg_col[h1] & tl_usable_reg_col[h2])
        active_both_tl = (tl_active_col[h1] & tl_active_col[h2] & tl_usable_reg_col[h1] & tl_usable_reg_col[h2])
        for h in (h1, h2):
            pooled_stats = spread_stats(df.loc[usable_both_tl, f"{h}_protein_log10"])
            active_stats = spread_stats(df.loc[active_both_tl, f"{h}_protein_log10"])
            rows.append({"pair": f"{h1}_{h2}", "readout": "translation", "host_measured": h,
                         "pooled_n": pooled_stats["n"], "pooled_iqr": pooled_stats["iqr"],
                         "pooled_std": pooled_stats["std"], "pooled_range": pooled_stats["range"],
                         "coactive_n": active_stats["n"], "coactive_iqr": active_stats["iqr"],
                         "coactive_std": active_stats["std"], "coactive_range": active_stats["range"],
                         "iqr_ratio_coactive_over_pooled": (active_stats["iqr"] / pooled_stats["iqr"])
                             if pooled_stats["iqr"] and pooled_stats["iqr"] > 0 else None})

    out_df = pd.DataFrame(rows)
    out_df.to_csv(RESULTS / "gate8_6_coactive_range_restriction_check.csv", index=False)
    pd.set_option("display.width", 200)
    print(out_df.round(3).to_string(index=False))
    print(f"\nWrote {RESULTS / 'gate8_6_coactive_range_restriction_check.csv'}")


if __name__ == "__main__":
    main()
