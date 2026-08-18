"""
GATE 2 - Task 0: Depth-matching audit on the Gate 1.5 physiology proxy vector.

Tests the hypothesis that shallow proteome coverage (V. natriegens: 665 proteins,
vs up to 5,034 for other hosts) mechanically inflates abundance-fraction metrics
for high-abundance protein classes (ribosomal proteins above all). Uses the same
classification rules as scripts/08_compute_physiology_proxies.py (imported
directly, not reimplemented, so there is no risk of drift between the two).
"""
import sys
import json
import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
proxymod = import_module("08_compute_physiology_proxies")

RAW = Path(__file__).resolve().parent.parent / "raw"
DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"

VN_DEPTH = 1032  # reprocessed pooled-81-sample depth (script 11), not the 665
                  # single-condition depth used in Gate 1.5 -- see script 11 docstring
TRUNCATION_DEPTHS = [1032, 1500, 2000, None]  # None = full depth

PAXDB_PRIMARY = {
    "EC": "paxdb_abundances_EC_integrated_3740039012.json",
    "BS": "paxdb_abundances_BS_integrated_4187424622.json",
    "PA": "paxdb_abundances_PA_integrated_1789377356.json",
    "SE": "paxdb_abundances_SE_integrated_1793422220.json",
    "CG": "paxdb_abundances_CG_integrated_3268883846.json",
}

CATS = list(proxymod.GENE_SYMBOL_RULES.keys())


def load_paxdb_df(fname):
    df, payload = proxymod.load_paxdb(fname)
    # rank by abundance descending (this is how PaxDb abundance/depth ordering works --
    # deeper coverage adds progressively lower-abundance proteins)
    df = df.sort_values("abundance_ppm", ascending=False).reset_index(drop=True)
    return df


def truncate_and_recompute(df, depth):
    sub = df if depth is None else df.head(depth)
    result, _ = proxymod.compute_fractions(sub, host="_", source_label=f"depth={depth}")
    return result


def load_vnatriegens_pooled():
    df = pd.read_parquet(DATA / "vnatriegens_proteome_pooled_81samples.parquet")
    import re
    gn_re = re.compile(r"GN=(\S+)")

    def extract_gene(row):
        m = gn_re.search(str(row["protein_names"]))
        return m.group(1) if m else ""

    df["gene_name"] = df.apply(extract_gene, axis=1)
    df["annotation_text"] = df["protein_names"]
    return df[["gene_name", "annotation_text", "abundance_ppm"]]


def main():
    vn_df = load_vnatriegens_pooled()
    print(f"V. natriegens reprocessed (pooled, 81-sample) depth: {len(vn_df)} proteins")

    # ------------------------------------------------------------------
    # 1. Depth-matched recomputation: truncate all 5 PaxDb hosts to top 665
    # ------------------------------------------------------------------
    print("\n" + "=" * 90)
    print("1. FULL-DEPTH vs DEPTH-MATCHED (top 665) -- all 6 metrics, all hosts")
    print("=" * 90)

    rows = []
    for host, fname in PAXDB_PRIMARY.items():
        df = load_paxdb_df(fname)
        full = truncate_and_recompute(df, None)
        matched = truncate_and_recompute(df, VN_DEPTH)
        row = {"host": host, "n_full": len(df), "n_matched": VN_DEPTH}
        for cat in CATS:
            row[f"{cat}_full"] = full[f"{cat}_mass_fraction"]
            row[f"{cat}_matched"] = matched[f"{cat}_mass_fraction"]
        rows.append(row)

    vn_full = truncate_and_recompute(vn_df, None)
    rows.append({"host": "VN", "n_full": len(vn_df), "n_matched": VN_DEPTH,
                 **{f"{cat}_full": vn_full[f"{cat}_mass_fraction"] for cat in CATS},
                 **{f"{cat}_matched": vn_full[f"{cat}_mass_fraction"] for cat in CATS}})

    depth_df = pd.DataFrame(rows)
    for cat in CATS:
        print(f"\n--- {cat} ---")
        print(depth_df[["host", "n_full", f"{cat}_full", f"{cat}_matched"]].to_string(index=False))

    depth_df.to_csv(DATA / "depth_matching_audit.csv", index=False)

    # ------------------------------------------------------------------
    # 2. Does cross-host ordering survive? Spearman rank correlation
    # ------------------------------------------------------------------
    print("\n" + "=" * 90)
    print("2. SPEARMAN RANK CORRELATION: full-depth ranking vs depth-matched ranking")
    print("=" * 90)
    spearman_results = {}
    for cat in CATS:
        full_vals = depth_df[f"{cat}_full"].values
        matched_vals = depth_df[f"{cat}_matched"].values
        rho, pval = spearmanr(full_vals, matched_vals)
        spearman_results[cat] = {"rho": float(rho), "pval": float(pval)}
        print(f"  {cat}: rho={rho:.3f} (p={pval:.3f})")
        full_rank = depth_df.assign(r=depth_df[f"{cat}_full"].rank(ascending=False))[["host", "r"]]
        matched_rank = depth_df.assign(r=depth_df[f"{cat}_matched"].rank(ascending=False))[["host", "r"]]
        top_full = full_rank.sort_values("r").iloc[0]["host"]
        top_matched = matched_rank.sort_values("r").iloc[0]["host"]
        print(f"    top host: full={top_full}, depth-matched={top_matched}")

    # ------------------------------------------------------------------
    # 3. Depth sensitivity curve for the 3 deepest hosts
    # ------------------------------------------------------------------
    print("\n" + "=" * 90)
    print("3. DEPTH SENSITIVITY CURVE (ribosomal_proteins fraction), 3 deepest hosts")
    print("=" * 90)
    deepest = depth_df.sort_values("n_full", ascending=False)["host"].head(3).tolist()
    print(f"3 deepest hosts (excl. VN): {deepest}")

    curve_rows = []
    for host in deepest:
        if host == "VN":
            continue
        fname = PAXDB_PRIMARY[host]
        df = load_paxdb_df(fname)
        for depth in TRUNCATION_DEPTHS:
            res = truncate_and_recompute(df, depth)
            n = len(df) if depth is None else min(depth, len(df))
            row = {"host": host, "truncation_depth": depth if depth else f"full({len(df)})", "n_actual": n}
            for cat in CATS:
                row[cat] = res[f"{cat}_mass_fraction"]
            curve_rows.append(row)
            print(f"  {host} depth={row['truncation_depth']:>10}: "
                  f"ribosomal={res['ribosomal_proteins_mass_fraction']*100:.2f}%  "
                  f"rnap={res['rnap_core_mass_fraction']*100:.3f}%  "
                  f"sigma={res['sigma_factors_mass_fraction']*100:.3f}%  "
                  f"chap={res['chaperones_mass_fraction']*100:.2f}%  "
                  f"ef={res['elongation_factors_mass_fraction']*100:.2f}%")

    curve_df = pd.DataFrame(curve_rows)
    curve_df.to_csv(DATA / "depth_sensitivity_curve.csv", index=False)

    # ------------------------------------------------------------------
    # Save everything for the memo / verdict step
    # ------------------------------------------------------------------
    with open(OUT / "depth_matching_results.json", "w") as f:
        json.dump({
            "depth_matched_table": depth_df.to_dict(orient="records"),
            "spearman": spearman_results,
            "depth_curve": curve_rows,
        }, f, indent=2, default=str)
    print(f"\nWrote {OUT / 'depth_matching_results.json'}")
    print(f"Wrote {DATA / 'depth_matching_audit.csv'}")
    print(f"Wrote {DATA / 'depth_sensitivity_curve.csv'}")


if __name__ == "__main__":
    main()
