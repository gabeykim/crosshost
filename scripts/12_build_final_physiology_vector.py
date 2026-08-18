"""
GATE 2 - Task 0 (final step): Freeze the corrected physiology vector.

Per the Task 0 depth-matching audit (scripts/10, /11), all 6 physiology metrics
survive as usable AFTER depth-matching correction (none are DEPTH-CONFOUNDED,
UNUSABLE by the Spearman-rank-stability / depth-curve evidence). This script
computes the final, frozen per-host values at the common depth-matched
truncation (1032 proteins -- V. natriegens's reprocessed native depth, see
script 11), for consumption by data/hosts_physiology.parquet (Task 1).

growth_rate_h is carried forward from Gate 1's literature review as-is (not
proteomics-depth-related, so not subject to this audit) -- see
out/task4_physiology_primary3.md / task4_physiology_secondary3.md for sourcing.
"""
import sys
import json
import pandas as pd
from pathlib import Path
from importlib import import_module

sys.path.insert(0, str(Path(__file__).resolve().parent))
proxymod = import_module("08_compute_physiology_proxies")
depthmod = import_module("10_depth_matching_audit")

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"

COMMON_DEPTH = 1032

PAXDB_PRIMARY = {
    "EC": ("paxdb_abundances_EC_integrated_3740039012.json", 3740039012, "PaxDb v6.0 Integrated"),
    "BS": ("paxdb_abundances_BS_integrated_4187424622.json", 4187424622, "PaxDb v6.0 Integrated"),
    "PA": ("paxdb_abundances_PA_integrated_1789377356.json", 1789377356, "PaxDb v6.0 Integrated"),
    "SE": ("paxdb_abundances_SE_integrated_1793422220.json", 1793422220, "PaxDb v6.0 Integrated"),
    "CG": ("paxdb_abundances_CG_integrated_3268883846.json", 3268883846, "PaxDb v6.0 Integrated"),
}

# [VERIFIED/COMPUTED in Gate 1 -- see out/task4_physiology_primary3.md, task4_physiology_secondary3.md]
# Best available published growth rate per host, with source condition noted.
# NOT condition-matched across hosts -- carried as metadata, not silently dropped.
#
# BS UPDATED IN GATE 4.5: previously None (Gate 1-4 imputed this downstream as
# the mean of the other 5 hosts' mu_h, 1.852 -- a real problem, since that made
# the feature informationless for BS by construction and mildly contaminated
# by the LOHO training hosts). Gate 4.5 re-attempted the search and downloaded
# the primary source directly (Zhu, Mori, Hwa & Dai 2025 PNAS, full-text PDF +
# SI Appendix PDF, via the PMC OA FTP mirror -- ftp://ftp.ncbi.nlm.nih.gov/
# pub/pmc/deprecated/oa_package/82/94/PMC12067254.tar.gz). The SI Appendix
# Tables S1/S2 (the exact source the paper's own text points to for numeric
# growth-rate values) were present in the downloaded PDF but do not contain a
# directly tabulated B. subtilis-in-LB-at-37C digit either (they are strain/
# medium reference tables, not a results table). The value WAS recoverable
# from Fig. 2B (main text) -- a bar chart of growth rate (1/h) for Ec/Bs/Vn
# across several media at 37C, LB being each species' tallest (leftmost) bar.
# Visually read: B. subtilis LB/37C bar height ~1.7 (h^-1), reasonably close
# to and consistent with the paper's own textual claim that "V. natriegens
# grows >50% faster than Ec and Bs...[implying]...E. coli and B. subtilis max
# growth rates in LB at 37C are much closer to each other than either is to
# V. natriegens" -- i.e. Bs and Ec's LB/37C rates should be similar, and they
# are (Bs~1.7 vs Ec's own already-adopted 1.7). This is READ FROM A FIGURE,
# not a tabulated digit -- flagged [INFERRED, medium-high confidence], the
# same confidence tier as EC's and PA's own values (both are also midpoints
# of a literature range, not a single verified digit either). Critically,
# this value is now independently sourced (not derived from any other host's
# value), resolving the "informationless by construction" and "contamination
# from training hosts" problems Gate 4.5 was tasked with fixing -- it is a
# genuine, condition-matched (LB, 37C, matching EC's and PA's own condition),
# B.-subtilis-specific measurement, just not extracted to full decimal
# precision from a table.
GROWTH_RATE = {
    "EC": {"mu_h": 1.7, "condition": "LB, 37C (midpoint of Sezonov et al. 2007 range 1.4-2.0)", "source": "Sezonov, Joseleau-Petit & D'Ari 2007, J Bacteriol", "confidence": "verified_range_midpoint"},
    "BS": {"mu_h": 1.7, "condition": "LB, 37C -- grown alongside E. coli in Zhu et al. 2025 PNAS; visually read from Fig. 2B bar chart (SI Appendix Tables S1/S2 do not tabulate this digit directly)", "source": "Zhu, Mori, Hwa & Dai 2025, PNAS 122:e2427091122, Fig. 2B", "confidence": "inferred_from_figure"},
    "PA": {"mu_h": 1.6, "condition": "LB, 37C (midpoint of Yang et al. 2008 range 1.54-1.73, PAO1/PA14)", "source": "Yang et al. 2008, J Bacteriol", "confidence": "verified_range_midpoint"},
    "SE": {"mu_h": 0.95, "condition": "M9+0.4% glucose, 37C, strain LT2", "source": "Mishra & Shashidhar 2022, J Bacteriol", "confidence": "verified_single_value"},
    "VN": {"mu_h": 4.43, "condition": "BHIN complex medium, 37C (fastest reported; 1.5-1.7 in minimal medium)", "source": "Hoffart et al. 2017, AEM", "confidence": "verified_single_value"},
    "CG": {"mu_h": 0.58, "condition": "CGXII minimal + glucose, 30C (midpoint of 0.5-0.65 range across 3 papers)", "source": "Baumgart et al./Unthan et al./Nat Commun 2023", "confidence": "verified_range_midpoint"},
}


def main():
    rows = []
    for host, (fname, dsid, source) in PAXDB_PRIMARY.items():
        df = depthmod.load_paxdb_df(fname)
        n_full = len(df)
        res = depthmod.truncate_and_recompute(df, COMMON_DEPTH)
        row = {
            "host": host,
            "physiology_source": source,
            "physiology_dataset_id": str(dsid),
            "proteome_depth_native": n_full,
            "proteome_depth_used": min(COMMON_DEPTH, n_full),
            "depth_matched": n_full > COMMON_DEPTH,
        }
        for cat in proxymod.GENE_SYMBOL_RULES:
            row[f"{cat}_fraction"] = res[f"{cat}_mass_fraction"]
            row[f"{cat}_n_proteins"] = res[f"{cat}_n_proteins"]
        rows.append(row)

    # V. natriegens: reprocessed pooled dataset, native depth 1032 == COMMON_DEPTH, no truncation needed
    vn_df = depthmod.load_vnatriegens_pooled()
    res = proxymod.compute_fractions(vn_df, "VN", "pooled_81samples")[0]
    rows.append({
        "host": "VN",
        "physiology_source": "PRIDE PXD027874 (Hervey et al., Naval Research Lab), pooled across 81 samples (3 timepoints x 3 NaCl levels x 3 temps x 3 reps)",
        "physiology_dataset_id": "PXD027874_pooled81",
        "proteome_depth_native": len(vn_df),
        "proteome_depth_used": len(vn_df),
        "depth_matched": False,
        **{f"{cat}_fraction": res[f"{cat}_mass_fraction"] for cat in proxymod.GENE_SYMBOL_RULES},
        **{f"{cat}_n_proteins": res[f"{cat}_n_proteins"] for cat in proxymod.GENE_SYMBOL_RULES},
    })

    df_final = pd.DataFrame(rows)

    # attach growth rate metadata
    df_final["growth_rate_mu_h"] = df_final["host"].map(lambda h: GROWTH_RATE[h]["mu_h"])
    df_final["growth_rate_condition"] = df_final["host"].map(lambda h: GROWTH_RATE[h]["condition"])
    df_final["growth_rate_source"] = df_final["host"].map(lambda h: GROWTH_RATE[h]["source"])
    df_final["growth_rate_confidence"] = df_final["host"].map(lambda h: GROWTH_RATE[h]["confidence"])
    # Gate 4.5: no host's growth_rate_mu_h is imputed-from-other-hosts anymore
    # (BS previously was; now independently sourced, see GROWTH_RATE comment
    # above). Column kept, not removed, so a future gate that re-introduces
    # an imputed value has an explicit place to flag it -- required by the
    # Gate 4.5 task spec ("an `imputed` flag column ... if any value remains
    # imputed").
    df_final["growth_rate_imputed"] = False

    print(df_final.to_string(index=False))

    out_path = DATA / "hosts_physiology_final.csv"
    df_final.to_csv(out_path, index=False)
    print(f"\nWrote {out_path}")

    with open(OUT / "task0_verdict.json", "w") as f:
        json.dump({
            "common_depth_used": COMMON_DEPTH,
            "vn_reprocessed_depth": len(vn_df),
            "vn_original_gate1_5_depth": 665,
            "per_metric_verdict": {
                "ribosomal_proteins": "DEPTH-SENSITIVE, CORRECTABLE -- Spearman rho=1.000 full-vs-matched ranking, smooth monotonic depth curve, use depth-matched value",
                "rnap_core": "DEPTH-SENSITIVE, CORRECTABLE -- Spearman rho=1.000, small absolute shifts, use depth-matched value",
                "chaperones": "DEPTH-SENSITIVE, CORRECTABLE -- Spearman rho=1.000, P. aeruginosa shows largest correction, use depth-matched value",
                "sigma_factors": "DEPTH-SENSITIVE, CORRECTABLE but LOW RELIABILITY -- Spearman rho=0.943, absolute values tiny (<0.2%), recommend reporting distinct-sigma-factor COUNT alongside fraction, treat fraction with caution",
                "elongation_factors": "DEPTH-SENSITIVE, CORRECTABLE but WEAKEST rank-stability -- Spearman rho=0.829, some mid-ranking reorder between full and matched depth",
                "growth_rate": "DEPTH-ROBUST -- not proteomics-derived, no depth confound applies (carries its own condition-matching caveats from Gate 1, unrelated to this audit)",
            },
            "n_metrics_surviving": 6,
            "decision_needed_triggered": False,
        }, f, indent=2)
    print(f"Wrote {OUT / 'task0_verdict.json'}")


if __name__ == "__main__":
    main()
