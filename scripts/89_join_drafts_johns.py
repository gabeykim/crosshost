"""
GATE 10 - Task 3 [LOAD-BEARING]: join DRAFTS to the existing Johns et al.
tables. Verified, not assumed, per the task's explicit instruction.

ID JOIN: DOES NOT WORK. DRAFTS's Oligo ID range (13097-14477) and this
project's three_host.parquet OLIGO ID range (14478-43726) are entirely
disjoint -- not a single value in common, and the ranges are contiguous
(DRAFTS ends exactly one below where three_host begins), strongly
suggesting the two libraries drew from the same overall synthesis/ID
numbering pipeline but non-overlapping ID batches. RS241's id range
(14772-48469) also shares zero values with DRAFTS's range by construction.

SEQUENCE-TEXT JOIN: works, and is the join key used throughout this
project's downstream analysis (as it is for source-genome-based fold
assignment, exact-duplicate detection, etc. elsewhere in this project).
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)


def main():
    drafts = pd.read_parquet(DATA / "drafts.parquet")
    johns = pd.read_parquet(DATA / "three_host.parquet")
    lib = pd.read_parquet(DATA / "three_host_library.parquet")
    rs241 = pd.read_parquet(DATA / "rs241.parquet")
    folds = pd.read_parquet(DATA / "splits" / "fold_assignment_FINAL.parquet")

    report = {}

    # --- ID join check ---
    drafts_ids = set(drafts["oligo_id"])
    johns_ids = set(johns["OLIGO ID"].astype(int))
    rs241_ids = set(rs241["id"].astype(int))
    report["id_join"] = {
        "drafts_id_range": [int(drafts["oligo_id"].min()), int(drafts["oligo_id"].max())],
        "johns_three_host_id_range": [int(johns["OLIGO ID"].astype(int).min()), int(johns["OLIGO ID"].astype(int).max())],
        "rs241_id_range": [int(rs241["id"].astype(int).min()), int(rs241["id"].astype(int).max())],
        "drafts_johns_id_overlap": len(drafts_ids & johns_ids),
        "drafts_rs241_id_overlap": len(drafts_ids & rs241_ids),
        "verdict": "IDs do NOT join -- ranges are disjoint. Join must use sequence text.",
    }
    print("=== ID JOIN ===")
    print(json.dumps(report["id_join"], indent=2))

    # --- Sequence-text join: DRAFTS vs three_host ---
    johns_seq_map = johns.set_index(johns["regulatory_sequence"].str.upper().str.strip())
    drafts_seq = drafts["sequence_165bp"]
    match_mask_three_host = drafts_seq.isin(set(johns_seq_map.index))
    n_match_three_host = int(match_mask_three_host.sum())

    # --- Sequence-text join: DRAFTS vs RS241 (via library ID lookup for RS241 sequence text) ---
    rs241_seq = lib.set_index("OLIGO ID").reindex(rs241["id"])["Regulatory Sequence"]
    n_rs241_seq_recoverable = int(rs241_seq.notna().sum())
    rs241_seq_set = set(rs241_seq.dropna().str.upper().str.strip())
    match_mask_rs241 = drafts_seq.isin(rs241_seq_set)
    n_match_rs241 = int(match_mask_rs241.sum())

    report["sequence_join"] = {
        "drafts_n_sequences": len(drafts),
        "drafts_vs_three_host_exact_match": n_match_three_host,
        "drafts_vs_three_host_match_pct_of_drafts": round(100 * n_match_three_host / len(drafts), 2),
        "rs241_sequences_recoverable_via_library_lookup": f"{n_rs241_seq_recoverable}/241",
        "drafts_vs_rs241_exact_match": n_match_rs241,
        "join_key": "sequence_165bp (uppercased, stripped) == three_host.regulatory_sequence / library-looked-up RS241 sequence",
    }
    print("\n=== SEQUENCE-TEXT JOIN ===")
    print(json.dumps(report["sequence_join"], indent=2))

    # --- Build the joined table (three_host overlap) ---
    # NOTE: three_host.parquet contains 181 groups of coincidentally-identical
    # sequence text across different genomic loci (Gate 2 finding) -- a
    # sequence-text join can hit >1 Johns row for the same DRAFTS sequence.
    # Keep the first match and report how many collisions occurred, rather
    # than silently dropping or duplicating rows.
    johns_seq_col = johns["regulatory_sequence"].str.upper().str.strip()
    dup_seq_count = johns_seq_col[johns_seq_col.isin(set(drafts["sequence_165bp"]))].duplicated().sum()
    johns_dedup = johns.assign(_seq=johns_seq_col).drop_duplicates(subset="_seq", keep="first").set_index("_seq")

    joined = drafts[match_mask_three_host].copy()
    matched_johns_rows = johns_dedup.reindex(joined["sequence_165bp"])
    joined["johns_oligo_id"] = matched_johns_rows["OLIGO ID"].values
    joined["johns_intgen_id"] = matched_johns_rows["intgen_id"].values if "intgen_id" in matched_johns_rows.columns else None
    joined["johns_source_genome"] = matched_johns_rows["genome_id"].values if "genome_id" in matched_johns_rows.columns else None
    if dup_seq_count > 0:
        print(f"  NOTE: {dup_seq_count} Johns-side sequence collisions among matched rows (kept first match each)")

    # --- Fold coverage check ---
    # NOTE: fold_assignment_FINAL.parquet stores "OLIGO ID" as string dtype,
    # not int -- a real dtype mismatch caught here (an initial int-keyed
    # reindex silently returned 0/112 matches; fixed by matching dtypes
    # explicitly rather than trusting pandas to coerce).
    fold_map = folds.assign(**{"OLIGO ID": folds["OLIGO ID"].astype(int)}).set_index("OLIGO ID")["fold"]
    joined["fold"] = fold_map.reindex(joined["johns_oligo_id"].astype(int)).values
    n_in_fold = int(joined["fold"].notna().sum())
    n_total_joined = len(joined)
    fold_counts = joined["fold"].value_counts().sort_index().to_dict()

    coverage_case = ("FULLY COVERED" if n_in_fold == n_total_joined else
                      "NOT COVERED" if n_in_fold == 0 else
                      f"PARTIALLY COVERED ({n_in_fold}/{n_total_joined} = {100*n_in_fold/n_total_joined:.1f}%)")

    report["fold_compatibility"] = {
        "n_drafts_sequences_matched_to_three_host": n_total_joined,
        "n_of_those_present_in_frozen_fold_assignment": n_in_fold,
        "case": coverage_case,
        "fold_distribution_of_matched_rows": {str(int(k)): int(v) for k, v in fold_counts.items() if pd.notna(k)},
        "note": ("The single fold=NaN row is a matched sequence whose OLIGO ID was excluded from the main-library "
                 "fold assignment for the same reason every RS241 ID was (see 'of_uncovered_rows_also_in_rs241' "
                 "below -- confirms it is literally an RS241 sequence, not a separate excluded category)."),
        "FROZEN SPLITS NOT MODIFIED THIS GATE": True,
    }
    # check whether the fold-uncovered matched rows are RS241 sequences specifically
    uncovered = joined[joined["fold"].isna()]
    uncovered_in_rs241 = uncovered["sequence_165bp"].isin(rs241_seq_set).sum()
    report["fold_compatibility"]["of_uncovered_rows_also_in_rs241"] = int(uncovered_in_rs241)

    print("\n=== FOLD COMPATIBILITY ===")
    print(json.dumps(report["fold_compatibility"], indent=2))

    # --- Per-species overlap with the three primary hosts this project actually models (EC, BS, PA) ---
    # CRITICAL: DRAFTS has no P. aeruginosa (its 'Pa' = P. agglomerans) -- only EC and BS overlap.
    overlap_ec = int((joined["Ec_usable"]).sum())
    overlap_bs = int((joined["Bs_usable"]).sum())
    report["primary_host_overlap"] = {
        "P_aeruginosa_in_DRAFTS": False,
        "P_aeruginosa_note": "DRAFTS's 'Pa' column is Pantoea agglomerans, confirmed against Appendix Table S1 and full-text species list -- NOT Pseudomonas aeruginosa. Zero occurrences of 'aeruginosa' anywhere in the paper.",
        "n_matched_rows_usable_in_Ec": overlap_ec,
        "n_matched_rows_usable_in_Bs": overlap_bs,
        "n_matched_rows_usable_in_both_Ec_and_Bs": int((joined["Ec_usable"] & joined["Bs_usable"]).sum()),
        "of_this_project_3_primary_hosts_EC_BS_PA_only_2_overlap_with_DRAFTS": ["EC", "BS"],
    }
    print("\n=== PRIMARY HOST OVERLAP (this project's EC/BS/PA vs DRAFTS) ===")
    print(json.dumps(report["primary_host_overlap"], indent=2))

    # --- Source-genome metadata carried natively ---
    report["source_genome_metadata"] = {
        "carried_natively_by_drafts": True,
        "columns": ["source_species", "source_genus", "source_family", "source_order", "source_class",
                    "source_phylum", "source_genome"],
        "note": "DRAFTS's Fig. 3 sheet carries its own full source-genome taxonomy per sequence, independent of "
                "Johns's metadata -- does not need to be joined from three_host_library for genome-blocked splitting, "
                "though cross-checking the two sources' agreement on the 112 overlapping sequences would be a natural "
                "sanity check for Gate 11, not required here.",
    }

    with open(RESULTS / "gate10_join_report.json", "w") as f:
        json.dump(report, f, indent=2, default=str)

    joined.to_parquet(DATA / "drafts_johns_overlap.parquet", index=False)
    print(f"\nWrote {RESULTS / 'gate10_join_report.json'}")
    print(f"Wrote {DATA / 'drafts_johns_overlap.parquet'} ({len(joined)} rows)")


if __name__ == "__main__":
    main()
