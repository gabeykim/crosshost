"""
GATE 1 - Task 1 (dataset shape) and Task 2 (RS241 integrity) analysis.

Definitions used (all justified in comments, cited to the manuscript body text
extracted from raw/nihms945382_bodytext.txt):

TRANSCRIPTION "usable value" (paper's own QC rule, "most analyses"):
    NOT (dna_count == 0) AND (rna_count + dna_count) >= 15
    Source: "We excluded constructs containing 0 DNA counts and also those whose
    RNA and DNA counts summed to less than 15 for most analyses."

TRANSCRIPTION "active" (paper's own headline stat, no total-count qualifier):
    rna_count > 0
    Source: "B. subtilis displayed the lowest number of measurably active RSs
    (18.9% with > 0 RNA reads)"

TRANSLATION "usable value":
    protein (log10) is not null. This is the FACS-seq-derived protein level;
    the paper does not state an explicit per-oligo count-based exclusion rule
    for this column analogous to the transcription one, so "usable" = present.

TRANSLATION "active" -- no directly reusable per-oligo boolean flag exists in
    the supplementary table; the paper's 3.3%/290-construct "active in all
    species" figure is a HEADLINE STAT computed by the authors with an
    unstated per-host threshold (likely FACS-seq bin above the no-GFP negative
    control). We report protein(log10) usable counts, and separately flag that
    an "active" (signal-above-baseline) definition could not be reconstructed
    per-oligo from the public table alone.

Part length: computed directly from Full Constructs sequences (RS + ATG + 12bp
    barcode) and cross-checked against the manuscript's stated "165 bp
    immediately upstream of annotated start codons."
"""
import pandas as pd
import numpy as np
from pathlib import Path
import json

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"
OUT.mkdir(exist_ok=True)

pd.set_option("display.width", 140)


def main():
    results = {}
    df = pd.read_parquet(DATA / "three_host_library.parquet")
    print(f"Loaded three_host_library.parquet: {df.shape}")
    print(f"Total unique OLIGO IDs: {df['OLIGO ID'].nunique()}")

    # ---------------------------------------------------------------
    # Q6 (do this first): column inventory
    # ---------------------------------------------------------------
    print("\n" + "=" * 70)
    print("Q6: COLUMN INVENTORY")
    print("=" * 70)
    for c in df.columns:
        dtype = df[c].dtype
        n_null = df[c].isna().sum()
        print(f"  {c!r:45s} dtype={str(dtype):10s} n_null={n_null}")

    # ---------------------------------------------------------------
    # Q3: part length
    # ---------------------------------------------------------------
    print("\n" + "=" * 70)
    print("Q3: PART (REGULATORY SEQUENCE) LENGTH")
    print("=" * 70)
    rs_len_from_meta = df["Regulatory Sequence"].dropna().str.len()
    print("From Metadata['Regulatory Sequence'] column directly:")
    print(f"  n={len(rs_len_from_meta)}  min={rs_len_from_meta.min()}  "
          f"median={rs_len_from_meta.median()}  max={rs_len_from_meta.max()}  "
          f"mean={rs_len_from_meta.mean():.2f}  std={rs_len_from_meta.std():.3f}")
    print("  value_counts (top 10):")
    print(rs_len_from_meta.value_counts().head(10).to_string())

    # cross-check via Full Constructs: seq = RS + ATG(3bp) + 12bp barcode
    fc_len = df["full_construct_seq"].dropna().str.len()
    implied_rs_len = fc_len - 3 - 12
    print("\nFrom Full Constructs length minus ATG(3bp) minus barcode(12bp):")
    print(f"  n={len(implied_rs_len)}  min={implied_rs_len.min()}  "
          f"median={implied_rs_len.median()}  max={implied_rs_len.max()}")
    print("  value_counts (top 10):")
    print(implied_rs_len.value_counts().head(10).to_string())

    results["part_length_bp_from_metadata_col"] = {
        "n": int(len(rs_len_from_meta)),
        "min": int(rs_len_from_meta.min()),
        "median": float(rs_len_from_meta.median()),
        "max": int(rs_len_from_meta.max()),
        "mean": float(rs_len_from_meta.mean()),
        "value_counts": rs_len_from_meta.value_counts().to_dict(),
    }
    results["part_length_variable"] = bool(rs_len_from_meta.nunique() > 1)

    # ---------------------------------------------------------------
    # Q2 + active fractions (Q4): per-host usable & active counts
    # ---------------------------------------------------------------
    print("\n" + "=" * 70)
    print("Q2 & Q4: PER-HOST TRANSCRIPTION/TRANSLATION USABLE & ACTIVE COUNTS")
    print("=" * 70)
    n_total = len(df)
    host_stats = {}
    for host in ["BS", "EC", "PA"]:
        rna = df[f"{host}_rna_count"]
        dna = df[f"{host}_dna_count"]
        protein = df[f"{host}_protein (log10)"]

        tx_usable = (~(dna == 0)) & ((rna + dna) >= 15) & rna.notna() & dna.notna()
        tx_active = rna > 0  # paper's headline definition, unqualified by total count
        tx_active_and_present = tx_active & rna.notna()

        tl_usable = protein.notna()

        print(f"\n--- Host: {host} ---")
        print(f"  transcription usable (paper QC: dna!=0 & rna+dna>=15): "
              f"{tx_usable.sum()} / {n_total} = {100*tx_usable.sum()/n_total:.2f}%")
        print(f"  transcription active (rna_count > 0, unqualified):     "
              f"{tx_active_and_present.sum()} / {rna.notna().sum()} = "
              f"{100*tx_active_and_present.sum()/rna.notna().sum():.2f}%  "
              f"[paper reports 18.9% BS / 52.0% EC / 83.8% PA]")
        print(f"  translation usable (protein(log10) not null):         "
              f"{tl_usable.sum()} / {n_total} = {100*tl_usable.sum()/n_total:.2f}%")
        print(f"  rna_count null count: {rna.isna().sum()}, dna_count null count: {dna.isna().sum()}")

        host_stats[host] = {
            "tx_usable_count": int(tx_usable.sum()),
            "tx_usable_pct": float(100 * tx_usable.sum() / n_total),
            "tx_active_gt0_count": int(tx_active_and_present.sum()),
            "tx_active_gt0_pct_of_nonnull": float(100 * tx_active_and_present.sum() / rna.notna().sum()),
            "tl_usable_count": int(tl_usable.sum()),
            "tl_usable_pct": float(100 * tl_usable.sum() / n_total),
        }
        # stash boolean masks for overlap calc
        df[f"_{host}_tx_usable"] = tx_usable
        df[f"_{host}_tl_usable"] = tl_usable

    results["per_host_counts"] = host_stats

    # ---------------------------------------------------------------
    # Q1: three-host overlap (RS-level, after collapsing barcode replicates)
    # ---------------------------------------------------------------
    print("\n" + "=" * 70)
    print("Q1: THREE-HOST OVERLAP (barcode/oligo level first, then RS-level)")
    print("=" * 70)

    tx_all3_oligo = df["_BS_tx_usable"] & df["_EC_tx_usable"] & df["_PA_tx_usable"]
    tl_all3_oligo = df["_BS_tl_usable"] & df["_EC_tl_usable"] & df["_PA_tl_usable"]
    print(f"OLIGO-level (barcode-level, includes any barcode replicates) three-host usable:")
    print(f"  transcription: {tx_all3_oligo.sum()}")
    print(f"  translation:   {tl_all3_oligo.sum()}")

    # RS-level: collapse by intgen_id (the mined-region identifier -- see barcode
    # analysis below for justification that this is the correct RS-level key)
    rs_key = "intgen_id"
    n_oligo_ids = df["OLIGO ID"].nunique()
    n_rs_ids = df[rs_key].nunique()
    print(f"\nOLIGO ID count: {n_oligo_ids}, unique {rs_key} (RS-level key) count: {n_rs_ids}")

    # An RS has a "usable" value in a host if ANY of its barcode replicates does
    rs_group = df.groupby(rs_key).agg(
        BS_tx_usable_any=("_BS_tx_usable", "any"),
        EC_tx_usable_any=("_EC_tx_usable", "any"),
        PA_tx_usable_any=("_PA_tx_usable", "any"),
        BS_tl_usable_any=("_BS_tl_usable", "any"),
        EC_tl_usable_any=("_EC_tl_usable", "any"),
        PA_tl_usable_any=("_PA_tl_usable", "any"),
        n_barcodes=("OLIGO ID", "count"),
    )
    tx_all3_rs = rs_group["BS_tx_usable_any"] & rs_group["EC_tx_usable_any"] & rs_group["PA_tx_usable_any"]
    tl_all3_rs = rs_group["BS_tl_usable_any"] & rs_group["EC_tl_usable_any"] & rs_group["PA_tl_usable_any"]

    print(f"\nRS-level (collapsed on {rs_key}, 'any barcode has usable value') three-host overlap:")
    print(f"  transcription: {tx_all3_rs.sum()} / {n_rs_ids} RSs")
    print(f"  translation:   {tl_all3_rs.sum()} / {n_rs_ids} RSs")

    # Stricter RS-level version: require ALL barcode replicates usable (conservative)
    rs_group_all = df.groupby(rs_key).agg(
        BS_tx_usable_all=("_BS_tx_usable", "all"),
        EC_tx_usable_all=("_EC_tx_usable", "all"),
        PA_tx_usable_all=("_PA_tx_usable", "all"),
        BS_tl_usable_all=("_BS_tl_usable", "all"),
        EC_tl_usable_all=("_EC_tl_usable", "all"),
        PA_tl_usable_all=("_PA_tl_usable", "all"),
    )
    tx_all3_rs_strict = rs_group_all["BS_tx_usable_all"] & rs_group_all["EC_tx_usable_all"] & rs_group_all["PA_tx_usable_all"]
    tl_all3_rs_strict = rs_group_all["BS_tl_usable_all"] & rs_group_all["EC_tl_usable_all"] & rs_group_all["PA_tl_usable_all"]
    print(f"\nRS-level STRICT (all barcode replicates usable) three-host overlap:")
    print(f"  transcription: {tx_all3_rs_strict.sum()} / {n_rs_ids} RSs")
    print(f"  translation:   {tl_all3_rs_strict.sum()} / {n_rs_ids} RSs")

    results["three_host_overlap"] = {
        "oligo_level": {
            "transcription": int(tx_all3_oligo.sum()),
            "translation": int(tl_all3_oligo.sum()),
            "n_total_oligos": int(n_oligo_ids),
        },
        "rs_level_any_barcode": {
            "transcription": int(tx_all3_rs.sum()),
            "translation": int(tl_all3_rs.sum()),
            "n_total_rs": int(n_rs_ids),
        },
        "rs_level_strict_all_barcodes": {
            "transcription": int(tx_all3_rs_strict.sum()),
            "translation": int(tl_all3_rs_strict.sum()),
            "n_total_rs": int(n_rs_ids),
        },
    }

    # cross-check against paper's stated "8,898 shared set" and "290 active in all 3" for translation
    print(f"\nSanity check against manuscript text:")
    print(f"  Paper states: 'shared set of 8,898 regulatory sequences' with translation data across 3 hosts")
    print(f"  Paper states: '290 constructs (3.3% of library)' expressing GFP (active) in all 3 species")
    print(f"  Our computed oligo-level translation-usable-in-all-3: {tl_all3_oligo.sum()}")

    # ---------------------------------------------------------------
    # Q5: barcode structure
    # ---------------------------------------------------------------
    print("\n" + "=" * 70)
    print("Q5: BARCODE STRUCTURE")
    print("=" * 70)
    barcode_counts = df.groupby(rs_key)["OLIGO ID"].count()
    multi_barcode = barcode_counts[barcode_counts > 1]
    print(f"Total unique {rs_key} (RS) values: {len(barcode_counts)}")
    print(f"RSs with exactly 1 barcode/oligo: {(barcode_counts == 1).sum()}")
    print(f"RSs with >1 barcode/oligo: {len(multi_barcode)}")
    print(f"Distribution of barcodes-per-RS:")
    print(barcode_counts.value_counts().sort_index().to_string())
    print(f"\nSum check: total oligos accounted for = {barcode_counts.sum()} (should equal {n_oligo_ids})")

    results["barcode_structure"] = {
        "n_unique_rs": int(len(barcode_counts)),
        "n_rs_single_barcode": int((barcode_counts == 1).sum()),
        "n_rs_multi_barcode": int(len(multi_barcode)),
        "barcodes_per_rs_distribution": {str(k): int(v) for k, v in barcode_counts.value_counts().sort_index().items()},
    }

    # Also verify: is "Regulatory Sequence" text identical for oligos sharing the same intgen_id?
    print("\nVerifying: do oligos sharing the same intgen_id have identical 'Regulatory Sequence' text?")
    sample_multi = multi_barcode.index[:50]
    mismatches = 0
    for rsid in sample_multi:
        seqs = df.loc[df[rs_key] == rsid, "Regulatory Sequence"].dropna().unique()
        if len(seqs) > 1:
            mismatches += 1
    print(f"  Checked {len(sample_multi)} multi-barcode RS groups; {mismatches} had >1 distinct RS sequence text")

    with open(OUT / "task1_intermediate_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nWrote {OUT / 'task1_intermediate_results.json'}")


if __name__ == "__main__":
    main()
