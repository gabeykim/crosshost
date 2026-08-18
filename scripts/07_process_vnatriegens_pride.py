"""
GATE 1.5 - Task B1: Process V. natriegens whole-proteome data from PRIDE PXD027874
(Hervey et al., "Vibrio natriegens Systems & Synthetic Biology: Proteome Profile
(Temperature & Salinity)"), since V. natriegens is ABSENT from PaxDb v6.0 (verified:
scanned the full 405-species list, only "Vibrio proteolyticus" -- a different
species -- is present).

PXD027874 is a MaxQuant-processed (LFQ intensity) factorial design: 3 timepoints
(TP1/TP2/TP3) x 3 added-NaCl levels (S0/S300/S540 mM) x 3 temperatures
(T20/T30/T37 degC) x 3 replicates = 81 samples, from raw/vnat_proteinGroups.txt
(MaxQuant proteinGroups.txt) and raw/vnat_expdesign.txt (sample->condition map).

Representative condition chosen: TP1 (first/earliest timepoint), S300 (300 mM
added NaCl -- closest available level to the ~256 mM/15 g/L NaCl used in Hoffart
et al. 2017's VN minimal medium, the best-characterized V. natriegens growth
condition from Gate 1 Task 4), T37 (37 degC, matching the temperature used for
most other hosts' best-available data). Replicates: Vnat85/86/87 -> LFQ.intensity.85/86/87.
This is a documented, explicit, single condition -- notably MORE precisely
specified than the "Integrated" PaxDb picks used for the other 5 hosts, which are
opaque weighted averages across many undocumented conditions. This asymmetry is
flagged in the memo.
"""
import pandas as pd
from pathlib import Path
import json

RAW = Path(__file__).resolve().parent.parent / "raw"

REP_COLS = ["LFQ.intensity.85", "LFQ.intensity.86", "LFQ.intensity.87"]
REP_LABELS = ["Vnat85", "Vnat86", "Vnat87"]
CONDITION = "TP1.S300.T37"


def main():
    print(f"Loading {RAW / 'vnat_proteinGroups.txt'} ...")
    usecols = ["Gene.names", "Majority protein IDs", "Protein.names",
               "Potential contaminant", "Reverse", "Only identified by site"] + REP_COLS
    df = pd.read_csv(RAW / "vnat_proteinGroups.txt", sep="\t", usecols=usecols, low_memory=False)
    print(f"Loaded {len(df)} protein groups")

    # Standard MaxQuant QC filters
    before = len(df)
    df = df[df["Potential contaminant"] != "+"]
    df = df[df["Reverse"] != "+"]
    df = df[df["Only identified by site"] != "+"]
    print(f"After removing contaminants/reverse-hits/site-only: {len(df)} (dropped {before - len(df)})")

    # Verify chosen condition/replicate mapping against the experimental design file
    expdesign = pd.read_csv(RAW / "vnat_expdesign.txt", sep="\t")
    check = expdesign[expdesign["label"].isin(REP_LABELS)]
    print("\nExperimental design cross-check for chosen replicates:")
    print(check.to_string(index=False))
    assert (check["condition"] == CONDITION).all(), "Condition label mismatch!"

    # Average LFQ intensity across the 3 replicates (0 = not quantified in that run;
    # treat as missing for the mean, matching standard MaxQuant practice)
    rep_df = df[REP_COLS].replace(0, pd.NA)
    df["lfq_mean"] = rep_df.mean(axis=1, skipna=True)
    n_zero_all3 = df["lfq_mean"].isna().sum()
    print(f"\nProteins with 0/NA LFQ in all 3 replicates of {CONDITION} (excluded): {n_zero_all3}")
    df = df.dropna(subset=["lfq_mean"])
    print(f"Proteins with a usable mean LFQ intensity in {CONDITION}: {len(df)}")

    # Convert to a ppm-like relative abundance (matches PaxDb's ppm convention:
    # each protein's share of total summed abundance, x 1e6) for cross-source comparability
    total = df["lfq_mean"].sum()
    df["abundance_ppm"] = df["lfq_mean"] / total * 1e6

    out = df[["Gene.names", "Majority protein IDs", "Protein.names", "lfq_mean", "abundance_ppm"]].copy()
    out.columns = ["gene_names", "protein_ids", "protein_names", "lfq_mean_intensity", "abundance_ppm"]
    out = out.sort_values("abundance_ppm", ascending=False).reset_index(drop=True)

    out_path = RAW.parent / "data" / "vnatriegens_proteome_TP1_S300_T37.parquet"
    out.to_parquet(out_path, index=False)
    print(f"\nWrote {out_path} ({len(out)} proteins)")
    print("\nTop 10 by abundance:")
    print(out.head(10).to_string())

    with open(RAW.parent / "out" / "vnat_processing_log.json", "w") as f:
        json.dump({
            "source": "PRIDE PXD027874",
            "condition": CONDITION,
            "replicate_labels": REP_LABELS,
            "replicate_columns": REP_COLS,
            "n_protein_groups_raw": int(before),
            "n_after_qc_filters": int(len(df) + n_zero_all3 - n_zero_all3),  # after contam/reverse/site filter
            "n_excluded_zero_lfq_all_reps": int(n_zero_all3),
            "n_final_usable": int(len(out)),
        }, f, indent=2)


if __name__ == "__main__":
    main()
