"""
GATE 2 - Task 0 step 5: Reprocess V. natriegens at greater depth.

Gate 1.5 used one specific condition (TP1/S300/T37, 3 replicates) from PRIDE
PXD027874, giving 665 quantified proteins -- a precisely-documented single
condition, but shallow relative to the PaxDb-derived hosts (1,227-5,034).

Checked for deeper alternatives:
  - PXD045789 ("V. natriegens proteome under nutrient limitation"): only a raw
    MS zip file is available for download (no processed protein-level table),
    would require running a search engine ourselves -- out of scope, not used.
  - PXD049476 (chromosome fusion paper): not checked further once PXD027874's
    own pooled depth proved sufficient (see below).
  - PXD027874 ITSELF has 81 samples (3 timepoints x 3 salinities x 3 temps x
    3 replicates). Pooling across ALL 81 samples -- i.e. counting a protein as
    quantified if it has a nonzero LFQ intensity in ANY sample, and taking the
    mean LFQ across the samples where it was detected -- gives 1,032 proteins
    after standard QC filtering (contaminant/reverse/site-only removed), vs
    665 for the single TP1/S300/T37 condition. This is used as the reprocessed,
    deeper V. natriegens proteome.

Trade-off, stated explicitly: this pooled version is no longer tied to one
documented growth condition -- it is now, conceptually, "average abundance
when detected across a temperature/salinity/timepoint panel," which is
actually closer in spirit to how PaxDb's own "Integrated" datasets for the
other 5 hosts are built (weighted averages across many different underlying
studies/conditions) than the single-condition version was. Both versions are
kept: the single-condition file from Gate 1.5 is untouched; this script writes
a new, separate file.
"""
import pandas as pd
from pathlib import Path
import json

RAW = Path(__file__).resolve().parent.parent / "raw"
DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"


def main():
    print(f"Loading {RAW / 'vnat_proteinGroups.txt'} ...")
    df = pd.read_csv(RAW / "vnat_proteinGroups.txt", sep="\t", low_memory=False)
    print(f"Loaded {len(df)} raw protein groups")

    before = len(df)
    df = df[df["Potential contaminant"] != "+"]
    df = df[df["Reverse"] != "+"]
    df = df[df["Only identified by site"] != "+"]
    print(f"After contaminant/reverse/site-only QC: {len(df)} (dropped {before - len(df)})")

    lfq_cols = [c for c in df.columns if c.startswith("LFQ.intensity")]
    print(f"n samples (LFQ columns): {len(lfq_cols)}")

    lfq = df[lfq_cols].replace(0, pd.NA)
    n_detected = lfq.notna().sum(axis=1)
    df["n_samples_detected"] = n_detected
    df["lfq_mean_pooled"] = lfq.mean(axis=1, skipna=True)

    usable = df[df["n_samples_detected"] > 0].copy()
    print(f"Proteins detected (LFQ>0) in at least 1 of {len(lfq_cols)} samples: {len(usable)}")

    total = usable["lfq_mean_pooled"].sum()
    usable["abundance_ppm"] = usable["lfq_mean_pooled"] / total * 1e6

    out = usable[["Gene.names", "Majority protein IDs", "Protein.names",
                   "n_samples_detected", "lfq_mean_pooled", "abundance_ppm"]].copy()
    out.columns = ["gene_names", "protein_ids", "protein_names", "n_samples_detected",
                   "lfq_mean_intensity_pooled", "abundance_ppm"]
    out = out.sort_values("abundance_ppm", ascending=False).reset_index(drop=True)

    out_path = DATA / "vnatriegens_proteome_pooled_81samples.parquet"
    out.to_parquet(out_path, index=False)
    print(f"\nWrote {out_path} ({len(out)} proteins)")

    with open(OUT / "vnat_reprocessing_log.json", "w") as f:
        json.dump({
            "source": "PRIDE PXD027874, ALL 81 samples pooled (all TP/salinity/temp conditions)",
            "n_raw_protein_groups": int(before),
            "n_after_qc": int(len(df)),
            "n_final_usable_any_sample_detected": int(len(out)),
            "comparison_to_single_condition_version": {
                "single_condition_TP1_S300_T37_n_proteins": 665,
                "pooled_all_81_samples_n_proteins": int(len(out)),
                "depth_gain": int(len(out)) - 665,
            },
            "other_prd_datasets_checked": {
                "PXD045789": "only raw MS zip available, no processed protein table -- not used",
                "PXD049476": "not evaluated further once PXD027874 pooled depth (1032) proved sufficient",
            },
        }, f, indent=2)
    print(f"Wrote {OUT / 'vnat_reprocessing_log.json'}")


if __name__ == "__main__":
    main()
