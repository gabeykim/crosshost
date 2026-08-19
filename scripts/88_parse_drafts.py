"""
GATE 10 - Task 2: parse the raw DRAFTS (Yim, Johns, et al. 2019, Mol Syst
Biol 15:e8875) source-data Excel files into clean parquet tables.

TEN-SPECIES ABBREVIATION KEY, verified against Appendix Table S1 and the
main-text species list (both independently, exact-string matches; zero
"aeruginosa" hits anywhere in the full text or appendix):
  Ec=E. coli, Ef=E. fergusonii, Se=S. enterica, Ko=K. oxytoca,
  Pa=P. AGGLOMERANS (NOT P. aeruginosa -- see below), Pp=P. putida,
  Vn=V. natriegens, Bs=B. subtilis, Cg=C. glutamicum, Ll=L. lactis.

**CRITICAL, LOAD-BEARING CORRECTION to a natural assumption:** DRAFTS's
"Pa" column is Pantoea agglomerans, a Proteobacterium isolated from
agricultural waste -- NOT Pseudomonas aeruginosa, one of this project's
three primary hosts. Confirmed three independent ways: (1) Appendix Table
S1 lists "P. agglomerans" at the exact column position "Pa" occupies in
the Fig. 3 sheet; (2) the full-text Results section names all ten species
explicitly by full binomial, "P. agglomerans" among them; (3) a full-text
and appendix-wide grep for "aeruginosa" returns zero matches. **This
project's *P. aeruginosa* is ABSENT from DRAFTS entirely.** Only *E. coli*
and *B. subtilis*, of this project's three primary hosts, appear in
DRAFTS's ten-species panel.

SOURCE FILES:
  raw/drafts/msb198875_SourceData_Fig3.xlsx, sheet "Fig. 3" -- the
    10-species x up-to-1047-sequence cell-free transcription matrix.
    Per-species columns: {sp}_rep{1,2}_{DNA,RNA} (per-replicate raw
    counts), {sp}_{DNA,RNA} (pooled), {sp}_tx (transcription level,
    linear-scale ratio), {sp}_zscore (log10(tx) z-scored, computed
    per-species independently), {sp}_zscore_union (z-score restricted to
    the 421-sequence "union set" usable in ALL 10 species -- this is the
    paper's own headline "421 sequences active in all ten" figure,
    independently reproduced exactly: 421/421 non-null for every species).
  raw/drafts/msb198875_SourceData_Fig2.xlsx, sheet "Appendix Fig. S5" --
    the RS234 in-vitro (DRAFTS cell-free) vs in-vivo comparison, 234
    sequences x 7 species (Ec, Se, Pp, Vn, Ko, Bs, Cg -- NOT Pa/Ef/Ll;
    matches the abstract's "seven species for in-vivo comparison" and,
    separately, still does not include P. aeruginosa).

USABILITY, per-species, per-sequence: the "_tx"/"_zscore" columns are NOT
NaN for unusable rows -- they contain STRING SENTINELS
('no_DNA_counts', 'no_RNA_counts', 'low_DNA_counts') in place of a numeric
value. pandas' default dtype inference silently treats these columns as
object/string dtype and `.notna()` returns True for the sentinel strings
too -- a real trap, caught here by explicit `pd.to_numeric(errors='coerce')`
before any usability check. Per-species "usable" = numeric value present
(not one of the three sentinel strings).
"""
import re
import numpy as np
import pandas as pd
from pathlib import Path

RAW = Path(__file__).resolve().parent.parent / "raw" / "drafts"
DATA = Path(__file__).resolve().parent.parent / "data"

SPECIES = ["Ec", "Ef", "Se", "Ko", "Pa", "Pp", "Vn", "Bs", "Cg", "Ll"]
SPECIES_FULL = {
    "Ec": "Escherichia coli", "Ef": "Escherichia fergusonii", "Se": "Salmonella enterica",
    "Ko": "Klebsiella oxytoca", "Pa": "Pantoea agglomerans", "Pp": "Pseudomonas putida",
    "Vn": "Vibrio natriegens", "Bs": "Bacillus subtilis", "Cg": "Corynebacterium glutamicum",
    "Ll": "Lactococcus lactis",
}
SPECIES_PHYLUM = {
    "Ec": "Proteobacteria", "Ef": "Proteobacteria", "Se": "Proteobacteria", "Ko": "Proteobacteria",
    "Pa": "Proteobacteria", "Pp": "Proteobacteria", "Vn": "Proteobacteria",
    "Bs": "Firmicutes", "Ll": "Firmicutes", "Cg": "Actinobacteria",
}
RS234_SPECIES = ["Ec", "Se", "Pp", "Vn", "Ko", "Bs", "Cg"]  # confirmed from the actual sheet columns


def to_numeric_with_sentinel_flag(series):
    """Returns (numeric_values, is_usable_bool, sentinel_reason). Sentinel
    strings ('no_DNA_counts','no_RNA_counts','low_DNA_counts') and any
    other non-numeric entries become NaN in the numeric series and False
    in usable; the original string is preserved as the reason where
    present."""
    numeric = pd.to_numeric(series, errors="coerce")
    is_string = series.apply(lambda x: isinstance(x, str))
    usable = numeric.notna()
    reason = series.where(is_string & ~usable, other=pd.NA)
    return numeric, usable, reason


def parse_main_matrix():
    df = pd.read_excel(RAW / "msb198875_SourceData_Fig3.xlsx", sheet_name="Fig. 3")
    assert len(df) == 1047, f"expected 1047 rows, got {len(df)}"

    out = pd.DataFrame({
        "oligo_id": df["Oligo ID"].astype(int),
        "sequence_165bp": df["Sequence_165bp"].str.upper().str.strip(),
        "gc_pct": df["GC_sequence_165bp"].astype(float),
        "source_species": df["species"], "source_genus": df["genus"],
        "source_family": df["family"], "source_order": df["order"],
        "source_class": df["class"], "source_phylum": df["phylum"],
        "source_genome": df["Genome"],
    })
    assert out["oligo_id"].duplicated().sum() == 0, "duplicate oligo_id in Fig. 3 sheet"
    assert out["sequence_165bp"].duplicated().sum() == 0, "duplicate sequence in Fig. 3 sheet"
    assert (out["sequence_165bp"].str.len() == 165).all(), "not all sequences are 165bp"

    for sp in SPECIES:
        dna, dna_usable, dna_reason = to_numeric_with_sentinel_flag(df[f"{sp}_DNA"])
        rna, rna_usable, _ = to_numeric_with_sentinel_flag(df[f"{sp}_RNA"])
        tx, tx_usable, tx_reason = to_numeric_with_sentinel_flag(df[f"{sp}_tx"])
        zscore, zscore_usable, _ = to_numeric_with_sentinel_flag(df[f"{sp}_zscore"])
        zscore_union = pd.to_numeric(df[f"{sp}_zscore_union"], errors="coerce")

        out[f"{sp}_DNA"] = dna
        out[f"{sp}_RNA"] = rna
        out[f"{sp}_tx"] = tx
        out[f"{sp}_zscore"] = zscore
        out[f"{sp}_zscore_union"] = zscore_union
        out[f"{sp}_usable"] = tx_usable  # a real numeric tx value exists
        out[f"{sp}_usable_union_set"] = zscore_union.notna()  # active-in-all-10-hosts subset
        out[f"{sp}_unusable_reason"] = tx_reason  # 'no_DNA_counts' / 'no_RNA_counts' / 'low_DNA_counts', else NaN

    n_union = out[[f"{sp}_usable_union_set" for sp in SPECIES]].all(axis=1).sum()
    assert n_union == 421, f"expected 421 sequences usable in all 10 species (paper's headline figure), recomputed {n_union}"
    print(f"[VERIFIED, recomputed independently] N usable in all 10 species (union set) = {n_union}")

    for sp in SPECIES:
        n_usable = out[f"{sp}_usable"].sum()
        print(f"  {sp} ({SPECIES_FULL[sp]}): {n_usable}/1047 usable "
              f"({out[f'{sp}_unusable_reason'].value_counts().to_dict()})")

    return out


def parse_rs234_invivo():
    df = pd.read_excel(RAW / "msb198875_SourceData_Fig2.xlsx", sheet_name="Appendix Fig. S5")
    assert len(df) == 234, f"expected 234 rows, got {len(df)}"

    out = pd.DataFrame({
        "oligo_id": df["Oligo ID"].astype(int),
        "sequence_full": df["Sequence+ATG+BC(rev)"],
    })
    for sp in RS234_SPECIES:
        for modality in ["invitro", "invivo"]:
            dna, dna_usable, _ = to_numeric_with_sentinel_flag(df[f"{sp}_{modality}_DNA"])
            rna, rna_usable, _ = to_numeric_with_sentinel_flag(df[f"{sp}_{modality}_RNA"])
            tx, tx_usable, _ = to_numeric_with_sentinel_flag(df[f"{sp}_{modality}_tx"])
            zscore, _, _ = to_numeric_with_sentinel_flag(df[f"{sp}_{modality}_zscore"])
            out[f"{sp}_{modality}_DNA"] = dna
            out[f"{sp}_{modality}_RNA"] = rna
            out[f"{sp}_{modality}_tx"] = tx
            out[f"{sp}_{modality}_zscore"] = zscore
            out[f"{sp}_{modality}_usable"] = tx_usable

    print(f"\nRS234 in-vivo/in-vitro comparison: {len(out)} sequences x {len(RS234_SPECIES)} species "
          f"({', '.join(RS234_SPECIES)}) x 2 modalities x 2 replicates")
    for sp in RS234_SPECIES:
        n_both = (out[f"{sp}_invitro_usable"] & out[f"{sp}_invivo_usable"]).sum()
        print(f"  {sp}: n usable in BOTH modalities = {n_both}/234")
    return out


def main():
    matrix = parse_main_matrix()
    matrix.to_parquet(DATA / "drafts.parquet", index=False)
    print(f"\nWrote {DATA / 'drafts.parquet'} ({matrix.shape[0]} rows, {matrix.shape[1]} cols)")

    rs234 = parse_rs234_invivo()
    rs234.to_parquet(DATA / "drafts_rs234_invivo.parquet", index=False)
    print(f"Wrote {DATA / 'drafts_rs234_invivo.parquet'} ({rs234.shape[0]} rows, {rs234.shape[1]} cols)")


if __name__ == "__main__":
    main()
