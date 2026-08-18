"""
GATE 2 - Task 2 (final step): Assemble data/hosts_genomic.parquet, the ~40-D
genomic host feature vector, from the outputs of scripts 14-17.

Composition (36 scalar features -- "~40-D" as specified):
  - genome_size_bp, gc_content, n_16s_rrna_copies                     (3)
  - n_trna_genes_total                                                (1)
  - tai_genome_wide_mean                                              (1)
  - anti_sd_rnaduplex_mfe_kcalmol                                     (1)
  - n_sigma_factors_total, n_sigma70_primary, n_sigma_ecf, n_sigma54  (4)
  - n_rnap_core_subunits                                              (1)
  - n_ribosomal_proteins, n_chaperones, n_heme_biosynthesis_genes     (3)
  - gc3_content, enc_nc (effective number of codons)                 (2)
  - tRNA gene copy number per amino-acid isotype (20 std aa + fMet)  (21)
  TOTAL: 3+1+1+1+4+1+3+2+21 = 37 dimensions

Sigma-factor family breakdown uses the HMM-name prefixes from script 17:
  primary  = hit includes Sigma70_r1_1 or Sigma70_r1_2 (region 1 is absent
             from ECF sigmas -- a well-established structural distinction)
  ecf      = hit includes Sigma70_ECF (the dedicated ECF-family HMM) and NOT primary
  sigma54  = hit includes Sigma54_DBD or Sigma54_CBD
"""
import json
import pandas as pd
import numpy as np
from pathlib import Path
from collections import Counter

OUT = Path(__file__).resolve().parent.parent / "out"
DATA = Path(__file__).resolve().parent.parent / "data"

HOSTS = ["EC", "BS", "PA", "SE", "CG", "VN"]

STANDARD_ISOTYPES = ["Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
                     "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "fMet"]


def effective_number_of_codons(codon_counts):
    """Nc (Wright 1990), computed the standard way: for each of the amino-acid
    synonymous codon families, compute homozygosity F, then combine per
    Wright's formula. Approximated here using the standard degeneracy-class
    grouping. STOP codons and Met/Trp (single-codon, degeneracy 1) excluded
    from the F-based classes per convention."""
    # standard genetic code codon->aa table (DNA alphabet)
    from Bio.Data import CodonTable
    table = CodonTable.unambiguous_dna_by_id[11].forward_table  # bacterial code
    aa_to_codons = {}
    for codon, aa in table.items():
        aa_to_codons.setdefault(aa, []).append(codon)

    degeneracy_classes = {2: [], 3: [], 4: [], 6: []}
    for aa, codons in aa_to_codons.items():
        deg = len(codons)
        if deg in degeneracy_classes:
            degeneracy_classes[deg].append(codons)
        # degeneracy 1 (Met, Trp) excluded per Wright 1990

    def homozygosity(codons):
        n = sum(codon_counts.get(c, 0) for c in codons)
        if n <= 1:
            return None
        freqs = [codon_counts.get(c, 0) / n for c in codons]
        F = (n * sum(f ** 2 for f in freqs) - 1) / (n - 1)
        return F

    F_by_class = {}
    for deg, groups in degeneracy_classes.items():
        Fs = [homozygosity(g) for g in groups]
        Fs = [f for f in Fs if f is not None and f > 0]
        F_by_class[deg] = np.mean(Fs) if Fs else None

    nc = 2  # the 2 single-codon families (Met, Trp) each contribute 1, summed as 2
    for deg in [2, 3, 4, 6]:
        F = F_by_class.get(deg)
        if F and F > 0:
            nc += deg / F
        else:
            nc += deg  # if F undefined (e.g. no synonymous variation captured), use max
    return min(nc, 61.0)


def main():
    basic = json.load(open(OUT / "genomic_features_basic.json"))
    antisd_trna = json.load(open(OUT / "antisd_trna_features.json"))
    tai = json.load(open(OUT / "tai_results.json"))
    sigma = json.load(open(OUT / "sigma_factor_hmmer_results.json"))

    rows = []
    for host in HOSTS:
        b = basic[host]
        at = antisd_trna[host]
        t = tai[host]
        s = sigma[host]

        # sigma factor family breakdown
        primary_names = {"Sigma70_r1_1", "Sigma70_r1_2"}
        ecf_names = {"Sigma70_ECF"}
        s54_names = {"Sigma54_DBD", "Sigma54_CBD"}
        n_primary = n_ecf = n_s54 = 0
        for pid, info in s["proteins"].items():
            hmms = set(info["matching_hmms"])
            if hmms & primary_names:
                n_primary += 1
            elif hmms & s54_names:
                n_s54 += 1
            elif hmms & ecf_names:
                n_ecf += 1

        codon_counts = b["codon_counts"]
        total_codons = sum(codon_counts.values())
        gc3 = sum(codon_counts.get(c, 0) for c in codon_counts if c[2] in "GC") / total_codons if total_codons else np.nan
        nc = effective_number_of_codons(codon_counts)

        row = {
            "host": host,
            "genome_size_bp": b["genome_size_bp"],
            "gc_content": b["gc_content"],
            "n_16s_rrna_copies": b["sixteen_s_copy_number"],
            "n_trna_genes_total": at["n_trna_genes_total"],
            "tai_genome_wide_mean": t["genome_wide_tai"],
            "anti_sd_rnaduplex_mfe_kcalmol": at["anti_sd_rnaduplex_mfe_kcalmol"],
            "anti_sd_nonamer": at["anti_sd_nonamer"],
            "n_sigma_factors_total": s["n_sigma_factors"],
            "n_sigma70_primary": n_primary,
            "n_sigma_ecf": n_ecf,
            "n_sigma54": n_s54,
            "n_rnap_core_subunits": len(b["gene_family_hits"]["rnap_core"]),
            "n_ribosomal_proteins": len(b["gene_family_hits"]["ribosomal_proteins"]),
            "n_chaperones": len(b["gene_family_hits"]["chaperones"]),
            "n_heme_biosynthesis_genes": len(b["gene_family_hits"]["heme_biosynthesis"]),
            "gc3_content": gc3,
            "enc_nc": nc,
        }
        for iso in STANDARD_ISOTYPES:
            row[f"trna_n_{iso}"] = at["trna_isotype_counts"].get(iso, 0)

        rows.append(row)

    df = pd.DataFrame(rows)
    print(f"hosts_genomic.parquet: {len(df)} rows, {len(df.columns) - 1} feature columns "
          f"(excluding 'host' identifier)")
    print(df[["host", "genome_size_bp", "gc_content", "n_16s_rrna_copies", "n_trna_genes_total",
              "tai_genome_wide_mean", "anti_sd_rnaduplex_mfe_kcalmol", "n_sigma_factors_total",
              "n_sigma70_primary", "n_sigma_ecf", "n_sigma54", "n_rnap_core_subunits",
              "n_ribosomal_proteins", "n_chaperones", "n_heme_biosynthesis_genes",
              "gc3_content", "enc_nc"]].to_string(index=False))

    out_path = DATA / "hosts_genomic.parquet"
    df.to_parquet(out_path, index=False)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
