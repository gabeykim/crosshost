"""
GATE 2.5 - Task A3: annotation-quality check on the keyword-derived genomic
features (Gate 2's disclosed deviation from KEGG KO / EggNOG-mapper).

Concern: annotation completeness varies across hosts (E. coli exhaustively
curated; V. natriegens, C. glutamicum less so), and keyword matching on
uneven annotation could produce systematically lower counts for
less-annotated organisms that masquerade as biological differences.
"""
import gzip
import re
import json
from pathlib import Path
from collections import Counter

from Bio import SeqIO

RAW = Path(__file__).resolve().parent.parent / "raw" / "genomes"
OUT = Path(__file__).resolve().parent.parent / "out"

HOSTS = ["EC", "BS", "PA", "SE", "CG", "VN"]

# "uninformative" annotation patterns: generic locus-tag-style names or
# explicit placeholder product descriptions
GENERIC_GENE_PATTERNS = re.compile(r"^(rs|jw|stm|bsu|pa|cg|b\d+)?\d+[a-z]?$", re.IGNORECASE)
HYPOTHETICAL_PRODUCT_PATTERNS = re.compile(
    r"hypothetical protein|uncharacterized protein|DUF\d+|putative protein\b|unknown function",
    re.IGNORECASE)


def parse_protein_faa_and_gff(host):
    hdir = RAW / host
    gene_map = {}
    with gzip.open(hdir / "genomic.gff.gz", "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attrs = dict(kv.split("=", 1) for kv in parts[8].split(";") if "=" in kv)
            protein_id = attrs.get("protein_id")
            gene = attrs.get("gene", "")
            if protein_id:
                gene_map[protein_id] = gene

    records = []
    with gzip.open(hdir / "protein.faa.gz", "rt") as f:
        for rec in SeqIO.parse(f, "fasta"):
            desc = rec.description
            desc_wo_id = desc[len(rec.id):].strip()
            desc_wo_org = re.sub(r"\s*\[[^\]]*\]\s*$", "", desc_wo_id)
            records.append({"protein_id": rec.id, "product": desc_wo_org,
                             "gene": gene_map.get(rec.id, "")})
    return records


def is_informative_gene_symbol(gene):
    if not gene:
        return False
    if GENERIC_GENE_PATTERNS.match(gene):
        return False
    return True


def is_hypothetical_product(product):
    return bool(HYPOTHETICAL_PRODUCT_PATTERNS.search(product))


def main():
    basic = json.load(open(OUT / "genomic_features_basic.json"))

    results = {}
    print(f"{'Host':6s} {'n_proteins':>10s} {'n_gene_symbol':>14s} {'gene_ratio':>10s} "
          f"{'n_hypothetical':>15s} {'hypo_ratio':>10s}")
    for host in HOSTS:
        records = parse_protein_faa_and_gff(host)
        n_total = len(records)
        n_gene = sum(1 for r in records if is_informative_gene_symbol(r["gene"]))
        n_hypo = sum(1 for r in records if is_hypothetical_product(r["product"]))
        gene_ratio = n_gene / n_total
        hypo_ratio = n_hypo / n_total
        print(f"{host:6s} {n_total:>10d} {n_gene:>14d} {gene_ratio:>10.3f} {n_hypo:>15d} {hypo_ratio:>10.3f}")
        results[host] = {"n_proteins": n_total, "n_with_gene_symbol": n_gene, "gene_symbol_ratio": gene_ratio,
                          "n_hypothetical_or_uncharacterized": n_hypo, "hypothetical_ratio": hypo_ratio}

    gene_ratios = [results[h]["gene_symbol_ratio"] for h in HOSTS]
    hypo_ratios = [results[h]["hypothetical_ratio"] for h in HOSTS]
    print(f"\nGene-symbol ratio range: {min(gene_ratios):.3f} - {max(gene_ratios):.3f} "
          f"(spread: {max(gene_ratios)-min(gene_ratios):.3f})")
    print(f"Hypothetical-product ratio range: {min(hypo_ratios):.3f} - {max(hypo_ratios):.3f} "
          f"(spread: {max(hypo_ratios)-min(hypo_ratios):.3f})")

    # correlation with feature counts
    print(f"\n{'Host':6s} {'gene_ratio':>10s} {'ribosomal':>10s} {'rnap':>6s} {'chaperone':>10s} {'heme':>6s}")
    feature_rows = []
    for host in HOSTS:
        b = basic[host]
        n_rib = len(b["gene_family_hits"]["ribosomal_proteins"])
        n_rnap = len(b["gene_family_hits"]["rnap_core"])
        n_chap = len(b["gene_family_hits"]["chaperones"])
        n_heme = len(b["gene_family_hits"]["heme_biosynthesis"])
        gr = results[host]["gene_symbol_ratio"]
        print(f"{host:6s} {gr:>10.3f} {n_rib:>10d} {n_rnap:>6d} {n_chap:>10d} {n_heme:>6d}")
        feature_rows.append({"host": host, "gene_ratio": gr, "n_ribosomal": n_rib, "n_rnap": n_rnap,
                              "n_chaperone": n_chap, "n_heme": n_heme})
        results[host]["feature_counts"] = {"ribosomal_proteins": n_rib, "rnap_core": n_rnap,
                                            "chaperones": n_chap, "heme_biosynthesis": n_heme}

    # simple correlation (Pearson) between gene_ratio and each feature count, n=6
    import numpy as np
    gr_arr = np.array([r["gene_ratio"] for r in feature_rows])
    for feat in ["n_ribosomal", "n_rnap", "n_chaperone", "n_heme"]:
        vals = np.array([r[feat] for r in feature_rows])
        if vals.std() > 0 and gr_arr.std() > 0:
            corr = np.corrcoef(gr_arr, vals)[0, 1]
        else:
            corr = float("nan")
        print(f"  Pearson r(gene_ratio, {feat}) = {corr:.3f} (n=6, eyeball only)")
        results[f"correlation_gene_ratio_vs_{feat}"] = float(corr) if not np.isnan(corr) else None

    with open(OUT / "annotation_quality_check.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nWrote {OUT / 'annotation_quality_check.json'}")


if __name__ == "__main__":
    main()
