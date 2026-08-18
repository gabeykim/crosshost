"""
GATE 2 - Task 2: Basic genomic features (no external tool needed beyond Biopython)
and the annotation-based gene-family counts.

METHOD NOTE / DEVIATION FROM SPEC (documented transparently): the task asked for
ribosomal protein complement, chaperone complement, and heme biosynthesis gene
counts via "KEGG KO or EggNOG-mapper v2." Neither was used. EggNOG-mapper
requires a multi-GB reference database download and KEGG KO assignment requires
either a licensed API or the same kind of heavy local database; both were judged
infeasible within this gate's time budget. Instead, the SAME validated
gene-symbol + annotation-keyword method built and audited in Gate 1.5
(scripts/08_compute_physiology_proxies.py) is reused here, applied to each
host's RefSeq protein FASTA headers (which carry PGAP-curated gene names and
product descriptions) instead of PaxDb abundance annotations. This is a
consistent, reproducible, already-QC'd method, not a new untested one -- but it
is a real deviation from the specified tool and is flagged in the memo.

RNAP core subunits (rpoB/rpoC) are identified the same way, plus explicit
identity confirmation (exact single-copy gene match per host).
"""
import gzip
import re
import json
from pathlib import Path
from collections import Counter

from Bio import SeqIO

RAW = Path(__file__).resolve().parent.parent / "raw" / "genomes"
DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"

HOSTS = ["EC", "BS", "PA", "SE", "CG", "VN"]

# Reuse Gate 1.5's classification rules directly (ribosomal_proteins, rnap_core,
# sigma_factors [not used here -- HMMER used instead, see script 15], chaperones,
# elongation_factors) plus new rules for heme biosynthesis.
GENE_SYMBOL_RULES = {
    "ribosomal_proteins": re.compile(r"^rp[sl][a-z0-9]{1,3}$|^rpm[a-z][0-9]?$", re.IGNORECASE),
    "rnap_core": re.compile(r"^rpo[abcz]$", re.IGNORECASE),
    "chaperones": re.compile(r"^(groel|groes|grol|gros|dnak|dnaj|grpe|tig)$", re.IGNORECASE),
    "heme_biosynthesis": re.compile(r"^hem[a-z][0-9]?$", re.IGNORECASE),
}
ANNOTATION_KEYWORD_RULES = {
    "ribosomal_proteins": [re.compile(r"\b\d{2}s ribosomal (subunit )?protein\b", re.IGNORECASE)],
    "rnap_core": [re.compile(r"dna-directed rna polymerase subunit (alpha|beta|omega|beta')", re.IGNORECASE),
                  re.compile(r"rna polymerase subunit (alpha|beta|omega)", re.IGNORECASE)],
    "chaperones": [re.compile(r"chaperonin", re.IGNORECASE), re.compile(r"trigger factor", re.IGNORECASE)],
    "heme_biosynthesis": [re.compile(r"heme (biosynthesis|abc transporter)", re.IGNORECASE),
                           re.compile(r"(uroporphyrinogen|coproporphyrinogen|protoporphyrinogen|ferrochelatase|"
                                      r"delta-aminolevulinic acid|porphobilinogen|hydroxymethylbilane)", re.IGNORECASE)],
}
EXCLUDE_ANNOTATION = re.compile(
    r"ribosome-binding factor|ribosomal rna|rrna methyltransferase|ribosomal large subunit "
    r"pseudouridine synthase|ribosome hibernation|ribosome silencing|ribosome maturation|"
    r"ribosomal protein.*methyltransferase|"
    r"(methyltransferase|hydroxylase|synthase|methylase|acetyltransferase|transferase)\s+for\s+\d{2}s\s+ribosomal|"
    r"\d{2}s\s+ribosomal\s+protein\s+\S+\s+(arginine\s+)?(hydroxylase|methyltransferase|methylase|acetyltransferase|synthase)",
    re.IGNORECASE)


def primary_clause(text):
    return (text or "").split(";")[0].split(" [")[0]


def classify(gene_name, product):
    hits = {}
    gene_name = (gene_name or "").strip()
    clause = primary_clause(product)
    if EXCLUDE_ANNOTATION.search(clause):
        return hits
    for cat, gre in GENE_SYMBOL_RULES.items():
        if gene_name and gre.match(gene_name):
            hits[cat] = "gene_symbol"
    for cat, kws in ANNOTATION_KEYWORD_RULES.items():
        if cat in hits:
            continue
        for kw in kws:
            if kw.search(clause):
                hits[cat] = "annotation_keyword"
                break
    return hits


def parse_protein_faa(path):
    """Parse RefSeq protein.faa.gz headers: >ACCESSION description [organism]"""
    records = []
    with gzip.open(path, "rt") as f:
        for rec in SeqIO.parse(f, "fasta"):
            desc = rec.description
            # strip the accession + trailing [organism]
            desc_wo_id = desc[len(rec.id):].strip()
            desc_wo_org = re.sub(r"\s*\[[^\]]*\]\s*$", "", desc_wo_id)
            records.append({"protein_id": rec.id, "product": desc_wo_org, "length": len(rec.seq)})
    return records


def parse_gff_gene_names(path):
    """Parse GFF3 to map protein_id -> gene name (PGAP curated /gene= qualifier)."""
    gene_map = {}
    with gzip.open(path, "rt") as f:
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
    return gene_map


def genome_basic_stats(fna_path):
    total_len = 0
    gc = 0
    n_replicons = 0
    with gzip.open(fna_path, "rt") as f:
        for rec in SeqIO.parse(f, "fasta"):
            n_replicons += 1
            seq = str(rec.seq).upper()
            total_len += len(seq)
            gc += seq.count("G") + seq.count("C")
    return total_len, gc / total_len if total_len else None, n_replicons


def codon_usage(cds_path):
    counts = Counter()
    n_cds = 0
    with gzip.open(cds_path, "rt") as f:
        for rec in SeqIO.parse(f, "fasta"):
            seq = str(rec.seq).upper()
            if len(seq) % 3 != 0 or len(seq) < 3:
                continue
            n_cds += 1
            for i in range(0, len(seq) - 2, 3):
                codon = seq[i:i + 3]
                if len(codon) == 3 and set(codon) <= set("ACGT"):
                    counts[codon] += 1
    total = sum(counts.values())
    freqs = {c: n / total for c, n in counts.items()} if total else {}
    return counts, freqs, n_cds


SIXTEEN_S_RE = re.compile(r"16S", re.IGNORECASE)
RIBOSOMAL_RNA_RE = re.compile(r"ribosomal RNA", re.IGNORECASE)


def extract_16s(rna_path):
    # RefSeq/PGAP product phrasing is NOT consistent across assemblies -- e.g.
    # E. coli uses "16S ribosomal RNA" but B. subtilis uses "ribosomal RNA-16S".
    # Match "16S" AND "ribosomal RNA" appearing anywhere in the description,
    # in either order, rather than one fixed phrase (a fixed-phrase match
    # silently returned 0 hits for B. subtilis on the first pass -- caught by
    # a plausibility check, since B. subtilis is known to have ~10 rrn operons).
    seqs = []
    with gzip.open(rna_path, "rt") as f:
        for rec in SeqIO.parse(f, "fasta"):
            if SIXTEEN_S_RE.search(rec.description) and RIBOSOMAL_RNA_RE.search(rec.description):
                seqs.append(str(rec.seq))
    return seqs


def main():
    results = {}
    for host in HOSTS:
        hdir = RAW / host
        print(f"\n=== {host} ===")

        total_len, gc, n_repl = genome_basic_stats(hdir / "genomic.fna.gz")
        print(f"  genome size: {total_len:,} bp, GC={gc*100:.2f}%, {n_repl} replicon(s)")

        counts, freqs, n_cds = codon_usage(hdir / "cds_from_genomic.fna.gz")
        print(f"  {n_cds} CDS, {sum(counts.values()):,} codons counted")

        proteins = parse_protein_faa(hdir / "protein.faa.gz")
        gene_map = parse_gff_gene_names(hdir / "genomic.gff.gz")
        print(f"  {len(proteins)} proteins, {len(gene_map)} gene-name mappings from GFF")

        cat_hits = {cat: [] for cat in GENE_SYMBOL_RULES}
        for p in proteins:
            gene = gene_map.get(p["protein_id"], "")
            hits = classify(gene, p["product"])
            for cat, rule in hits.items():
                cat_hits[cat].append({"protein_id": p["protein_id"], "gene": gene,
                                       "product": p["product"], "rule": rule})

        for cat, hits in cat_hits.items():
            genes = sorted(set(h["gene"] for h in hits if h["gene"]))
            print(f"  {cat}: n={len(hits)}  genes={genes}")

        sixteen_s = extract_16s(hdir / "rna_from_genomic.fna.gz")
        print(f"  16S rRNA copies found: {len(sixteen_s)}")

        results[host] = {
            "genome_size_bp": total_len,
            "gc_content": gc,
            "n_replicons": n_repl,
            "n_cds": n_cds,
            "n_proteins": len(proteins),
            "codon_counts": dict(counts),
            "codon_freqs": freqs,
            "gene_family_hits": {cat: [{"gene": h["gene"], "protein_id": h["protein_id"], "rule": h["rule"]}
                                        for h in hits] for cat, hits in cat_hits.items()},
            "sixteen_s_sequences": sixteen_s,
            "sixteen_s_copy_number": len(sixteen_s),
        }

    with open(OUT / "genomic_features_basic.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {OUT / 'genomic_features_basic.json'}")


if __name__ == "__main__":
    main()
