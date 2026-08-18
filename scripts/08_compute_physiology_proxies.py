"""
GATE 1.5 - Task B2: Compute uniform omics-derived physiology proxies for all 6 hosts.

IDENTIFICATION METHOD (applied identically to every host -- this uniformity is the
entire point of this approach):

For each protein, we have a gene/preferred name and a free-text functional
annotation (PaxDb: 'preferred_name' + 'annotation'; V. natriegens/PRIDE:
gene symbol parsed from the UniProt-style 'GN=' tag in 'protein_names' +
the protein_names description text itself).

A protein is assigned to a category if EITHER:
  (a) its gene name exactly matches (case-insensitive) a fixed, curated symbol
      list for that category (conserved bacterial nomenclature -- see CATEGORY
      definitions below), OR
  (b) its annotation/description text contains one of a fixed set of keyword
      phrases for that category (case-insensitive substring match).
Matching via (a) OR (b) -- not AND -- because gene-symbol nomenclature is
NOT uniform across our 6 hosts (e.g. sigma factors are rpoD/rpoS/rpoH/rpoE/rpoN
in Gammaproteobacteria but sigA/sigB/sigC... in Firmicutes/Actinobacteria), while
annotation text ("sigma factor") IS uniform. Gene-symbol matching is kept for the
categories where bacterial nomenclature genuinely is near-universal (ribosomal
proteins rps*/rpl*/rpm*, RNAP core rpoA/B/C/Z, chaperones, EF-Tu/EF-G) as a
higher-precision primary signal; annotation-text matching is what makes sigma
factors (and anything gene-symbol matching misses) work uniformly.

Every match is logged (gene name + which rule fired) to out/physiology_proxy_matches.json
for audit -- nothing here is a black box.

CATEGORIES:
  ribosomal_proteins : structural 30S/50S ribosomal proteins (NOT rRNA-modifying
                        enzymes, NOT ribosome assembly/hibernation factors)
  rnap_core           : DNA-directed RNA polymerase core subunits alpha/beta/beta'/omega
                        (rpoA/rpoB/rpoC/rpoZ) -- explicitly EXCLUDES sigma factors
  sigma_factors        : any RNA polymerase sigma factor (primary or alternative/ECF)
  chaperones           : GroEL/GroES (HSP60 system) + DnaK/DnaJ/GrpE (HSP70 system) +
                        trigger factor
  elongation_factors  : EF-Tu (tuf/tufA/tufB) and EF-G (fusA)
"""
import json
import re
import pandas as pd
from pathlib import Path

RAW = Path(__file__).resolve().parent.parent / "raw"
DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"

# (dataset_id_or_label, is_primary)
PAXDB_SOURCES = {
    "EC": [("integrated", "paxdb_abundances_EC_integrated_3740039012.json", True)],
    "BS": [("integrated", "paxdb_abundances_BS_integrated_4187424622.json", True),
           ("abele2025", "paxdb_abundances_BS_abele2025_945251378.json", False)],
    "PA": [("integrated", "paxdb_abundances_PA_integrated_1789377356.json", True),
           ("abele2025", "paxdb_abundances_PA_abele2025_938070019.json", False)],
    "SE": [("integrated", "paxdb_abundances_SE_integrated_1793422220.json", True),
           ("abele2025", "paxdb_abundances_SE_abele2025_1071551919.json", False)],
    "CG": [("integrated", "paxdb_abundances_CG_integrated_3268883846.json", True)],
}

# ---- category definitions ----
GENE_SYMBOL_RULES = {
    "ribosomal_proteins": re.compile(r"^rp[sl][a-z0-9]{1,3}$|^rpm[a-z][0-9]?$", re.IGNORECASE),
    "rnap_core": re.compile(r"^rpo[abcz]$", re.IGNORECASE),
    "sigma_factors": re.compile(r"^(rpod|rpos|rpoh|rpoe|rpon|rpof|fliA|siga|sigb|sigc|sigd|sige|sigf|sigg|sigh|sigi|sigj|sigk|sigl|sigm|sign|sigo|sigp|sigv|sigw|sigx|sigy|sigz)$", re.IGNORECASE),
    "chaperones": re.compile(r"^(groel|groes|grol|gros|dnak|dnaj|grpe|tig)$", re.IGNORECASE),
    "elongation_factors": re.compile(r"^(tuf|tufa|tufb|fusa)$", re.IGNORECASE),
}

ANNOTATION_KEYWORD_RULES = {
    "ribosomal_proteins": [re.compile(r"\b\d{2}s ribosomal (subunit )?protein\b", re.IGNORECASE)],
    "rnap_core": [re.compile(r"dna-directed rna polymerase subunit (alpha|beta|omega|beta')", re.IGNORECASE),
                  re.compile(r"rna polymerase subunit (alpha|beta|omega)", re.IGNORECASE)],
    "sigma_factors": [re.compile(r"sigma factor", re.IGNORECASE),
                       re.compile(r"sigma-70", re.IGNORECASE),
                       re.compile(r"rna polymerase sigma", re.IGNORECASE)],
    "chaperones": [re.compile(r"chaperonin", re.IGNORECASE),
                    re.compile(r"trigger factor", re.IGNORECASE)],
    "elongation_factors": [re.compile(r"elongation factor tu\b", re.IGNORECASE),
                            re.compile(r"elongation factor g\b", re.IGNORECASE)],
}

# explicit exclusions: things that would false-positive on gene-symbol or keyword
# regexes above but are NOT the structural protein/enzyme we want.
# Found by manual audit of scripts/08 output round 1 (see out/physiology_proxy_matches.json
# history / memo comparability section for the specific false positives this fixes):
#   - "anti.{0,20}sigma" catches anti-sigma-factor REGULATORS (e.g. B. subtilis ylaD
#     "anti-YlaC sigma factor", yxlD "sigma-Y antisigma factor component") being
#     wrongly counted as sigma factors themselves.
#   - the enzyme-suffix pattern catches enzymes that MODIFY a ribosomal protein
#     (e.g. E. coli prmA "Methyltransferase for 50S ribosomal subunit protein L11",
#     roxA "50S ribosomal protein L16 arginine hydroxylase") being wrongly counted
#     as the ribosomal protein itself.
EXCLUDE_ANNOTATION = re.compile(
    r"ribosome-binding factor|ribosomal rna|rrna methyltransferase|ribosomal large subunit "
    r"pseudouridine synthase|ribosome hibernation|ribosome silencing|ribosome maturation|"
    r"ribosomal protein.*methyltransferase|sigma factor regulator|"
    r"anti.{0,20}sigma|"
    r"(methyltransferase|hydroxylase|synthase|methylase|acetyltransferase|"
    r"transferase)\s+for\s+\d{2}s\s+ribosomal|"
    r"\d{2}s\s+ribosomal\s+protein\s+\S+\s+(arginine\s+)?"
    r"(hydroxylase|methyltransferase|methylase|acetyltransferase|synthase)",
    re.IGNORECASE)


def primary_clause(annotation_text):
    """PaxDb/UniProt-style annotations consistently put the protein's primary
    name/function as the first clause before the first semicolon, e.g.
    '50S ribosomal subunit protein L23; One of the early assembly proteins,
    it binds 23S rRNA...(mentions trigger factor much later)'. Restricting
    keyword matching to this first clause -- rather than the full free-text
    description -- was added after finding real false positives where a
    protein's EXTENDED description mentioned another category's keyword in
    passing (e.g. ribosomal protein L23's full annotation mentions "trigger
    factor" because L23 physically contacts it at the ribosome exit tunnel,
    which wrongly classified L23 as a chaperone before this fix)."""
    return annotation_text.split(";")[0]


def classify(gene_name, annotation_text):
    """Return the set of categories this protein matches, plus which rule fired."""
    hits = {}
    gene_name = (gene_name or "").strip()
    annotation_text = annotation_text or ""
    clause = primary_clause(annotation_text)

    if EXCLUDE_ANNOTATION.search(clause):
        return hits

    for cat, gene_re in GENE_SYMBOL_RULES.items():
        if gene_name and gene_re.match(gene_name):
            hits[cat] = "gene_symbol"

    for cat, kw_list in ANNOTATION_KEYWORD_RULES.items():
        if cat in hits:
            continue
        for kw_re in kw_list:
            if kw_re.search(clause):
                hits[cat] = "annotation_keyword"
                break
    return hits


def load_paxdb(path):
    with open(RAW / path) as f:
        payload = json.load(f)
    rows = payload["data"]
    df = pd.DataFrame(rows)
    df = df.rename(columns={"preferred_name": "gene_name", "annotation": "annotation_text",
                             "abundance": "abundance_ppm"})
    return df[["gene_name", "annotation_text", "abundance_ppm"]], payload


def load_vnatriegens():
    df = pd.read_parquet(DATA / "vnatriegens_proteome_TP1_S300_T37.parquet")
    gn_re = re.compile(r"GN=(\S+)")

    def extract_gene(row):
        m = gn_re.search(str(row["protein_names"]))
        return m.group(1) if m else ""

    df["gene_name"] = df.apply(extract_gene, axis=1)
    df["annotation_text"] = df["protein_names"]
    return df[["gene_name", "annotation_text", "abundance_ppm"]]


def compute_fractions(df, host, source_label):
    total = df["abundance_ppm"].sum()
    matches = {cat: [] for cat in GENE_SYMBOL_RULES}
    for _, row in df.iterrows():
        hits = classify(row["gene_name"], row["annotation_text"])
        for cat, rule in hits.items():
            matches[cat].append({"gene_name": row["gene_name"], "abundance_ppm": float(row["abundance_ppm"]),
                                  "rule": rule, "annotation": row["annotation_text"][:100]})

    result = {"host": host, "source": source_label, "n_proteins_total": int(len(df)),
              "total_abundance_ppm": float(total)}
    for cat, hit_list in matches.items():
        cat_abundance = sum(h["abundance_ppm"] for h in hit_list)
        result[f"{cat}_n_proteins"] = len(hit_list)
        result[f"{cat}_mass_fraction"] = cat_abundance / total if total > 0 else None
        result[f"{cat}_matched_genes"] = sorted(set(h["gene_name"] for h in hit_list if h["gene_name"]))
    return result, matches


def main():
    all_results = []
    all_matches_log = {}

    for host, sources in PAXDB_SOURCES.items():
        for label, fname, is_primary in sources:
            df, payload = load_paxdb(fname)
            print(f"\n=== {host} / {label} (n={len(df)} proteins, dataset_id={payload['dataset_id']}) ===")
            result, matches = compute_fractions(df, host, label)
            result["is_primary"] = is_primary
            result["dataset_id"] = payload["dataset_id"]
            all_results.append(result)
            all_matches_log[f"{host}_{label}"] = matches
            for cat in GENE_SYMBOL_RULES:
                print(f"  {cat}: n={result[f'{cat}_n_proteins']}, "
                      f"fraction={result[f'{cat}_mass_fraction']:.4f}" if result[f'{cat}_mass_fraction'] is not None else "N/A")

    # V. natriegens
    vn_df = load_vnatriegens()
    print(f"\n=== VN / pride_pxd027874 (n={len(vn_df)} proteins) ===")
    result, matches = compute_fractions(vn_df, "VN", "pride_pxd027874")
    result["is_primary"] = True
    result["dataset_id"] = "PXD027874_TP1_S300_T37"
    all_results.append(result)
    all_matches_log["VN_pride_pxd027874"] = matches
    for cat in GENE_SYMBOL_RULES:
        print(f"  {cat}: n={result[f'{cat}_n_proteins']}, "
              f"fraction={result[f'{cat}_mass_fraction']:.4f}" if result[f'{cat}_mass_fraction'] is not None else "N/A")

    # Save full results table
    results_df = pd.DataFrame(all_results)
    results_df.to_csv(OUT.parent / "data" / "physiology_proxy_results.csv", index=False)
    print(f"\nWrote {DATA / 'physiology_proxy_results.csv'}")

    # Save match audit log (trimmed -- just gene names + rule, not full abundance dumps)
    audit = {}
    for key, matches in all_matches_log.items():
        audit[key] = {cat: [{"gene_name": h["gene_name"], "rule": h["rule"]} for h in hits]
                       for cat, hits in matches.items()}
    with open(OUT / "physiology_proxy_matches.json", "w") as f:
        json.dump(audit, f, indent=2)
    print(f"Wrote {OUT / 'physiology_proxy_matches.json'}")

    # Print final summary matrix
    print("\n" + "=" * 100)
    print("SUMMARY MATRIX (primary source per host)")
    print("=" * 100)
    primary = results_df[results_df["is_primary"]]
    cats = list(GENE_SYMBOL_RULES.keys())
    header = f"{'host':6s}" + "".join(f"{c[:18]:>20s}" for c in cats)
    print(header)
    for _, row in primary.iterrows():
        line = f"{row['host']:6s}"
        for c in cats:
            v = row[f"{c}_mass_fraction"]
            line += f"{v*100:>19.2f}%" if v is not None else f"{'N/A':>20s}"
        print(line)


if __name__ == "__main__":
    main()
