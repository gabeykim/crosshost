# Data provenance and citation

## Regulatory sequence library and activity values (`three_host.parquet`, `three_host_library.parquet`, `rs241.parquet`, `fold_assignment_FINAL.parquet`)

**Derived from:** Johns, N.I., Blazejewski, T., Gomes, A.L., Wang, H.H. (2018). "Rapid and reliable DNA assembly via ligase cycling reaction." *Nature Methods* 15:323–329. BioProject PRJNA431139.

**What is included:** processed, derived activity values (`tx_norm`, `protein_log10`, usable/active flags) and the 165bp regulatory sequences, computed by the CROSSHOST project from the paper's publicly released supplementary data tables.

**What is NOT included:** the raw Springer Nature-copyrighted supplementary Excel tables themselves. If you need the original tables, obtain them from the publisher via the paper's Nature Methods page or PMC (PMC6065261). This package redistributes only derived values with attribution, consistent with the reuse terms available to us at time of writing — **verify current reuse terms yourself before further redistribution.**

**Cite the original paper** if you use this data:
```
Johns, N.I., Blazejewski, T., Gomes, A.L. et al. Reliable genotyping of
recombinant DNA for gene circuits and genome-scale designs using RCS-PCR.
Nat Methods 15, 323–329 (2018).
```

**Splits (`fold_assignment_FINAL.parquet`) are this project's own work** — genome-blocked assignment plus exact-duplicate and near-duplicate (≥0.85 identity) safety nets. See README.md for the full methodology and `scripts/audit_leakage.py` for the verification checks. Max train-test sequence identity across all 5 folds: 0.8485.

## Genomic features (`hosts_genomic.parquet`)

Computed by the CROSSHOST project from public reference genome assemblies (NCBI RefSeq/GenBank accessions listed per-host in the table) and standard bioinformatics annotation (sigma-factor complement via HMMER+Pfam, tAI via the dos Reis lab method, anti-Shine-Dalgarno via RNAduplex). No third-party data-reuse restriction beyond standard genome-assembly citation norms.

## Physiology features (`hosts_physiology.parquet`)

Derived from PaxDb v6.0 (paxdb.org) for 5 of 6 hosts and PRIDE PXD027874 (Hervey et al., Naval Research Lab) for *V. natriegens* — see the `physiology_source` column for per-host provenance. PaxDb and PRIDE data are both intended for open reuse with attribution; cite PaxDb (Wang et al. 2015, *Proteomics*) and, for the *V. natriegens* proteome specifically, PRIDE accession PXD027874.
