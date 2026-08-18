# CROSSHOST — Gate 1 Reality Check Memo

**Date:** 2026-08-03
**Prepared for:** Gabriel
**Purpose:** Verify five load-bearing numbers before any modeling code is written. No models were built, no training pipeline was written. Everything below is either [VERIFIED] against a primary source, [COMPUTED] from verified inputs with arithmetic shown, [INFERRED] with confidence stated, or [NOT FOUND] with the search trail documented.

**Bottom line up front: PASS.** All three Gate 1 pass conditions are met and none of the four kill criteria trigger. The one item that does NOT cleanly resolve — physiology-data availability for the genome-vs-physiology scientific comparison — is not a Gate 1 kill criterion per the charter, but it is materially weaker than the charter implies and needs your decision before Gate 2. See Task 4 and "Decision needed" at the end.

---

## How the data was obtained

The Johns et al. 2018 *Nature Methods* paper (DOI 10.1038/nmeth.4633, PMID 30052624) is behind Nature's paywall, but the NIH-mandated author manuscript is deposited in PubMed Central as **PMC6065261** (NIHMS945382). The PMC Open Access API (`https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=PMC6065261`) resolved to a tarball at `ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/deprecated/oa_package/1c/27/PMC6065261.tar.gz` (NCBI is mid-migration off the old FTP layout as of this writing; the `deprecated/` prefix was required). That tarball contains the full manuscript text (`nihms945382.nxml`), the PDF, and **all five supplementary data tables as XLSX**, plus the supplementary methods/figures DOCX. This is a complete, non-paywalled, primary-source copy of everything needed for Task 1 and Task 2. Nothing was estimated or taken from a review/secondary summary.

All raw files are in `raw/`. All processing scripts are in `scripts/` and are re-runnable (`source .venv/bin/activate && python scripts/0N_*.py`).

---

## TASK 1 — Dataset shape

### Q6 — Column inventory (answered first; needed context for everything else)

Five supplementary XLSX files, each with a "Read Me" sheet the authors wrote themselves. [VERIFIED — read directly from `raw/NIHMS945382-supplement-{3,4,5,6,7}.xlsx`]

**Supplementary Data Table 1** (`supplement-3.xlsx`) — Oligo Library Info, 29,249 rows:
| Sheet | Columns | Notes |
|---|---|---|
| Metadata | OLIGO ID (int), intgen_id (int, intergenic-region ID), genome_id (int, IMG genome ID), [nc_id] (str, NCBI accession), intgen_lp/intgen_rp (int, mined-region genomic coordinates), intgen_td (str, strand), gene_id (int), gene_lrp (str, gene coordinates), gene_td (str), gene_cogid (str, COG ID), gene_cogcat (str, COG functional category), core count (int, # genomes the homolog is found in), [gene_cogdesc] (str), phylum/class/family/genus/order (str), [strain_name] (str), Regulatory Sequence (str, the raw 165 bp sequence) | 21 columns, 0 nulls except structurally-bracketed optional fields, all populated |
| Full Constructs | Oligo ID (int), Regulatory Sequence + ATG + 12bp barcode (str) | 29,248 rows — one OLIGO ID (43726) present in Metadata is missing here |
| Phylums | Kingdom, Phylum, # species per phylum | 26 rows |

**Supplementary Data Table 2** (`supplement-4.xlsx`) — BS/EC/PA expression, 29,249 rows each:
| Sheet | Columns |
|---|---|
| BS, EC, PA (one sheet per host) | OLIGO ID (int), rna_count (int), dna_count (int), primary TSS (int), fraction_around_primary_TSS (float), tx_raw (float, raw transcription = rna/dna), tx_norm (float, host-normalized transcription), protein (log10) (float, FACS-seq-derived protein level), delta_G (float, UTR folding energy), best_sigma70_match_score (float) |
| Phylum Fraction Active | Phylum, # Organisms, # RS in "union set", Avg. GC, Fraction on BS/EC/PA (>0 RNA reads) | 26 rows — this is the source of the paper's headline active-fraction percentages, see Q4 |

**Supplementary Data Table 3** (`supplement-5.xlsx`) — *E. coli*-only condition-dependent transcription (LB, NaCl-stress, Fe-starvation, LB-stationary, M9-minimal; 29,249 rows each) plus a 100-row "Robust RSs" sheet. Not needed for the 3-host or RS241 questions; noted for completeness since it's part of the same paper.

**Supplementary Data Table 4** (`supplement-6.xlsx`) — **RS241**, 241 rows: `id` (Oligo ID) + Log2 Transcription and Log10 Translation values for BS, CG, PA, VN, SE, EC. See Task 2.

**Supplementary Data Table 5** (`supplement-7.xlsx`) — cloning oligos, strain list (with growth temperatures), gene sequences, and 10 "SsGC" synthetic genetic circuit constructs. Not part of the core prediction dataset; the Strains sheet's growth-temperature column was cross-checked in Task 4.

### Q3 — Part (regulatory sequence) length

**[VERIFIED, three independent ways, all in agreement.]**

1. Direct measurement of the `Regulatory Sequence` column in Metadata: **165 bp exactly, for all 29,249 rows, zero variance** (`min=165, median=165, max=165, std=0.000`).
2. Cross-check via `Full Constructs` (RS + ATG[3bp] + 12bp barcode), sequence length minus 15: **165 bp exactly, for all 29,248 rows.**
3. The manuscript's own Methods text: *"we...extracted the 165 bp immediately upstream of annotated start codons. These sequences will be referred to as RSs for convenience."* (`raw/nihms945382_bodytext.txt`)

**The charter's assumed value of 165 bp is CORRECT.** Part length is uniform, not variable — no distribution/histogram is needed beyond "single spike at 165." CNN kernel widths sized to 165 bp require no change.

### Q1 — Three-host overlap (the load-bearing number)

**Definition used**, taken directly from the paper's own Methods text (not invented): *"We excluded constructs containing 0 DNA counts and also those whose RNA and DNA counts summed to less than 15 for most analyses."* This is the paper's own QC rule for a usable transcription value. Applied per-host, per-oligo:

```
tx_usable = (dna_count != 0) AND (rna_count + dna_count >= 15)
```

For translation, no equivalent per-oligo exclusion rule is stated in the text or reconstructable from the released columns, so "usable" = `protein (log10)` is not null (i.e., the FACS-seq pipeline produced a value for that construct in that host).

| Readout | Definition | RSs usable in **all three** hosts (BS+EC+PA) | vs. ≥5,000 pass bar | vs. ~3,000 kill bar |
|---|---|---|---|---|
| **Transcription** | `dna≠0 AND rna+dna≥15`, all 3 hosts | **11,276 / 29,249** | ✅ clears by 2.3× | ✅ clears by 3.8× |
| **Translation** | `protein(log10)` not null, all 3 hosts | **7,887 / 29,249** | ✅ clears by 1.6× | ✅ clears by 2.6× |

[COMPUTED — full arithmetic and code in `scripts/03_task1_task2_analysis.py`, run output archived in `out/task1_intermediate_results.json`]

Because RS-level deduplication turned out to be a non-issue for this table (see Q5 below — every OLIGO ID already maps 1:1 to a unique regulatory sequence), the oligo-level and RS-level numbers are **identical**: 11,276 and 7,887 either way.

**Cross-validation against the paper's own aggregate statistics:**
- The paper's per-phylum "Phylum Fraction Active" sheet reports a "# Regulatory Sequence in union set" total of **11,265**. A weighted average of that sheet's per-phylum active fractions reproduces the paper's headline 18.9%/52.0%/83.8% figures almost exactly (see Q4). The 11,265 figure is 11 short of our independently-computed 11,276 — [INFERRED, high confidence] the gap is explained by the "Host_associated" phylum category (45 RSs) which is present in the raw Metadata but absent from the 26-row Phylum Fraction Active table; a handful of those 45 evidently pass the 3-host QC filter. **This is strong independent confirmation that our 11,276 transcription figure and the paper's own "union set" are the same underlying quantity**, computed two different ways.
- The paper states translation data form *"a shared set of 8,898 regulatory sequences"* across the three hosts. Our directly-computed 7,887 is 89% of that. [COMPUTED, with an acknowledged unresolved gap] We could not find an additional per-oligo QC column that would close this gap exactly — it's plausible the paper applied an unstated additional filter to the protein data (e.g., also requiring `dna_count>0`, which alone removes several thousand null-protein rows per host, see script output) that isn't fully specified in the released Methods text or columns. **We report our own number (7,887) as the reproducible, code-derived figure and flag the paper's 8,898 as a close but not exactly reconciled aggregate.** Either number clears the pass bar by a wide margin, so this does not affect the Gate 1 decision.

### Q2 — Per-host usable and active counts

[COMPUTED — `scripts/03_task1_task2_analysis.py`]

| Host | Transcription usable (paper's QC rule) | Transcription "active" (rna_count>0, full-library denominator) | Translation usable (protein not null) |
|---|---|---|---|
| *B. subtilis* | 15,848 / 29,249 (54.2%) | 5,223 (17.9% of 29,249) | 11,689 / 29,249 (39.9%) |
| *E. coli* | 24,613 / 29,249 (84.1%) | 16,575 (56.7% of 29,249) | 28,226 / 29,249 (96.5%) |
| *P. aeruginosa* | 21,473 / 29,249 (73.4%) | 18,564 (63.5% of 29,249) | 19,735 / 29,249 (67.5%) |

### Q4 — *B. subtilis* active fraction (and the other two, for comparison)

**Important nuance the charter's ~19% number does not convey on its own, but which does not overturn it:** the paper's headline active-fraction statistic — *"B. subtilis displayed the lowest number of measurably active RSs (18.9% with >0 RNA reads), while E. coli and P. aeruginosa had substantially higher fractions... (52.0% and 83.8% respectively)"* [VERIFIED, quoted directly from manuscript body text] — is computed over the **11,265-RS "union set"** (i.e., specifically over the subset of RSs that already have usable transcription data in all three hosts — the same population as our Q1 answer), **not over the full 29,249-member library.**

When we compute the same "active" definition (`rna_count>0`) over the **full library** instead, the numbers shift meaningfully for two of the three hosts:

| Host | Paper's headline (union-set denominator, n=11,265) | Our full-library recomputation (n=29,249) |
|---|---|---|
| *B. subtilis* | **18.9%** [VERIFIED, manuscript text] | 17.9% [COMPUTED] — close, consistent |
| *E. coli* | 52.0% [VERIFIED] | 56.7% [COMPUTED] — same order |
| *P. aeruginosa* | 83.8% [VERIFIED] | 63.5% [COMPUTED] — **20 points lower** |

[INFERRED, high confidence] The union-set denominator is enriched for RSs that already produced usable signal in all three hosts simultaneously, which mechanically biases it toward higher activity — this effect is largest for *P. aeruginosa* because its full-library active fraction (63.5%) is furthest from its union-set fraction (83.8%).

**For the charter's specific question — is *B. subtilis* usable as a held-out host? — the number that matters is essentially unchanged (17.9–18.9%) regardless of which denominator you use, so the charter's ~19% figure is CONFIRMED, and the practical conclusion (B. subtilis is real but the hardest host) holds either way.** But anyone reporting the EC/PA numbers in a paper needs to be explicit about which denominator is in play — we recommend always stating "N usable transcription/translation values" alongside any activity percentage in Gate 2+ outputs.

### Q5 — Barcode structure — **contradicts a charter assumption**

**[VERIFIED, and this is the most consequential correction in this memo for Gate 2 planning.]**

The charter (Part III, Gate 2; Part V rule 1) assumes: *"the library contains double-barcoded replicates; barcode-level splitting would leak."* This assumption traces to the manuscript text, which does describe a barcode-replicate QC analysis: *"we randomly selected a subset of 4,778 RSs from the total library to encode a different set of 12 bp barcodes as an internal control"* and reports *"duplicate RSs with alternate barcodes (Pearson r = 0.88 and 0.86 respectively)"* for transcription/translation concordance, with a supplementary figure caption (extracted from the DOCX) specifying *"n = 2,273"* barcode pairs compared in *E. coli*.

**But in the publicly released Supplementary Data Table 1 (the 29,249-oligo Metadata/Full Constructs sheets actually used for the 3-host expression data in Table 2), this barcode-replicate structure is not present:**
- Every one of the 29,249 `OLIGO ID`s maps to a **globally unique** `intgen_id` (verified: `n unique intgen_id == n unique (genome_id, intgen_id) == 29,249`, i.e. no coordinate-based collisions).
- Every `OLIGO ID` maps to a Regulatory Sequence string that is unique in **98.6%** of cases. The exceptions: **181 groups of exactly-identical 165 bp sequence text, covering 413 oligos total (1.4% of the library)**, spread across *different* `intgen_id`/genomic loci — [INFERRED, medium-high confidence] this reads as biological coincidence (the same short regulatory motif independently mined from closely related strains/genomes in the 184-genome mining set), not the deliberate barcode-duplicate design described in the text.
- The magnitude doesn't match either: 181 groups / 413 oligos is far short of the "4,778 RSs" / "n=2,273 pairs" the text describes.

**[INFERRED, medium confidence] Most likely explanation:** the barcode-replicate QC analysis (Suppl. Fig. S3, Pearson r=0.88/0.86, n=2,273 pairs in *E. coli*) was performed on a separate, earlier/rawer dataset — plausibly closer to the initial 34,027-oligo synthesized pool mentioned in the text — that was not itself included in these five released XLSX supplementary files. The **released** 29,249-member library used for the main BS/EC/PA/RS241 expression data is, empirically, already RS:barcode 1:1.

**Practical consequence for Gate 2:** the specific barcode-deduplication step the charter calls for (dedup before splitting to prevent leakage) is not actually needed to prevent *barcode*-replicate leakage in this table, because there is essentially no barcode-replicate structure to leak from. It would still be good practice to dedup on exact sequence text (catching the 181/413 coincidental-match group) before building folds, for the ordinary reason that near-identical training/test sequences are a leakage risk regardless of *why* they're identical — but this is a much smaller and different concern than the charter describes, and MMseqs2 clustering (already planned for Gate 2) would catch it anyway via sequence identity, independent of any barcode/ID-based dedup logic.

### Manuscript-verified QC/pseudo-value details (for Gate 2 ingestion spec)

Since Q6 asked for enough detail to write an ingestion spec, these methodological facts are recorded here too, all [VERIFIED] from manuscript body text:
- Constructs with 0 DNA reads and >15 RNA reads (135 BS / 373 EC / 172 PA constructs) were given a pseudo-value equal to the **maximum** value in the shown range for visualization purposes — interpreted by the authors as likely high-expression constructs that dropped out due to fitness cost, not true zeros.
- Constructs with 0 RNA reads and >15 DNA reads were given a pseudo-value equal to the **minimum** value in the shown range (transcriptionally inactive).
- These pseudo-values were used **only for the paper's own visualizations** ("for Figures 2A, 3A, and 4B"); the excluded/threshold rule above (`dna≠0 AND rna+dna≥15`) is what the paper calls out as governing "most analyses." We used the latter, not the pseudo-value scheme, for all counts in this memo — recommend the same for Gate 2 (pseudo-values are a plotting convenience, not real measurements, and should not enter a training set as if they were).
- Translation efficiency (a distinct derived quantity from raw protein level) is defined as *"the ratio of the measured transcription rate by the GFP protein levels"* and was only computed for RSs meeting: single-TSS classification (>80% of reads within ±5bp of median TSS), ≥1 read of each of RNA and DNA, and total (RNA+DNA) reads >15 — i.e., a stricter version of the transcription QC rule. This governs the *derived* translation-efficiency metric, not the raw `protein (log10)` column used above.
- Anti-Shine-Dalgarno sequences used for the RBS Calculator baseline (relevant to Gate 3's biophysical baseline): `ACCTCCTTA` for *E. coli*/*P. aeruginosa*, `ACCTCCTTT` for *B. subtilis*.

---

## TASK 2 — RS241 integrity

**[VERIFIED] 241 sequences confirmed** — `raw/NIHMS945382-supplement-6.xlsx`, "Log2 Transcription" and "Log10 Translation" sheets, both 241 data rows, author-written Read Me sheet confirms: *"This file contains two sheets containing expression data for transcription and translation for our small library RS241 in 6 bacterial species."* Manuscript text confirms the design: *"we selected 241 library members (RS241 library), cloned and introduced them into additional industrially useful hosts Salmonella enterica, Vibrio natriegens... and Corynebacterium glutamicum."*

**Coverage is genuinely ragged — not a clean 241×6 matrix.** "Usable" = the cell is not null (these are already paper-processed log-scale values; no further threshold is reconstructable from this table).

| Host | Transcription usable | Translation usable |
|---|---|---|
| *E. coli* | 219 / 241 (90.9%) | 230 / 241 (95.4%) |
| ***B. subtilis*** | **149 / 241 (61.8%)** | **149 / 241 (61.8%)** |
| *P. aeruginosa* | 240 / 241 (99.6%) | 235 / 241 (97.5%) |
| *S. enterica* | 237 / 241 (98.3%) | 237 / 241 (98.3%) |
| *V. natriegens* | 166 / 241 (68.9%) | 229 / 241 (95.0%) |
| *C. glutamicum* | 211 / 241 (87.6%) | 218 / 241 (90.5%) |
| **All 6 hosts simultaneously usable** | **111 / 241 (46.1%)** | **126 / 241 (52.3%)** |

[COMPUTED — `scripts/04_task2_rs241_analysis.py`, output archived in `out/task2_rs241_results.json`]

**Key, non-obvious finding: the limiting host in RS241 is *B. subtilis* — one of the original three primary hosts — not any of the three newly-added hosts.** *S. enterica*, *V. natriegens*, and *C. glutamicum* (the hosts the charter specifically worried about) are all reasonably well covered (68.9–98.3%). This is consistent with the Task 1 finding that *B. subtilis* is the hardest host to get signal in generally, and it means the charter's stated worry ("a subset where only two of the extra three hosts have usable data") does not describe the actual failure mode — the real constraint is B. subtilis coverage, in both the main library and RS241.

**Practical implication for the charter's H-MAIN held-out-*B.-subtilis* evaluation:** if you require all 6 hosts present per sequence for a clean apples-to-apples RS241 evaluation set, you have **111 sequences** (transcription) or **126** (translation) to evaluate on, not 241 — still workable for a held-out-host test, but smaller than "241" suggests at face value, and this should be stated explicitly in Gate 5's evaluation methodology.

Distribution of "how many of the 6 hosts have a usable value" per RS (transcription): 2 hosts→11 RSs, 3→16, 4→29, 5→74, 6→111. Not a small tail — a third of RS241 (85/241) has data in fewer than 5 of the 6 hosts.

**ID-space note (not a data-quality problem, just documented):** 34 of the 241 RS241 `id`s fall outside the main library's OLIGO ID range (14478–43726); RS241 IDs run up to 48469. [INFERRED, high confidence] Consistent with "selected 241 library members" — these were evidently assigned IDs from a later/separate synthesis batch, not re-using the original 29,249-oligo numbering exactly. This does not affect usability; RS241 sequences can still be matched to the main library by exact sequence text where needed.

---

## TASK 3 — "From Context to Code" disposition

Full report: `out/task3_context_to_code.md`. Summary:

- **Identity [VERIFIED, bioRxiv official metadata API]:** "From Context to Code: Rational De Novo DNA Design and Predicting Cross-Species DNA Functionality Using Deep Learning Transformer Models," Dahiya, Bakken, Fages-Lartaud, Lale. bioRxiv 10.1101/2023.10.15.562386, posted 2023-10-15, **v1 only, never updated in ~3 years.**
- **Affiliation [VERIFIED]:** Norwegian University of Science and Technology (NTNU) + **Syngens AS**, an NTNU spinout. **This is explicitly commercial** — the paper's own competing-interests statement discloses that two authors (Dahiya, Lale) are Syngens co-founders and the other two are Syngens employees. Syngens markets an "AI-powered DNA Design Platform" ("Imagine ChatGPT, but for designing DNA sequences") and cites this preprint as its supporting science.
- **Publication status [VERIFIED]:** Not published in any peer-reviewed venue as of 2026-08-03 — bioRxiv's own Crossref-linked "published" field returns `"NA"`, confirmed by Google Scholar and the lab's own publications page.
- **Code/data/model release [VERIFIED]:** None found. Syngens' GitHub org has zero public repositories. No Zenodo/HuggingFace artifacts located.
- **Uses Johns et al. 2018? [VERIFIED for citation / NOT FOUND for data use]:** Cited in the reference list as general background on cross-species regulatory-sequence variability. No mention of PRJNA431139 or any indication the dataset was used for training/testing.
- **The critical question — does it test transfer to an unseen host? [INFERRED, medium-high confidence: NO.]** Systematic full-text search for the vocabulary a leave-one-host-out study would need — "generalize," "unseen host," "held-out," "leave-one-out," "zero-shot," "cross-host," "transfer" — returned **zero matches for all of them.** The paper's framing throughout is "host-specific" design across a fixed roster of ~6 hosts (*B. subtilis, C. glutamicum, E. coli, P. putida, V. natriegens, S. venezuelae*) that all appear in both training/platform-description and validation roles. Its one explicit "novel data" holdout test (Figure 7) is explicitly defined by the authors as novel *sequences*, not novel *hosts*. Confidence is not maximal because the retrieved Methods section (biorxiv.org direct fetch was blocked with HTTP 403; a proxy renderer was used instead) is thin and doesn't give an explicit per-host train/test breakdown — but no evidence for unseen-host testing was found despite deliberate, repeated searching.

**Verbatim abstract:**
> Synthetic biology currently operates under a framework dominated by trial-and-error approaches, which hinders the effective engineering of organisms and the expansion of large-scale biomanufacturing. Motivated by the success of computational designs in areas like architecture and aeronautics, we aspire to transition to a more efficient and predictive methodology in synthetic biology. In this study, we report a DNA Design Platform that relies on the predictive power of Transformer-based deep learning architectures. The platform transforms the conventional paradigms in synthetic biology by enabling the context-sensitive and host-specific engineering of 5′ regulatory elements—promoters and 5′ untranslated regions (UTRs) along with an array of codon-optimised coding sequence (CDS) variants. This allows us to generate context-sensitive 5′ regulatory sequences and CDSs, achieving an unparalleled level of specificity and adaptability in different target hosts. With context-aware design, we significantly broaden the range of possible gene expression profiles and phenotypic outcomes, substantially reducing the need for laborious high-throughput screening efforts. Our context-aware, AI-driven design strategy marks a significant advancement in synthetic biology, offering a scalable and refined approach for gene expression optimisation across a diverse range of expression hosts. In summary, this study represents a substantial leap forward in the field, utilising deep learning models to transform the conventional design, build, test, learn-cycle into a more efficient and predictive framework.

**Assessment: adjacent work, not a direct competitor.** It targets host-specific design within a fixed, trained-on host roster and is not shown to test transfer to a genuinely unseen host — the opposite of this project's specific claim. It is also a commercial platform with no public artifacts, which further limits direct overlap with a public benchmark. The kill criterion ("From Context to Code already published this benchmark") is **NOT TRIGGERED.**

---

## TASK 4 — Physiology data availability (the scientific core — the weakest finding in this memo)

Full reports: `out/task4_physiology_primary3.md` (E. coli, B. subtilis, P. aeruginosa) and `out/task4_physiology_secondary3.md` (S. enterica, V. natriegens, C. glutamicum).

### First, the Bernstein-lab-specific question

**[VERIFIED] The Bernstein lab's own chassis-effect papers do not cover this project's host panel at all.** Two papers were read in full:
- Chan, Baldwin & Bernstein, *BioDesign Research* 2023 (PMC10432152) — hosts: *E. coli* + *Halopseudomonas aestusnigri*, *H. oceani*, *Pseudomonas deceptionensis*, *P. fluorescens*, *P. putida*. All grown in LB, 30°C.
- Chan et al., *mSystems* 2024 (PMC11406997) — hosts: six *Stutzerimonas* species, with *E. coli* used only as a cloning reference, not an experimental chassis.

Both papers do publish physiology tables (growth rate, carrying capacity, GC content, plasmid copy number, codon adaptation index) — but **their host panel never includes *B. subtilis* or *P. aeruginosa*, so no Bernstein-lab table is directly reusable for this project's three (or six) target hosts.** The charter's suggestion to check whether Bernstein's own papers "tabulate exactly what is needed" is answered: no, physiology data for our specific hosts has to be assembled independently from the broader literature.

### The 6-host × 6-metric matrix

Legend: ●●● = growth-rate-resolved, directly measured, primary source read in full. ●● = measured but not growth-rate-resolved, or existence-verified with numeric table not extractable (paywall). ● = only a related/proxy quantity found. ○ = NOT FOUND.

| Metric | *E. coli* | *B. subtilis* | *P. aeruginosa* | *S. enterica* | *V. natriegens* | *C. glutamicum* |
|---|:---:|:---:|:---:|:---:|:---:|:---:|
| Growth rate / doubling time | ●●● | ●●● | ●●● | ●●● | ●●● | ●●● |
| Ribosomes/cell | ●●● | ● | ○ | ○ | ●●● | ●●● |
| RNA polymerase/cell | ●●● | ○ | ○ | ○ | ○ | ○ |
| RNA/protein ratio | ●●● | ●●● | ● (different metric) | ●● (paywalled) | ●● | ●●● |
| Proteome allocation to translation (φ_R) | ●●● | ●●● | ○ | ○ | ●●● | ●● (exists, not extracted) |
| Other growth-law params | ●●● | ●●● | ● | ●● | ●●● | ●●● |

**Every host has at least a verified growth rate. RNA polymerase abundance is essentially unpublished for every host except *E. coli* — recommend dropping it from the feature set entirely rather than imputing it.** Beyond growth rate, data quality splits sharply into two tiers:

- **Strong (3+ metrics, growth-rate-resolved):** *E. coli*, *B. subtilis*, *V. natriegens*, and (with one metric's numbers not yet extracted) *C. glutamicum*.
- **Thin (growth rate only, everything else NOT FOUND or paywalled/unverifiable):** ***P. aeruginosa* and *S. enterica*** — exactly two of this project's six hosts, and notably **one of the three primary hosts** (*P. aeruginosa*) is in the thin tier.

**A genuine bright spot:** Zhu, Mori, Hwa & Dai, *PNAS* 2025 (122:e2427091122) grew *E. coli*, *B. subtilis*, **and** *V. natriegens* side-by-side in **identical media and temperatures** (LB and multiple minimal media, 37°C primary, plus 30/43°C) and measured RNA/protein ratio and proteome-sector allocation directly comparably across all three. This is the single strongest, most directly comparable multi-host physiology dataset found anywhere in this search — it happens to cover 3 of the 6 target hosts with genuinely matched conditions. Its central finding (E. coli and B. subtilis converge on nearly identical translation-proteome allocation, ~40% in LB, despite being phylogenetically distant) is itself in the spirit of the Bernstein claim this project is testing.

### Comparability warning

**[VERIFIED, both physiology sub-reports independently flag the same problem]** Even where numbers exist for multiple hosts, growth conditions are often not matched:
- *C. glutamicum*'s best (and essentially only) quantitative growth-law dataset (Nat. Commun. 2023, PMC10497606 — ribosome counts, R/P ratio, translation elongation rate, all growth-rate-resolved) is measured at **30°C**, its physiological optimum. Every other host's best data is at **37°C**. There is no 37°C growth-law dataset for *C. glutamicum* and no 30°C dataset for the others in the sources found.
- Medium composition varies within and across hosts (rich/complex vs. minimal-defined), and "M9 minimal" recipes differ in exact composition between the *E. coli* and *B. subtilis* literatures even when both are nominally "M9-based, 37°C."
- *S. enterica*'s foundational datasets (Schaechter/Maaløe/Kjeldgaard 1958; Rosset/Julien/Monier 1966 — the papers that *originated* the RNA/protein growth-law concept) exist and were confirmed to exist, but sit behind Microbiology Society/Elsevier paywalls that could not be penetrated in this pass — their numeric tables are **[NOT FOUND]**, not "don't exist." A follow-up institutional-access attempt could plausibly recover them.

### Verdict against the charter's own kill bar

The charter's explicit kill criterion for this concern lives at **Gate 2**, not Gate 1: *"Physiology features are unrecoverable for ≥3 of 6 hosts → the central comparison is dead."* On the evidence gathered: every host has at least one usable, verified metric (growth rate), so by a strict "zero data" reading, **0 of 6 hosts are fully unrecoverable — this specific Gate 2 kill bar is NOT TRIGGERED.** But the charter's own framing in this Gate 1 prompt says physiology "is now the scientific core," and the honest picture is: a rich, multi-metric physiology vector is buildable for 4 of 6 hosts (E. coli, B. subtilis, V. natriegens, C. glutamicum), and only a thin, growth-rate-only vector is available for the other 2 (P. aeruginosa, S. enterica) — one of which is a primary host the H-MAIN gate depends on for the "easy proteobacterial hop" partial-pass path. **This does not kill Gate 1, but it is materially weaker than "physiology data is recoverable" and it changes what the Gate 2 physiology feature vector can actually contain.** See "Decision needed," below.

---

## TASK 5 — Funder sweep

Full report: `out/task5_funder_sweep.md`. Five sources queried with all 8 specified terms (~500 core records manually scanned), plus a direct DARPA BAA check.

| Source | Method | Records scanned | KILL_CANDIDATE | Notes |
|---|---|---|---|---|
| NIH RePORTER | Full `advanced_text_search` (title+terms+abstract), FY2024–26 filter — **worked on first attempt, no fallback needed** | 91 | 0 | 3 ADJACENT (broad-host-range/biocontainment engineering, single-host) |
| NSF Awards API | `keyword=` search | 200 (+250 diagnostic) | 0 | **[Documented limitation]** multi-word keyword queries don't reliably filter; single-word fallback used, flagged as weakened |
| DOE / Agile BioFoundry | No award-index API exists — **[VERIFIED, NOT FOUND — no fetchable index]**; fell back to OSTI.gov (publications, not awards) | 2 (+25 supplementary) | 0 | 2 ADJACENT, incl. a 2025 Agile BioFoundry "Tier System" chassis-standardization framework paper |
| CORDIS (EU) | Native project-search API | 43 | 0 | 4 ADJACENT (single-host chassis engineering) |
| UKRI Gateway to Research | `q=` keyword search | 160 | 0 | **[Documented limitation]** same multi-word degradation as NSF; 1 close-but-ADJACENT match |

**No KILL_CANDIDATE found anywhere.** The single closest conceptual match across all five sources: an active EPSRC-funded PhD studentship at Imperial College London ("Graph learning methods for engineering mammalian promoters in bioproduction," 2023–2027) building a GNN to predict *context-dependent* promoter activity — but the "context" is different mammalian cell lines/bioproduction conditions, not cross-bacterial-species prediction. Classified **ADJACENT**, not a competitor, because the domain (mammalian) doesn't overlap this project's bacterial scope.

**DARPA BTO AIxBio, BAA HR001126S0003 [VERIFIED, three independent sources — primary BAA PDF, grants.gov, darpa.mil]:** Still open. **Deadline: September 30, 2026, 4:00 PM ET**, rolling abstracts and full proposals, posted October 1, 2025. The charter's working assumption is correct. Note: the BAA's formal title is just "Biological Technologies" — "AIxBio" is informal shorthand for its AI+biology topic area, not its official name.

**Funder-competitor kill criterion: NOT TRIGGERED.**

---

## Pass conditions and kill criteria — evaluated

### PASS requires all three:

| Condition | Status | Evidence |
|---|---|---|
| ≥5,000 RSs with usable three-host values | ✅ **MET** | Transcription: 11,276. Translation: 7,887. Both readouts individually clear 5,000. |
| *B. subtilis* has enough active sequences to serve as a held-out host | ✅ **MET** | 15,848 usable transcription values (54.2% of library), 5,223 active (rna>0). In RS241, 149/241 (61.8%) usable in both readouts. Not the strongest host, but far from unscoreable. |
| No funded direct competitor | ✅ **MET** | Zero KILL_CANDIDATEs across NIH, NSF, DOE, CORDIS, UKRI (~500 records scanned). |

**→ PASS condition fully met.**

### KILL if any of these:

| Criterion | Status | Evidence |
|---|---|---|
| Three-host overlap below ~3,000 sequences | ❌ **NOT TRIGGERED** | 11,276 (tx) / 7,887 (tl), both 2.6–3.8× the kill bar |
| *B. subtilis* active fraction too low to score | ❌ **NOT TRIGGERED** | 17.9–18.9% active (paper-verified and independently recomputed), 5,223 active sequences absolute — same order of magnitude as the other two hosts, not a degenerate near-zero case |
| A funded competitor with institutional standing is building this | ❌ **NOT TRIGGERED** | Confirmed via 5-source funder sweep |
| "From Context to Code" already published this benchmark | ❌ **NOT TRIGGERED** | Commercial, host-specific-within-trained-hosts design tool; no evidence of unseen-host transfer testing; no public artifacts |

**→ No kill criterion triggers. GATE 1 OUTCOME: PASS.**

---

## Everything else worth knowing that no single task question asked for

- **Translation efficiency** (a *derived* quantity — GFP normalized by transcription level) is a distinct column/concept from the raw `protein (log10)` value used throughout this memo, and it has its own, stricter QC rule (single-TSS classification + ≥1 read of each of RNA/DNA + total reads >15). If Gate 4's translation head is trained on "translation efficiency" rather than raw protein level, the usable-count arithmetic above needs to be redone under that stricter rule — the numbers in this memo are for the raw protein/GFP level, which the charter's Gate 4 spec (two output heads: "transcription and translation... not a ratio") actually calls for, so this is very likely a non-issue, just flagged for precision.
- Full-construct sequences confirm the vector chemistry described in the text: BamHI/PstI cut sites + start codon (ATG) + 12bp Levenshtein-distance->2 barcode, appended to each 165bp RS.
- Anti-Shine-Dalgarno sequences (needed for Gate 2's ~40-D genomic feature vector and Gate 3's RBS Calculator baseline) are given explicitly in the Methods: `ACCTCCTTA` (*E. coli*, *P. aeruginosa*), `ACCTCCTTT` (*B. subtilis*).
- *P. aeruginosa* PAO1 test strain carries a Δ*psy2* deletion specifically to remove pyocin S2 autofluorescence — a strain-construction detail worth carrying into any Gate 2 strain-metadata table.
- The paper's own "four general groups" finding — universally active (16.9%), differentially active in 2/3 species (33.3%), specifically active in 1 species (37.4%), inactive in all 3 (12.4%) — is stated in the text but its exact denominator wasn't independently re-derived in this pass; worth reproducing in Gate 2 if this categorical breakdown is useful for the paper's framing.

---

## WHAT I COULD NOT DO

- **Reconcile the translation three-host-overlap number exactly.** Our own reproducible computation (7,887, `protein(log10)` not null in all 3 hosts) is 89% of the paper's stated aggregate ("shared set of 8,898"). The gap is unresolved — plausibly an unstated additional QC filter (e.g. also requiring `dna_count>0`) that isn't fully specified in the available Methods text or the released columns. Does not affect the Gate 1 decision (both numbers clear the pass bar), but should be resolved before Gate 2 locks in a canonical "usable translation value" definition, ideally by re-reading the full Methods PDF/DOCX more carefully or contacting the authors.
- **Extract the exact B. subtilis LB/37°C growth rate digit** from Zhu et al. 2025 PNAS — the paper's main text confirms the experiment (B. subtilis grown alongside E. coli, same media/temps) and shows the comparison in Figure 2, but the specific numeric value is in a Supplementary Appendix table that WebFetch could not access in this pass.
- **Recover the original Schaechter/Maaløe/Kjeldgaard 1958 and Rosset/Julien/Monier 1966 numeric tables for S. enterica** — both papers are confirmed to exist and are exactly on-topic (they originated the RNA/protein growth-law concept, using this exact organism), but sit behind Microbiology Society and Elsevier paywalls that returned HTTP 403 to every fetch attempt. This is the single most promising avenue for closing the *S. enterica* physiology gap and would be worth a manual/institutional-access pass before concluding *S. enterica* physiology is unrecoverable.
- **Confirm the exact growth medium for the Aiyar et al. 2002 *V. natriegens* ribosome count (115,000/cell)** — the value itself is well-corroborated via two independent secondary citations, but the original paper is behind an ASM Journals login wall, so the specific medium is [INFERRED, low-medium confidence] rather than confirmed word-for-word.
- **Access biorxiv.org directly for "From Context to Code."** Direct WebFetch returned HTTP 403 for both the abstract page and the full PDF; all full-text claims about that paper were obtained via the official bioRxiv metadata API (for structured fields) plus a third-party read-through proxy (`r.jina.ai`) for prose sections. This means the "zero mentions of unseen-host-transfer vocabulary" finding is bounded by what that proxy rendering surfaced and could in principle miss content in a part of the PDF the proxy didn't render (e.g., deep supplementary methods) — hence that finding is tagged [INFERRED, medium-high confidence] rather than [VERIFIED].
- **Query a dedicated DOE Office of Science or Agile BioFoundry award index** — no such public API/index exists; fell back to OSTI.gov, which indexes publications/technical reports rather than the underlying awards. This is a weaker proxy than the other four funder sources and is flagged as such in Task 5.
- **Get fully clean, non-degraded multi-word search behavior from the NSF Awards API and UKRI Gateway to Research** — both empirically only filter reliably on single distinctive words; multi-word phrase queries fall through to near-unfiltered, date-sorted result sets. Documented with control-query evidence in `out/task5_funder_sweep.md`; single-word fallback queries were run as a partial mitigation and are flagged everywhere they were used.
- **Access sam.gov's DARPA BAA listing directly** — it's JavaScript-rendered and did not return usable content via WebFetch. Not a gap in the actual finding (three other independent sources — the primary BAA PDF, grants.gov, and darpa.mil — all agree on the September 30, 2026 deadline), just a documented access failure for one of the four sources checked.

---

## CONTRADICTIONS WITH THE CHARTER

1. **Barcode-replicate deduplication is not the leakage risk the charter describes, for this specific released table.** The charter assumes the library "contains double-barcoded replicates; barcode-level splitting would leak." In the actual released Supplementary Data Table 1/2 (29,249 oligos), every OLIGO ID maps 1:1 to a unique genomic locus and — in 98.6% of cases — a unique sequence. The manuscript's described barcode-duplicate QC set (4,778 RSs, n=2,273 pairs in E. coli, Pearson r=0.88/0.86) is not reconstructable from the released tables and most likely lived in a separate, unreleased raw dataset. **Practical effect: smaller than expected, not larger** — this makes Gate 2's data-prep step easier than planned, not harder, but the specific mechanism described in the charter needs correcting before it's written into code (dedup on exact sequence text to catch a small 181-group/413-oligo coincidental-match set, rather than on a nonexistent barcode-pair key).

2. **Physiology data is recoverable, but unevenly, and the unevenness lands partly on a primary host.** The charter frames physiology-data availability as a binary Gate 2 kill check ("≥3 of 6 hosts unrecoverable"). The actual picture is not binary: every host has *some* usable data (so the literal Gate 2 kill bar is not triggered), but *P. aeruginosa* — one of the three **primary** hosts, not one of the three "extra" RS241 hosts — turns out to be one of the two weakest-characterized hosts physiologically (growth rate only; no ribosome, RNAP, RNA/protein, or φ_R data found). The charter's Task 4 prompt anticipated this risk only in the abstract ("this determines whether the comparison is runnable at all") without flagging that the weak link might be a primary host rather than an RS241-only host. This doesn't kill the project, but it does mean the genome-vs-physiology comparison (the paper's stated scientific core) will likely need to be scoped to a 4-host subset (E. coli, B. subtilis, V. natriegens, C. glutamicum) for the metrics beyond growth rate, with P. aeruginosa and S. enterica carried only on growth rate — a real, if partial, renegotiation of the headline as the charter itself anticipated might be necessary.

3. **The Bernstein lab's own physiology tables are not reusable for this project's hosts**, despite the charter's Task 4 prompt suggesting this be checked as a likely shortcut ("they may have tabulated exactly what is needed"). Their published host panel (E. coli + various Pseudomonads/Halopseudomonads/Stutzerimonas) never includes B. subtilis or P. aeruginosa. This is a clarification/narrowing of a charter assumption, not a contradiction of a stated fact, but it forecloses a shortcut the charter held open.

4. **Two charter numbers previously flagged as "unverified" are now independently confirmed correct**, which is worth stating plainly since the charter explicitly asked for this: **165 bp part length** (three independent confirmations: direct sequence measurement, construct-length arithmetic, and the manuscript's own stated methodology) and **~19% B. subtilis active fraction** (18.9% in the paper's own union-set framing, 17.9% in our full-library recomputation — stable either way). These were not contradicted; they held up.

---

## FILES WRITTEN

- `out/GATE1_MEMO.md` — this memo
- `out/state.json` — machine-readable carry-forward state
- `out/task1_intermediate_results.json` — raw computed values backing Task 1
- `out/task2_rs241_results.json` — raw computed values backing Task 2
- `out/task3_context_to_code.md` — full "From Context to Code" research report
- `out/task4_physiology_primary3.md` — full physiology report, E. coli/B. subtilis/P. aeruginosa
- `out/task4_physiology_secondary3.md` — full physiology report, S. enterica/V. natriegens/C. glutamicum
- `out/task5_funder_sweep.md` — full funder-sweep report, all 5 sources + DARPA BAA check
- `data/three_host_library.parquet` — 29,249-row merged core library table (Metadata + Full Constructs + BS/EC/PA expression), built from raw XLSX
- `data/rs241.parquet` — 241-row RS241 six-host table
- `scripts/01_inspect_supplements.py` — sheet/column inspector for the raw XLSX files
- `scripts/02_build_core_tables.py` — builds the two Parquet tables above from raw XLSX
- `scripts/03_task1_task2_analysis.py` — Task 1 computations (overlap, part length, active fractions, barcode structure)
- `scripts/04_task2_rs241_analysis.py` — Task 2 (RS241) computations
- `scripts/funder_sweep.sh` — runnable funder-sweep API queries (written by the Task 5 agent)
- `raw/` — all primary-source downloads: 5 supplementary XLSX files, manuscript PDF/NXML/DOCX, PMC tarball, funder-sweep raw JSON/PDF responses, DARPA BAA PDF
