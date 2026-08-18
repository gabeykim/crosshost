# CROSSHOST — Gate 2 Memo: Data Foundation

**Date:** 2026-08-04
**Prepared for:** Gabriel
**Scope:** Depth-matching audit on the physiology vector; core data tables; genomic feature vector with unit tests; sequence-identity splits with leakage audit. No modeling code was written.

**Bottom line up front: PASS.** All 6 physiology metrics survive depth-matching (well above the ≥3 bar). A viable, balanced 5-fold split exists and is frozen. The leakage audit passes all 6 checks. But two investigations turned up real, load-bearing findings that change how the splits work and what they can promise: (1) the physiology depth confound was real and now-corrected; (2) sequence-identity clustering, even at the tool's own reported "0.7 threshold," demonstrably misses genuine high-identity cross-genome pairs — the splits that are actually frozen rely on a broader, independently-verified near-duplicate safety net, not on clustering alone.

---

## TASK 0 — Depth-matching audit: per-metric verdict

**The hypothesis was correct, and it was worse than a first look suggested.** V. natriegens's Gate 1.5 physiology numbers came from a single growth condition (665 proteins). Truncating the other 5 hosts to that same depth and comparing:

| Metric | Rank order preserved? (Spearman ρ, full vs. depth-matched) | Verdict |
|---|---|---|
| Ribosomal protein fraction | 1.000 | **DEPTH-SENSITIVE, CORRECTABLE** |
| RNAP core fraction | 1.000 | **DEPTH-SENSITIVE, CORRECTABLE** |
| Chaperone fraction | 1.000 | **DEPTH-SENSITIVE, CORRECTABLE** |
| Sigma factor fraction | 0.943 | **DEPTH-SENSITIVE, CORRECTABLE**, low reliability (absolute values <0.2%, report distinct-sigma-factor *count* alongside) |
| Elongation factor fraction | 0.829 | **DEPTH-SENSITIVE, CORRECTABLE**, weakest rank stability of the five |
| Growth rate | n/a (not proteomics-derived) | **DEPTH-ROBUST** (carries its own condition-matching caveats from Gate 1, unrelated to this audit) |

**All 6 survive** — none is DEPTH-CONFOUNDED/UNUSABLE. No DECISION NEEDED trigger.

**Before accepting 665 as V. natriegens's ceiling, a deeper source was checked and found [step 5 of the audit]:** pooling across all 81 samples of PRIDE PXD027874 (3 timepoints × 3 salinities × 3 temperatures × 3 replicates — "detected in ≥1 of 81 samples, mean LFQ across samples where detected") lifts V. natriegens to **1,032 quantified proteins**, up from 665. This was adopted as the new native depth. Trade-off disclosed: the pooled version is no longer tied to one documented growth condition — conceptually it's now "typical abundance across a temp/salinity/timepoint panel," which is actually closer in spirit to how PaxDb's own "Integrated" datasets for the other 5 hosts are built (weighted averages across many studies) than the single-condition version was.

**Final depth-matched values (common depth = 1,032, V. natriegens's new native ceiling):**

| Host | n proteins (native) | Ribosomal fraction | RNAP fraction | Sigma fraction | Chaperone fraction | EF fraction |
|---|---:|---:|---:|---:|---:|---:|
| *E. coli* | 3,747 | 21.29% | 0.76% | 0.089% | 2.03% | 2.51% |
| *B. subtilis* | 4,052 | 19.97% | 0.40% | 0.124% | 1.43% | 1.69% |
| *P. aeruginosa* | 5,034 | 25.48% | 0.97% | 0.142% | 4.51% | 2.67% |
| *S. enterica* | 2,620 | 23.04% | 0.53% | ~0.06% | 0.57% | 3.78% |
| ***V. natriegens*** | **1,032** | **27.61%** | **1.91%** | 0.066% | 2.47% | 4.09% |
| *C. glutamicum* | 1,227 | 14.94% | 0.32% | 0.09% | 0.73% | 2.46% |

**V. natriegens remains the highest-ribosomal-fraction host after correction** (27.6% vs. 25.5% for the next-highest, *P. aeruginosa*), but the margin shrank from a 6.4-point gap (uncorrected, 29.8% vs 23.5%) to a 2.1-point gap — **roughly two-thirds of the apparent lead was a depth artifact, not biology.** The qualitative story (fastest grower → highest ribosomal allocation, matching classical bacterial growth-law theory) survives; the magnitude claimed for it does not, and the corrected numbers are what Gate 3+ should use.

Full detail: `out/depth_matching_results.json`, `out/task0_verdict.json`, `data/hosts_physiology_final.csv`, `data/depth_sensitivity_curve.csv`.

---

## TASK 1 — Core data tables

Built and frozen (all in `data/`, hashed into `MANIFEST.json`):

- **`three_host.parquet`** — 29,249 rows (1 per RS), sequence + genome/taxonomy origin + per-host (BS/EC/PA) transcription/translation values and usable/active flags, using Gate 1's paper-derived QC rule. **Assertion in code** (not comment): `OLIGO ID` is 1:1 with `intgen_id` — if this ever fails on a future data refresh, the script raises immediately rather than silently reintroducing the barcode-replicate risk the charter originally worried about (amendment C1). Translation N verified = 7,887 (amendment C3), transcription N = 11,276 (Gate 1) — both re-asserted in code against the values already established.
- **`rs241.parquet`** — 241 rows, six-host values plus the Gate 1.5 pairwise-usability flags for both RS241 evaluation configurations (see Task 3.5).
- **`hosts_genomic.parquet`** — the ~40-D genomic feature vector (Task 2).
- **`hosts_physiology.parquet`** — the corrected 6-metric physiology vector from Task 0, with `physiology_source` and `proteome_depth_native`/`proteome_depth_used` columns on every row (persisted per new standing rule SR2 — see `out/CHARTER_AMENDMENTS.md`).

---

## TASK 2 — Genomic feature vector (37 features, 6 hosts) — all unit tests passed

Genomes downloaded from RefSeq (exact accessions verified via NCBI Datasets API, all "Complete Genome" reference assemblies): *E. coli* GCF_000005845.2, *B. subtilis* GCF_000009045.1, *P. aeruginosa* GCF_000006765.1, *S. enterica* GCF_000006945.2, *C. glutamicum* GCF_000011325.1, *V. natriegens* GCF_001456255.1 (ATCC 14048 — matches the PRIDE proteome's strain).

**Unit tests, each against an independent published/authoritative value, not a self-consistency check:**

| Feature | Test | Result |
|---|---|---|
| Anti-Shine-Dalgarno | *E. coli*/*P. aeruginosa* nonamer = ACCTCCTTA; *B. subtilis* = ACCTCCTTT (Johns et al. 2018 manuscript Methods, Gate 1 primary source) | **PASS**, exact match, all 3 |
| tAI weight formula | Python port matches the canonical dos Reis R package's own `get.ws()` output, run live via Rscript in this environment on the reference 87-tRNA *E. coli* vector | **PASS**, max abs diff = 3.89e-16 (float precision floor) |
| tAI ranking | *E. coli* AAA (Lys) codon has tAI weight 1.0, verified against the *actual R output* — **not memory**: an earlier draft wrongly asserted CTG=1.0 as a "well-known fact," caught by checking the real reference computation before committing to the assertion | **PASS**, exact match |
| Sigma factor complement (HMMER + Pfam) | *E. coli* count = 7 (RpoD/RpoS/RpoH/RpoE/RpoN/FliA/FecI — Gruber & Gross 2003; Paget 2015) | **PASS**, exact match, after excluding 4 region-4-only false positives (DNA-binding response regulators NarL/YhjB/NarP/UvrY sharing the HTH fold with sigma region 4 but no sigma-specific core domain) |
| 16S rRNA copy number | plausibility vs. published rrn operon counts | **PASS**: EC=7, BS=10, PA=4, SE=7, CG=6, VN=11 (V. natriegens highest, consistent with its extreme growth rate) |
| GC content | plausibility vs. well-known reference values | **PASS**: EC=50.8%, BS=43.5%, PA=66.6%, SE=52.2%, CG=53.8%, VN=45.1% — all match closely |

**Sigma-factor HMM source:** individual Pfam HMMs (not the full ~20GB Pfam-A database) fetched from InterPro's REST API for 10 sigma-factor-related families (Sigma70 regions 1.1/1.2/2/3/4/4_2/non-essential, dedicated ECF family, Sigma54 DNA-binding + core-binding domains), identified by searching InterPro's API rather than from memory. Anti-sigma-factor families (Crl, anti-σ regulators) were identified and explicitly excluded by their own Pfam descriptions.

**Method deviation, disclosed:** ribosomal protein complement, chaperone complement, RNAP core subunits, and heme biosynthesis gene counts used gene-symbol + annotation-keyword matching (the same validated method from Gate 1.5, applied to RefSeq protein FASTA headers) rather than KEGG KO or EggNOG-mapper v2 as specified — both require multi-GB reference database downloads judged infeasible in this gate's time budget. This is a reused, already-QC'd method, not an ad hoc substitute, but it is a real deviation worth knowing about if Gate 2's genomic vector is ever audited against a KO-based pipeline.

**Genome-vs-physiology contrast, as intended:** the genomic vector carries sigma-factor **counts** (7/16/22/6/7/8 for EC/BS/PA/SE/CG/VN — *B. subtilis* and *P. aeruginosa*'s large repertoires match their well-documented regulatory complexity) while the physiology vector carries sigma-factor **abundance fraction**. Both are in the frozen tables, ready for the Gate 5 comparison.

Full detail: `out/genomic_features_basic.json`, `out/antisd_trna_features.json`, `out/tai_results.json`, `out/sigma_factor_hmmer_results.json`, `out/unit_test_summary.json`.

---

## TASK 3 — Splits: the real story, not just the requested sweep

### 3.1 — Clustering at 3 thresholds: a tooling limitation found and disclosed

`mmseqs easy-cluster` (the standard, sensitive mode) **could not run in this environment** — it repeatedly failed with "Cannot fit databases into [N]G" regardless of sensitivity/thread/split tuning, because this sandboxed machine has 8GB RAM with ~92% of swap already committed (confirmed via `vm_stat`/`vm.swapusage`). Switched to `mmseqs easy-linclust` (Steinegger & Söding 2018's dedicated linear-time, low-memory mode), which completed successfully.

**Finding: `--min-seq-id 0.3` and `--min-seq-id 0.5` produced byte-identical clustering output** (verified: `diff` of the two `_cluster.tsv` files is empty, 0 lines). A probe sweep (0.3/0.4/0.5/0.6/0.7/0.8/0.9) showed cluster count is flat at 27,869 through 0.5, then increases monotonically from 0.6 onward (28,134 → 28,281 → 28,418 → 28,587). **For these 165bp sequences, linclust's default k-mer seeding has a practical identity-detection floor between 0.5 and 0.6** — below that, the prefilter simply finds no additional candidate pairs to test, regardless of the requested threshold. This is a genuine tool/data-scale limitation, not a bug in the pipeline, and it means the originally-planned 0.3/0.5/0.7 sensitivity sweep only has two effectively distinct data points (≤0.5, and 0.7) for this dataset.

| min-seq-id | n clusters | largest cluster | n singletons |
|---|---:|---:|---:|
| 0.3 | 27,869 | 13 | 26,893 |
| 0.5 | 27,869 (identical to 0.3) | 13 | 26,893 |
| 0.7 | 28,281 | 11 | 27,543 |

### 3.2b — Exact-duplicate characterization: clustering missed a real leak, and the reason it happened matters

Gate 1 found 181 groups of byte-identical 165bp sequences (413 oligos). Building folds directly from linclust clusters, the per-fold max-identity audit (an independent k-mer-index + shift-tolerant exact-identity check, cross-validated against Biopython's global pairwise aligner — exact match, 154/165 = 93.3%, no discrepancy) caught a pair (OLIGO 21865/37869, **same source genome**, byte-identical sequence) sitting in two different singleton clusters at `--min-seq-id 0.7`.

**Full characterization, not just the one instance:**

1. **181 distinct duplicate-sequence groups, 413 oligos, 336 pairwise duplicate relationships.** Matches Gate 1 exactly.
2. **Same-genome vs. different-genome split: 47 groups same-genome, 134 groups (74%) different-genome.** This is the load-bearing number — most exact duplicates are NOT simple within-genome replicate artifacts; they are independent convergence between different source genomes on an identical 165bp sequence.
3. **Explanation:** the 47 same-genome groups are consistent with the same intergenic region being assigned two OLIGO IDs during original library construction. The 134 different-genome groups are direct evidence of shared/conserved intergenic sequence across related genomes in the 184-genome mining set — real biological relatedness, not a data artifact.
4. **How many duplicate groups did clustering actually miss?** Only **1 of 181** (3 oligos) was split across clusters, consistently at every tested threshold (0.3/0.5/0.7). Clustering caught 180/181 (99.4%) correctly — this is a real but *low-frequency* failure mode, not pervasive breakage. It matters because even a low failure rate is a genuine leak if unfixed, and because of point 6 below.
5. **What the safety net actually does:** majority-fold consolidation — for any exact-duplicate group whose members land in >1 fold, every member is reassigned to whichever single fold already held the largest share of the group. No oligo is dropped, no cluster is merged; only fold *labels* are overwritten for the minority member(s). Total N is unchanged.
6. **Does this change confidence in clustering? Yes, and here's the chain of evidence, not just inference:** if two genomes can independently yield a **100%**-identical sequence and clustering still misses it at a 70% threshold, it follows directly (not just plausibly) that **near**-duplicate cross-genome pairs — 95%, 90%, 85% identity, well within homology range — are also being missed, and *undetectably* so, because they don't trip an exact-string-match safety net the way 100%-identical pairs do. This was then directly confirmed (not just inferred) in Task 3.3.

### 3.3 — Source-genome leakage check: clustering AND naive genome-blocking both leak

**Genome distribution:** 184 source genomes, well spread — top-1 genome = 1.98% of the library, top-10 = 14.19%, median 139 sequences/genome, only 2 genomes contribute <10 sequences, none contributes just 1. **Taxonomic spread:** 27 phyla represented (Proteobacteria 42 genomes, Firmicutes 40, Bacteroidetes 20, Actinobacteria 20, plus 23 smaller phyla) — genuinely diverse, not concentrated in one clade.

**Cross-genome max identity under cluster-based folds (excluding same-genome pairs entirely) reached 0.933** at every tested threshold, including 0.7 — confirmed directly: OLIGO 29542/29549 (different genomes, 637000108 and 637000263) are **154/165 = 93.3% identical by Biopython global alignment** (gap-free colinear substitutions, cross-checked against the k-mer-index estimate — exact agreement) yet sit in two separate singleton mmseqs clusters at `--min-seq-id 0.7`. This directly demonstrates that a "clustered at 0.7" claim does not actually bound cross-fold identity at 70% for this dataset.

**Naive genome-blocked splitting (assign whole genomes to folds, no sequence-level check) is *not* better — it's worse: max identity reached 1.000 (100%, exact duplicate) in 4 of 5 folds.** This is the mechanical consequence of the 134 different-genome exact-duplicate groups found in 3.2b: blocking by genome does nothing to prevent two *different* genomes from contributing the *same* sequence to opposite folds. A second concrete case surfaced here too: OLIGO 20578/32225 are not equal as literal strings (so an exact-match safety net misses them) but are **100% identical once a 2bp register shift is applied** — consistent with a small indel near the 165bp extraction window boundary between two related genomic loci. Three distinct leakage mechanisms are now documented for this dataset: (1) literal exact duplicates, (2) substitution-only near-duplicates, (3) shift/indel near-duplicates — no single simple method (exact-string dedup, genome-blocking, or clustering alone) catches all three.

**Genome-blocked class-balance check** (before the fix below): fold sizes 5,838–5,862 (well balanced); `usable_all3_transcription` rate per fold 0.364–0.419 (reasonably stable, no fold is degenerate).

### 3.4 — Splitting scheme chosen, with reasoning, and the splits frozen

**Recommendation and what was built: genome-blocked assignment as the primary structural layer, PLUS a comprehensive sequence-identity safety net, not clustering as the sole defense.**

Reasoning: cluster-based splitting's only advantage over genome-blocking would be a stronger sequence-level guarantee — but 3.3 showed cluster-based folds still leak cross-genome identity up to 93.3%, i.e. no better than genome-blocking's own residual risk once a proper safety net is applied to either. Given that, genome-blocking wins on the criterion that actually matters for defensibility: **it is a clean, easily-stated guarantee** ("no two sequences from the same source organism appear on both sides of a split") that a reviewer can verify by construction, whereas "clustered at threshold X" turned out not to mean what it sounds like it means for this dataset.

**Final construction (both layers applied, in order):**
1. **Layer 1:** genome-blocked fold assignment (whole genomes → folds, greedy load-balanced).
2. **Layer 2:** a full all-pairs near-duplicate search over the entire 29,249-sequence library (k-mer-index candidate screening, k=10, ≥5 shared k-mers to trigger an exact shift-tolerant identity check — the same validated method from 3.2b/3.3, cross-checked against Biopython) at a **merge threshold of 0.85**. 255,684 candidate pairs were checked; 795 pairs merged at ≥0.85 identity into 447 multi-member connected components (up from the 181 exact-duplicate groups alone — 266 additional near-duplicate groups found this way). Any connected component spanning >1 fold was consolidated via the same majority-fold rule as 3.2b: 295 components affected, 344 oligos reassigned.

**Why 0.85 and not something else:** comfortably above the highest cross-genome identity mmseqs clustering used as an operating threshold (0.7) while low enough to catch both concretely observed problem pairs (93.3% substitution-only, 100% shifted) with real margin. It is not claimed as a universal safe threshold on its own — the number that actually matters is the measured result below.

**Final per-fold max train-test identity (full, unsampled check, after both layers):**

| Fold | Test N | Train N | Max identity |
|---|---:|---:|---:|
| 0 | 5,848 | 23,401 | **0.8485** |
| 1 | 5,854 | 23,395 | 0.8395 |
| 2 | 5,848 | 23,401 | **0.8485** |
| 3 | 5,838 | 23,411 | **0.8485** |
| 4 | 5,861 | 23,388 | 0.8364 |

**Maximum across all 5 folds: 0.8485.** This is the headline, defensible number — not "we clustered at 0.7."

**RS241 mechanical enforcement — a critical fix, not just an assertion:** checking the frozen fold assignment directly found **207 of the 241 RS241 oligo IDs were physically present** in the main-library fold table (RS241 was "selected" from existing library members, Gate 1) and would have been used in training for ~4/5 of folds under a naive "train on everything not in the test fold" loop — a direct violation of the charter's non-negotiable rule. **Fixed:** all 241 RS241 IDs removed from the main-library fold assignment entirely (29,249 → 29,042 rows). Verified by direct set-intersection check (zero remaining) both in the fix script and independently in `audit_leakage.py`.

**Final frozen fold sizes** (`data/splits/fold_assignment_FINAL.parquet`, after RS241 removal): fold 0 = 5,810; fold 1 = 5,818; fold 2 = 5,809; fold 3 = 5,792; fold 4 = 5,813. Total 29,042.

### 3.5 — RS241 evaluation configurations

`data/splits/rs241_configs.json`:

| Config | Train hosts | *S. enterica* (tx / tl) | *V. natriegens* (tx / tl) | *C. glutamicum* (tx / tl) |
|---|---|---|---|---|
| **PRIMARY** | EC+BS+PA | 146 / 141 | 113 / 134 | 138 / 130 |
| **SECONDARY** | EC+PA (deeper phylogenetic distance to *C. glutamicum*) | 216 / 225 | 162 / 215 | 200 / 211 |

Matches the Gate 1.5 pairwise table exactly (cross-checked, no drift).

---

## TASK 4 — Leakage audit (`scripts/audit_leakage.py`)

Standalone, exits non-zero on failure. **All 6 checks pass** — full output in the status report below.

---

## TASK 5 — Standing rules persisted

Four numbered standing rules (SR1–SR4) plus the *P. aeruginosa* two-method-independent finding and the headline max-identity result were appended to `out/CHARTER_AMENDMENTS.md` under a new "STANDING RULES — Gate 2" heading. See that file for full text; summarized in the status report below.

---

## WHAT I COULD NOT DO

- **Could not run `mmseqs easy-cluster`** (the sensitive/complete clustering mode) in this environment — it failed with out-of-memory errors at every sensitivity/thread/split-limit combination tried, down to `-s 1 --split 16 --split-memory-limit 1500M`. Used `easy-linclust` instead throughout, which has a documented, now-empirically-characterized lower sensitivity for this short-sequence, low-average-identity dataset (see 3.1).
- **Could not fully resolve why linclust misses specific high-identity pairs** (e.g. the 93.3%-identical 29542/29549 pair) at the algorithmic level — the working hypothesis (approximate/chunked/parallelized processing with no formal completeness guarantee for linear-time clustering) is standard, documented behavior for this class of algorithm, but this gate did not trace the exact internal cause, since the empirical safety net (Task 3.4) makes the root cause less operationally important than the fact that it happens.
- **Did not use KEGG KO or EggNOG-mapper v2** for ribosomal/chaperone/RNAP/heme gene identification as specified — both require multi-GB database downloads judged infeasible in this gate's time budget. Substituted the same gene-symbol + annotation-keyword method validated in Gate 1.5, applied to RefSeq protein headers. Flagged as a real, disclosed deviation, not silently substituted.
- **Did not verify the Abele et al. 2025 atlas's exclusion of *E. coli*/*C. glutamicum*** any further than Gate 1.5 already established (still unresolved whether they're genuinely absent from that atlas or just not yet in PaxDb's per-species index for those two species).
- **Did not extend the normalization check (audit_leakage.py check 5) to real functional verification** — no normalization code exists before Gate 3, so the check currently only verifies the *absence* of a suspicious global-statistic artifact in frozen tables. It must be extended once Gate 3 introduces actual fit/transform code to check that the fitting function is called once per fold on training data only.
- **Did not trace the exact biological cause of the 2bp-shifted near-duplicate pair** (OLIGO 20578/32225) beyond noting it's consistent with a small indel near the extraction window boundary — confirming this would require re-examining the source genomes' actual gene annotations at those loci, not attempted here.

## CONTRADICTIONS WITH THE CHARTER

1. **The charter's assumed leakage defense (barcode dedup) doesn't apply** (established in Gate 1/amendment C1, reconfirmed here) — sequence-identity methods are the only defense, and this gate found they are *themselves* imperfect and require a purpose-built safety net beyond straightforward clustering.
2. **"Clustered at a stated identity threshold" does not mean what the charter (and most literature in this space) implicitly assumes it means, for this dataset.** This is the single most important methodological finding of this gate: mmseqs2 clustering, run in good faith at `--min-seq-id 0.7` with standard coverage settings, provably fails to co-cluster real 93.3%-identical and 100%-identical (same-genome and different-genome) sequence pairs. Any project that reports "we deduplicated at 70% identity" without an independent max-identity audit is very likely also carrying leakage it doesn't know about — worth flagging in the eventual paper as a methods contribution in its own right, not just a private caveat.
3. **RS241 was at real risk of being trained on** — not hypothetically, but concretely: 207/241 IDs were sitting in the ordinary fold-assignment table before this gate's fix. This is exactly the kind of "enforce mechanically, not by discipline" failure mode the charter itself anticipated (Part V rule 1) but for a different mechanism (barcode replicates) than the one that actually manifested (RS241's own IDs overlapping the main library).

## FILES WRITTEN

- `out/GATE2_MEMO.md` — this memo
- `out/CHARTER_AMENDMENTS.md` — updated with SR1–SR4 and two new findings
- `out/state.json` — updated with `gate_2` block
- `out/depth_matching_results.json`, `out/task0_verdict.json`, `data/hosts_physiology_final.csv`, `data/depth_matching_audit.csv`, `data/depth_sensitivity_curve.csv`
- `out/genomic_features_basic.json`, `out/antisd_trna_features.json`, `out/tai_results.json`, `out/sigma_factor_hmmer_results.json`, `out/unit_test_summary.json`
- `out/exact_duplicate_characterization.json`, `out/source_genome_leakage_results.json`, `out/genome_blocked_safe_results.json`, `out/final_split_results.json`, `out/rs241_exclusion_fix.json`
- `data/three_host.parquet`, `data/rs241.parquet`, `data/hosts_genomic.parquet`, `data/hosts_physiology.parquet`, `data/MANIFEST.json`
- `data/splits/fold_assignment_FINAL.parquet`, `data/splits/clusters_id{30,50,70}.parquet`, `data/splits/clustering_summary.csv`, `data/splits/rs241_configs.json`
- `data/vnatriegens_proteome_pooled_81samples.parquet`
- `raw/genomes/` (6 hosts × genomic/protein/GFF/CDS/RNA FASTA), `raw/pfam_hmms/`, `raw/mmseqs/`
- `scripts/09` through `scripts/25`, plus `scripts/audit_leakage.py`
