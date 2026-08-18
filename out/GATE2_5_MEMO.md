# CROSSHOST — Gate 2.5 Memo: Three Open Items Resolved

**Date:** 2026-08-04
**Scope:** A1 (clustering finding tool-specificity), A2 (splitting scheme documentation — see `out/SPLITTING_SCHEME.md`), A3 (annotation-quality check on genomic features).

---

## A1 — The clustering-miss finding is tool-specific to linclust, not general

**Method:** built a tractable 3,415-sequence subset (413 exact-duplicate oligos + the known 93.3%-identity pair OLIGO 29542/29549 + 3,000 random background sequences), computed ground truth via the same k-mer-index + shift-tolerant exact-identity method already validated against Biopython (407 true pairs at ≥0.7 identity), then ran three tools against that ground truth.

| Tool | Ran successfully? | Precision | Recall | Co-clusters known 93.3% pair? |
|---|---|---:|---:|---|
| `mmseqs easy-cluster` (sensitive/full mode) | **Yes**, on this subset (still fails on the full 29,242-sequence library — see below) | 0.967 | **1.000** | **YES** |
| `cd-hit-est` (independent, greedy-exhaustive tool) | Yes, at c=0.8 — see note | 0.961 | **1.000** | **YES** |
| `mmseqs easy-linclust` (the mode actually used for the full library) | Yes | 0.981 | 0.993 | **NO** |

**Note on CD-HIT:** this build (4.8.1) hard-rejects any `-c` below 0.8 for nucleotide clustering ("invalid clstr threshold, should >=0.8"), at every word size tried (3–8), confirmed empirically rather than assumed — CD-HIT's nucleotide mode is not designed to run below 80% identity for short sequences at all (its own documentation points to the separate `psi-cd-hit` workflow for that regime, not attempted here). Since the known problem pair is at 93.3% — comfortably above 0.8 — running CD-HIT at 0.8 still directly tests whether it catches this specific pair. It does, with zero false negatives (fn=0) against its own 0.8-threshold ground truth.

**Verdict: the failure is specific to `easy-linclust`'s linear-time approximation, not a general property of sequence clustering.** Both a more sensitive mode of the *same* tool (mmseqs `easy-cluster`) and a completely independent implementation (CD-HIT) correctly identify and co-cluster the pair that motivated the whole investigation, at 100% recall on this subset.

**What this changes and what it doesn't:**
- **Changes:** the paper's claim narrows from "clustering-based leakage defenses are broadly unreliable for this class of data" to "the fast/approximate clustering mode we were forced to use in a memory-constrained environment has a real, measurable gap; sensitive clustering does not show it at the tested scale." This is a narrower, more defensible, and frankly more mundane claim — worth being honest about rather than overselling the finding.
- **Doesn't change:** `easy-cluster` *still cannot run on the full 29,242-sequence library* in this environment (confirmed repeatedly in Gate 2, OOM at every tuning attempt down to `-s 1 --split-memory-limit 1500M`). The practical problem — no trustworthy clustering tool runs at full scale here — persists regardless of this finding. **The frozen splitting scheme (genome-blocked + independent all-pairs near-duplicate safety net, `out/SPLITTING_SCHEME.md`) is unchanged and still the right choice**, because it was already designed to not depend on trusting any clustering tool's output — it only got a better, narrower justification.

Full detail: `out/clustering_tool_comparison.json`, `scripts/26_clustering_tool_comparison.py`.

---

## A2 — Splitting scheme documentation

Written as a standalone document: **`out/SPLITTING_SCHEME.md`**. Covers every step in execution order with exact parameters, the safety-net mechanics (detection method, action taken, tie-breaking), where RS241 exclusion happens, the three-way comparison table justifying the chosen scheme (0.933 cluster-blocked-only, 1.000 genome-blocked-only, **0.8485 the adopted scheme**), determinism, and four explicitly-stated known limitations. Written for a skeptical reviewer, per the instruction.

---

## A3 — Annotation-quality check: features are usable as-is, with a stated caveat

**Gene-symbol informativeness varies enormously across the six hosts** — from *E. coli* at 100% (4,300/4,300 proteins carry an informative gene symbol) down to *C. glutamicum* at 20.9% (614/2,941):

| Host | n proteins | % with gene symbol | % hypothetical/uncharacterized |
|---|---:|---:|---:|
| *E. coli* | 4,300 | 100.0% | 7.8% |
| *B. subtilis* | 4,237 | 98.6% | 17.2% |
| *S. enterica* | 4,554 | 72.6% | 0.8% |
| *V. natriegens* | 4,504 | 32.5% | 12.2% |
| *P. aeruginosa* | 5,572 | 31.5% | 40.5% |
| *C. glutamicum* | 2,941 | 20.9% | 21.2% |

This confirms the concern was worth checking — the spread is real and large (79 percentage points on gene-symbol ratio).

**But the specific feature counts used in `hosts_genomic.parquet` do not track this at all:**

| Host | gene-symbol ratio | ribosomal | RNAP | chaperone | heme |
|---|---:|---:|---:|---:|---:|
| EC | 1.000 | 57 | 4 | 7 | 16 |
| BS | 0.986 | 51 | 4 | 6 | 14 |
| SE | 0.726 | 61 | 4 | 4 | 14 |
| VN | 0.325 | 58 | 4 | 6 | 20 |
| PA | 0.315 | 56 | 4 | 6 | 12 |
| CG | 0.209 | 55 | 4 | 8 | 11 |

Pearson correlation (n=6, eyeball-level only, as instructed) between gene-symbol ratio and feature count: ribosomal r=**−0.18**, chaperone r=**−0.27**, heme r=**+0.16**, RNAP r=**undefined** (zero variance — every host has exactly 4 RNAP core subunits found, the theoretically correct count for all of them). **No positive correlation exists; if anything the sign is slightly negative** for two of the four features (worse-annotated genomes show marginally *higher*, not lower, counts) — the opposite of what an annotation-completeness artifact would predict.

**Why:** ribosomal proteins, RNAP subunits, chaperones (GroEL/DnaK family), and heme biosynthesis genes are among the most highly conserved, essential, long-studied gene families in bacterial genomics — annotation pipelines (including the automated PGAP pipeline behind all six RefSeq assemblies here) handle these specific families reliably even in genomes with low overall curation depth, because they're recognized via strong homology to extremely well-characterized reference sequences, not via case-by-case manual curation.

**Verdict: the four count-based genomic features (ribosomal protein complement, RNAP core subunits, chaperone complement, heme biosynthesis genes) are usable as-is.** No normalization or annotation-quality covariate is needed for these specific features. The caveat the paper must carry: this robustness is a property of *these particular* gene families (universal, essential, highly conserved) and should not be assumed to generalize to any *other* keyword-matched feature added later — Gate 2's standing rule SR1 (depth-matching check before use) has an analogous spirit here and a similar check (correlate any new keyword-derived feature against overall annotation completeness before trusting it) should be applied to future additions. Flagged for Gate 5's ablation to test directly rather than assumed settled.

Full detail: `out/annotation_quality_check.json`, `scripts/27_annotation_quality_check.py`.
