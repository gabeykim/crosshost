# CROSSHOST — Splitting Scheme Specification

**Purpose:** precise enough that an independent implementer reproduces the frozen fold assignment in `data/splits/fold_assignment_FINAL.parquet` exactly. This is the paper's methods section on splitting.

**Update (Gate 2.5, see `out/GATE2_5_MEMO.md` Task A1):** the original justification for this scheme ("clustering is generally unreliable for this dataset") has been narrowed. Sensitive clustering (`mmseqs easy-cluster`) and an independent tool (`cd-hit-est`) both correctly co-cluster the specific 93.3%-identity pair that motivated this scheme, at 100% recall on a tractable subset. The failure is specific to `easy-linclust`, the fast/approximate mode used because `easy-cluster` could not run at full library scale (29,242 sequences) in this environment's 8GB RAM. **The scheme below is unchanged** — it was already designed to not depend on trusting any single clustering method's output, and it remains the right choice on practical grounds (it must work at full scale, and it is independently verified regardless of which clustering method, if any, underlies it) — but the reasoning is now: "we could not run a fully-trustworthy clusterer at full scale in this environment, so the split does not rely on clustering at all; it relies on an independently-verified all-pairs identity check," not "clustering doesn't work for this kind of data."

---

## Inputs

- `data/three_host.parquet` — 29,249 rows, one per RS (`OLIGO ID`, `regulatory_sequence`, `genome_id`, plus per-host transcription/translation usability flags). Source: Gate 1/Gate 2 processing of Johns et al. 2018 Supplementary Data Tables 1–2.
- `data/rs241.parquet` — 241 rows, the RS241 six-host evaluation subset. **Reserved entirely for held-out-host evaluation; must never appear in the output of this pipeline.**

## Step-by-step, in execution order

### Step 1 — Genome-blocked base assignment

1. Count sequences per `genome_id` across all 29,249 rows of `three_host.parquet`. There are 184 distinct source genomes.
2. Sort genomes descending by sequence count.
3. Greedy load-balance: initialize 5 fold counters to 0. For each genome (largest first), assign its **entire** sequence set to whichever fold currently has the smallest total sequence count, then add that genome's count to the chosen fold's running total.
4. Result: every `OLIGO ID` is assigned to exactly one of folds {0,1,2,3,4} via its `genome_id`. No two sequences from the same source genome are ever in different folds at this step.

**Determinism:** fully deterministic given the sorted genome-count order — no random number generator is used in this step. (`np.argmin` ties, if any, resolve to the lowest fold index; with 184 genomes of varying size this has not produced an observed tie in practice.)

### Step 2 — All-pairs near-duplicate detection (independent of Step 1 and of any clustering tool)

This step does **not** use mmseqs2, cd-hit, or any external clustering tool. It is a self-contained k-mer-index + exact-identity computation, run once over the full 29,249-sequence set:

1. **Candidate screening:** build an inverted index mapping every 10-mer (`k=10`) to the list of sequences containing it. For each sequence, look up all other sequences sharing **≥5** 10-mers with it (`MIN_SHARED_KMERS = 5`) — these are the only pairs subjected to the next step. (This candidate filter is what makes the method tractable: random unrelated 165bp sequences essentially never share 5+ ten-mers by chance, given 4¹⁰ ≈ 1M possible 10-mers against ~156 target k-mers per sequence.)
2. **Exact identity computation, for candidates only:** for each candidate pair `(a, b)`, compute identity as the best ungapped match fraction over a small window of relative offsets, `offset ∈ [-5, +5]`: for each offset, align `a[i]` against `b[i+offset]` for all valid `i`, count matches, divide by the number of overlapping positions; take the max over all 11 offsets. This tolerates small indel-sized register shifts (found necessary — see `out/GATE2_MEMO.md` Task 3.3, the OLIGO 20578/32225 case, a pair 100% identical only after a 2bp shift) without requiring a full gapped aligner. **Validated** against Biopython's `Align.PairwiseAligner` (global mode, match=1/mismatch=0/gap_open=-2/gap_extend=-0.5) on a specific pair (OLIGO 29542/29549): both methods agree exactly at 93.3% (154/165).
3. **Merge threshold:** any candidate pair with identity **≥ 0.85** is unioned into the same group via union-find (disjoint-set). `IDENTITY_MERGE_THRESHOLD = 0.85`.
4. Result on the full library (this run): 255,684 candidate pairs checked, 795 pairs merged, forming 447 connected components with ≥2 members (up from the 181 exact-duplicate-only groups found in Gate 1 — 266 additional near-duplicate groups found by this broader check).

**Why 0.85:** chosen to sit comfortably above the highest cross-genome identity that clustering had been assumed to bound (0.7) while catching both concretely observed problem cases (93.3% substitution-only, 100% shift-only) with margin. It is not claimed to be a theoretically optimal cutoff — the number that matters is the *measured* final max identity (Step 5), which is reported regardless of what this threshold was set to.

**Determinism:** fully deterministic — no randomness in index construction, candidate screening, or identity computation. Iteration order over Python dict/set structures is insertion-order-stable in this implementation (CPython ≥3.7), so results are reproducible run-to-run on the same input.

### Step 3 — Consolidate near-duplicate groups that span multiple folds

1. For every connected component found in Step 2 with ≥2 members, check how many distinct folds (from Step 1) its members currently occupy.
2. If a component spans >1 fold: **majority-fold consolidation** — count how many of the component's members are in each fold; move every member of the component to whichever single fold already holds the largest share. Ties broken by `pandas.Series.value_counts().idxmax()`, which returns the first-encountered value among ties (stable but not documented as a formal tie-break rule — no observed tie has driven a consequential outcome in this dataset).
3. This run: 295 components (of 447 multi-member components) spanned >1 fold before consolidation; 344 oligos were reassigned. **No oligo is dropped and no two components are merged into one** — only fold *labels* are overwritten for the minority member(s) of each affected component. Total sequence count is unchanged by this step.

### Step 4 — RS241 exclusion (mechanical, not by discipline)

1. Take the set of all 241 `id` values in `data/rs241.parquet`.
2. Remove **every** row of the Step 3 output whose `OLIGO ID` is in that set — regardless of whether removal is strictly necessary for a given ID (defense in depth: all 241 are removed, not just the 207 that happened to physically overlap the main library at the time of the Gate 2 fix).
3. **Assertion, enforced in code, not just checked after the fact:** `assert len(set(output["OLIGO ID"]) & rs241_ids) == 0` — the pipeline raises immediately if this fails.
4. Result: 29,249 → 29,042 rows (207 removed; the other 34 RS241 IDs were never in the main library's ID space to begin with, per Gate 1).

### Step 5 — Verification (not part of fold construction, but mandatory before freezing)

For each fold `f` in {0,1,2,3,4}: let `test = {oligos with fold==f}`, `train = {all other oligos}`. Using the same k-mer-index + exact-identity method as Step 2 (an independent re-check, not a re-use of Step 2's cached candidate list), compute `max(identity(t, tr) for t in test for tr in train)` restricted to k-mer-index candidates. This is the number reported as the "final max train-test identity" and is what `scripts/audit_leakage.py` re-verifies on every subsequent gate.

**Result, this run:** fold 0 = 0.8485, fold 1 = 0.8395, fold 2 = 0.8485, fold 3 = 0.8485, fold 4 = 0.8364. **Maximum across all folds: 0.8485.**

## Output

`data/splits/fold_assignment_FINAL.parquet` — columns `OLIGO ID` (str), `fold` (int, 0–4). 29,042 rows. Every downstream gate should treat "train for held-out fold `f`" as "all rows with `fold != f`" and "test for held-out fold `f`" as "all rows with `fold == f`" — RS241 rows are absent from this table entirely and must be sourced separately from `data/rs241.parquet` for the RS241 evaluation track (see `data/splits/rs241_configs.json`).

---

## Why this scheme over the alternatives — with the measured numbers

| Scheme | Max cross-fold identity measured | Verdict |
|---|---:|---|
| Cluster-blocked only (mmseqs easy-linclust, min-seq-id 0.7, whole clusters → folds) | 0.933 | Rejected: clustering (in the mode that could run at full scale) provably misses real high-identity pairs |
| Genome-blocked only, no sequence-level check | **1.000** | Rejected: *worse* than cluster-blocked — does nothing to prevent different genomes independently yielding identical sequences (134 of 181 exact-duplicate groups are cross-genome) |
| **Genome-blocked + exact-duplicate safety net + near-duplicate safety net (this scheme)** | **0.8485** | **Adopted** |

Genome-blocking is retained as the Step-1 structural layer (rather than dropped in favor of pure identity-based clustering via Step 2 alone) because it gives a clean, reviewer-verifiable-by-construction guarantee ("no two sequences from the same source organism are ever on both sides of a split") in addition to the measured sequence-identity bound — belt and suspenders, not redundant.

## Known limitations

1. **The 0.8485 residual is a real, disclosed, non-zero number, not a guarantee of zero leakage.** Pairs below the 0.85 merge threshold but above whatever a reader considers "too similar" (e.g. 0.80–0.8485) do exist across folds by construction. This should be stated plainly in the paper, not rounded away.
2. **The k-mer-index candidate screen (Steps 2 and 5) could in principle miss a true high-identity pair that shares fewer than 5 ten-mers** — this is only possible for pairs with identity concentrated in a way that avoids exact 10-mer matches (e.g. many short, evenly-spaced single-base differences). No such case has been found in this dataset, but it is not mathematically excluded. A shorter k (more sensitive, slower) or a formal minimum-identity-implies-minimum-shared-kmers bound would close this gap rigorously; not done here for time.
3. **The merge threshold (0.85) and candidate-screen parameters (k=10, min_shared_kmers=5) were chosen by the same people who ran the audit** — they were not tuned against an external ground truth beyond the two concretely observed problem pairs. Gate 2.5's tool comparison (Task A1) used the same parameters for its own ground-truth computation, which is internally consistent but means the "ground truth" in that comparison is itself defined by this method, not by a fully independent oracle.
4. **Class balance across folds was checked only for `usable_all3_transcription` rate** (range 0.364–0.419 under the genome-blocked-only precursor scheme; not re-verified after Steps 2–4's reassignments in this document, though the reassignment volume — 344 + 207 oligos out of 29,249, ~1.9% — is small enough that material balance drift is unlikely). Should be re-checked explicitly in Gate 3 if fold-level performance looks anomalous for any host.
