# crosshost — Project Charter and Execution Plan

**Version:** 1.0
**Date:** 2026-08-03
**Owner:** Gabriel
**Working directory:** `crosshost/`
**Status:** Approved to execute. Gate 1 not yet run.

---

## PART I — WHAT THIS IS

### The one-line version

A public benchmark measuring whether anyone — including large genomic foundation models — can predict how a regulatory DNA sequence will behave in a bacterial host that their model has never seen.

### The problem, stated as the field states it

A promoter characterized in *E. coli* frequently does not behave the same way in *B. subtilis* or *P. aeruginosa*. The field calls this the **chassis effect** or **context dependence**. Multiple independent groups have published, between 2023 and 2025, that:

- The chassis effect hinders prediction of circuit function from part composition alone, causes costly repetitions of the design-build-test cycle, and discourages use of non-model organisms.
- A predictive understanding of which biological properties underpin chassis effects is **a major knowledge gap left unanswered**.
- It remains a major challenge to identify which genome-encoded determinants govern chassis effects.

This is voiced demand from practitioners in primary literature. It is the reason this project exists and it is the thing that distinguishes it from the previously-killed crispr-eval project, where no such demand could be found.

### The published objection this project is built to test

The group most invested in the chassis effect (Bernstein lab and collaborators, *mSystems* / bioRxiv, Stutzerimonas work) has published a direct claim:

> It is impractical to postulate that the observable chassis-effect between a given set of hosts can be explained by a single or even a set of predictable genome-encoded functions without experimental insight.

And offered a competing hypothesis with supporting evidence:

> Hosts exhibiting more similar metrics of growth and molecular physiology also exhibit more similar performance of the genetic inverter, indicating that specific bacterial physiology underpins measurable chassis effects.

**This project's central scientific question is whether that claim is correct.** Can host *genomic* features predict cross-host regulatory activity, or does it require host *physiology* — the thing you have to measure in a lab?

That question is falsifiable, answerable from public data, and publishable whichever way it resolves. It does not require the project to succeed at prediction in order to succeed as a contribution.

### What is being built, precisely

**Primary asset — the benchmark.** A versioned, citable evaluation suite for cross-host regulatory activity prediction in bacteria. Fixed data splits, held-out hosts, withheld test labels, defined metrics, one-command reproduction, permanent DOI. No such benchmark exists — every major regulatory-activity benchmark (BEND, Genomic Benchmarks, DART-Eval, DNALONGBENCH, the Nucleotide Transformer suite) is human, animal, or plant.

**Secondary asset — the demonstration model.** A small host-conditioned CNN, built to (a) prove the benchmark is answerable, (b) set a baseline others must beat, (c) test the Bernstein claim by comparing genomic against physiological host representations, and (d) build the ML engineering depth that is a stated primary goal of this work.

### What is explicitly NOT being built

- **A commercial predictor or SaaS product.** Five prior research passes established that value in this industry accrues to wet-lab capability, proprietary data, and molecule ownership — not standalone software. The founder has no lab by design.
- **Level 2 (whole-enzyme expression prediction).** CUT. Reasons: (1) the required negative-outcome data does not exist publicly — established across five independent research threads; (2) TargetTrack failure causes are unparseable at acceptable precision, conceded in the original spec; (3) MPB-EXP already covers heterologous expression prediction across 88 species using a pretrained protein transformer with transfer learning at ~0.78 average accuracy. Do not reopen without new evidence on all three points.
- **A hosted live-submission leaderboard**, at least initially. Ship a versioned evaluation package. Stand up a board only on demonstrated demand.

### Success criteria, in priority order

1. **A published benchmark with a DOI that another group cites or uses.** This is the durable asset.
2. **A preprint answering the genome-vs-physiology question**, with an honest negative result being fully acceptable.
3. **Demonstrated ML engineering depth** — the founder can point to a FiLM-conditioned model, deep ensembles, conformal prediction, and a transfer-learning evaluation as work he did.
4. **Optionality** — a credible artifact for a genomics/ML role, a PhD application, or a DARPA-style proposal.

Explicitly **not** a success criterion: revenue, users, or a startup. If those emerge, they are upside.

---

## PART II — CONSTRAINTS AND KNOWN RISKS

### Hard constraints

- **No wet lab, ever.** Any step requiring newly generated experimental data is a scope violation, not a caveat.
- **Solo builder**, part-time-to-full-time, Python-strong, ML-developing.
- **Compute budget under $200 total.** The entire plan runs in well under 50 GPU-hours. Do not fine-tune a 40B-parameter model.
- **Public data only**, with licensing respected per Part V.

### The three known structural weaknesses

These are not solvable. They must be stated plainly in every output.

**1. n_hosts ≤ 6.** Three bacteria with full library coverage (*E. coli*, *B. subtilis*, *P. aeruginosa*), plus three more (*S. enterica*, *V. natriegens*, *C. glutamicum*) at ~241 sequences via the RS241 subset. The unit of generalization is the **host**, not the sequence. Any uncertainty estimate that bootstraps over sequence clusters while holding the host fixed measures the wrong thing and will produce misleadingly tight intervals. **Host-level uncertainty must be reported, and the n≤6 ceiling stated in the abstract.**

**2. The data is from 2018 and is the only dataset of its kind.** Johns et al. 2018 (*Nat. Methods* 15:323–329, BioProject PRJNA431139). Eight years old. The Tsinghua group behind DeepCROSS used it and then generated their own new MPRA data — which is evidence that they judged the 2018 data insufficient alone. Treat that as a signal about the ceiling, not a coincidence.

**3. Part of the field routes around the problem rather than predicting it.** Orthogonal-machinery approaches (T7 RNAP-based "universal" expression modules) partly neutralize host variation; one group reports striking linear relationships in circuit activity across non-model microbes using such modules. Users of orthogonal machinery are not in the addressable audience. This caps relevance but does not touch the scientific question.

### Competitive watchlist and flip conditions

| Entity | What they have | Flip condition | Response |
|---|---|---|---|
| **Wang lab, Tsinghua** (DeepCROSS, *Nat Commun* 2025) | Cross-species RS design, 2 bacteria, own MPRA capability, already uses the Johns data | Posts a preprint framed on unseen-host transfer or a benchmark | Narrow to the calibration-protocol novelty; add hosts they lack; do not compete head-on — they can measure, you can only wait |
| **Bernstein lab** (chassis-effect physiology work) | The published objection this project tests; broad-host-range experimental capability | Publishes a genome-vs-physiology predictive comparison | The benchmark becomes more valuable, not less; cite and incorporate |
| **de Boer lab, UBC** (GAME) | Benchmark framework with a species-mapping layer, human/mammalian only | Announces a microbial module | Approach as a contribution channel rather than a rival — this is the preferred outcome |
| **"From Context to Code"** (bioRxiv 2023) | Transformer models for cross-species DNA functionality, host-specific 5′ element design | Already published or commercialized | **UNRESOLVED — must be checked in Gate 1** |
| **Evo 2 / NTv3 / Borzoi successors** | Scale | Convincing cross-host quantitative activity results beating supervised baselines | Model contribution obsolete; benchmark becomes the evaluation standard — keep it, drop the model ambition |

Current empirical state, for context: foundation models **lose** to lightweight supervised baselines on regulatory tasks (DART-Eval, NeurIPS 2024; "Specialized Foundation Models Struggle to Beat Supervised Baselines," ICLR 2025). No major benchmark suite contains a cross-host quantitative regulatory-activity task. That gap is the opening.

### Time-bounded external item

**DARPA BTO AIxBio, BAA HR001126S0003**, open through **September 2026**, explicitly names sparse-data predictive models. Roughly seven weeks out from this charter. If funding matters to this track, the window closes before a 12-week build completes. Decide early whether to pursue.

---

## PART III — THE PLAN

Nine gates. Each has a deliverable, a binary pass condition, a kill criterion, and a fallback. **No phase may fail into nothing.**

Time estimates assume roughly 20–25 focused hours per week. Adjust proportionally.

---

### GATE 1 — Reality check (2–3 days, $0)

**This runs before any modeling code exists.** Everything downstream depends on numbers not yet personally verified.

**Deliverable:** A memo answering five questions from primary sources.

1. **Exact dataset shape.** Open the Johns et al. 2018 supplementary tables. Count rows. How many regulatory sequences have real measured values in *all three* primary hosts, separately for transcription and translation? What is the actual part length? What is the *B. subtilis* active fraction? (All three of these were previously taken from paywalled text or secondary summary and are unverified.)
2. **RS241 integrity.** Is the six-host subset genuinely ~241 sequences across six hosts with usable values, or is coverage patchier?
3. **"From Context to Code" disposition.** bioRxiv 2023, transformer models for cross-species DNA functionality. Published? Abandoned? Commercial? Who are the authors and where are they now?
4. **Physiology data availability.** For the six hosts, are physiological metrics recoverable from published literature — growth rate, doubling time, ribosome content, RNAP abundance, proteome allocation? This determines whether the genome-vs-physiology comparison, the core scientific question, is runnable at all.
5. **Funder check.** NIH RePORTER (via the v2 POST API), NSF, DOE, ERC, UKRI: is anyone funded to build cross-host regulatory prediction or a benchmark for it? The prior negative result is provisional.

**Pass condition:** ≥5,000 sequences with usable three-host values AND *B. subtilis* has enough active sequences to serve as a held-out host AND no funded direct competitor.

**Kill criteria:**
- Three-host overlap below ~3,000 sequences → the benchmark is too small to be credible.
- *B. subtilis* active fraction so low it cannot be scored → the stringent test is gone; only the easy proteobacterial hop remains, which is not a real cross-host claim.
- A funded competitor with standing is building this.
- "From Context to Code" turns out to have already published this benchmark.

**Fallback if killed:** A short "state of cross-host regulatory prediction data" note documenting why the benchmark cannot be built from public data. Genuinely useful to the field and cheap to write.

---

### GATE 2 — Data foundation (1 week)

**Deliverable:** `data/` containing Parquet tables — one for the three-host library, one for RS241, one host-feature table, one host-physiology table. Plus MMseqs2 cluster assignments and frozen split files.

Requirements:
- Sequences deduplicated to RS level **before** splitting (the library contains double-barcoded replicates; barcode-level splitting would leak).
- MMseqs2 `easy-cluster --min-seq-id 0.5 -c 0.8`, whole clusters assigned to single folds.
- **Two host representations built in parallel:** (a) ~40-D genomic features — sigma-factor complement, anti-Shine-Dalgarno free energy, tRNA copy numbers and tAI, codon usage, GC, chaperone/heme/RNAP annotations; (b) a physiology feature vector from Gate 1's literature findings. **The comparison between these two is the project's scientific core.**
- Every host feature unit-tested against a known value (e.g. the *E. coli* anti-SD sequence, published *E. coli* tAI ranks).
- Within-host quantile normalization fit **on training folds only**, inside each fold.
- **RS241 reserved entirely for held-out-host evaluation. Never trained on. Enforce in code, not by discipline.**

**Pass condition:** Splits frozen and hashed; every host feature unit-tested; a leakage audit script confirms zero RS overlap across folds and zero RS241 leakage into training.

**Kill criterion:** Physiology features are unrecoverable for ≥3 of 6 hosts → the central comparison is dead. Continue with genomic features only, but the paper's headline changes and must be renegotiated.

**Fallback:** The curated dataset and splits are themselves publishable as a data descriptor even with no model.

---

### GATE 3 — Baselines before the model (4–5 days)

**Deliverable:** Four baselines, evaluated under the frozen splits with leave-one-host-out.

1. **Mean/majority predictor** — the floor.
2. **Per-host model trained only on N calibration examples**, N ∈ {0, 10, 30, 100, 300, 1000, 3000}. **This is the real competitor.** A practitioner can simply measure N parts in their host. The project must beat that.
3. **Free per-host lookup embedding.** Note carefully: this cannot run on an unseen host by construction. It must therefore be given an explicit, documented fallback (mean of training-host embeddings, or nearest phylogenetic neighbor) or the comparison is undefined. **Specify the fallback before running.**
4. **Salis Promoter Calculator / RBS Calculator** where applicable — the biophysical baseline.

**Pass condition:** All four run end-to-end on frozen splits and produce host-level results with host-level uncertainty.

**Note on a corrected error:** an earlier version of this plan made "beat the free per-host embedding" the primary kill gate. That test cannot be failed, because the baseline cannot produce a number on an unseen host. **The real gate is beating the per-host-trained-on-N baseline (Gate 5).** Baseline 3 is a diagnostic, not a gate.

---

### GATE 4 — The model (1 week)

**Deliverable:** A FiLM-conditioned CNN, trained and evaluated under the frozen splits.

- 4 convolutional layers, kernel widths 15/9/5/3, channels 128/128/64/64, ~200–300k parameters. **Kernel widths must be re-sized to the part length verified in Gate 1, not the assumed 165 bp.**
- Host vector enters via FiLM at layers 2 and 3 — multiplicative conditioning, because host identity should modulate *which* sequence motifs matter rather than shift the output.
- Two output heads: transcription and translation as separate tasks. **Not a ratio** — a ratio destroys the translation-initiation signal that distinguishes hosts.
- Two-stage: active/inactive classifier plus a strength regressor on actives. The classifier is the product-relevant output ("will this fire at all in your host").
- **Trained twice: once with genomic host features, once with physiology host features.** This is the experiment.

**Pass condition:** Both variants train, converge, and produce leave-one-host-out predictions with host-level uncertainty.

---

### GATE 5 — THE KILL GATE (3–4 days)

**The decision point. Nothing beyond this gate is authorized until it clears.**

**Deliverable:** Calibration-cost curves and the pre-registered hypothesis tests.

**Pre-registered hypotheses — write these down before looking at results:**

- **H-MAIN (the gate):** The cross-host model at N=100 calibration examples matches or beats the per-host-only baseline at N ≥ 3,000, on *B. subtilis* held out, with host-level uncertainty reported. That is a 30× data-efficiency advantage. (100× is a stretch target only; quote 30× until data says otherwise.)
- **H-SCIENCE (the paper's core):** Genomic host features and physiology host features differ measurably in leave-one-host-out performance. **Direction not predicted.** If physiology wins, the Bernstein claim is supported and that is a strong, citable result. If genomics wins, the Bernstein claim is challenged and that is a stronger one.
- **H-DIAGNOSTIC:** The host-biology representation outperforms the free per-host embedding under its documented fallback — evidence the model learned host biology rather than host identity.

**Pass condition:** H-MAIN holds on *B. subtilis* with honest host-level uncertainty.

**Partial pass:** H-MAIN holds on *P. aeruginosa* (the easy proteobacterial hop) but not *B. subtilis* → proceed, but the claim narrows to "host-dependent" and the abstract must say so.

**Kill criterion:** H-MAIN fails on both hosts.

**Fallback if killed — and this is the important part:** The benchmark ships anyway, with a negative result. "We built the first cross-host regulatory activity benchmark, evaluated genomic and physiological host representations plus foundation models, and none beat simply measuring 100 parts in your target host." **That paper is more useful to the field than a marginal positive**, it directly answers a published question, and it is a legitimate NeurIPS D&B or *Scientific Data* contribution. The model work is not wasted; it becomes the baseline suite.

There is no outcome at Gate 5 that produces nothing.

---

### GATE 6 — Foundation model head-to-head (1 week)

**Deliverable:** Evo 2 and NTv3 (and any other available genomic FM) evaluated on the benchmark, frozen-embedding or likelihood-based, no fine-tuning.

This is the single most quotable result available: *here is a benchmark, and here is a 250k-parameter host-conditioned model against a 40B-parameter foundation model.* Existing evidence suggests the small model wins on regulatory tasks. If it does, that is the headline. If it does not, that is a bigger headline and the benchmark becomes the standard by which it was shown.

**Pass condition:** At least two FMs evaluated with documented methodology and honest handling of the fact that neither was designed for this task.

---

### GATE 7 — Uncertainty and ablations (1 week)

**Deliverable:** 5-model deep ensembles, split-conformal prediction intervals, reliability diagrams, ECE, conformal coverage at 80/90%. Host-feature ablations showing which features carry the transfer signal.

**Why this matters:** the product-relevant output is a risk flag, and a risk flag with no calibration is worthless. Also the section that most distinguishes competent ML work from a demo.

**Pass condition:** Calibration reported honestly, including where it fails.

---

### GATE 8 — Benchmark packaging (1.5 weeks)

**Deliverable:** The durable asset.

- pip-installable package, permissive license on new work
- Hugging Face Datasets mirror plus **Zenodo DOI**
- **Withheld test labels** with server-side or documented held-out evaluation (following ProteinGym and the DREAM promoter challenge)
- Fixed folds so out-of-fold predictions are the evaluation
- Baseline scores for everything from Gates 3, 4, and 6
- Single `make reproduce` regenerating every reported number from raw tables
- Pinned environment lockfile, all seeds set
- **Licensing quarantine, enforced in the directory structure:** DeepCROSS-derived data (MIT academic-only, commercial use requires author contact) and PromoGen2 (CC-BY-NC-4.0) in separately labeled partitions so the core benchmark stays cleanly licensed. Johns supplementary tables are Springer Nature copyright — distribute derived activity values with attribution and verify reuse terms; do not re-host raw tables without checking.

**Pass condition:** A clean clone reproduces every number with one command.

---

### GATE 9 — Publication and distribution (2 weeks)

**Deliverable:** Preprint plus outreach.

- bioRxiv preprint. **The abstract must state n_hosts ≤ 6 and that this caps generality.**
- Primary venue: NeurIPS Datasets & Benchmarks. Archival: *Scientific Data* descriptor. Timing: LMRL/MLCB workshop.
- **Approach the de Boer lab about contributing a microbial module to GAME** rather than launching a rival leaderboard. This converts the largest benchmark-side competitor into a distribution channel and borrows standing the founder does not have. Highest-leverage single action in the plan.
- Notify the Bernstein group directly: their published claim was tested. Whichever way it resolved, they are the most likely citers.
- Contact the Wang lab (Tsinghua) — generous framing, offer inclusion.

**Pass condition:** Preprint posted, DOI minted, at least three groups contacted.

---

## PART IV — TOTAL SHAPE

| Gate | Content | Time | Cumulative |
|---|---|---|---|
| 1 | Reality check | 3 days | 3 days |
| 2 | Data foundation | 1 wk | ~1.5 wk |
| 3 | Baselines | 5 days | ~2.5 wk |
| 4 | Model | 1 wk | ~3.5 wk |
| 5 | **KILL GATE** | 4 days | ~4 wk |
| 6 | FM head-to-head | 1 wk | ~5 wk |
| 7 | Uncertainty + ablations | 1 wk | ~6 wk |
| 8 | Packaging | 1.5 wk | ~7.5 wk |
| 9 | Publication | 2 wk | ~9.5 wk |

**~10 weeks** at 20–25 hours/week. Compute under $200. The kill gate lands at week 4 — before the expensive packaging and writing work.

---

## PART V — STANDING RULES

These apply to every gate and to every prompt written for the Claude Code session.

1. **Never train on RS241.** Enforced in code with an assertion, not by memory.
2. **Fit all normalization and clustering on training folds only**, inside the fold loop.
3. **Report host-level uncertainty, not just cluster-level.** The generalization unit is the host and n ≤ 6.
4. **Distinguish evidence from inference in every output.** Label estimates as estimates. Cite the file and record ID for anything factual.
5. **A negative result is a deliverable.** Never tune toward a positive.
6. **Verify before trusting a stated number** — including numbers in this charter. Several were inherited from secondary sources and are flagged as such.
7. **No new experimental data, ever.** If a step needs it, the step is out of scope.
8. **Licensing quarantine is structural**, not a note in a README.
9. **Anti-retrofit:** if evidence points somewhere other than this plan, say so. Do not reshape findings to fit the charter.

---

## PART VI — WORKFLOW

Execution runs in a Claude Code session in the `crosshost/` folder. Prompts are written per gate and pasted in. Each prompt ends with a mandatory **STATUS REPORT** whose format is fixed, so that the report alone is sufficient to decide the next move without a follow-up question.

Required status report structure:

```
=== CROSSHOST STATUS REPORT ===
GATE: [number and name]
OUTCOME: PASS | PARTIAL PASS | FAIL | BLOCKED

VERIFIED FACTS
  [each with the file/record it came from]

NUMBERS PRODUCED
  [table; state which are final and which provisional]

PASS CONDITION
  [restated verbatim] → MET / NOT MET, with the evidence

KILL CRITERIA
  [each restated] → TRIGGERED / NOT TRIGGERED

WHAT I COULD NOT DO
  [failed endpoints, blocked sources, weaker fallbacks used,
   assumptions made where data was missing — mandatory section,
   never omit, never soften]

CONTRADICTIONS WITH THE CHARTER
  [anything found that contradicts a stated premise — the most
   valuable thing to report]

FILES WRITTEN
  [paths]

STATE
  [contents of state.json]

DECISION NEEDED FROM GABRIEL
  [specific question, or "none — ready for next gate"]
=== END ===
```

The agent reports and proposes. It does not decide whether to proceed, and it does not assess whether the project is worthwhile.
