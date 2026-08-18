# GATE 8.5 MEMO — Closing the Open Objections

> **UPDATE — Gate 8.6 retracted Task 1B's headline finding.** The
> regression-to-the-mean control found that 10 of 11 "winning" cells below
> are explained by the reference value alone, with no sequence
> contribution. See `out/GATE8_6_MEMO.md` for the full control and the
> corrected picture. The rest of this memo (1A, 1C, 1D, and all of Task 2)
> is unaffected and stands as originally written.

**STOP-AND-REPORT NOTICE, per the task's own instruction ("If any of 1A–1D
produced a positive result, restructure around it and tell me"): TWO of the
four tasks produced real, statistically supported positive results — 1A
and 1B — and 1B's directly involves *B. subtilis*.**

- **Task 1A**: two alternative conditioning mechanisms beat both
  sequence-only and FiLM distinguishably, but only at **EC-transcription**.
  *B. subtilis* is unaffected.
- **Task 1B**: predicting the cross-host **shift** (given a reference
  host's measured value) recovers real, sequence-specific signal beyond a
  constant host-to-host offset for ***B. subtilis* transcription, in both
  directions (EC↔BS, PA↔BS), in 7 of 8 tested configurations** — a
  materially different picture from the traditional absolute-level
  framing, and one that **does** touch the host the central claim is about.
  *B. subtilis* **translation** shows no such signal (0 of 8), and neither
  does one of the two transcription directions for the EC-PA pair.

**Neither finding overturns H-MAIN or reverses the central claim** (see the
per-task verdicts for exactly what does and does not change), but both
require the central claim's language to be corrected from an unqualified
"conditioning does not help" / "cross-host signal is unrecoverable" to a
more precise, narrower statement. Task 3's reframing is held until this
memo section is reviewed, per the task's explicit instruction.

## Verdicts on 1A–1D — did any of the four recover cross-host signal the current framing says is unrecoverable?

| Task | Verdict | Status |
|---|---|---|
| 1A — alternative conditioning mechanisms | **YES, narrowly — 2 of 18 comparisons (both EC-transcription) beat sequence-only AND FiLM distinguishably; *B. subtilis* unaffected** | COMPLETE |
| 1B — predict the shift, not the level | ~~YES for *B. subtilis* transcription (7/8 configs beat the mean-shift baseline); NO for *B. subtilis* translation (0/8)~~ **RETRACTED by Gate 8.6 — regression to the mean, see update above** | COMPLETE, then retracted |
| 1C — restrict to the co-active subset | **NO — the hypothesis is contradicted for the host that matters, confirmed for only the pair that already didn't need it** | COMPLETE |
| 1D — separate the two stages for transfer | **NO signal recovered, but a genuine stage-specific asymmetry found** | COMPLETE |

---

## TASK 1A — alternative conditioning mechanisms [COMPLETE — LOAD-BEARING, POSITIVE RESULT]

Script: `scripts/78_conditioning_mechanisms.py` (training, 30 from-scratch
fits: 2 mechanisms x 3 held-out hosts x 5 folds), `scripts/82_conditioning_mechanism_comparison.py`
(bootstrap-CI'd comparison). Results: `out/results/gate8_5_conditioning_mechanisms.csv`,
`out/results/gate8_5_conditioning_mechanisms_verdict.csv`. Identical protocol
to Gate 4's own FiLM LOHO run: same 5 frozen folds, same seeds (=fold
index), same conv trunk, same epochs/patience/batch size/optimizer, same
masked joint 4-head loss, same `evaluate()` function, same 90% bootstrap CI
utility. Only the conditioning mechanism differs. Genomic host-feature
variant only (disclosed scope choice, see script docstring). Hypernetwork
(the third, optional mechanism) was not run — disclosed compute-bounded
choice, see script docstring.

**Full 18-cell comparison (3 hosts x 2 readouts x 3 mechanisms vs. sequence-only and FiLM-genomic):**

| Host | Readout | seq-only | FiLM-genomic | concat | perhost-avg | perhost-nearest |
|---|---|---|---|---|---|---|
| **EC** | **transcription** | 0.367 [0.348,0.387] | 0.371 [0.215,0.516] | **0.555** [0.533,0.580] | **0.615** [0.577,0.654] | 0.233 [0.171,0.299] |
| EC | translation | 0.562 [0.548,0.574] | 0.155 [-0.059,0.346] | 0.576 [0.555,0.599] | 0.524 [0.506,0.543] | 0.423 [0.377,0.466] |
| BS | transcription | 0.263 [0.246,0.282] | 0.131 [0.008,0.238] | 0.271 [0.252,0.291] | 0.277 [0.257,0.297] | 0.279 [0.259,0.299] |
| BS | translation | 0.291 [0.260,0.327] | 0.151 [0.111,0.191] | 0.298 [0.263,0.331] | 0.296 [0.267,0.326] | 0.303 [0.274,0.333] |
| PA | transcription | 0.479 [0.456,0.500] | 0.338 [0.237,0.434] | 0.470 [0.433,0.506] | 0.309 [0.282,0.337] | 0.481 [0.455,0.506] |
| PA | translation | 0.390 [0.366,0.417] | 0.131 [0.050,0.209] | 0.415 [0.390,0.448] | 0.376 [0.351,0.402] | 0.397 [0.370,0.425] |

Bold = the two cells where a mechanism beats both sequence-only and FiLM
with non-overlapping 90% CIs (`beats_both_distinguishably=True` in the
verdict table). Of 18 (host, readout, mechanism) comparisons, **2 win — and
both are EC-transcription.**

**Why this is very unlikely to be a bug (checked directly, not assumed):**
1. **Per-fold consistency.** Concat's 5 individual EC-transcription fold
   values are [0.514, 0.569, 0.536, 0.551, 0.606] — tight, all in a narrow
   band. Per-host-heads-avg: [0.559, 0.691, 0.563, 0.623, 0.641] — same
   pattern. A leakage bug would typically produce an implausibly high or
   unstable result, not a modest, cross-fold-consistent one.
2. **FiLM's own EC-transcription numbers on the SAME 5 folds are highly
   unstable by contrast**: [0.315, 0.553, 0.045, 0.606, 0.337] — ranging
   from near-zero to 0.61 fold-to-fold. This is the exact signature Task
   1A's own objection predicted: a FiLM γ/β generator fit from only 2
   training hosts is underdetermined and fold-sensitive, while concat's
   single unconstrained linear combination at the final layer is not.
3. **Code review**: both mechanisms reuse the identical trunk, identical
   `build_pooled_arrays` fold-filtering (train on the 2 non-held-out hosts'
   `fold != test_fold` rows, evaluate on the held-out host's `fold ==
   test_fold` rows — exactly Gate 4's own split), and the held-out host's
   genomic feature vector is the only "new" information reaching the model
   at eval time — the same public covariate FiLM was already allowed to
   see. No test-set information enters training in either mechanism.

**The averaging rule matters, and the nearest-host rule actively misleads.**
`perhost_heads_nearest` for EC selects only the single closest training
host's private head by genomic-feature Euclidean distance — which turns out
to be ***B. subtilis*** (distance 5.42, vs. *P. aeruginosa* at 8.26,
despite EC and PA being the two Gammaproteobacteria and BS the outgroup
Firmicute — genomic-feature-space proximity does not track phylogeny here).
Using only BS's head scores 0.233 — *worse* than sequence-only and FiLM
both — while averaging BS's and PA's heads scores 0.615. **Feature-space
nearness is a poor proxy for "which training host's head will transfer
well"; averaging is the safer default**, confirmed empirically, not just
assumed as originally stated in the task.

**What does NOT change: *B. subtilis*.** Every BS cell (both readouts) is
statistically indistinguishable across sequence-only, concat, and
per-host-heads (all three overlap heavily) — concat/perhost are only
mildly, non-distinguishably above sequence-only, and only clearly ahead of
FiLM specifically (FiLM is the worst performer in every BS cell, consistent
with the underdetermined-generator hypothesis, but sequence-only was
already competitive with or ahead of FiLM there too, per the original
Gate 5.5 finding). **No mechanism recovers cross-host signal for the host
the central claim is actually about.**

**PA shows an asymmetric failure mode worth flagging honestly:**
`perhost_heads_avg` at PA-transcription scores 0.309 [0.282, 0.337] — below
sequence-only (0.479) AND below FiLM (0.338), making it the *worst* of all
five systems in that one cell. Averaging two private heads is not a
universally safe choice; here it underperforms even the single-shared-model
baseline it was meant to improve on. Not explained further here — reported
as an observed asymmetry, not built on.

**VERDICT: does any mechanism beat sequence-only?** Yes, at exactly one
(host, readout) cell — **EC transcription** — where **two independent,
architecturally distinct mechanisms (concatenation and per-host-heads-
averaging) both win, distinguishably, against both sequence-only and
FiLM.** This is not "FiLM specifically is weak" (Task 1A's opening
objection) generalizing to "conditioning never helps" (Gate 5's original
framing) — it is a **third, more precise statement**: conditioning CAN
recover a real, reproducible signal, but only for the training/host
configuration where it did in this test (EC held out, trained on
BS+PA) — and does not for the configuration ( *B. subtilis* held out) the
project's central claim is actually about. **The paper's claim must be
corrected from "conditioning fails under multiple mechanisms" to
"conditioning fails for *B. subtilis* under three distinct mechanisms;
for *E. coli* specifically, two of three mechanisms show a real,
reproducible improvement over both sequence-only and FiLM."** This is a
narrowing of the central claim's generality, not a reversal of it — see
Task 3 (held pending review) for exactly how this is written into
`PAPER_FRAMING.md`.

---

## TASK 1D — separate the two stages for transfer [COMPLETE]

Reused the already-computed, already-saved Gate 4/5.5 zero-shot LOHO results
(no retraining) and split scripts/42's `evaluate()` output — which already
computes classifier (MCC/AUC) and regressor (Spearman rho) metrics
independently — into its own bootstrap-CI'd comparison. Script:
`scripts/79_two_stage_transfer.py`. Full table:
`out/results/gate8_5_two_stage_transfer.csv`; verdict table:
`out/results/gate8_5_two_stage_transfer_verdict.csv`.

**Classification stage (AUC):** of 12 (host x readout x FiLM-variant)
comparisons, FiLM beats sequence-only distinguishably (90% CI, no overlap)
in **1 of 12** — *E. coli* translation, physiology variant (seq-only AUC
0.551 [n/a CI shown as point since n_folds draws] vs. physiology-FiLM AUC
0.666, distinguishable).

**Regression stage (Spearman rho):** of 12 comparisons, FiLM beats
sequence-only distinguishably in **0 of 12.** Sequence-only wins
distinguishably in 9 of 12; the remaining 3 are statistically
indistinguishable but favor sequence-only in point estimate.

**Reading:** conditioning shows a real, if narrow, foothold on the binary
"will this fire at all" question in exactly one cell, and no foothold
whatsoever on "how strong will it be" anywhere tested. This refines the
project's blanket "conditioning does not help" claim to a **stage-specific**
one: the strength regressor is where the negative result is uniform and
strong; the classifier has one documented exception. Framed honestly, not
dismissed — see Task 3 evidence-hierarchy placement below.

---

## TASK 1B — predict the shift, not the level [COMPLETE — POSITIVE RESULT FOR B. SUBTILIS TRANSCRIPTION]

> **CORRECTION — Gate 8.6, dated after this section was written.** The
> regression-to-the-mean control (Task 1 of Gate 8.6) found that this
> section's headline finding does not survive scrutiny. 10 of the 11
> "winning" cells below are explained, as well or better, by a
> reference-value-only baseline with no sequence input at all — including
> every *B. subtilis*-involving transcription cell this section's verdict
> was built on. **This section's conclusion is retracted.** See
> `out/GATE8_6_MEMO.md` Task 1 for the full control and the one narrow
> exception that does survive (PA→EC translation only, unrelated to the
> claim below). Left in place, uncorrected in its own body text below, for
> the historical record of what was found and how it was later caught —
> matching this project's standing practice of disclosing corrections in
> place rather than silently rewriting (`out/GATE7_MEMO.md`,
> `out/GATE5_5_MEMO.md` dated corrections).

Scripts: `scripts/81_shift_prediction.py` (training, 120 fits: 6 ordered
host pairs x 2 readouts x 2 conditioning variants x 5 folds), analysis done
inline (results: `out/results/gate8_5_shift_prediction.json`,
`out/results/gate8_5_shift_prediction_summary.csv`). Target =
target-host-strength minus reference-host-strength, on the same
log1p(tx_norm) / protein_log10 scale used everywhere else in this project,
restricted to sequences with a usable strength value in both hosts. Same 5
frozen folds (standard train/test split within-fold, not LOHO — the
reference host's true value is handed to the model at both train and test
time, so there is no unseen-host generalization question here, only a
held-out-*sequence* one). **DISCLOSED COMPUTE-BOUNDED DEVIATION** from
"identical protocol": epochs=25/patience=6 vs. Gate 4's 35/8, justified by
much smaller per-pair training pools (294–9,741 rows vs. Gate 4's
46,000–58,000) — see script docstring. Two metrics: Spearman rho of the
model's own predictions (meaningless for the constant baselines, which have
zero prediction variance) and **fraction of shift-variance explained**
relative to each baseline (`1 - MSE(model)/MSE(baseline)`), which is the
metric that actually answers "can the model beat this baseline."

**Beats "assume it transfers unchanged" (zero-shift): 21 of 24
(pair, readout, variant) cells, distinguishably (90% CI lower bound > 0).**
The shift is very much learnable in the loose sense — knowing the sequence
plus the reference host's value beats assuming no change almost everywhere.

**Beats "predict the mean shift" (a per-pair constant offset, computed from
TRAIN-fold data only): only 11 of 24 cells, distinguishably.** This is the
metric that actually isolates *sequence-specific* shift information from a
trivial systematic host-to-host bias (e.g., "PA reads on average ~0.3
log-units higher than BS" would pass the zero-shift test for free without
containing a single bit of sequence information). Beating the mean-shift
baseline is the real bar, and it is cleared in fewer than half of all
cells — **but the 11 successes are not randomly distributed.**

**The pattern, broken out by whether *B. subtilis* is involved:**

| Group | Transcription: beats mean-shift | Translation: beats mean-shift |
|---|---|---|
| ***B. subtilis*-involving pairs** (EC↔BS, PA↔BS, 4 directions x 2 variants = 8 cells each) | **7 / 8** (all except BS→PA, no-conditioning variant — which the WITH-conditioning variant then does clear) | **0 / 8** |
| **EC-PA pair** (non-*B. subtilis*, 2 directions x 2 variants = 4 cells each) | 1 / 4 (PA→EC, with-conditioning only; EC→PA transcription no-condition is *distinguishably worse* than mean-shift, CI entirely negative) | 3 / 4 |

**This is real and directly relevant to the central claim.** Every
*B. subtilis*-involving transcription direction and variant except one
clears the mean-shift bar — a materially different result from every
absolute-level framing tested in this project (Gate 5's H-MAIN, Gate 5.5's
four-way comparison, Task 1A above), all of which found *B. subtilis*
essentially unpredictable from sequence alone. **Handed the reference
host's measured value, the model recovers real, sequence-specific
information about how a promoter's transcription strength will shift when
moving into or out of *B. subtilis* — something no absolute-level model in
this project has shown.** *B. subtilis* **translation shows the opposite
result: zero configurations beat mean-shift, several distinguishably
worse** (e.g. BS→PA translation, no-conditioning: fve_mean = −0.507
[−0.895, −0.157] — actively worse than a constant prediction). This is
consistent with this project's own repeated finding that *B. subtilis*
translation is the noisiest, most floor-dominated signal in the whole
dataset (`out/KNOWN_ISSUES.md` item 2) — the shift reframing does not
rescue a signal that may not be there to find.

**Does conditioning help, and is this the place it helps most (as the task
predicted)?** Inconsistent, not a clean yes. It clearly helps in two
places — BS→PA translation (fve_zero flips from −0.281 to +0.112, a real
qualitative change) and PA→EC transcription (fve_mean 0.144→0.235, and only
the with-conditioning variant clears the mean-shift bar distinguishably) —
and clearly hurts in two others (PA→BS transcription, rho drops
0.464→0.413; PA→BS translation, fve_mean drops from −0.102 to a
distinguishably-worse −0.270). No consistent direction. The task's
prediction that conditioning would help "most likely" in the shift framing
is only partially borne out.

**VERDICT: for the readout and host this project's central claim rests on
most (*B. subtilis*, transcription), shift-prediction recovers a real
signal absolute-level prediction does not. For the readout the claim's own
prior evidence already flagged as noisiest (*B. subtilis* translation), it
does not.** This narrows, but does not reverse, the central claim, and adds
a genuinely new, positive, actionable finding: a practitioner who already
has one host's measurement can predict the *direction and rough magnitude*
of a *B. subtilis* transcription shift from sequence, even though
predicting *B. subtilis*'s absolute transcription level from sequence alone
remains near-impossible. See Task 3 (held pending review) for how this
changes the evidence hierarchy.

---

## TASK 1C — restrict to the co-active subset [COMPLETE]

Script: `scripts/80_coactive_subset.py`. Results:
`out/results/gate8_5_coactive_correlation.csv` (Part 1),
`out/results/gate8_5_coactive_model_rerun.csv` +
`out/results/gate8_5_coactive_model_rerun_summary.csv` (Part 2, bootstrap
CI'd).

**A correction to the task's own premise, found while implementing it:**
the task asks to compute correlation "restricted to sequences active in
both hosts... alongside the pooled figures," implying this restriction was
new. It is not — `scripts/57` (Gate 5.5), whose numbers PAPER_FRAMING.md's
finding #1 already leads with, already restricts to active-in-both-hosts.
That number **is** the co-active-restricted figure, not a pooled one. The
genuinely new comparison (Part 1) is co-active vs. the previously-uncomputed
**pooled** figure (usable in both hosts, activity status ignored).

**Part 1 result — raw measurement correlation, pooled vs. co-active:**

| Pair | Readout | n pooled | ρ pooled | n co-active | ρ co-active | Direction |
|---|---|---|---|---|---|---|
| EC-BS | transcription | 14,088 | 0.655 | 3,668 | **0.258** | co-active WORSE, sharply |
| EC-BS | translation | 866 | 0.161 | 866 | 0.161 | identical (masks coincide for this pair) |
| **EC-PA** | transcription | 19,643 | 0.621 | 9,741 | **0.754** | co-active better |
| **EC-PA** | translation | 4,513 | 0.724 | 3,826 | **0.742** | co-active slightly better |
| BS-PA | transcription | 11,969 | 0.508 | 2,099 | **0.257** | co-active WORSE, sharply |
| BS-PA | translation | 352 | 0.297 | 314 | 0.263 | co-active slightly worse |

**The pattern is exactly the opposite of what the task's hypothesis
predicts, for every pair that includes *B. subtilis*.** Restricting to
sequences that manage to activate in *both* hosts of a pair does not reveal
hidden conservation — for EC-BS and BS-PA transcription specifically, it
roughly **halves** the correlation (0.655→0.258, 0.508→0.257). The only
pair where co-active restriction helps is EC-PA — already this project's
strongest, most-Gammaproteobacteria-close pair, the one place the
hypothesis was least needed.

**Part 2 result — does the model's predictability improve on the co-active
subset?** Re-scored the ALREADY-TRAINED, already-saved sequence-only and
genomic-FiLM LOHO checkpoints (no retraining) on the co-active-restricted
test subset per (held-out host, training-partner host, readout, fold).
Bootstrap-CI'd summary, sequence-only system:

| Held-out vs. partner | Readout | Full-test ρ [90% CI] | Co-active-subset ρ [90% CI] |
|---|---|---|---|
| EC vs. BS | transcription | 0.367 [0.348, 0.387] | **0.079** [0.041, 0.115] |
| EC vs. BS | translation | 0.562 [0.548, 0.574] | **0.237** [0.178, 0.292] |
| **EC vs. PA** | transcription | 0.367 [0.348, 0.387] | **0.526** [0.501, 0.551] |
| **EC vs. PA** | translation | 0.562 [0.548, 0.574] | 0.475 [0.435, 0.515] (slightly lower) |
| BS vs. EC | transcription | 0.263 [0.246, 0.282] | 0.283 [0.270, 0.297] (small overlap-adjacent gain) |
| BS vs. PA | transcription | 0.263 [0.246, 0.282] | 0.252 [0.211, 0.292] (flat) |
| PA vs. BS | transcription | 0.479 [0.456, 0.500] | **0.391** [0.354, 0.428] (lower) |
| PA vs. EC | transcription | 0.479 [0.456, 0.500] | 0.504 [0.478, 0.530] (small gain) |

The same EC-PA-only pattern holds for model predictability, not just raw
correlation: co-active restriction distinguishably **helps** prediction
only where EC and PA are both host and partner, and distinguishably
**hurts** it for the EC-BS pair specifically (both directions, both
readouts) — the pair this entire project's central claim is about.

**HONESTY CHECK (required by the task):** this is exactly the
outcome-conditioning the task warned about, so read Part 2 as a statement
about where cross-host conservation lives, not an improved predictor — a
real deployment could not know in advance which of a new host's sequences
would be "co-active" without already having measured them.

**Verdict: the co-active-subset hypothesis ("activity is host-specific,
strength above threshold is conserved") is NOT supported — it is actively
**contradicted** for *B. subtilis*, the host the central claim depends on.**
A plausible (not verified, stated as speculation) mechanistic reading: the
minority of sequences that manage to activate in the low-activity-rate host
(*B. subtilis*, 17.9% active) may do so via a different, less
title-transferable route than the majority-active hosts, making their exact
strength level *less* predictable from another host's measurement, not
more. No further claim is built on this speculation.

---

## TASK 2 — reconcile with the external adversarial review

Per-claim verification against primary sources, not the second-hand summary.

### 2a. Finding #1 (cross-host measurement correlation) — **CONFIRMED, replication with quantification**

Source: `raw/nihms945382_bodytext.txt` (Johns et al. 2018 main text, verified
against `raw/nihms945382.pdf`, 25pp, title *"Metagenomic mining of
regulatory elements enables programmable species-selective gene
expression"* — the paper's own title already centers species-selectivity)
and `raw/NIHMS945382-supplement-2.docx` (supplementary figure legends,
extracted directly, not summarized).

**Exact main-text wording** (body text, ~position 14,679): *"We
additionally assessed the transcription activity and translation efficiency
of 212 RSs that contained both active transcription and translation data
across all species (Suppl. Figure S13a). Translation efficiency of each RS
was determined by normalizing its GFP level to its transcription level.
Interestingly, we find that between recipients, only E. coli and P.
aeruginosa showed significant correlations between regulatory sequences in
terms of transcription levels and translation efficiencies."*

**Exact Supplementary Figure S13 caption:** *"Cross species and in silico
comparisons of gene expression levels. (a) Correlation of regulatory
sequence activity in terms of transcription level and translation
efficiency (calculated as the ratio of GFP protein levels and transcription
levels) between recipient species. Each point corresponds to a single
regulatory sequence that has measurable transcription and translation data.
Pearson correlation coefficient (r) and statistical significance values (p)
are shown for each subplot (n=212 for all six panels)."* Six panels = 3
host pairs x 2 metrics (transcription, translation efficiency).

**Exact Supplementary Figure S15 caption:** *"Cross-species transcription
and translation level correlations. (a) Pairwise Pearson correlation of
transcription (blue triangle) and translation (green triangle) activity
profiles of the RS241 library across six host species... Numbers in each
box correspond to the Pearson correlation coefficients (n = 241)."*

**Verdict: CONFIRMED**, with the exact figure numbers and N the second-hand
report claimed (S13/n=212/three hosts; S15/n=241/RS241 six hosts), and the
qualitative claim ("only E. coli and P. aeruginosa showed significant
cross-host correlations") is Johns et al.'s own stated result, in their own
words, not an inference.

**What this project's finding #1 adds, precisely (three real, disclosed
differences from Johns's S13, not just a re-statement):**
1. **Different translation operationalization.** Johns's S13 correlates
   *translation efficiency* (GFP/transcription ratio) across hosts. This
   project's finding #1 correlates *raw translation level* (protein_log10)
   — a more directly practitioner-relevant quantity ("how much protein do I
   get"), not a derived ratio, and NOT the same thing Johns measured.
2. **Different statistic and restriction.** Johns used Pearson correlation
   on a three-way-intersected n=212 (active in transcription AND
   translation, in ALL THREE hosts simultaneously). This project uses
   Spearman (more robust to the heavy-tailed activity distributions
   documented throughout this project) on a pairwise-only co-active
   restriction (active in both hosts of a given PAIR, not all three) —
   yielding far larger N per pair (this project's EC-PA transcription:
   n=9,741 vs. Johns's n=212 for all six S13 panels combined).
3. **Genuinely new: the disattenuation analysis (Gate 8)** connecting
   measurement reliability to the achievable ceiling on model performance —
   this does not exist in Johns et al. at all; they report the correlation
   descriptively and do not connect it to what a downstream predictive
   model could achieve.

**WHAT I COULD NOT DO:** raw/ contains the S13/S15 figure *legends*
(via the supplementary docx) and the body-text summary sentence, but not
the S13/S15 figure *images* themselves (only main-text Figures 1-5 have
image files in raw/; no full supplementary-figures PDF is present). I could
not independently read the exact per-cell Pearson r values Johns reported.
Given the qualitative claim is confirmed unambiguously via exact quoted
text (both the main-text sentence and the caption), and this project's own
independently-computed numbers already show the same qualitative pattern
via a different method, I judged fetching the external PMC/journal
supplementary-figure PDF to extract exact per-cell Pearson r values as not
required to close this verification — disclosed, not silently skipped.

**REQUIRED ACTION (Task 3):** demote finding #1 from "our discovery" to
"replication with quantification, credited to Johns et al." — done below.

---

### 2b. Finding #3 (sequence-only beats two genomic FMs) — **PARTIALLY CONFIRMED; one sub-claim false, but the core objection holds and is more serious than the second-hand report implied**

Verified by an independent research agent against the actual PromoGen2
paper (Xia et al., *"Design prokaryotic cis-regulatory elements using
language model,"* Nucleic Acids Research 54(4):gkag122, 2026; PMC12907563 —
confirmed as the paper behind the `jinyuan22/promogen2-base` HuggingFace
model this project actually used, via its Data Availability section citing
the exact HF URLs and GitHub repo).

**Sub-claim "NT and DNABERT-2 both failed": NOT CONFIRMED for DNABERT-2.**
DNABERT-2 was **never evaluated** on PromoGen2's zero-shot promoter-strength
task at all. Its only appearance in the paper is in an unrelated
sequence-*generation* benchmark, with the explicit statement: *"We did not
include Nucleotide Transformer or DNABERT in this benchmark because these
models do not have the ability to directly generate sequences and therefore
cannot be evaluated under the same generative setting."* That sentence is
about generation capability, not zero-shot scoring performance — "not
included" is not "failed." Nucleotide Transformer *was* tested on the
zero-shot scoring task (Fig. 2D/2F, via pseudo-perplexity), but its exact
value appears only as an unlabeled figure bar, not a quoted number.

**Sub-claim "avg Spearman ~0.475–0.481": NOT CONFIRMED.** No such figure
appears anywhere in the published paper (full-text search: zero hits). The
paper's actual stated headline number is *"improving the average Spearman
correlation from 0.27 to 0.50 compared to the best baseline [MegaDNA]"* —
0.50, not 0.475–0.481. (0.68+0.52+0.30)/3 = 0.50 exactly, consistent with
the per-host figures below, so 0.50 is very likely the correct number and
the second-hand report's 0.475-0.481 is simply wrong.

**Sub-claim "per-host EC/PA/BS ≈ 0.68/0.52/0.30, same dataset": CONFIRMED
verbatim.** Fig. 2D text: *"particularly in P. aeruginosa (Spearman rho =
0.52) and E. coli. (Spearman rho = 0.68)... In B. subtilis, PromoGen2 also
achieved better performance (0.3)."* The task dataset is explicitly cited
as *"promoter sequences with measured strength... from Johns et al.,
Metagenomic mining of regulatory elements enables programmable
species-selective gene expression, Nat Methods 2018;15:323–9"* — the
**same** dataset and species this project uses. Directly comparable.

**THE COMPARISON THAT ACTUALLY MATTERS, computed directly (not in the
second-hand report, found by pulling this project's own numbers once the
above was confirmed):**

| Host | Sequence-only (this project, zero-shot LOHO) | PromoGen2, this project's frozen-embed+head protocol (Gate 6) | PromoGen2, its own native published protocol |
|---|---|---|---|
| EC | 0.367 [0.348, 0.387] | 0.495 [0.448, 0.540] | **0.68** |
| BS | 0.263 [0.246, 0.282] | 0.243 [0.202, 0.277] | **0.30** |
| PA | 0.479 [0.456, 0.500] | 0.490 [0.463, 0.513] | **0.52** |

(Sequence-only and this project's own PromoGen2 numbers: transcription,
zero-shot LOHO, 90% bootstrap CI. PromoGen2 native: single published point
estimate, no CI available.)

This project's own Gate 6 comparison already shows PromoGen2 (via this
project's frozen-embedding+shallow-head protocol) beating sequence-only
distinguishably at EC-transcription-zero-shot — one of only 2 of 37 cells
where an FM wins (verified directly against `out/results/gate6_fm_verdict_table.csv`;
both winning cells are EC-transcription, DNABERT-2 and PromoGen2). **But
PromoGen2's own native protocol scores even higher at EC (0.68 vs. this
project's measured 0.495)**, and is directionally ahead of sequence-only at
**all three hosts**, not just EC.

**Verdict: the "DNABERT-2 failed" claim is false and should not be
repeated. But the underlying objection — that this project's FM evaluation
protocol may understate what these models can do, making the
"sequence-only wins nearly everywhere" framing less robust than presented —
is CONFIRMED, and by a wider margin than the second-hand report's own
(wrong) numbers suggested.** This project evaluated both FMs uniformly via
frozen-embedding+shallow-head (Gate 6's own defensible choice, for a fair
head-to-head between architecturally different FMs) rather than giving
PromoGen2 its own best/native zero-shot protocol (direct
likelihood/perplexity-based scoring, which is what the model was actually
trained to do). That protocol choice measurably cost PromoGen2 performance,
most visibly at EC.

**REQUIRED ACTION (Task 3):** demote finding #3 to an appendix, drop any
"model was too small" capacity-limit framing (irrelevant here — the issue
is protocol fairness, not model size), and state the native-vs-measured gap
explicitly. Translation readout and RS241 were not checked against
PromoGen2's paper (its own benchmark is transcription/"promoter strength"
only; a translation-efficiency comparison isn't well-defined against their
numbers) — scope limit disclosed, not silently extended.

---

### 2c. Finding #8 (feature ablation vs. Bernstein claim) — **CONFIRMED, with an important two-paper scoping correction**

Verified by an independent research agent, which located and read both
papers directly (PMC full text + journal page, not abstracts alone):

- **Chan, Baldwin & Bernstein 2023** — bioRxiv 2023.02.27.529268 → published
  as *"Revealing the Host-Dependent Nature of an Engineered Genetic Inverter
  in Concordance with Physiology,"* BioDesign Research 5:0016 (Aug 2023),
  DOI 10.34133/bdr.0016. Host panel: **6 hosts, 3 genera, zero
  Stutzerimonas** — *E. coli*, *Halopseudomonas aestusnigri*,
  *Halopseudomonas oceani*, *Pseudomonas deceptionensis*, *Pseudomonas
  fluorescens*, *Pseudomonas putida*. Physiology: direct wet-lab (plate
  reader growth/fluorescence, qPCR plasmid copy number, computed codon
  adaptation index) — **no transcriptomes**, no reference-proteome
  proxying. Contains the charter's "genetic inverter" quote, in the
  Abstract, verbatim.
- **Chan & Bernstein 2024** — bioRxiv 2024.02.15.580380 → published as
  *"Pangenomic landscapes shape performances of a synthetic genetic circuit
  across Stutzerimonas species,"* mSystems 9(9):e00849-24 (Aug 2024), DOI
  10.1128/msystems.00849-24. Host panel: **6 hosts, ALL Stutzerimonas** —
  *S. chloritidismutans*, *S. perfectomarina*, *S. degradans*, *S.* sp.
  pgs16, *S.* sp. pgs17, *S. stutzeri* (*E. coli* used only as a cloning
  reference, not a performance-panel host). Physiology: direct wet-lab
  (plate-reader growth curves **and** RNA-seq transcriptomes, Illumina
  NovaSeq 6000, DESeq2). Contains the charter's "impractical to postulate"
  quote, in the Discussion, verbatim.

**Verdict: CONFIRMED, but the claim must be attributed to the 2024 paper
specifically, not "the Bernstein lab" generically.** The charter's Part I
quote (both sentences presented together) actually splices language from
two different papers with two different host panels — the 2023 paper's
6-host panel spans three genera and does not use transcriptomes; the 2024
paper's 6-host panel is exactly what the objection describes (closely
related Stutzerimonas, growth curves + transcriptomes, direct
measurement). This composite-citation issue was not previously flagged
anywhere in this project's records and is itself worth disclosing.

**Why this collides with finding #8 (the feature-group ablation):** Gate
7's ablation zeroed out groups within the 37-D **genomic** feature vector
(sigma-factor, anti-SD, tAI/codon, RNAP, chaperone/heme — all computed from
**reference genome/proteome annotations**), not the 6-D physiology proxy
vector, and found no group carrying detectable signal. But (a) this
ablation never touched the physiology vector at all — the vector that is
actually analogous to what Bernstein's papers measure — and (b) even the
physiology vector this project built is a coarse 6-feature
**reference-proteome-derived proxy** (PaxDb-style aggregate abundance
fractions), categorically different from either Bernstein paper's **direct
per-experiment measurement** (growth curves + qPCR/transcriptomes measured
by that lab, on that host, in that experiment). Comparing a null result
from ablating reference-genome-derived features against a positive result
from direct wet-lab physiology measurement is not a fair test of the same
hypothesis, regardless of host range.

**REQUIRED ACTION (Task 3):** drop finding #8 as any kind of general test or
refutation of the Bernstein claim. Retain strictly scoped: "within our 6
coarse, reference-proteome-derived proxy features and a 3-primary-host
regime spanning two phyla, no genomic feature group carried a detectable
rank-correlation signal." Retain the phylogenetic scoping observation
(Bernstein's tested regime is close relatives, exactly where this project's
own data shows cross-host correlation is highest) as an observation, now
correctly attributed to the 2024 Stutzerimonas paper specifically — not a
counterargument to either paper.

---

### 2d. Reliability coverage limitation — **CONFIRMED already substantially disclosed; strengthened here**

`out/KNOWN_ISSUES.md` item 3 already states plainly: *"What we could NOT
establish: a direct reliability estimate for B. subtilis or P. aeruginosa
specifically, or for translation in any host."* `out/PAPER_FRAMING.md`
finding #1 already states the correction is *"anchored by an
empirically-measured E. coli transcription reliability of 0.91."* Both were
accurate before this gate. Strengthened per Task 3 below with one additional
explicit sentence tying the *B. subtilis* reliability gap directly to the
fact that *B. subtilis* is the host the central claim depends on most.

---

*(Tasks 1A, 1B, 1C Part 2, Task 3 evidence-hierarchy rewrite, Task 4
git/Zenodo, and the Part VI status report follow once the in-flight
background runs complete.)*
