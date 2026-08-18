**Task 1 verdict: the shift-prediction result does NOT survive the regression-to-the-mean control. 10 of 11 winning cells fail to beat a reference-value-only baseline; the one cell that survives is not among the *B. subtilis*-transcription cells the original headline was built on. Task 1B is retracted as a positive finding. The framing reverts to the negative, with Task 1A (EC-transcription only, narrow) as the sole surviving positive result.**

# GATE 8.6 MEMO — Validate 1B, then Reframe

Per Task 3's own instruction ("only proceed if 1B survives Task 1... if it
collapses, stop and report"), **Task 3 (the `PAPER_FRAMING.md` reframing
built around 1B) was NOT run.** `PAPER_FRAMING.md` is unchanged by this
gate. Everything below is Task 1 (the control) and Task 2 (characterization
of the two secondary findings, corrected for Task 1's outcome).

---

## TASK 1 — the regression-to-the-mean control [LOAD-BEARING, COMPLETE]

Script: `scripts/84_shift_regression_to_mean_control.py`. Results:
`out/results/gate8_6_reference_only_baseline.csv`,
`out/results/gate8_6_shift_vs_reference_correlation.csv`,
`out/results/gate8_6_shuffled_sequence_control.csv`,
`out/results/gate8_6_shift_controls_decisive_comparison.csv`. Figure:
`out/figures/gate8_6_shift_regression_to_mean_control.png`.

**Check 3 first (the mechanism), because it explains everything else.**
Spearman(shift, reference value), full co-usable population, no fold split:

| Pair | Readout | ρ(shift, reference) |
|---|---|---|
| EC→BS | transcription | **−0.707** |
| EC→BS | translation | −0.668 |
| EC→PA | transcription | +0.114 |
| EC→PA | translation | **−0.705** |
| BS→EC | transcription | −0.419 |
| BS→EC | translation | −0.564 |
| BS→PA | transcription | −0.270 |
| BS→PA | translation | −0.643 |
| PA→EC | transcription | **−0.649** |
| PA→EC | translation | +0.148 |
| PA→BS | transcription | **−0.812** |
| PA→BS | translation | −0.502 |

**Strong, pervasive regression to the mean, exactly as the task predicted.**
10 of 12 pairs show |ρ| > 0.4, several above 0.7. A sequence with a high
reference-host value mechanically has more room to fall than to rise, and
this dataset shows that mechanism dominating almost every pair.

**Checks 1+2 — reference-only baseline vs. the decisive comparison.** OLS
on the reference value alone (1 feature, fit per train-fold, evaluated on
the matching test-fold — same 5 folds as scripts/81) against the original
sequence+reference model, for exactly the 11 (pair, readout, variant) cells
scripts/81 reported as distinguishably beating the mean-shift baseline:

| Pair | Readout | Variant | Original (seq+ref) ρ | Reference-only ρ | Shuffled-seq ρ | Beats ref-only? |
|---|---|---|---|---|---|---|
| EC→BS | tx | no_cond | 0.432 | **0.709** | 0.557 | NO |
| EC→BS | tx | with_cond | 0.405 | **0.709** | 0.494 | NO |
| EC→PA | tl | no_cond | 0.486 | **0.703** | 0.371 | NO |
| EC→PA | tl | with_cond | 0.447 | **0.703** | 0.344 | NO |
| BS→EC | tx | no_cond | 0.299 | **0.413** | 0.249 | NO |
| BS→EC | tx | with_cond | 0.320 | **0.413** | 0.234 | NO |
| BS→PA | tx | with_cond | 0.185 | **0.276** | 0.077 | NO |
| PA→EC | tx | with_cond | 0.552 | **0.650** | 0.521 | NO |
| **PA→EC** | **tl** | **with_cond** | **0.341** | 0.142 | −0.015 | **YES** |
| PA→BS | tx | no_cond | 0.464 | **0.811** | 0.641 | NO |
| PA→BS | tx | with_cond | 0.413 | **0.811** | 0.565 | NO |

**In 10 of 11 cells, the reference-value-only baseline — no sequence at
all — scores AS HIGH OR HIGHER than the original sequence+reference model.**
The effect sizes (orig − ref-only) are large and consistently negative:
−0.09 to −0.40 in rho. **Every single *B. subtilis*-involving transcription
cell — the entire basis of Task 1B's headline claim — fails this control.**
Handing the model a sequence adds nothing (or actively hurts, likely by
giving the optimizer a noisier, higher-variance path to the same
reference-value-dominated answer) beyond what the reference value alone
already provides.

**Check 4 (shuffled-sequence) confirms this from a second, independent
angle.** Re-running the exact same architecture with sequences randomly
permuted relative to their true (reference, shift) pairing — so the
sequence pathway carries pure noise — produces performance close to (and in
several cells above) the ORIGINAL model's own performance (e.g. EC→BS
transcription: original 0.432, shuffled-sequence 0.557). **A model handed
uninformative sequences does about as well as the model handed real ones.**
This is only possible if the real sequences were not contributing much
information the reference value didn't already supply.

**The one survivor.** PA→EC translation, with-conditioning: original
ρ=0.341 [0.309,0.373] clearly exceeds reference-only ρ=0.142 [0.111,0.174]
(non-overlapping 90% CI) and clearly exceeds shuffled-sequence ρ=−0.015
[−0.034,0.003]. This is a real, controlled, sequence-specific signal — but
it is **translation, not transcription**, and **PA→EC, not a *B. subtilis*
pair**. It shares no host or readout with the finding Task 1B's headline
was built on. Reported honestly as a small positive footnote, not
generalized.

**Effect size, for the one survivor (per the task's Check 5):** the
sequence-carrying model explains rho=0.341 vs. the reference-only model's
0.142 — an absolute gain of +0.199 in Spearman rho, or in
variance-explained terms, `fve_mean` (fraction of shift-variance explained
beyond a constant mean-shift baseline) goes from a value indistinguishable
from zero for reference-only to 0.073 [0.026, 0.115] for the full model
(from `out/results/gate8_5_shift_prediction_summary.csv`) — a real but
modest effect, not a large one. Not built into anything further.

**VERDICT: Task 1B, as headlined in `out/GATE8_5_MEMO.md` ("shift-prediction
recovers *B. subtilis* transcription signal, 7/8 configurations"), does NOT
survive.** It was overwhelmingly a regression-to-the-mean artifact. **This
finding is retracted.** `out/GATE8_5_MEMO.md`'s Task 1B section should be
read historically, superseded by this memo — a correction notice is added
there (see below). The corrected picture across all four Gate 8.5 tasks:

| Task | Corrected verdict |
|---|---|
| 1A — alternative conditioning mechanisms | **Real, narrow positive: EC-transcription only. Stands.** |
| 1B — predict the shift, not the level | **RETRACTED. Regression to the mean, not sequence signal, for every *B. subtilis* cell. One unrelated cell (PA→EC translation) survives as a minor, separately-scoped finding.** |
| 1C — restrict to the co-active subset | Negative, hypothesis contradicted for *B. subtilis* pairs. Stands. |
| 1D — separate the two stages for transfer | Negative, one narrow classifier exception. Stands. |

---

## Correction added to `out/GATE8_5_MEMO.md`

The following blockquote was inserted at the top of that file's Task 1B
section (full text, not paraphrased):

> **CORRECTION — Gate 8.6, dated after this section was written.** The
> regression-to-the-mean control (Task 1 of Gate 8.6) found that this
> section's headline finding does not survive scrutiny. 10 of the 11
> "winning" cells below are explained, as well or better, by a
> reference-value-only baseline with no sequence input at all — including
> every *B. subtilis*-involving transcription cell this section's verdict
> was built on. **This section's conclusion is retracted.** See
> `out/GATE8_6_MEMO.md` Task 1 for the full control and the one narrow
> exception that does survive (PA→EC translation only, unrelated to the
> claim below). Left in place, uncorrected in its own body text, for the
> historical record of what was found and how it was later caught — see
> this project's own standing practice of disclosing corrections in place
> rather than silently rewriting (`out/GATE7_MEMO.md`, `out/GATE5_5_MEMO.md`
> dated corrections).

---

## TASK 2 — characterizing the two secondary findings [COMPLETE]

### 2a. FiLM's fold instability, elevated and generalized

Script: `scripts/86_mechanism_fold_variance.py`. Results:
`out/results/gate8_6_mechanism_fold_variance.csv`,
`out/results/gate8_6_mechanism_fold_variance_summary.json`. Figure:
`out/figures/gate8_6_mechanism_fold_variance.png`.

**Fold-to-fold standard deviation of zero-shot Spearman rho, all 4 systems,
all 6 (host, readout) cells (from the ALREADY-COMPUTED Gate 4/8.5 LOHO
results — no retraining):**

| Host | Readout | sequence-only | **FiLM** | concat | per-host-heads |
|---|---|---|---|---|---|
| EC | transcription | 0.026 | **0.200** | 0.031 | 0.050 |
| EC | translation | 0.018 | **0.275** | 0.030 | 0.025 |
| BS | transcription | 0.023 | **0.152** | 0.026 | 0.027 |
| BS | translation | 0.048 | **0.056** | 0.045 | 0.038 |
| PA | transcription | 0.030 | **0.130** | 0.049 | 0.038 |
| PA | translation | 0.036 | **0.110** | 0.039 | 0.035 |

**FiLM has the highest fold-to-fold variance in every single cell, usually
by a wide margin** — at EC transcription specifically, FiLM's std (0.200)
is **6.4× concat's, 4.0× per-host-heads', and 7.8× sequence-only's.** This
is not an EC-specific artifact; it holds across all three primary hosts and
both readouts. **This is a methodological finding independent of Task 1A's
own positive result: a γ/β generator fit from two training-host feature
vectors is measurably, universally unstable in this project's few-domain
regime, regardless of which host is held out.** Useful to anyone building a
host-conditioned model on a handful of domains — FiLM's specific
parameterization (2 conv-layer insertion points, a small MLP generator) is
a bad default choice when the number of training domains is this small,
independent of whether conditioning helps at all.

**DISCLOSED SCOPE LIMIT — draw variance not computed.** Gate 3.5/Gate 4's
fold-vs-draw variance ratio required a repeated-draw dimension (Gate 4's
N=100 calibration curves fine-tuned the same base model 10 times per fold
at fixed data). Task 1A's zero-shot LOHO protocol trains exactly ONE model
per (host, fold) for each mechanism — there is no draw dimension in the
already-computed data. Computing a true draw-level (weight-init) variance
would require re-training every mechanism multiple times per fold at fixed
data, which was not part of Task 1A's protocol and was not additionally run
here — a disclosed scope limit, not a silent substitution. The fold-to-fold
std reported above is itself the primary quantity Gate 3.5's own headline
finding (SR5) was about, and is sufecient to support the claim being made
(mechanism-dependent instability, not draw-noise-dependent).

**1A's positive result, re-scoped precisely per this task's instruction:**
concatenation and per-host-heads beat sequence-only and FiLM **only at
*E. coli* transcription** — the one primary host with a measured (not
sensitivity-bounded) reliability estimate (0.912, Gate 8) and the largest,
cleanest usable-N of any primary host/readout combination. **Frame as
"conditioning can help where measurement is reliable and training data is
plentiful," not as a general property of these mechanisms.** *B. subtilis*
is unaffected by any of the three mechanisms tested — every BS cell across
concat/per-host-heads is statistically indistinguishable from sequence-only
(see `out/GATE8_5_MEMO.md` Task 1A's full 18-cell table).

### 2b. 1C's refuted hypothesis — direction, coherence, and the range-restriction check

**Reframed as a refuted hypothesis, direction named, per the task's
instruction:** the hypothesis under test was "activity is host-specific;
strength above threshold is conserved." The data says the **opposite** for
*B. subtilis* pairs: restricting to sequences co-active in both hosts
**roughly halves** the cross-host correlation (EC-BS transcription:
0.655→0.258; BS-PA transcription: 0.508→0.257). **Implication, stated
directly: the shared cross-host signal for *B. subtilis* pairs lives
predominantly in whether a sequence fires at all, not in how strongly it
fires once it does.** Once both hosts fire, magnitude agreement drops
sharply rather than strengthening — the opposite of "conserved strength
above threshold."

**Range-restriction check (new this gate).** Script:
`scripts/85_coactive_range_restriction_check.py`. Results:
`out/results/gate8_6_coactive_range_restriction_check.csv`. Range
restriction (a classical mechanical cause of lower observed correlation,
independent of any real relationship change) requires the restricted
subsample to have a **narrower** distribution than the pooled one. Checked
directly: for EC-BS and BS-PA transcription (the pairs that halved), the
co-active-restricted subset's IQR is **5×–53× WIDER**, not narrower, than
the pooled subset's (e.g. *B. subtilis*'s own tx_norm IQR: 0.058 pooled →
3.082 co-active for EC-BS, a 53× widening — because the pooled set is
dominated by a dense cluster of near-zero "usable but inactive" values that
active-restriction removes). **Range restriction is ruled out as the
explanation — if anything, a wider range should preserve or increase
observed correlation, making the actual drop more notable, not less.** N
does shrink substantially (14,088→3,668 for EC-BS transcription), but the
point-estimate shift (0.40 in rho) is roughly 20× the standard error at
that N — far beyond what sampling noise from the smaller N alone would
produce.

**Coherence with Task 1's finding on 1B, revised:** the task asked whether
1B (shift-prediction working for transcription) and 1C (co-activity
carrying the signal) point toward the same "threshold behavior is
conserved" conclusion. **Given Task 1B's *B. subtilis*-transcription
finding did not survive the regression-to-the-mean control, this specific
coherence question is now largely moot** — 1B no longer contributes an
independent transcription-specific data point to compare against 1C. What
remains: 1C's own finding (co-activity does not carry a conserved-magnitude
signal, but the DROP itself shows the shared information is concentrated
in the fire/no-fire decision, not the analog level) stands on its own,
unsupported and uncontradicted by 1B's now-retracted result. No inference
is drawn from a coincidence that no longer exists.

### 2c. Transcription/translation asymmetry — corrected to 3 converging lines, not 4

**One of the task's four proposed converging lines does not survive.**
Reassessed honestly:

1. **Sigma-70 motif signal exists for transcription, not translation**
   (Gate 3, `scripts/33`: ρ=0.24–0.43 transcription vs. |ρ|<0.02
   translation-ΔG-alone). Stands, independently verified, unaffected by
   this gate.
2. **The floor artifact hits translation hardest** (89.9% of *B. subtilis*
   translation rows pinned at a detection-floor pseudo-value, Gate 3).
   Stands.
3. **Usable *B. subtilis* translation N is an order of magnitude smaller
   than the raw column implies** (1,101 floor-corrected vs. 11,564 raw,
   `out/KNOWN_ISSUES.md` item 2). Stands.
4. ~~Shift-prediction works 7/8 for transcription, 0/8 for translation~~ —
   **DOES NOT SURVIVE.** The transcription side of this specific comparison
   was retracted by Task 1 above (regression-to-the-mean). **Removed from
   the convergent-evidence list.** Citing Task 1B's original framing here
   would be citing a result this same gate found does not hold.

**Restated as a result on 3 independent lines, not 4:** transcription
carries recoverable cross-host structure (a real biophysical feature and,
per Task 1A, recoverable conditioning signal at least at *E. coli*);
translation, in this dataset, shows neither — no biophysical signal, the
most severe floor artifact, and (for *B. subtilis* specifically) the
smallest real N of any primary host/readout combination.

**The honest caveat the task asked for, labeled as inference:** three
lines of evidence is enough to state the asymmetry as a *pattern in this
dataset*, but not enough to rule out "the FACS-seq translation readout is
too noisy to support this class of analysis" as a full explanation
independent of any real biological difference in cross-host
conservation. **I would bet on a mix, weighted toward the noise
explanation for *B. subtilis* specifically** (n=1,101 floor-corrected,
smallest of any primary host/readout cell in this project, `out/KNOWN_ISSUES.md`
item 2) but NOT for *E. coli*/*P. aeruginosa* (whose translation N is
9,146/17,630 respectively — plenty of data, and Gate 5.5's four-way
comparison still found no biophysical or conditioning signal there
either, which a pure-noise explanation does not by itself account for).
**This is inference, not a measured fact** — labeled as such, not
presented as established. **The gap that would settle it: no measurement
reliability estimate exists for translation in any host** (`out/KNOWN_ISSUES.md`
item 3, unchanged by this gate) — the same disattenuation approach Gate 8
used for transcription (EC reliability 0.912, from 5 growth-condition
replicates) has no translation-readout analog in the released data. Until
that gap closes, "translation is fundamentally less cross-host-conserved"
and "translation is too noisy to tell" remain empirically indistinguishable
for *B. subtilis* specifically.

### Task 2 (from the prior prompt) — literature reconciliation, reconfirmed unchanged

2a (Johns et al. 2018), 2b (PromoGen2), 2c (Bernstein-lab Stutzerimonas
scope), and 2d (reliability coverage) were fully verified against primary
sources in Gate 8.5 (`out/GATE8_5_MEMO.md` Task 2) and are unaffected by
this gate's findings — none of Gate 8.6's results touch those claims.
Re-read and confirmed still accurate; not re-verified from scratch since
nothing in this gate bears on them. Full detail remains in
`out/GATE8_5_MEMO.md`.

---

## TASK 3 — NOT RUN

Per the task's own contingency ("only proceed if 1B survives Task 1"), and
Task 1's verdict above (1B does not survive), **`PAPER_FRAMING.md` was not
touched.** The evidence hierarchy proposed in the Gate 8.6 prompt (led by
shift-prediction) is not applicable. The framing to carry forward is:

- **1A stands**: a narrow, real positive (EC-transcription, two mechanisms)
  that requires "conditioning fails under multiple mechanisms" to be
  corrected to "conditioning fails for *B. subtilis* under three
  mechanisms; helps only at *E. coli* transcription, where measurement
  reliability is highest and data is most plentiful" — this correction
  (already identified in Gate 8.5, unaffected by 1B's retraction) is still
  owed to `PAPER_FRAMING.md` and should be made in a future gate.
- **1B is negative** (retracted positive), contributing nothing to the
  evidence hierarchy beyond the methodological lesson (regression-to-the-
  mean is a serious, real confound for any future shift/delta-prediction
  work on this dataset — worth a `KNOWN_ISSUES.md` entry).
- **1C stands** as a refuted-hypothesis, direction-named negative finding,
  now with a range-restriction check ruling out the obvious artifact
  explanation.
- **1D stands** unchanged.

`PAPER_FRAMING.md`'s actual reframing (incorporating all of Gate 8.5's
still-valid findings — 2a/2b/2c/2d's literature corrections, 1A's
E.-coli-only scoping, 1C's refuted-hypothesis framing, 1D's stage-specific
note — none of which depended on 1B) remains outstanding and should be
the next gate's Task 1, done cleanly without the shift-prediction claim
in it.

---

## TEMPTATIONS TO ADJUST

1. **Reporting only the raw "beats zero-shift" number (21/24) and not the
   much stricter "beats mean-shift" number (11/24) for 1B**, which would
   have made the original finding look stronger than it was even before
   Task 1's control. Not done — both were reported in Gate 8.5, which is
   exactly what made this gate's control possible to specify precisely.
2. **Treating the one surviving cell (PA→EC translation) as partial
   vindication of the *B. subtilis*-transcription headline.** Not done —
   it shares no host or readout with the original claim and is reported as
   an unrelated, minor, separately-scoped footnote.
3. **Softening "RETRACTED" to "weakened" or "requires more evidence."**
   Not done — 10 of 11 cells is not an ambiguous result; the word used is
   RETRACTED, stated plainly, first line of this memo.
4. **Quietly editing `out/GATE8_5_MEMO.md`'s Task 1B section to remove the
   original (wrong) claim**, rather than leaving it with a dated correction
   notice. Not done — matches this project's established practice
   (Gate 5.5, Gate 7) of correcting in place with the error visible, not
   erasing it.
5. **Running Task 3's reframing anyway, on the theory that 1A's positive
   result alone justifies touching the framing this gate.** Not done — the
   task's contingency was specific to Task 1B's outcome and explicit
   ("only proceed if 1B survives"); 1A's status was already known before
   this gate started and doesn't retroactively authorize Task 3 here.
