# CROSSHOST: A Benchmark for Predicting Bacterial Regulatory DNA Activity Across Host Species, and Why Host Conditioning Does Not Help

**Status: preprint-ready draft. Not submitted anywhere.**

---

## Abstract

A regulatory DNA sequence characterized in one bacterial species frequently
behaves differently in another — the "chassis effect." Multiple groups have
published this as a major open problem: chassis effects hinder prediction
of circuit function from part composition, cause costly design-build-test
repetition, and discourage engineering of non-model organisms, and a
predictive understanding of which biological properties underpin them
remains an unanswered gap. One published hypothesis holds that host
physiology, not genome sequence, explains chassis effects, and that genome-encoded
features alone cannot. We built a benchmark to test this directly: three
bacterial hosts with dense regulatory-activity coverage (*Escherichia coli*,
*Bacillus subtilis*, *Pseudomonas aeruginosa*; **n_hosts ≤ 6** including
three additional hosts at greatly reduced N), frozen genome-blocked splits,
a pre-registered hypothesis test, and a small (~214K-parameter)
FiLM-conditioned CNN trained to predict transcription and translation
activity from sequence plus host features. **The pre-registered kill-gate
hypothesis failed on all 8 primary comparisons** — host-conditioned
prediction at N=100 calibration examples did not beat a simple per-host
baseline trained on N=3,000 examples, not even on the closest call
(model ρ=0.213 vs. baseline ρ=0.218). We then tested whether this was a
failure of the specific conditioning mechanism, the host-feature
representation, the model's capacity, or the prediction target — three
conditioning mechanisms, two host-feature vectors, two foundation models,
and two target reframings, systematically. Genomic and physiological host
features never distinguishably beat a free identity tag. Alternative
conditioning mechanisms helped in exactly one of 18 tested cells — the one
host/readout pair with the highest measurement reliability and largest
sample size — and left *B. subtilis*, the host the central question depends
on, unaffected under every mechanism tested. A promising-looking
"shift-prediction" reframe recovering *B. subtilis* signal was retracted
after a control found it was regression to the mean, not sequence
information. What does explain the pattern: raw cross-host measurement
agreement itself — a replication, at larger scale and with a
measurement-error correction, of a result Johns et al. (2018) already
reported — bounds what any model can achieve, and it varies sharply with
host phylogenetic distance, not with any model's sophistication. The
sequence-only model is also poorly calibrated off-distribution (Expected
Calibration Error 4–9× worse on non-*E. coli* hosts), a caveat that applies
regardless of the negative conditioning result. We release the benchmark,
splits, evaluation code, and full baseline suite — including this negative
result and both documented failed rescue attempts — as a public,
versioned resource.

---

## 1. Introduction

Synthetic biology relies on characterized genetic parts to build predictable
circuits. A large fraction of those parts are characterized in one or a
handful of model organisms — overwhelmingly *E. coli* — and their behavior
in other hosts is often assumed rather than measured. This assumption
regularly fails. The literature calls the resulting mismatch the **chassis
effect** or **context dependence**, and describes it in terms that leave
little ambiguity about its cost to the field: it *"hinders prediction of
circuit function from part composition alone, causes costly repetitions of
the design-build-test cycle, and discourages use of non-model
organisms,"* and *"a predictive understanding of which biological
properties underpin chassis effects"* is described as *"a major knowledge
gap left unanswered."* Identifying which genome-encoded determinants
govern chassis effects is stated as remaining a major open challenge. This
is voiced demand from practitioners in the primary literature, not a
problem this project invented to have something to work on.

Part of the field routes around the question rather than answering it.
Orthogonal-machinery approaches — T7 RNA-polymerase-based "universal"
expression modules and similar systems — partially neutralize host
variation for the subset of users willing to adopt an orthogonal
transcription/translation system. This is a real, useful mitigation, but it
caps relevance rather than resolving the underlying predictive question: it
helps users who can redesign their circuit around orthogonal machinery, not
users who need to predict how an *existing*, non-orthogonal part will
behave in a new host.

One direct, falsifiable hypothesis about the chassis effect's cause comes
from work on engineered genetic circuits across related *Pseudomonas* and
*Halopseudomonas* hosts and, in follow-up work, across six closely related
*Stutzerimonas* strains. That work reports that *"hosts exhibiting more
similar metrics of growth and molecular physiology also exhibit more
similar performance of the genetic inverter, indicating that specific
bacterial physiology underpins measurable chassis effects"* (Chan, Baldwin
& Bernstein, 2023), and separately argues that *"it is impractical to
postulate that the observable chassis-effect between a given set of hosts
can be explained by a single or even a set of predictable genome-encoded
functions without experimental insight"* (Chan & Bernstein, 2024) — i.e.,
that host *physiology*, measured directly, explains chassis effects, and
that genome sequence alone should not be expected to.

This is a testable claim, and testing it does not require a new wet-lab
program: Johns et al. (2018) already generated the raw material. Their
study measured transcription and translation activity for tens of
thousands of regulatory sequences mined from 184 bacterial genomes, across
three primary recipient hosts (*E. coli*, *B. subtilis*, *P. aeruginosa*)
and, for a 241-sequence subset, three additional hosts. This is, to our
knowledge, the only publicly available dataset that measures the *same*
regulatory sequences' activity across multiple bacterial host species at
scale. Existing genomic-sequence benchmarks (BEND, Genomic Benchmarks,
DART-Eval, DNALONGBENCH, the Nucleotide Transformer suite among them) are,
to our knowledge, predominantly built on human, animal, or plant genomes;
**we have not independently re-verified the full task composition of each
of these suites against their own primary sources** and state this as our
working understanding, not a checked fact — a gap disclosed here rather
than written around (Section 6). What we can state with confidence, having
searched directly: no existing benchmark suite poses a cross-*host*
(same sequence, multiple bacterial recipient species) prediction task,
which is the specific gap this benchmark fills regardless of how the
broader landscape claim above resolves.

We built CROSSHOST: a frozen, versioned benchmark on this dataset, with a
small demonstration model used to test whether genome-encoded or
proteome-derived physiological host features let a model transfer
regulatory-activity prediction to a bacterial host it did not train on. We
pre-registered the central hypothesis before training the model. It failed.
We then spent four independent gates trying to determine whether that
failure was about the conditioning mechanism, the host representation, model
capacity, or the prediction target — because a negative result that has
survived several genuine attempts to overturn it is a stronger contribution
than one that was accepted at face value. Two of those attempts initially
looked positive; one survived scrutiny in a single narrow cell, the other
did not survive a control we ran specifically because it was the single
most likely thing to invalidate it. We report both, including how the
second was caught, because a negative result whose own author tried hardest
to break it is more trustworthy than one presented without a fight, and
because a documented failed rescue attempt is itself useful information to
the next person who has this same idea.

**n_hosts ≤ 6.** Every claim in this paper is bounded by that number, and
the abstract states it explicitly for a reason: the unit of generalization
in this study is the host, not the sequence, and six is not many hosts.

---

## 2. Methods

### 2.1 Data

Johns et al. (2018), *Nature Methods* 15:323–329 (BioProject PRJNA431139).
A library of regulatory sequences (RSs), 165 bp each, mined from 184
bacterial genomes, transformed into three primary recipient hosts (*E.
coli*, *B. subtilis*, *P. aeruginosa*) and measured for transcription
(RNA-seq-derived, normalized DNA/RNA read ratio) and translation (FACS-seq
GFP fluorescence) activity. The released core library totals 29,249 RS/oligo
records; per-host usable counts (a measurement exists) are *E. coli*
24,613 / *B. subtilis* 15,848 / *P. aeruginosa* 21,473 for transcription.
**Two different, both-correct active-fraction denominators are in use in
this literature and must be distinguished explicitly, as recommended by
our own Gate 1 verification:** computed over each host's own usable
population (the denominator this project uses operationally throughout,
since per-host baselines and folds are host-specific), active fractions
are *E. coli* 61.1%, *B. subtilis* 28.0%, *P. aeruginosa* 83.3%. Computed
instead over the three-host "union set" (RSs with usable transcription
data in all three hosts simultaneously, n=11,265–11,276) — the denominator
Johns et al. use in their own headline statistic — active fractions are
*E. coli* 52.0%, *B. subtilis* 18.9% (Johns et al.'s own reported figure,
independently reproduced at 17.9% by this project on the same
union-set population), *P. aeruginosa* 83.8%. **The direction (*B.
subtilis* lowest by a wide margin, *P. aeruginosa* highest) is identical
under either denominator; the absolute percentages are not directly
comparable across denominators and should not be quoted without stating
which is meant.** A
241-sequence extension panel ("RS241") adds three further hosts at reduced
N (*S. enterica*, *V. natriegens*, *C. glutamicum*); 34 of the 241 nominal
RS241 oligo IDs have no recoverable sequence text in the released tables,
leaving 207 usable.

**Translation floor artifact.** `protein_log10` (the translation-strength
regression target) is pinned at a per-host floor/pseudo-value for a large
fraction of nominally usable rows — 66.5% (*E. coli*), 89.9% (*B.
subtilis*), 10.1% (*P. aeruginosa*) — almost certainly Johns et al.'s own
below-detection-limit convention, not real quantitative measurements. We
exclude floor-pinned rows from the regression target throughout;
floor-corrected usable translation-regression N is 9,146 (EC) / 1,101 (BS)
/ 17,630 (PA) — an order of magnitude smaller than the raw column implies
for *B. subtilis* specifically.

### 2.2 Host features

Two independent host-representation vectors, built to test the
genome-vs-physiology question directly:

- **Genomic (37-D):** sigma-factor complement, anti-Shine-Dalgarno
  hybridization free energy, tRNA gene copy number and codon adaptation
  index (tAI), codon usage, genomic GC content, and chaperone/heme/RNAP
  subunit-count annotations, unit-tested against published reference
  values (e.g. the *E. coli* anti-SD sequence, published tAI ranks). Fit
  once per fold on training-host data only where any statistic is derived.
- **Physiology (6-D):** ribosomal-protein fraction, RNAP core fraction,
  sigma-factor fraction, chaperone fraction, elongation-factor fraction,
  and growth rate, derived from public reference proteomes (primarily
  PaxDb) depth-matched across hosts to correct for a confirmed 7-fold
  variation in native proteome-coverage depth across the six hosts. This
  is a **coarse, reference-database-derived proxy** for physiology, not a
  direct wet-lab measurement — a distinction load-bearing for how this
  project's results relate to prior physiology-based claims (Discussion).

Both vectors are z-scored across all six hosts (not fold-locally): host
identity/covariate information is treated as known at deployment time for
a genuinely new host, unlike sequence-level label statistics, which are
always fit fold-locally.

### 2.3 Architecture

A four-layer 1D convolutional trunk (kernel widths 15/9/5/3, channels
128/128/64/64, matching the verified 165bp part length), BatchNorm, a
residual connection wrapping the final conv layer, global average pooling,
and two jointly-trained two-stage output heads (transcription,
translation; each an active/inactive classifier plus a strength regressor
on actives) — 213,956 parameters with no conditioning pathway
(`SequenceOnlyCNN`, the primary ablation baseline throughout this paper).
The primary conditioning mechanism, pre-registered and used for H-MAIN, is
**FiLM** (Feature-wise Linear Modulation): a small MLP maps the host
feature vector to per-channel (γ, β) applied multiplicatively at
convolutional layers 2 and 3 (~227K genomic-variant / ~226K
physiology-variant parameters). Two further conditioning mechanisms —
**concatenation** (host vector appended to the pooled representation before
the final layer) and **per-host output heads** (a private head set per
training host, combined at inference by prediction-space averaging) — were
added later specifically to test whether FiLM's own weakness, not
conditioning per se, explained the negative result (Section 3.3).

### 2.4 Splitting scheme

Sequences are deduplicated to unique regulatory-sequence text before
splitting. The primary leakage defense is sequence-homology clustering
(MMseqs2, `--min-seq-id 0.5 -c 0.8`), genome-blocked, with an explicit
exact- and near-duplicate safety net (a k-mer-index-based all-pairs check,
force-merging any connected component spanning multiple folds) added after
clustering alone was found to miss byte-identical duplicate pairs at every
tested identity threshold. **Maximum train-test sequence identity across
all 5 folds: 0.8485** (fold-by-fold: 0.8485, 0.8395, 0.8485, 0.8485,
0.8364) — reported as a number, not a claim that clustering was applied,
because approximate clustering tools carry no formal completeness
guarantee and this dataset's 184 source genomes include phylogenetically
close pairs where that gap matters. RS241 is never trained on, enforced in
code (an assertion, checked by a standing leakage audit run before and
after every gate) and independently verified: 0 of 241 RS241 IDs appear in
any main-library fold. Within-host quantile normalization of sequence-level
statistics is fit on training-fold data only, inside the fold loop, verified
by the same standing audit.

### 2.5 Pre-registration

`PREREGISTRATION.md`, committed 2026-08-04, before any host-conditioned
model training began, and frozen at commit time. Three hypotheses:

- **H-MAIN** (the kill gate, dual-arm per readout, respecified before Gate
  4 training on data-availability grounds — transcription and translation
  were found mechanistically distinct in Gate 3.5): the cross-host model
  at N=100 host-specific calibration examples matches or beats a
  per-host-only baseline at N=3,000 (transcription) or its own measured
  saturation ceiling, N≈300 (translation), on *B. subtilis* held out.
- **H-SCIENCE**: genomic and physiology host features differ measurably in
  leave-one-host-out performance, direction not predicted in advance —
  either outcome (physiology wins, supporting the Bernstein-lab hypothesis;
  genomics wins, challenging it) is a publishable result under the
  pre-registration.
- **H-DIAGNOSTIC**: the host-biology representation outperforms a free
  per-host lookup embedding (identity only, zero biological content) under
  a pre-specified fallback (mean of training-host embeddings) — evidence
  the model learned host biology rather than merely host identity.

**Amendment 1** (2026-08-05, before H-MAIN evaluated): report both the
pre-registered frozen-trunk transfer mechanism and a supplementary
unfrozen-conv4 mechanism at N=100, after Gate 4's own validation work found
the frozen-trunk mechanism interacts with a small-N degradation. The
frozen-trunk mechanism remained primary; the supplementary mechanism was
never substituted into the headline result.

**Amendment 2** (2026-08-05, before H-MAIN evaluated): widen the decision
interval from an 80% Student's-t interval to a **90% percentile bootstrap**
(10,000 resamples, hierarchical: fold identities resampled, then
within-fold draws resampled), after measuring the model's fold-to-fold
variance exceeding within-fold draw-to-draw variance by a factor
inconsistent with the t-interval's symmetric-unimodal assumption. This
change is strictly conservative — it makes MET *harder* to declare, not
easier — and was committed before any H-MAIN comparison was computed under
either rule.

### 2.6 Decision rule

H-MAIN is MET only if the model's 90% bootstrap CI lower bound exceeds the
baseline's 90% CI upper bound; overlapping intervals, regardless of which
mean is nominally higher, are NOT MET. This rule, and the two amendments
above, were fixed before Gate 5 evaluation and were not revisited after
seeing results.

---

## 3. Results

### 3.1 The pre-registered kill gate: H-MAIN failed on all 8 primary comparisons

Neither readout, neither host-feature variant, on either primary held-out
host (*B. subtilis*, the official test; *P. aeruginosa*, the charter's
partial-pass check), met the pre-registered bar. The closest result in the
entire gate was *B. subtilis* transcription with genomic features: model
ρ=0.213 [0.154, 0.260] vs. baseline ρ=0.218 [0.160, 0.251] — an overlap of
0.091, "NOT MET, point estimate favors baseline," and this cell was not
revisited under any alternative mechanism for the headline result. Every
*P. aeruginosa* and *E. coli* primary-mechanism cell also failed, several
by wide margins (e.g. *P. aeruginosa* genomic transcription: model 0.104
vs. baseline 0.447). A supplementary (non-pre-registered) unfrozen-conv4
mechanism showed several *P. aeruginosa*/*E. coli* cells with the model's
point estimate slightly exceeding the baseline's, but with overlapping
intervals (NOT MET) and, per the pre-registration, this mechanism was never
eligible to supply the headline verdict.

### 3.2 The model is poorly calibrated off-distribution

The sequence-only model's zero-shot active/inactive classifier is
well-calibrated on *E. coli* (Expected Calibration Error 0.05–0.10) and
severely miscalibrated on *B. subtilis* and *P. aeruginosa* (ECE
0.34–0.47) — a 4–9× degradation. Split-conformal interval coverage for the
strength regressor stayed close to nominal (77–94% empirical against
80/90% targets) regardless of this classifier miscalibration, since
conformal intervals guarantee marginal coverage independent of the
underlying model's calibration — the ECE failure is a point-probability
problem specifically, not something interval-coverage numbers alone would
reveal. **Practical implication: a model that wins on rank correlation is
not automatically trustworthy in absolute probability terms cross-host**,
and any deployment of a model from this suite on a host meaningfully
different from its training hosts should recalibrate before trusting raw
predicted probabilities.

### 3.3 Testing whether the mechanism was the problem: three conditioning mechanisms

To address the objection that FiLM specifically — fit from only two
training-host feature vectors per leave-one-host-out fold — is a
known-weak choice in a low-data regime, we implemented two architecturally
distinct alternatives (concatenation; per-host output heads, combined by
prediction-space averaging) under an identical protocol to FiLM's own
zero-shot evaluation (same 5 folds, same seeds, same trunk, same training
budget, same bootstrap procedure), genomic host-feature variant.

**Result: 2 of 18 (host, readout, mechanism) comparisons show a mechanism
beating both sequence-only and FiLM distinguishably, and both are *E. coli*
transcription.** Concatenation (ρ=0.555) and per-host-heads-averaging
(ρ=0.615) both clear sequence-only (0.367) and FiLM (0.371) with
non-overlapping 90% intervals; per-fold values are tight and consistent
(concatenation: 0.514–0.606 across 5 folds), unlike FiLM's own highly
unstable per-fold trajectory on the identical folds (0.045–0.606). **Every
*B. subtilis* cell, both readouts, is statistically indistinguishable
across all three mechanisms** — no mechanism recovers signal for the host
the central question is actually about. One asymmetric failure mode is
worth reporting rather than smoothing over: per-host-heads-averaging at
*P. aeruginosa* transcription scores 0.309, below *both* sequence-only
(0.479) and FiLM (0.338) — the worst of five systems in that one cell,
demonstrating that averaging private per-host heads is not a universally
safe combination rule.

A secondary combination rule tested for per-host heads — selecting only
the single *nearest* training host's head by genomic-feature Euclidean
distance, rather than averaging both — performed markedly worse for *E.
coli* (ρ=0.233, below every baseline). The nearest training host for *E.
coli* by this metric turned out to be *B. subtilis* (distance 5.42) rather
than the phylogenetically closer *P. aeruginosa* (distance 8.26) —
feature-space proximity in this 37-D genomic vector does not track
phylogeny, and averaging both training hosts' heads is the empirically
safer default.

**FiLM's own instability is a separate, generalizable finding.** FiLM's
fold-to-fold standard deviation of zero-shot Spearman ρ exceeds every
alternative mechanism's — sequence-only, concatenation, per-host heads —
in all 6 (host, readout) cells tested, not only where a positive result
was found: 4.0–7.8× higher at *E. coli* transcription specifically (FiLM
std=0.200 vs. concatenation 0.031, per-host-heads 0.050, sequence-only
0.026). A γ/β generator fit from two training-host feature vectors is a
measurably poor default in this few-domain regime, independent of whether
conditioning helps at all — useful to anyone conditioning a model on a
handful of domains, regardless of this paper's negative result.

### 3.4 Testing whether the host representation was the problem: a free identity tag

H-DIAGNOSTIC: a free per-host lookup embedding carrying zero biological
content matched or beat both host-biology feature vectors at zero-shot in
all 12 tested (host, readout, variant) cells — 0 clean wins for the
biological features, 10 cells where the free embedding's point estimate
won outright. **Neither the 37-D genomic nor the 6-D physiology vector
carries detectable content beyond an arbitrary host tag**, at this n_hosts.

### 3.5 Testing whether the readout mattered: co-activity, not magnitude, carries what signal exists

Restricting cross-host measurement correlation to sequences active in
*both* hosts of a pair — testing the hypothesis that activity is
host-specific while strength above threshold is conserved — produces the
**opposite** result for *B. subtilis* pairs: correlation roughly halves
(*E. coli*–*B. subtilis* transcription: 0.655→0.258; *B. subtilis*–*P.
aeruginosa*: 0.508→0.257). Range restriction, a classical mechanical cause
of lower observed correlation independent of any real effect, was checked
directly and ruled out: the restricted subset's interquartile range is
5×–53× *wider*, not narrower, than the pooled set's — if anything this
should preserve or inflate the correlation, making the actual drop more
notable. The one pair where restriction helps (*E. coli*–*P. aeruginosa*:
0.621→0.754) is the pair with the highest correlation to begin with, the
one place the hypothesis was least needed. **What shared cross-host signal
exists for *B. subtilis* pairs appears concentrated in whether a sequence
fires at all, not in how strongly it fires once it does.**

Separating the model's own classification stage (does it fire) from its
regression stage (how strongly) on the same already-trained checkpoints
found a comparable pattern: FiLM beats sequence-only distinguishably on
the classification (AUC) metric in 1 of 12 cells, and on the regression
(Spearman ρ) metric in 0 of 12 — the strength regressor is where the
negative result is uniform; the classifier has one narrow, isolated
exception.

### 3.6 A retracted rescue attempt: shift-prediction and the regression-to-the-mean control

A fourth attempt reframed the prediction target: instead of absolute
activity level, predict the *shift* from a reference host's measured value
to a target host, given sequence and the reference value as inputs. This
initially appeared to succeed specifically for *B. subtilis* transcription
— 7 of 8 tested configurations distinguishably beat a mean-shift constant
baseline, a result no absolute-level framing in this project had produced
for that host.

**It did not survive a control we ran specifically because it was the
single most likely thing to invalidate it.** Spearman correlation between
the shift and the reference value itself is strongly negative for most
host pairs (as low as −0.812, *P. aeruginosa*→*B. subtilis* transcription)
— textbook regression to the mean: a sequence with a high reference-host
value mechanically has more room to fall than to rise. A reference-value-only
baseline (ordinary least squares, one feature, no sequence input at all)
matched or beat the original sequence-plus-reference model in **10 of 11**
cells originally reported as wins, with effect sizes of −0.09 to −0.40 in
Spearman ρ. A second, independent control — retraining the identical
architecture with sequences randomly permuted relative to their targets,
reference values kept intact — reproduced most of the original model's
performance despite the sequence input now carrying zero real information
(*E. coli*→*B. subtilis* transcription: original ρ=0.432, shuffled-sequence
ρ=0.557, *higher* with no sequence signal at all). **We retract this
finding.** The one cell that does survive both controls (*P.
aeruginosa*→*E. coli* translation, with conditioning: original ρ=0.341 vs.
reference-only ρ=0.142 vs. shuffled ρ=−0.015, non-overlapping) shares no
host or readout with the retracted claim and is reported as an unrelated,
minor, separately-scoped result, not partial vindication.

### 3.7 What does explain the pattern: raw cross-host measurement agreement

*E. coli* and *P. aeruginosa* activity measurements correlate at ρ≈0.75
(transcription, n=9,741 co-active pairs); any *B. subtilis* pair correlates
at only ρ≈0.16–0.26. **This is a replication of Johns et al.'s own reported
result** — their Supplementary Fig. S13 (n=212, three primary hosts) and
Fig. S15 (n=241, six-host RS241 panel) report exactly this pattern, and
their main text states plainly that *"between recipients, only E. coli and
P. aeruginosa showed significant correlations."* We verified this directly
against their primary-source text and figure captions rather than assuming
a second-hand summary. This project's addition is not the raw contrast
itself but two methodological extensions: a different, more directly
practitioner-relevant translation operationalization (raw protein level
rather than Johns et al.'s translation-efficiency ratio) on substantially
larger per-pair N via a pairwise rather than three-way-intersected
restriction, and a disattenuation analysis absent from Johns et al.
entirely. Anchored by an empirically-measured *E. coli* transcription
measurement reliability of 0.912 (five independent growth-condition
replicates, cross-checked to 0.929 by an independent method — the only
host×readout combination in the released data with usable replicate
structure), a classical measurement-error correction leaves the
*B. subtilis*-vs-*E. coli*/*P. aeruginosa* gap ratio **unchanged**: 0.341
observed, 0.341 corrected at reliability 0.9, and only 0.515 even under a
deliberately pessimistic reliability of 0.5. **This gap is not a
measurement-noise artifact at any tested reliability level, and it directly
bounds what any model — however conditioned — could achieve on this
benchmark.**

### 3.8 Appendix-level results: foundation models and feature ablation

Two genomic foundation models (DNABERT-2, 117M parameters; PromoGen2, 148M
parameters), evaluated via frozen-embedding-plus-shallow-head, uniformly
with each other, lost to the 213,956-parameter sequence-only model in 35 of
37 statistically distinguishable comparisons across primary hosts and
RS241. This comparison protocol, however, did not give either model its
best shot: PromoGen2's own published native zero-shot protocol (direct
likelihood scoring, the identical Johns et al. dataset) reports
0.68/0.52/0.30 (*E. coli*/*P. aeruginosa*/*B. subtilis* transcription) —
higher than both this project's own measured PromoGen2-embedding numbers
(0.495/0.490/0.243) and sequence-only itself, at every host. We demote this
comparison to an appendix-level result and drop any "the model was too
small" capacity-limit framing — the issue is protocol fairness, not model
size. A five-group ablation of the genomic feature vector (sigma-factor,
anti-Shine-Dalgarno, tAI/codon usage, RNAP subunit count, chaperone/heme)
found no group's removal changing rank correlation beyond fold-to-fold
noise (max |Δ|=0.089) — scoped strictly to this project's own 6-feature,
reference-genome-derived setup, and not framed as testing the Bernstein-lab
physiology hypothesis, since the physiology vector itself was never
ablated and is a coarse database-derived proxy, categorically different
from either Bernstein-lab paper's direct wet-lab measurement.

---

## 4. Discussion

**A negative result that survived four independent, controlled attempts to
overturn it is harder to dismiss than one accepted at face value.** Three
conditioning mechanisms, a target reframing, and a subset-restriction
hypothesis were tried specifically to find the conditions under which host
information helps cross-host prediction on this benchmark. Two initially
looked promising. One (alternative conditioning mechanisms) survives, but
only in the single cell with the best data and the only measured
reliability estimate — a scope-narrowing result, not a general win. The
other (shift-prediction) was retracted after a control designed
specifically to catch its most likely failure mode caught exactly that
failure mode. What remains standing across every test is the raw
cross-host measurement-agreement pattern (Section 3.7), which requires no
model at all and which every model's own performance tracks.

**For host-conditioned modeling in general:** this project's own
methodological finding (Section 3.3) — that FiLM's fold-to-fold instability
is 4–8× worse than two simple alternatives in a few-training-domain
regime, independent of whether conditioning helps at all — is a caveat
worth carrying into any future host-conditioned model built on a handful of
domains, this dataset or otherwise. A negative result obtained with an
unstable conditioning mechanism is a weaker negative result than one
obtained after checking that the mechanism itself was not the confound;
this project checked, and the instability, once quantified, retroactively
strengthens rather than undermines the original finding, since it shows the
negative result was not an artifact of an unexamined unstable default.

**The measurement-agreement bound is, we think, the paper's most durable
contribution**, independent of any model's architecture: cross-host
prediction on this dataset is bounded below by how well the hosts' own
measurements agree with each other, and that bound (i) is directly
computable without training anything, (ii) survives a measurement-noise
correction at every tested reliability level, and (iii) tracks host
phylogenetic distance closely enough that the one place any conditioning
mechanism helped (*E. coli* transcription) and the one place restricting to
co-active sequences helped (the *E. coli*–*P. aeruginosa* pair) are the
same underlying phenomenon — the Gammaproteobacteria pair, the closest
relatives in this project's host panel, and not coincidentally also the
pair with the only measured (not sensitivity-bounded) reliability estimate.

**What this does and does not say about the Bernstein-lab physiology
hypothesis:** nothing decisive, in either direction. This project's
physiology feature vector is a coarse, six-dimension, reference-proteome-derived
proxy, not the direct wet-lab measurement (growth curves, transcriptomes)
either Bernstein-lab paper used, and the ablation that found no genomic
feature group carrying signal never touched the physiology vector at all.
What is worth stating precisely: the regime where the 2024 Chan & Bernstein
paper found physiology predictive — six closely related *Stutzerimonas*
strains — is exactly the kind of phylogenetically close comparison this
project's own data independently shows has the highest cross-host
measurement agreement (Section 3.7) and the only regime where any
conditioning mechanism in this project showed a benefit (Section 3.3).
These are not competing claims tested in the same regime; they are
compatible claims tested in non-overlapping regimes of host similarity,
and the most useful next step is a physiology-vs-genomics comparison
conducted with the Bernstein lab's own direct measurement methodology
extended to more phylogenetically distant hosts — precisely what neither
group has yet done and what this project cannot do without a wet lab.

**What new data would settle the open questions:** (1) a direct
measurement-reliability estimate for *B. subtilis* and *P. aeruginosa*
specifically, and for translation in any host — currently bounded only by
sensitivity analysis (Section 3.7); (2) a physiology feature vector built
from direct per-host measurement rather than a reference-database proxy,
tested on this project's own three-host panel; (3) a genuinely large-scale
genomic foundation model (Evo 2 or a successor) evaluated on this
benchmark under its own native protocol — not attempted here for a
documented infrastructure reason (Appendix); (4) new multi-host regulatory
MPRA data on hosts phylogenetically intermediate between *E. coli*/*P.
aeruginosa* and *B. subtilis*, which would let a future study distinguish
"conditioning helps only in the Gammaproteobacteria-close regime" from
"conditioning helps only where N and reliability happen to be highest" —
these two explanations are confounded in the current three-primary-host
data and cannot be separated with the data this project has.

---

## 5. Limitations

*(Reproduced in full from `out/KNOWN_ISSUES.md`, the project's user-facing
limitations document, not abbreviated.)*

1. **n_hosts ≤ 6 — the generalization unit is the host, not the sequence.**
   Three primary hosts with dense coverage, three more at ~207 recoverable
   sequences via RS241. Any claim about "genome-encoded functions" or
   "bacterial regulatory prediction" as a general class is not supported by
   6 data points. Every interval in this project's own results bootstraps
   over folds/hosts, not sequences, for exactly this reason.
2. **The translation floor-value artifact.** `protein_log10` is pinned at a
   per-host floor/pseudo-value for a large fraction of nominally usable
   rows (66.5% EC, 89.9% BS, 10.1% PA). Corrected usable-for-regression N
   after excluding floor-pinned rows: EC 9,146 / BS 1,101 / PA 17,630.
3. ***B. subtilis* measurement noise, and what we could establish about
   it.** The EC–PA vs. *B. subtilis*-pair correlation gap survives a
   disattenuation correction essentially unchanged at the one
   empirically-grounded reliability level (0.91, EC transcription), and
   corrected *B. subtilis*-pair correlations never exceed ~53% of EC–PA's
   even under a deliberately pessimistic reliability of 0.5. What we could
   NOT establish: a direct reliability estimate for *B. subtilis* or *P.
   aeruginosa* specifically, or for translation in any host. **The one
   reliability number this project actually measured is for the host and
   readout the central claim depends on least; *B. subtilis*'s own
   reliability is bounded by a sensitivity grid, not measured** — a
   sensitivity grid is real evidence, but a different strength of evidence
   than a measurement.
4. **Fold variance exceeds draw variance, and grows with N.** Fold-to-fold
   standard deviation exceeded within-fold draw-to-draw standard deviation
   in 15 of 23 baseline combinations checked, growing sharply with N (e.g.
   *E. coli* transcription at N=3,000: 91× ratio). A single fixed test fold
   can look precise while being an unrepresentative sample of true
   host-level performance; use the full 5-fold rotation, not a single
   fold, when extending this benchmark.
5. **The batching cliff — a real ~300× performance bug, now fixed.**
   Unbatched inference over a large pooled array showed a severe
   non-linear cost cliff (measured: 4× the data at ~300× the cost) on both
   CPU and Apple Silicon MPS. Now fixed (chunked at 1,024 rows by default,
   numerically verified identical to the unbatched path); any new
   inference code against this benchmark should chunk large arrays.
6. **Evo 2 was not evaluated — a specific technical reason.** The official
   `evo2` package requires CUDA/Flash Attention/Transformer Engine on a
   Hopper GPU (verified against the ArcInstitute repository and GitHub
   issue #67); this project's environment is Apple Silicon, MPS-only. Two
   smaller genomic foundation models (DNABERT-2, PromoGen2) were evaluated
   instead and also did not beat sequence-only, which weakens but does not
   fully close the "the model was too small" objection.
7. **The retired ceiling metric.** An earlier "percent of cross-host
   measurement-correlation ceiling" metric went through two rounds of
   correction (a wrong "fix" of an originally-correct number, later itself
   corrected) before being retired outright, because the sequence-only
   model genuinely exceeds the (correctly-computed) ceiling for *B.
   subtilis* on both readouts — a real, explainable result (measurement
   noise attenuates the raw correlation in a way a model trained on
   thousands of examples is not), not a bug, but not a metric worth
   shipping either. Does not appear anywhere in this paper's baseline
   tables.
8. **Cross-host calibration failure.** See Section 3.2. ECE 0.05–0.10 (EC)
   vs. 0.34–0.47 (BS, PA) — a 4–9× degradation. Recalibrate before trusting
   raw predicted probabilities cross-host.
9. **Regression to the mean is a serious confound for any shift/delta-prediction
   extension of this dataset.** See Section 3.6. If you build a
   shift-prediction model on this or a similarly small-host dataset, run a
   reference-value-only baseline and a shuffled-sequence control before
   trusting the result — this project did not, initially, and the result
   was wrong.
10. **The held-out evaluation split is not cryptographically enforced.**
    The public model-development table still contains the withheld fold's
    real labels, since it is the same table used for development elsewhere
    in the project. Like most small academic benchmarks, the protocol
    relies on the convention being respected, not a technical barrier.

---

## 6. What I could not do

*(Kept as its own section, not folded into Limitations, per explicit
instruction — this is the honest account of blocked or partial work,
distinct from properties of the released data.)*

- **Evo 2 was not run**, for the documented infrastructure reason in
  Limitations item 6. A free hosted-API path exists but requires external
  account creation this project's tooling cannot perform without browser
  automation it does not have.
- **PromoGen2 and DNABERT-2 were evaluated via one protocol
  (frozen-embedding + shallow head), not each model's own best/native
  protocol.** This was a defensible choice for a fair head-to-head between
  two architecturally different foundation models, but it measurably cost
  PromoGen2 performance relative to its own published native zero-shot
  numbers, and we did not re-run the comparison under PromoGen2's native
  protocol specifically — this is disclosed as an open gap, not resolved.
- **No direct reliability estimate exists for *B. subtilis* or *P.
  aeruginosa* specifically, or for translation in any host** — no
  replicate/condition-series data exists for these in the released tables.
  Everything downstream of this gap is a sensitivity analysis, not a
  measurement, and is reported as such throughout.
- **We could not separate two explanations for the transcription/translation
  asymmetry** — that translation is fundamentally less cross-host-conserved,
  versus that the FACS-seq translation readout in this dataset is too
  noisy to support this class of analysis — and say so plainly in Section
  3.5/Discussion rather than picking one. We stated a bet (a mix, weighted
  toward noise for *B. subtilis* specifically, not for EC/PA) and labeled
  it as inference, not fact.
- **We did not re-verify the Johns et al. 2018 per-cell Pearson correlation
  values** from their Supplementary Figures S13/S15 directly — only the
  figure captions (via the supplementary text document) and the main-text
  summary sentence were available in this project's source materials; the
  images themselves were not. The qualitative claim is confirmed
  unambiguously from exact quoted text; the exact per-cell r values were
  not independently re-derived.
- **`make reproduce-full`** (the complete from-raw-data pipeline,
  encompassing all model training) **was not re-executed end-to-end** in
  this project's packaging pass — it would cost the cumulative 40+ hours
  of CNN/foundation-model training this project's gates already spent.
  Every individual script in the dependency graph has been run and
  produced its output at least once; the full chain was not re-run in one
  sitting from raw data alone. `make reproduce` (the fast path, regenerating
  figures/tables from already-computed intermediate results) was tested
  and verified working from a genuine clean git clone.
- **We did not independently re-verify the task composition of the five
  named comparison benchmarks in the Introduction** (BEND, Genomic
  Benchmarks, DART-Eval, DNALONGBENCH, the Nucleotide Transformer suite)
  against their own primary sources. The claim that they are "predominantly
  human, animal, or plant" is carried from this project's own planning-stage
  research (the charter), which itself flags several of its numbers as
  inherited from secondary sources — this specific claim was never put
  through this project's own primary-source verification pass the way the
  Johns et al., PromoGen2, and Bernstein-lab claims were. Stated as working
  understanding in the Introduction, not fact, and should be checked before
  submission if a reviewer would reasonably expect it verified.
- **We did not run a fifth attempt to recover a positive cross-host
  signal.** Four independent, pre-specified attempts (three conditioning
  mechanisms, a target reframing, a subset-restriction hypothesis) is the
  number we judged sufficient to call the negative result robust rather
  than under-tested; running further attempts after two already failed
  under scrutiny would risk exactly the anti-retrofit failure — reshaping
  the search until something sticks — this project's own standing rules
  warn against.

---

## References

**Verified directly against primary sources within this project** (exact
DOI/venue/date confirmed by reading the actual paper, not carried from
training-data memory — see `out/GATE8_5_MEMO.md` Task 2 for the
verification record):

- Johns, N.I., Gomes, A.L.C., Yim, S.S., et al. (2018). Metagenomic mining
  of regulatory elements enables programmable species-selective gene
  expression. *Nature Methods* 15, 323–329.
- Chan, K., Baldwin, G.S. & Bernstein, H.C. (2023). Revealing the
  Host-Dependent Nature of an Engineered Genetic Inverter in Concordance
  with Physiology. *BioDesign Research* 5, 0016. DOI 10.34133/bdr.0016.
- Chan, K. & Bernstein, H.C. (2024). Pangenomic landscapes shape
  performances of a synthetic genetic circuit across *Stutzerimonas*
  species. *mSystems* 9(9), e00849-24. DOI 10.1128/msystems.00849-24.
- Xia, Y. et al. (2026). Design prokaryotic cis-regulatory elements using
  language model. *Nucleic Acids Research* 54(4), gkag122. [PromoGen2;
  PMC12907563.]

**NOT independently verified within this project — carried from the
charter's planning-stage research or general knowledge, and flagged here
rather than presented with false confidence. Check against primary sources
before submission:**

- DNABERT-2 (Zhou, Z. et al., "DNABERT-2: Efficient Foundation Model and
  Benchmark for Multi-Species Genome") — citation details (venue, exact
  arXiv ID) not re-verified against the primary source in this project;
  the model itself was used directly via its HuggingFace weights
  (`zhihan1996/DNABERT-2-117M`), which is independently confirmed correct.
- The DART-Eval (NeurIPS 2024 D&B) citation and the ICLR 2025
  "Specialized Foundation Models Struggle to Beat Supervised Baselines"
  citation, both referenced in Section 1/Discussion as prior evidence that
  foundation models underperform supervised baselines on regulatory
  tasks — inherited from the charter's own text, not independently
  re-verified.
- BEND, Genomic Benchmarks, DNALONGBENCH, and the Nucleotide Transformer
  suite — named in Section 1 as comparison benchmarks; see Section 6 for
  the explicit disclosure that their task composition was not
  independently re-verified.

None of the above affects this paper's own experimental results, all of
which are independently computed and audited within this project; the gap
is confined to background/related-work citations.
