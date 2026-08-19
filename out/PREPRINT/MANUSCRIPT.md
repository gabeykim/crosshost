# Cell-free systems do not reproduce the bacterial chassis effect: host specificity resides in cellular context rather than transcription machinery

**Status: preprint-ready draft. Not submitted anywhere.**

**Title options considered** (Gate 12, Task 5 — presented for review, the first is used throughout this draft):

1. *Cell-free systems do not reproduce the bacterial chassis effect: host specificity resides in cellular context rather than transcription machinery* (used below)
2. *Host specificity in bacterial regulatory activity is contextual, not mechanistic: a cell-free comparison and a pre-registered test of genome-encoded features*
3. *The bacterial chassis effect vanishes in cell-free lysates: evidence that host specificity is cellular context, not transcriptional machinery*

---

## Abstract

A regulatory DNA sequence characterized in one bacterial species often
behaves differently in another — the chassis effect, a barrier
practitioners describe as unresolved for engineering non-model hosts.
Comparing cell-free (DRAFTS, Yim et al. 2019) against in-vivo (Johns et
al. 2018) measurements on the same 165bp library, cross-host correlations
are uniformly high in cell-free lysates
(**0.623–0.911** across 45 species pairs) but collapse in vivo for the
same species pair (**0.597 cell-free vs. 0.258 in vivo**, *E.
coli*–*B. subtilis*) — a comparison DRAFTS's own authors did not make.
**We read this, as a hypothesis consistent with the data rather than a
demonstrated causal claim, as evidence that host specificity resides
substantially in cellular context** (membrane, supercoiling, resource
competition, growth physiology) **rather than in transcriptional
machinery**, which cell-free lysates retain. This predicts host-descriptor
features encoding machinery composition should carry no signal; we tested
that directly: **a pre-registered hypothesis test failed on all 8 primary
comparisons**; an arbitrary host-identity tag matched or beat both a
37-feature genomic and a 6-feature physiology vector in all 12 cells; a
feature-group ablation found no group carrying signal beyond fold-to-fold
noise; and sequence-only was not distinguishably beaten in 16 of 18 tests
across three conditioning mechanisms. The same model is also poorly
calibrated off its training host (Expected Calibration Error 4–9× worse on
non-*E. coli* hosts). **n_hosts ≤ 6** throughout: three
primary hosts, three more at reduced N. We release the benchmark, splits,
evaluation code, and baseline suite, including the negative result and a
retracted rescue attempt.

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
behave in a new host. **We name our audience up front, not in a limitations
section at the end: this paper is about the second group** — practitioners
working with native, non-orthogonal regulatory parts, for whom the
predictive question is still open.

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

This project tests a related but sharper version of that hypothesis, motivated
by a comparison across measurement modalities rather than by assumption.
Johns et al. (2018) measured transcription and translation activity for
tens of thousands of regulatory sequences mined from 184 bacterial genomes,
across three primary recipient hosts (*E. coli*, *B. subtilis*, *P.
aeruginosa*) and, for a 241-sequence subset, three additional hosts — to
our knowledge the only publicly available dataset measuring the *same*
regulatory sequences' activity across multiple bacterial host species at
scale. A second dataset from the same laboratory, DRAFTS (Yim et al.,
2019), measured cell-free (TXTL) transcription of the same 165bp library
across ten species. Comparing the two (Section 3.2) shows cross-host
differences that dominate the in-vivo data are largely absent in the
cell-free data — cell-free lysates retain transcriptional machinery but not
cellular context, so this comparison motivates a specific, testable
prediction: if host specificity is substantially contextual rather than
mechanistic, then host-descriptor features built from genome- and
proteome-encoded machinery composition should carry no signal for
predicting cross-host activity. We pre-registered a test of this prediction
before training any model. It held on every primary comparison. We then
spent four independent attempts trying to find a case where it did not —
because a prediction that survives genuine attempts to overturn it is a
stronger contribution than one accepted at face value. Two of those
attempts initially looked like counterexamples; one survives narrowly, in a
single cell explainable by data quality rather than biological content, the
other did not survive a control we ran specifically because it was the
single most likely thing to invalidate it. We report both, including how
the second was caught, because a result whose own author tried hardest to
break it is more trustworthy than one presented without a fight.

**n_hosts ≤ 6.** Every claim in this paper is bounded by that number, and
the abstract states it explicitly for a reason: the unit of generalization
in this study is the host, not the sequence, and six is not many hosts.

**The benchmark suites we examined — BEND (Marin et al., 2024), Genomic
Benchmarks (Grešová et al., 2023), DART-Eval (Patel et al., 2024),
DNALONGBENCH (Cheng et al., 2025), and the Nucleotide Transformer suite
(Dalla-Torre et al., 2025) — define tasks on human, mouse, *C. elegans*,
*D. melanogaster*, and yeast genomes; we found no bacterial cross-host
regulatory activity task among them.** This is a bounded claim about five
named suites, each read directly rather than assumed — not a claim about
every benchmark suite that exists, which no finite check could support.
Full per-suite coverage, with sources: `out/results/gate13_benchmark_landscape.md`.
We built CROSSHOST — a frozen, versioned benchmark on the Johns et al.
data, with a small demonstration model — to carry out this test and to
leave a durable, reusable artifact behind it.

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

**Secondary data source: DRAFTS.** Yim, Johns et al. (2019), *Molecular
Systems Biology* 15:e8875, from the same laboratory, measured cell-free
(TXTL) transcription of the same 165 bp regulatory-sequence library across
ten bacterial species (1,047 sequences at full coverage, 421 at the
ten-way-usable intersection). Acquired via NCBI's PMC Open Access bulk
package after primary download routes (EMBO/Springer Link, PMC per-file
download, bioRxiv) were blocked by bot-detection or rate-limiting; every
file hash-verified against `data/MANIFEST.json`. DRAFTS's oligo-ID
numbering is entirely disjoint from this project's own library despite
sharing its 184-genome source-mining effort, so we join the two datasets on
exact sequence text, not ID, recovering 112 of DRAFTS's 1,047 sequences as
also present in the three-host library (Section 3.2). DRAFTS is used
descriptively, as a source of independent evidence about cross-host
structure in a different modality — never for model training, and never in
contact with this project's own frozen splits. **We verified directly
against the primary source, before treating this dataset as load-bearing,
that DRAFTS does not itself make the cross-modality comparison this paper
makes** — see Section 3.2 and `out/results/gate12_drafts_verification.md`
for the full, exhaustive check.

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
  variation in native proteome-coverage depth across the six hosts. **This
  is a coarse, reference-database-derived composition proxy, not a
  measurement of dynamic physiological state** — the abundance fractions
  it encodes are closer in kind to the genomic vector's machinery
  annotations than to the membrane potential, supercoiling, resource
  competition, or growth-phase transitions the mechanistic hypothesis in
  Section 3.3 is about. This distinction is load-bearing for how this
  project's results relate to prior physiology-based claims (Discussion)
  and for why this vector's failure (Section 3.4) does not contradict that
  hypothesis.

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
conditioning per se, explained the negative result (Section 3.4).

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

Both amendments are dated, both predate any result they touch, and neither
was revisited after seeing results.

### 2.6 Decision rule

H-MAIN is MET only if the model's 90% bootstrap CI lower bound exceeds the
baseline's 90% CI upper bound; overlapping intervals, regardless of which
mean is nominally higher, are NOT MET. This rule, and the two amendments
above, were fixed before Gate 5 evaluation and were not revisited after
seeing results.

---

## 3. Results

### 3.1 Cross-host measurement structure (Johns et al.'s finding; we confirm and quantify it) (Figures 1-2)

*E. coli* and *P. aeruginosa* activity measurements correlate at ρ≈0.75
(transcription, n=9,741 co-active pairs); any *B. subtilis* pair correlates
at only ρ≈0.16–0.26. **This is Johns et al.'s own reported finding, not a
new observation** — their Supplementary Fig. S13 (n=212, three primary
hosts) and Fig. S15 (n=241, six-host RS241 panel) report exactly this
pattern, and their main text states plainly that *"between recipients,
only E. coli and P. aeruginosa showed significant correlations."* We
verified this directly against their primary-source text and figure
captions rather than assuming a second-hand summary. **Our contribution is
not the raw contrast itself but a corrected, larger-scale quantification of
it and two robustness checks Johns et al. did not run:** a different, more
directly practitioner-relevant translation operationalization (raw protein
level rather than Johns et al.'s translation-efficiency ratio) on
substantially larger per-pair N via a pairwise rather than
three-way-intersected restriction, and a measurement-noise correction plus
a source-composition-confound correction, described next.

**We confirm the gap survives measurement noise.** Anchored by an
empirically-measured *E. coli* transcription measurement reliability of
0.912 (five independent growth-condition replicates, cross-checked to
0.929 by an independent method — the only host×readout combination in the
released Johns et al. data with usable replicate structure; *B. subtilis*
and *P. aeruginosa* have no equivalent replicate structure and are not
directly measured, see Limitations), a classical measurement-error
(disattenuation) correction leaves the *B. subtilis*-vs-*E. coli*/*P.
aeruginosa* gap ratio **unchanged**: 0.341 observed, 0.341 corrected at
reliability 0.9, and only 0.515 even under a deliberately pessimistic
reliability of 0.5.

**We confirm the gap survives a real, independently-surfaced composition
confound.** DRAFTS (Section 3.2) revealed a universal source-genome
GC-activity confound in its own cell-free data (ρ −0.49 to −0.74 across
all ten hosts). The same relationship holds, more weakly, in our own
in-vivo data (ρ −0.20 to −0.61 across the three primary hosts, both
readouts) — raising the question of whether the EC–PA-vs-*B. subtilis*
contrast is itself partly a GC-composition artifact rather than a host
effect. It is not, for either readout, though the two readouts differ in
degree and **we report them separately rather than averaging them**:
partial-correlation control (source-genome GC% held constant, confirmed to
three decimal places by two independent methods — a closed-form
partial-correlation formula and rank-residual regression) leaves the
**transcription** gap ratio at 0.308 (from 0.341 raw, a 9.7% relative
change — **survives intact**) and the **translation** gap ratio at 0.240
(from 0.286 raw, a 15.9% relative change — **survives, but attenuates
appreciably more than transcription does**). The pair carrying the central
claim's positive evidence is nearly unaffected by GC control on its own:
EC–PA moves only −4.7% (transcription) and −0.6% (translation), while the
*B. subtilis* pairs, already close to zero, shrink by a larger relative
amount from a smaller base — expected mechanically, not evidence of a
different underlying effect. A second, assumption-free check — stratifying
by source phylum instead of modeling GC as a continuous covariate — agrees:
within Proteobacteria alone, EC–PA transcription ρ=0.634 and *B.
subtilis*-pairs remain near zero or negative; within Firmicutes alone,
EC–PA ρ=0.812, if anything higher, not lower. Neither correction, alone or
combined, brings the ratio anywhere near parity (1.0). **This gap is not a
measurement-noise or source-composition artifact at any tested level, in
either readout, and it directly bounds what any model — however
conditioned — could achieve on this benchmark.**

### 3.2 The modality contrast (Figure 3)

To test whether the cross-host structure above reflects transcriptional
machinery or broader cellular context, we compared it against DRAFTS (Yim,
Johns et al., 2019, *Molecular Systems Biology* 15:e8875), an independent
dataset from the same laboratory that characterized transcription from the
same 165 bp regulatory-sequence library using cell-free lysates across ten
bacterial species. Cell-free transcription-translation (TXTL) systems
retain core transcriptional machinery — RNA polymerase, sigma factors,
ribonucleotides — but lack an intact membrane, native chromosomal
supercoiling, macromolecular resource competition, and growth-phase-dependent
physiology.

**It does not persist.** Cross-host transcription correlations in DRAFTS
span **0.623–0.911** across all 45 pairwise species comparisons, including
comparisons that cross phylum boundaries (Proteobacteria vs. Firmicutes
vs. Actinobacteria; same-phylum mean 0.852, cross-phylum mean 0.769). This
is a strikingly narrow, uniformly high band. For the one pair directly
comparable between the two studies — *E. coli* and *B. subtilis*, the only
two of our three primary hosts present in DRAFTS (*P. aeruginosa* was not
tested; DRAFTS's own "Pa" abbreviation denotes *Pantoea agglomerans*, an
unrelated species — see Limitations) — the cell-free cross-host
correlation is **0.597** (n=82), more than double the in-vivo figure of
**0.258** (n=3,668) computed on the identical library. *B. subtilis*
itself, a severe outlier in vivo, is unremarkable in cell-free: its mean
correlation with the other nine DRAFTS species (0.783) sits near the
middle of the ten-species distribution.

**DRAFTS's own authors do not make this comparison, and we checked this
directly rather than trusting a second-hand summary before treating it as
load-bearing** — a prior second-hand species list produced a real error in
this project (DRAFTS's "Pa" mistaken for *P. aeruginosa* rather than *P.
agglomerans*, caught in Gate 10), and this claim carries far more weight
than that one did. We searched DRAFTS's full text, all 6 Appendix tables,
all 14 Appendix figures, all 5 source-data spreadsheets' sheet and column
structure, and the peer-review correspondence for any comparison of
cross-species correlation between the two modalities. **None exists.**
DRAFTS reports two separate, unconnected things: within-species in-vitro/
in-vivo fidelity (Pearson r 0.71–0.90 across seven species, their own
Appendix Fig. S5; 94.2% of *E. coli* activities within 1-log variation,
their Fig. 1D), and within-cell-free cross-species structure (the 45-pair
matrix above, their Fig. 3). These two findings sit next to each other in
their Discussion but are never combined into the comparison this paper
makes — that comparison, and the reading that follows from it, is this
project's own contribution. Full verification, checked surface by surface:
`out/results/gate12_drafts_verification.md`.

**This population-level comparison — 82 cell-free sequences against 3,668
in-vivo sequences, not a within-sequence paired test — is the form the
evidence takes, and we state this here, in the paragraph making the claim,
not only in Limitations.** The two libraries share only 112 of DRAFTS's
1,047 sequences by exact sequence-text match (their oligo-ID numbering is
disjoint despite both being drawn from the same 184-genome mining effort by
the same laboratory); of those 112, only 15 have usable in-vivo data in
both *E. coli* and *B. subtilis* simultaneously — far too few to support a
same-sequence paired estimate on its own (computed for completeness:
ρ=0.147, n=15, not used as evidence). The 0.597-versus-0.258 contrast
instead compares two populations measured on the same library by the same
lab under each modality's own full usable set. This is a real limitation on
precision, not a hidden one, and the comparison remains informative: both
figures are large-N, host-modality-consistent estimates of the same
underlying quantity, not noise from a small sample dressed up as a paired
result.

This is not a general failure of the cell-free system to reproduce a given
host's own biology. DRAFTS separately measured 234 sequences in both
cell-free and in-vivo formats for seven species (RS234), allowing a direct
same-host, cross-modality comparison. **Within-species agreement is
good**: Spearman correlations between the cell-free and in-vivo
measurement of the same species range from **0.69 to 0.90** across all
seven species, including *B. subtilis* (0.69) — recomputed directly from
DRAFTS's own released source data, since no clean extractable per-species
table exists in the paper beyond scatter-plot annotations (the paper's own
stated range, 0.71–0.90, is close but not identical; our recomputation is
labeled as such, not presented as a literal quote). Cell-free measurement
is therefore a reasonably faithful proxy for a given host's own
transcriptional behavior; what it does not reproduce is the *difference
between hosts*.

**We are not the first to note cell-free and in-vivo measurements can
diverge, though we believe we are the first to make this specific
cross-host contrast.** Pandi et al. (2022), evaluating an active-learning
optimization workflow (METIS) on a cell-free transcription-translation
system, report that *"the cell-free and in vivo yields for the 20
combinations showed a relatively low correlation of 0.41,"* concluding that
*"although cell-free systems offer rapid prototyping solutions, the
optimal candidates are not necessarily directly transferable in vivo."*
That finding is about a different quantity — absolute yield agreement for
a small set of optimized constructs within a single host, not cross-host
correlation structure across ten species — but it independently supports
the general caution this section's specific finding sharpens: cell-free
fidelity to in-vivo behavior is not something to assume, and here we show
precisely which axis of fidelity holds (within-host ranking) and which does
not (cross-host difference).

**Implications for cell-free prototyping.** Cell-free expression systems
are an active area of synthetic-biology tooling investment, in part
because they promise faster design-build-test cycles than transforming and
growing a living host. This comparison suggests a specific, actionable
boundary on that promise: cell-free measurement appears reliable for
predicting how a regulatory part will behave *within* a given host
(within-species agreement 0.69–0.90), but unreliable for predicting how
that behavior will *differ across hosts* — a practitioner characterizing a
part in one host's cell-free system should not expect the cross-host
comparison to transfer to living cells, where the two hosts may look far
more similar cell-free than they behave in vivo.

**A host-count sweep on DRAFTS was considered and explicitly not run.**
DRAFTS's ten-host panel raised the possibility of training a
host-conditioned model across varying numbers of DRAFTS hosts, to test
whether this project's conditioning failure (Section 3.4) was an artifact
of having only two training hosts per fold. We did not run it: this
section's own finding shows cell-free cross-host correlations are
uniformly high and *B. subtilis* is not an outlier in that modality, so a
sweep conducted there would be well-powered to answer a real question about
cell-free systems specifically — whether training-host count matters for
cell-free cross-species transfer — but could not resolve whether more
training hosts would fix conditioning on the *in-vivo* problem this
project's central claim is about, since the two modalities measure
substantially different phenomena.

### 3.3 The mechanistic hypothesis

Cell-free lysates retain the enzymatic core of transcription — RNA
polymerase and sigma factors — but remove the membrane, native chromosomal
architecture, growth-phase physiology, and resource-competition context a
living cell provides, and cross-host differences largely disappear when
that context is removed (Section 3.2). This is consistent with a broader
pattern in the literature connecting each of those removed elements to
transcriptional output individually: DNA supercoiling, a genome-architecture
property absent in linear cell-free reactions, functions as a global
transcriptional regulator whose effect on a promoter depends quantitatively
on the GC content of that promoter's discriminator sequence (El Houdaigui
et al., 2019); resource competition between a gene of interest and the
host's own expression machinery is directly measurable from lysate
capacity assays and predicts in-vivo expression burden well (in
vitro/in vivo capacity R²=0.74, Borkowski et al., 2018); and growth-phase-dependent
physiology is, by construction, a property only a living, dividing cell
has. **This is a mechanistic hypothesis consistent with the data, not a
demonstrated causal claim.** The individual mechanisms above are
established in the literature; the link from them to this project's own
modeling failure (Section 3.4) is this paper's inference, not a fact we
measured directly. If host specificity resides substantially in this kind
of context rather than in fixed differences in the transcriptional
machinery itself, then host-descriptor features built from machinery and
composition annotations — genomic or proteomic — should carry no signal for
predicting cross-host activity, regardless of how well-constructed they
are. Confirming the mechanism directly would require an experiment this
project cannot run: titrating physiological context back into a cell-free
system, or measuring the same library in vivo across a matched
growth-condition series in multiple hosts.

### 3.4 Testing the prediction (Figures 4-5)

Four independent lines of evidence test the prediction above, and each
comes out the way the hypothesis in Section 3.3 says it should.

**An arbitrary host tag matches or beats both feature vectors, in every
cell tested (pre-registered, H-DIAGNOSTIC).** A free per-host lookup
embedding carrying zero biological content matched or beat both the 37-D
genomic and the 6-D physiology-proxy host-feature vectors at zero-shot in
all **12 of 12** tested (host, readout, variant) cells — 0 clean wins for
the biological features, 10 cells where the free embedding's point
estimate won outright. Neither vector carries detectable content beyond an
arbitrary host label, at this n_hosts — consistent with both vectors
describing machinery/composition rather than the physiological context the
hypothesis says matters (Section 2.2).

**A feature-group ablation finds no group carrying signal.** A five-group
ablation of the genomic feature vector (sigma-factor, anti-Shine-Dalgarno,
tAI/codon usage, RNAP subunit count, chaperone/heme) found no group's
removal changing rank correlation beyond fold-to-fold noise (**max
|Δ|=0.089**) — scoped strictly to this project's own 6-feature,
reference-genome-derived setup, and not framed as testing the Bernstein-lab
physiology hypothesis, since the physiology vector itself was never
ablated and is, as established in Section 2.2, a coarse database-derived
proxy categorically different from either Bernstein-lab paper's direct
wet-lab measurement.

**Sequence-only prediction is not distinguishably beaten under three
conditioning mechanisms, in 16 of 18 tests.** To address the objection
that FiLM specifically — fit from only two training-host feature vectors
per leave-one-host-out fold — is a known-weak choice in a low-data regime,
we implemented two architecturally distinct alternatives (concatenation;
per-host output heads, combined by prediction-space averaging, and a
secondary nearest-host-only variant) under an identical protocol to FiLM's
own zero-shot evaluation. Across 3 mechanism-variants × 6 (host, readout)
cells = 18 comparisons, an alternative mechanism beats sequence-only
distinguishably in exactly **2 of 18** — both at *E. coli* transcription,
the one primary host/readout combination with the largest usable N and the
only measured (not sensitivity-bounded) reliability estimate (Section
3.1). **Every *B. subtilis* cell, both readouts, is statistically
indistinguishable across sequence-only and all three mechanisms** — no
mechanism recovers signal for the host the central question is actually
about. This is a scope-narrowing exception, not a general counterexample:
conditioning helps only where measurement is most reliable and training
data most plentiful, exactly the regime this project's own data (Section
3.1) shows carries the least of the *B. subtilis* uncertainty the
hypothesis is about.

**The pre-registered kill gate failed on all 8 primary comparisons.**
Neither readout, neither host-feature variant, on either primary held-out
host (*B. subtilis*, the official test; *P. aeruginosa*, the charter's
partial-pass check), met the pre-registered bar (Sections 2.5–2.6). The
closest result in the entire gate was *B. subtilis* transcription with
genomic features: model ρ=0.213 [0.154, 0.260] vs. baseline ρ=0.218 [0.160,
0.251] — an overlap of 0.091, "NOT MET, point estimate favors baseline,"
and this cell was not revisited under any alternative mechanism for the
headline result. Every *P. aeruginosa* and *E. coli* primary-mechanism cell
also failed, several by wide margins (e.g. *P. aeruginosa* genomic
transcription: model 0.104 vs. baseline 0.447). This hypothesis was
committed in `PREREGISTRATION.md` before any host-conditioned model was
trained; both amendments to it (Section 2.5) are dated and were made before
the results they touch were seen, and both make the bar harder to clear,
not easier. **The pre-registration, not any individual result, is this
project's principal credibility asset for this claim** — the prediction
was stated and the bar was fixed before either could be shaped by what the
data turned out to say.

### 3.5 Calibration (Figure 6)

The sequence-only model's zero-shot active/inactive classifier is
well-calibrated on *E. coli* (Expected Calibration Error 0.05–0.10) and
severely miscalibrated on *B. subtilis* and *P. aeruginosa* (ECE
0.34–0.47) — a 4–9× degradation. Split-conformal interval coverage for the
strength regressor stayed close to nominal (77–94% empirical against
80/90% targets) regardless of this classifier miscalibration, since
conformal intervals guarantee marginal coverage independent of the
underlying model's calibration — the ECE failure is a point-probability
problem specifically, not something interval-coverage numbers alone would
reveal.

**We quantify the magnitude of this decay; we do not claim it is a novel
warning.** Practitioners deploying a model off its training distribution
already validate at the bench rather than trust raw probabilities, and
published cross-host/cross-dataset predictors show comparable
off-distribution decay elsewhere: LaFleur, Hossain & Salis (2022), a
promoter-design model trained and validated on their own data, report R²
= 0.80 on that data but R² = 0.45–0.65 (Spearman ρ 0.67–0.70) across three
independent external in-vivo datasets they did not use for training — a
real, if less extreme, version of the same off-distribution pattern this
section measures directly for cross-host transfer specifically. **Practical
implication, unchanged from a novel-warning framing: a model that wins on
rank correlation is not automatically trustworthy in absolute probability
terms cross-host**, and any deployment of a model from this suite on a
host meaningfully different from its training hosts should recalibrate
before trusting raw predicted probabilities.

### 3.6 Secondary findings (Figure 7)

**Co-activity, not magnitude, carries what signal exists.** Restricting
cross-host measurement correlation to sequences active in *both* hosts of
a pair — testing the hypothesis that activity is host-specific while
strength above threshold is conserved — produces the **opposite** result
for *B. subtilis* pairs: correlation roughly **halves** (*E.
coli*–*B. subtilis* transcription: 0.655→0.258; *B. subtilis*–*P.
aeruginosa*: 0.508→0.257). Range restriction, a classical mechanical cause
of lower observed correlation independent of any real effect, was checked
directly and ruled out: the restricted subset's interquartile range is
**5×–53× wider**, not narrower, than the pooled set's — if anything this
should preserve or inflate the correlation, making the observed drop more
notable, not an artifact. The one pair where restriction helps (*E.
coli*–*P. aeruginosa*: 0.621→0.754) is the pair with the highest
correlation to begin with, the one place the hypothesis was least needed.

**FiLM is measurably unstable, independent of whether conditioning
helps.** FiLM's fold-to-fold standard deviation of zero-shot Spearman ρ
exceeds every alternative mechanism's — sequence-only, concatenation,
per-host heads — in all 6 (host, readout) cells tested, **4.0–7.8× higher**
at *E. coli* transcription specifically (FiLM std=0.200 vs. concatenation
0.031, per-host-heads 0.050, sequence-only 0.026). A γ/β generator fit from
only two training-host feature vectors is a measurably poor default in
this few-domain regime — useful to anyone conditioning a model on a
handful of domains, regardless of this paper's own negative result, and,
because the instability was checked and quantified rather than assumed
absent, it retroactively strengthens rather than undermines Section 3.4's
finding that FiLM specifically shows no signal.

### 3.7 Two retractions, and why they strengthen the rest (Figure 8)

Two findings that initially looked like results were checked hard enough
to fail, and we report both with their mechanisms in the results, not as
footnotes.

**The ceiling metric.** An earlier "percent of cross-host
measurement-correlation ceiling" metric was proposed, "corrected" in a
later gate, and that correction was itself found wrong (it compared the
new number against the wrong baseline model) in the gate after that. It
was then retired outright — not because either correction was
irreparable, but because the sequence-only model genuinely exceeds the
correctly-computed ceiling for *B. subtilis* on both readouts (102–111%),
a real and explainable consequence of disattenuation (measurement noise
depresses a raw pairwise correlation in a way a model trained on thousands
of examples is not depressed), not a bug and not a metric worth shipping.
It does not appear anywhere in this paper's results.

**Shift-prediction.** A fourth attempt reframed the prediction target:
instead of absolute activity level, predict the *shift* from a reference
host's measured value to a target host, given sequence and the reference
value as inputs. This initially appeared to succeed specifically for *B.
subtilis* transcription — 7 of 8 tested configurations distinguishably
beat a mean-shift constant baseline, a result no absolute-level framing in
this project had produced for that host. **It did not survive a control
we ran specifically because it was the single most likely thing to
invalidate it.** Spearman correlation between the shift and the reference
value itself is strongly negative for most host pairs (as low as −0.812,
*P. aeruginosa*→*B. subtilis* transcription) — textbook regression to the
mean: a sequence with a high reference-host value mechanically has more
room to fall than to rise. A reference-value-only baseline (ordinary least
squares, one feature, no sequence input at all) matched or beat the
original sequence-plus-reference model in **10 of 11** cells originally
reported as wins. A second, independent control — retraining the identical
architecture with sequences randomly permuted relative to their targets,
reference values kept intact — reproduced most of the original model's
performance despite the sequence input now carrying zero real information
(*E. coli*→*B. subtilis* transcription: original ρ=0.432, shuffled-sequence
ρ=0.557, **higher** with no sequence signal at all). **We retract this
finding.** The one cell that does survive both controls (*P.
aeruginosa*→*E. coli* translation, with conditioning: original ρ=0.341 vs.
reference-only ρ=0.142 vs. shuffled ρ=−0.015, non-overlapping) shares no
host or readout with the retracted claim and is reported as an unrelated,
minor, separately-scoped result, not partial vindication.

**We report both retractions in the results, with their mechanisms, as
evidence the controls work, not as embarrassments to minimize.** A reader
who sees a project retract its own most exciting finding, once, after
running the control specifically designed to catch it, has a concrete
reason to trust the findings it did not retract — this is not an assertion
about our own trustworthiness in general, only about this specific,
checkable instance of it.

### 3.8 Appendix: foundation models

Two genomic foundation models (DNABERT-2, 117M parameters; PromoGen2, 148M
parameters), evaluated via frozen-embedding-plus-shallow-head, uniformly
with each other, lost to the 213,956-parameter sequence-only model in 35 of
37 statistically distinguishable comparisons across primary hosts and
RS241. **We make no capacity claim from this comparison; the protocol, not
either model's inherent capability, explains the result.** PromoGen2's own
published native zero-shot protocol (direct likelihood scoring, the
identical Johns et al. dataset) reports 0.68/0.52/0.30 (*E.
coli*/*P. aeruginosa*/*B. subtilis* transcription) — **higher than both
this project's own measured PromoGen2-embedding numbers (0.495/0.490/0.243)
and sequence-only itself, at every host.** PromoGen2's authors did not
report these models failing at this task; under their own protocol,
PromoGen2 succeeds at it. It is this project's frozen-embedding-plus-shallow-head
comparison protocol — chosen for a fair head-to-head between two
architecturally different foundation models — that did not give either
model its best shot. We demote this comparison to an appendix-level result
accordingly and drop any "the model was too small" framing.

---

## 4. Discussion

**What this means for practice.** Cell-free prototyping is a genuinely
useful shortcut for predicting how a regulatory part will behave *within*
a target host (Section 3.2, within-species agreement 0.69–0.90) — but it is
not a reliable shortcut for predicting how that part's behavior will
*differ* if moved to a different host. A practitioner screening constructs
for a single, fixed chassis via cell-free lysate should expect results
broadly consistent with in-vivo behavior. A practitioner trying to rank a
part's suitability *across* candidate hosts using cell-free measurements
alone should not trust that ranking to hold in living cells — the hosts may
look far more similar to each other cell-free than they behave in vivo.
This is the sentence in this paper a working synthetic biologist can act on
today, without waiting on any further modeling result.

**Robustness, not proof: two independent attacks on the central contrast,
both survived.** The *E. coli*–*P. aeruginosa*-vs-*B. subtilis* contrast
(Section 3.1) has now survived two attacks aimed specifically at
explaining it away. It is not attributable to measurement noise: a
classical disattenuation correction, anchored by an empirically-measured
*E. coli* transcription reliability of 0.912, leaves the gap ratio
essentially unchanged (0.341 observed, 0.341 corrected) at that reliability
level. It is also not attributable to source-genome GC composition: a
real, universal GC-activity confound, surfaced independently by the DRAFTS
comparison (Section 3.2), leaves the transcription gap ratio at 0.308 and
the translation gap ratio at 0.240 under partial-correlation control,
confirmed a second, assumption-free way by phylum stratification. **We
call this robustness, deliberately, not proof.** Neither correction was
pre-registered — both were motivated after the fact, by an external
adversarial-review process and by a finding in an independent dataset,
respectively — and surviving two checks does not rule out a third, unchecked
one. What can be said is narrower and, we think, still meaningful: the two
confounds most readily proposed by a careful reader — "your measurements
are just noisy" and "you're measuring source-genome composition, not host
biology" — have each been tested directly on this project's own data and
neither survives contact with it.

**A negative modeling result that confirms a stated prediction is a
different kind of evidence than a negative result accepted at face
value.** Section 3.4's four lines of evidence were not run to see what
would happen; they were run because Sections 3.2–3.3 predicted, in
advance of any model being trained, that they should come out this way.
That two of the four independent attempts to find an exception (Section
3.7) initially looked like counterexamples, and were checked hard enough
that one survives only in a single explainable cell and the other did not
survive at all, is further evidence the prediction is not an artifact of
under-testing.

**For host-conditioned modeling in general:** this project's own
methodological finding (Section 3.6) — that FiLM's fold-to-fold instability
is 4–8× worse than two simple alternatives in a few-training-domain
regime, independent of whether conditioning helps at all — is a caveat
worth carrying into any future host-conditioned model built on a handful of
domains, this dataset or otherwise.

**What this does and does not say about the Bernstein-lab physiology
hypothesis:** nothing decisive, in either direction. This project's
physiology feature vector is a coarse, six-dimension, reference-proteome-derived
composition proxy, not the direct wet-lab measurement (growth curves,
transcriptomes) either Bernstein-lab paper used, and the ablation that
found no genomic feature group carrying signal never touched the
physiology vector at all. What is worth stating precisely: the regime
where the 2024 Chan & Bernstein paper found physiology predictive — six
closely related *Stutzerimonas* strains — is exactly the kind of
phylogenetically close comparison this project's own data independently
shows has the highest cross-host measurement agreement (Section 3.1) and
the only regime where any conditioning mechanism in this project showed a
benefit (Section 3.4). These are not competing claims tested in the same
regime; they are compatible claims tested in non-overlapping regimes of
host similarity, and the most useful next step is a physiology-vs-genomics
comparison conducted with the Bernstein lab's own direct measurement
methodology extended to more phylogenetically distant hosts — precisely
what neither group has yet done and what this project cannot do without a
wet lab.

**What would settle the open questions — and this requires a lab, not more
analysis of this dataset.** Two things, stated plainly as an invitation
rather than a hedge: (1) a direct measurement-reliability estimate for *B.
subtilis* and *P. aeruginosa* specifically, and for translation in any
host — currently bounded only by sensitivity analysis (Section 3.1), and
the single gap most likely to change how confidently this paper's central
finding can be stated; (2) new multi-host regulatory-activity data
measured in vivo, ideally on hosts phylogenetically intermediate between
*E. coli*/*P. aeruginosa* and *B. subtilis*, which would let a future study
distinguish "conditioning helps only in the Gammaproteobacteria-close
regime" from "conditioning helps only where N and reliability happen to be
highest" — two explanations this project's own three-primary-host data
cannot separate. Neither of these is a modeling problem; both require wet-lab
capability this project does not have, which is exactly why they are
stated here rather than left implicit — the labs positioned to run either
experiment will read this as what it is.

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
   instead and also did not beat sequence-only under this project's own
   protocol, though PromoGen2's own native protocol beats it (Section
   3.8) — the capacity objection is weakened, not closed, by either point.
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
8. **Cross-host calibration failure.** See Section 3.5. ECE 0.05–0.10 (EC)
   vs. 0.34–0.47 (BS, PA) — a 4–9× degradation. Recalibrate before trusting
   raw predicted probabilities cross-host.
9. **Regression to the mean is a serious confound for any shift/delta-prediction
   extension of this dataset.** See Section 3.7. If you build a
   shift-prediction model on this or a similarly small-host dataset, run a
   reference-value-only baseline and a shuffled-sequence control before
   trusting the result — this project did not, initially, and the result
   was wrong.
10. **The held-out evaluation split is not cryptographically enforced.**
    The public model-development table still contains the withheld fold's
    real labels, since it is the same table used for development elsewhere
    in the project. Like most small academic benchmarks, the protocol
    relies on the convention being respected, not a technical barrier.
11. **"Pa" means two different organisms depending on which dataset you're
    in.** In this project's own tables, `PA` always means *Pseudomonas
    aeruginosa*, one of the three primary hosts. In DRAFTS's tables (Section
    3.2), `Pa` means *Pantoea agglomerans*, an unrelated Proteobacterium —
    *P. aeruginosa* does not appear anywhere in DRAFTS. A join or comparison
    script matching on the bare two-letter code instead of the full species
    name will silently merge two unrelated organisms' data.
12. **DRAFTS and this project's own tables do not share an ID space.**
    DRAFTS's Oligo IDs (range 13097–14477) and this project's `three_host.parquet`
    OLIGO IDs (range 14478–43726) are entirely disjoint, despite both
    libraries being built from the same 184-genome mining effort by the same
    laboratory. Any join between the two datasets must use sequence text, not
    ID; sequence-text join recovers 112 of DRAFTS's 1,047 sequences (10.7%)
    as also present in the three-host library (Section 3.2).

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
- **We could not run a true within-sequence paired in-vivo/cell-free
  comparison at adequate N.** The DRAFTS-Johns sequence overlap usable in
  both *E. coli* and *B. subtilis* simultaneously is only 15 sequences
  (Section 3.2); the headline 0.597-vs-0.258 modality comparison is
  necessarily between two populations sharing the same library and
  laboratory, not a single paired sample. DRAFTS also measures transcription
  only — no cell-free translation comparison exists, so Section 3.2's
  conclusions are scoped to transcriptional host-specificity and say
  nothing about whether the same pattern would hold for translation. And
  *P. aeruginosa*, one of this project's three primary hosts, is entirely
  absent from DRAFTS's ten-species panel — the modality bridge rests on two
  shared primary hosts, not three.
- **We could not verify two citations we were given for prior evidence of
  cell-free/in-vivo divergence, and dropped them rather than cite them
  unverified.** Two candidate citations (a 2013 Chappell et al. paper and a
  2013 Sun et al. paper) were checked directly against primary sources
  before use, per this project's standing citation-verification practice.
  Chappell, Jensen & Freemont (2013, *Nucleic Acids Research* 41(5):3471)
  in fact reports the opposite of divergence — strong cell-free/in-vivo
  agreement for both promoters (R²=0.946) and RBS elements (R²=0.968) — and
  is not cited here for that claim. No 2013 paper by a first author named
  Sun could be confirmed to report a promoter-ranking divergence finding;
  the closest identifiable match is a methods/protocol paper with no such
  comparison as its subject. Section 3.2 instead cites Pandi et al. (2022),
  independently verified to report exactly the claim needed (cell-free/in-vivo
  yield correlation of 0.41, "not necessarily directly transferable in
  vivo"). Reported here because silently substituting one citation for
  another, without disclosing that the originally-suggested ones did not
  hold up, would understate how close this paper came to citing a source
  for the opposite of what it actually found.
- **We could not separate two explanations for the transcription/translation
  asymmetry** — that translation is fundamentally less cross-host-conserved,
  versus that the FACS-seq translation readout in this dataset is too
  noisy to support this class of analysis — and say so plainly rather
  than picking one. We stated a bet (a mix, weighted toward noise for
  *B. subtilis* specifically, not for EC/PA) and labeled it as inference,
  not fact.
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
  and verified working from a genuine clean git clone (Gate 8.5), and
  re-run in place, successfully, after being extended to cover Gates
  8.5–10.5 (Gate 11).
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

**Every citation below is verified directly against a canonical source**
(CrossRef for DOI-bearing works, the arXiv API for preprints, cross-checked
against the peer-reviewed venue by direct search where no DOI exists) —
mechanically, via `scripts/97_verify_citations.py`, and by hand where the
script flagged something for manual review. Full comparison output:
`out/results/gate13_citation_verification.json`/`.csv`. No citation in this
manuscript remains unverified as of Gate 13.

- Johns, N.I., Gomes, A.L.C., Yim, S.S., et al. (2018). Metagenomic mining
  of regulatory elements enables programmable species-selective gene
  expression. *Nature Methods* 15, 323–329. DOI 10.1038/nmeth.4633.
- Yim, S.S., Johns, N.I., et al. (2019). Multiplex transcriptional
  characterizations across diverse bacterial species using cell-free
  systems. *Molecular Systems Biology* 15, e8875. DOI 10.15252/msb.20198875.
  [DRAFTS; PMC6692573; acquired and verified directly, Gate 10 —
  `out/GATE10_MEMO.md` Task 1; the specific absence of a cross-modality
  cross-species comparison independently re-verified, Gate 12 —
  `out/results/gate12_drafts_verification.md`.]
- Chan, K., Baldwin, G.S. & Bernstein, H.C. (2023). Revealing the
  Host-Dependent Nature of an Engineered Genetic Inverter in Concordance
  with Physiology. *BioDesign Research* 5, 0016. DOI 10.34133/bdr.0016.
- Chan, K. & Bernstein, H.C. (2024). Pangenomic landscapes shape
  performances of a synthetic genetic circuit across *Stutzerimonas*
  species. *mSystems* 9(9), e00849-24. DOI 10.1128/msystems.00849-24.
- Xia, Y. et al. (2026). Design prokaryotic cis-regulatory elements using
  language model. *Nucleic Acids Research* 54(4), gkag122. [PromoGen2;
  PMC12907563.]
- El Houdaigui, B., Forquet, R., Hindré, T., Schneider, D., Nasser, W.,
  Reverchon, S. & Meyer, S. (2019). Bacterial genome architecture shapes
  global transcriptional regulation by DNA supercoiling. *Nucleic Acids
  Research* 47(11), 5648–5657. DOI 10.1093/nar/gkz300.
- Borkowski, O., Bricio, C., Murgiano, M., Rothschild-Mancinelli, B.,
  Stan, G.-B. & Ellis, T. (2018). Cell-free prediction of protein
  expression costs for growing cells. *Nature Communications* 9, 1457.
  DOI 10.1038/s41467-018-03970-x.
- Pandi, A. et al. (2022). A versatile active learning workflow for
  optimization of genetic and metabolic networks. *Nature Communications*
  13, 3876. DOI 10.1038/s41467-022-31245-z.
- LaFleur, T.L., Hossain, A. & Salis, H.M. (2022). Automated
  model-predictive design of synthetic promoters to control
  transcriptional profiles in bacteria. *Nature Communications* 13, 5159.
  DOI 10.1038/s41467-022-32829-5.
- Zhou, Z., Ji, Y., Li, W., Dutta, P., Davuluri, R.V. & Liu, H. (2024).
  DNABERT-2: Efficient Foundation Model and Benchmark for Multi-Species
  Genome. *ICLR 2024*. arXiv:2306.15006. [The model itself was used
  directly via its HuggingFace weights, `zhihan1996/DNABERT-2-117M`, which
  is separately confirmed correct. No DOI exists for ICLR papers; the
  arXiv preprint (2023) predates the ICLR 2024 acceptance — both dates are
  real, ICLR 2024 is the peer-reviewed venue.]
- Marin, F.I., Teufel, F., Horlacher, M., Madsen, D., Pultz, D., Winther,
  O. & Boomsma, W. (2024). BEND: Benchmarking DNA Language Models on
  Biologically Meaningful Tasks. *ICLR 2024*. arXiv:2311.12570. [Named in
  Section 1's benchmark-landscape comparison — tasks defined on the human
  genome only, `out/results/gate13_benchmark_landscape.md`.]
- Grešová, K., Martinek, V., Čechák, D., Šimeček, P. & Alexiou, P. (2023).
  Genomic benchmarks: a collection of datasets for genomic sequence
  classification. *BMC Genomic Data* 24, 25. DOI 10.1186/s12863-023-01123-8.
  [Named in Section 1 — tasks on human, mouse, *C. elegans*, and
  *D. melanogaster*, `out/results/gate13_benchmark_landscape.md`.]
- Patel, A., Singhal, A., Wang, A., Pampari, A., Kasowski, M. & Kundaje, A.
  (2024). DART-Eval: A Comprehensive DNA Language Model Evaluation
  Benchmark on Regulatory DNA. *Advances in Neural Information Processing
  Systems* 37 (NeurIPS 2024 Datasets and Benchmarks Track). DOI
  10.52202/079017-1981. arXiv:2412.05430. [Named in Section 1 — tasks
  entirely on human ENCODE cis-regulatory elements,
  `out/results/gate13_benchmark_landscape.md`.]
- Cheng, W., Song, Z., Zhang, Y., Wang, S., Wang, D., Yang, M., Li, L. &
  Ma, J. (2025). DNALONGBENCH: a benchmark suite for long-range DNA
  prediction tasks. *Nature Communications*. DOI 10.1038/s41467-025-65077-4.
  [Named in Section 1 — tasks on human and mouse only,
  `out/results/gate13_benchmark_landscape.md`.]
- Dalla-Torre, H., Gonzalez, L., Mendoza-Revilla, J. et al. (2025).
  Nucleotide Transformer: building and evaluating robust foundation models
  for human genomics. *Nature Methods* 22(2), 287–297. DOI
  10.1038/s41592-024-02523-z. [Named in Section 1 — 18-task benchmark on
  human, mouse, and yeast, `out/results/gate13_benchmark_landscape.md`.]

**Removed in Gate 13, not carried forward:** a citation to "Specialized
Foundation Models Struggle to Beat Supervised Baselines" (Xu, Gupta, Cheng,
Shen, Shen, Talwalkar & Khodak, ICLR 2025, arXiv:2411.02796 — full author
list confirmed via the arXiv API, not assumed) previously sat in this
section's disclaimer, described as "referenced in Section 1/Discussion."
Checked directly: it is not actually cited anywhere in the current
manuscript body and supports no claim this paper makes. The paper itself
is real and independently verified (arXiv:2411.02796, ICLR 2025 acceptance
confirmed via OpenReview) — it was simply never load-bearing here, and
keeping an unused reference around understates how thin the connection
was. Removed rather than retrofitted into the text to justify keeping it.
