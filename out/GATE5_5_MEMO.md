# CROSSHOST — Gate 5.5 Memo: The Sequence-Only Ablation

**Date:** 2026-08-06
**Scope:** test the post-hoc hypothesis that the model learned largely host-invariant sequence grammar, and that host conditioning — biological or arbitrary — contributes little in either direction. **No Gate 5 verdict is revisited anywhere in this memo.** Where a Gate 5.5 finding bears on a Gate 5 result, the observation is noted and the Gate 5 verdict is left exactly as it stands.

**This entire gate is post-hoc** — generated after seeing Gate 5's results, not pre-registered. Every claim below is labeled as such and must be labeled as such in the paper.

---

## TASK 1 — The sequence-only ablation: **verdict**

**Sequence-only matches or exceeds both conditioned models, and where a difference is statistically detectable, sequence-only usually wins.** This is a stronger, more decisive result than the pre-specified "INDISTINGUISHABLE" category anticipated — the data doesn't just fail to show conditioning helping, it shows conditioning *actively underperforming* an unconditioned model in the majority of cases where a confident comparison is possible.

**Full grid (36 comparisons: 3 hosts × 2 readouts × 2 eval points [N=0 zero-shot, N=100] × mechanism where applicable):**

| | Count |
|---|---:|
| Distinguishable (90% CI, no overlap) | 18 / 36 |
| — of which sequence-only wins | **17** |
| — of which a conditioned model (genomic) wins | **1** |
| Not distinguishable (intervals overlap) | 18 / 36 |

**The one exception:** *B. subtilis*-transcription, N=100, top-conv mechanism (genomic 0.170 vs. sequence-only 0.081). This is a narrow, likely idiosyncratic result — the sequence-only model's own N=100 top-conv value (0.081) is actually *lower* than its own N=100 frozen-trunk value (0.234) on the identical host/readout, an internal inconsistency suggesting overfitting when conv4 unfreezes at small N for this specific cell, not a general pattern of conditioning being beneficial. It is reported honestly as the one case that goes the other way, not smoothed over.

**RS241 (genuinely unseen hosts) confirms the same pattern**: sequence-only is the top performer in 10 of 12 (config × host × readout) cells, often by a wide margin (PRIMARY *C. glutamicum*-transcription: sequence-only 0.695 vs. genomic 0.596 vs. physiology 0.097). Full table: `out/results/gate5_5_rs241_fourway.csv`.

Full four-way data: `out/results/gate5_5_fourway_comparison.{json,csv}`, `gate5_5_seqonly_verdict_table.csv`. Figures: `out/figures/gate5_5_fourway_comparison.png`, `gate5_5_rs241_fourway.png`.

### Architecture note, honoring the explicit distinction requested

The sequence-only model (`scripts/58_sequence_only_model.py`, `SequenceOnlyCNN`) has **no host-vector parameter anywhere in its `forward()` signature and no `FiLMGenerator` import** — this is categorically different code from Gate 4's FiLM-disabled sanity-check control, which retained a `host_vec` argument (ignored) and a live `use_film` branch. Confirmed parameter count: **213,956 — identical** to the Gate 4 control's count, which in retrospect *proves* that control already carried zero conditioning parameters (its FiLM generator was never instantiated when `use_film=False`); the two models are numerically equivalent at inference but only this gate's model is architecturally incapable of receiving host information, which is what "the conditioning pathway absent" requires.

---

## TASK 2 — Where the transferable signal lives

### Task 2.1 (priority): cross-host measurement correlation — computed directly from the raw data

**This is the single most explanatory number in this entire investigation, and it should have been computed at Gate 2.**

| Host pair | Transcription ρ (n) | Translation ρ (n) |
|---|---:|---:|
| *E. coli* – *P. aeruginosa* | **0.754** (n=9,741) | **0.742** (n=3,826) |
| *E. coli* – *B. subtilis* | 0.258 (n=3,668) | 0.161 (n=866) |
| *B. subtilis* – *P. aeruginosa* | 0.257 (n=2,099) | 0.263 (n=314) |

Computed on floor-corrected actives measured in **both** hosts of each pair (Gate 3's translation floor-correction applied identically). Full data: `out/results/gate5_5_crosshost_measurement_correlation.json`.

**This single table explains most of Gate 5's puzzle.** *E. coli* and *P. aeruginosa* measurements correlate at ρ≈0.75 — a sequence that's active in one is very likely active in the other, **independent of any model, feature, or host-conditioning mechanism whatsoever**. *B. subtilis* correlates with both other hosts at only ρ≈0.16–0.26 — its regulatory response to the same sequences is genuinely, measurably more different. This is a property of the 2018 Johns et al. data itself, not of anything this project built.

### Task 2.2: per-host ceiling vs. actual model performance

| Host | Readout | Ceiling (max pairwise ρ) | Actual zero-shot (genomic) | % of ceiling |
|---|---|---:|---:|---:|
| *E. coli* | transcription | 0.754 | 0.371 | 49.2% |
| *B. subtilis* | transcription | 0.258 | 0.131 | 50.7% |
| *P. aeruginosa* | transcription | 0.754 | 0.338 | 44.8% |
| *E. coli* | translation | 0.742 | 0.155 | 20.9% |
| *B. subtilis* | translation | 0.263 | 0.151 | 57.4% |
| *P. aeruginosa* | translation | 0.742 | 0.131 | 17.6% |

Full data: `out/results/gate5_5_per_host_ceiling.json`. Figure: `gate5_5_performance_vs_ceiling.png`.

**Transcription reaches a remarkably consistent ~45–51% of its implied ceiling across all three hosts, regardless of whether that ceiling is high (*E. coli*/*P. aeruginosa*, 0.75) or low (*B. subtilis*, 0.26).** This reframes *B. subtilis*'s "hardest host" status: it isn't that the model does something categorically worse for *B. subtilis* — it reaches essentially the *same fraction* of a ceiling that happens to be much lower because *B. subtilis*'s biology is genuinely less predictable from *E. coli*/*P. aeruginosa* data. **Translation is a different, honest story**: *E. coli* and *P. aeruginosa* reach only 18–21% of their (high) translation ceiling, while *B. subtilis* reaches 57% of its (low) one — translation transfer is systematically weaker relative to what the data would allow, for the two hosts with the most headroom to lose.

> **DATED CORRECTION — 2026-08-08, Gate 7.** This table's "% of ceiling" numbers are for the **genomic** host-conditioned CNN (as the table's own column header says) and remain **correct and reproducible** — re-verified exactly, to 6 decimal places, against `out/gate4_loho_results.json` in Gate 7. They are **not retracted**. Two things did need correcting:
> 1. **A Gate 6 "correction" to this table was itself wrong.** Gate 6 diffed this file against the *sequence-only* model's zero-shot rho (a different model), found a mismatch, and wrongly concluded this file was buggy. It was comparing the wrong baseline, not catching a real error. That Gate 6 correction is superseded — see `out/GATE6_MEMO.md`'s Gate 7 addendum.
> 2. **The percent-of-ceiling *metric itself* is retired as of Gate 7**, for a real, separate reason found while investigating the above: the sequence-only model (Gate 5.5's stronger model) **exceeds** this same ceiling for *B. subtilis*, both readouts (102–111%) — a genuine result, not a bug, and not fixable by correcting an input. A ceiling a valid model can exceed is not functioning as a ceiling; see `out/GATE7_MEMO.md` Task 1 for the full derivation (raw cross-host correlation is attenuated by measurement noise in *both* hosts, while a trained model's correlation with one host's measurements is attenuated by only that host's noise — so there is no general guarantee a model cannot exceed the raw pairwise figure, and *B. subtilis* is exactly the host with the most measurement noise: smallest N, a 89.9%-floor-pinned translation column, 17.9%/10.1% activity rates). **The genomic numbers above are retained as a correct historical finding about the genomic model specifically; the paper no longer uses "percent of ceiling" as a general framing for any model.** Raw cross-host measurement correlation (Task 2.1 above) carries the interpretive weight going forward on its own.

### Task 2.3: sigma-70 asymmetry revisited — **does not track**

Gate 3 found sigma-70 motif score carries transcription signal (ρ=0.24–0.43) but essentially none for translation (|ρ|<0.02) via a single hand-crafted biophysical feature. The sequence-only CNN's own zero-shot transfer does **not** reproduce this asymmetry:

| | Sequence-only zero-shot ρ (mean across 3 hosts) |
|---|---:|
| Transcription | 0.370 |
| Translation | **0.415** |

**Translation is, if anything, slightly higher.** This is a genuine, honestly-reported non-match, not forced into the expected pattern. The most defensible reading: the CNN learns transferable translation-relevant sequence structure (plausibly RBS spacing/secondary-structure patterns a convolutional receptive field can capture) that a single PWM-style sigma-70-analog feature cannot — Gate 3's ΔG-based translation proxy specifically was shown to carry no signal, but that was one hand-built feature, not a ceiling on what any sequence-based model can extract. Full data: `out/results/gate5_5_sigma70_revisit.json`.

---

## Resolving the three Gate 5 tensions — post-hoc, not a revision of any Gate 5 verdict

1. **H-DIAGNOSTIC failed 12/12** (free embedding matches/beats biological features): now explained by a broader pattern — conditioning of *any* kind (biological or arbitrary identity) adds little over pure sequence learning. H-DIAGNOSTIC's specific comparison (identity vs. biology) was a special case of a more general finding (conditioning vs. none) that this gate makes visible for the first time.
2. **RS241 zero-shot transfer is strong** (ρ 0.4–0.65): explained directly by Task 2.1 — cross-host measurement correlation is often substantial, and a model that learns the shared sequence grammar captures most of that transferable signal without needing correct host-specific calibration.
3. ***V. natriegens* breaks the extrapolation pattern** (extreme distance, unremarkable performance drop): now explained almost tautologically — a model with no host-conditioning pathway has no notion of "distance from the training hosts' feature vectors" to break down under. The extrapolation-distance covariate (Gate 5, Task 2.3) is a property of *conditioned* models specifically; it was never going to predict a sequence-only model's behavior, and *V. natriegens*'s RS241 sequence-only zero-shot performance (Task 1 above) is unremarkable in exactly the way this predicts.

**None of this changes any Gate 5 number.** H-MAIN, H-SCIENCE, and H-DIAGNOSTIC's verdicts stand exactly as reported. What changes is the explanation for *why* they landed where they did — and that explanation is now testable, quantified, and, per Task 2.1, mostly traceable to a property of the underlying 2018 dataset rather than to anything about this project's modeling choices.

---

## TEMPTATIONS TO ADJUST

Retained as a standing section per instruction.

1. **The one cell where conditioning "won" (*B. subtilis*-transcription-top_conv) was tempting to investigate and potentially exclude or re-run** as a likely fluke, given its internal inconsistency (worse than the same model's own N=100-frozen result). I did not re-run it, exclude it, or adjust the bootstrap procedure to smooth it over — it is reported as the one exception it is, with the specific reason for suspecting it's noise stated plainly, not acted upon.
2. **The strength of the sequence-only result (17/18 wins, not just "indistinguishable") was tempting to undersell** — softening "sequence-only usually wins" into the more cautious "conditioning contributes nothing detectable" would have matched the pre-specified verdict categories more comfortably but would understate what the bootstrap intervals actually show. I reported the stronger, more precise finding instead.
3. **Task 2.1's cross-host correlation numbers were tempting to fold into a claim about Gate 5's H-MAIN result being "explained away"** — i.e., to imply the negative H-MAIN result doesn't really count because the ceiling was always low for *B. subtilis*. I did not make this move: Gate 5's H-MAIN comparison was against a *per-host* baseline (data-efficiency), not a *cross-host* ceiling — the two are different questions, and a low cross-host ceiling does not excuse or explain away a data-efficiency loss against a same-host baseline. The two findings are complementary context, not a retroactive defense.
4. **The sigma-70 non-match (Task 2.3) was tempting to explain away by pointing to CNN capacity** without stating clearly that this is speculation. The memo states the CNN-capacity explanation as the "most defensible reading," explicitly hedged, not as an established fact.

**No Gate 5 verdict, threshold, interval, or metric was changed as a result of any of the above.**

---

## WHAT I COULD NOT DO

- **Did not determine why the one BS-transcription-top_conv cell favors conditioning** — flagged as likely noise given its internal inconsistency, not root-caused further (would need a multi-seed sweep of that specific cell, out of scope for a single ablation gate).
- **Did not extend the cross-host measurement correlation analysis to RS241 hosts** — would require the same floor-correction and pairwise-usability logic applied to SE/VN/CG, a natural follow-up not done here given the primary-host analysis was sufficient to establish Task 2.1's point.
- **Did not investigate why translation transfer is systematically weaker relative to its ceiling for *E. coli*/*P. aeruginosa* specifically (18-21% vs. *B. subtilis*'s 57%)** — noted as an honest open question in Task 2.2, not resolved.
- **Did not test whether an even smaller/simpler sequence-only architecture would do just as well** — this gate ran the identical architecture (minus conditioning) to isolate the conditioning variable specifically, not to find the minimal sufficient model.
- **Did not attempt to identify which convolutional features drive the translation transfer signal that the sigma-70/ΔG proxies miss (Task 2.3)** — a natural mechanistic follow-up, out of scope for this gate.

## CONTRADICTIONS WITH THE CHARTER

1. **The charter's premise for building a host-conditioned model — that host identity, biological or otherwise, would materially improve cross-host prediction — is not supported by the evidence gathered across Gates 4, 5, and 5.5.** This is the project's most consequential finding to date: extensive, carefully-controlled experimentation (FiLM conditioning validated as functionally real in Gate 4, two independent feature vectors tested in Gate 5, a clean ablation in this gate) converges on host conditioning contributing little to negative marginal value for this architecture, this dataset, and this task.
2. **Task 2.1's finding — that cross-host measurement correlation should have been computed at Gate 2 — is a real, if minor, process gap**, not a contradiction of a specific charter claim, but worth recording: a cheap, purely-data-driven analysis that would have materially informed the project's expectations going into Gate 4 was not run until forced to by Gate 5's puzzling results. Future projects in this vein should compute cross-condition measurement correlation as a Gate 2-stage task, not a post-hoc one.
3. **The Bernstein claim interpretation (H-SCIENCE, Gate 5) is sharpened, not contradicted, by this gate**: if neither genomic nor physiology conditioning materially helps, and a model with no host information at all does comparably or better, this is evidence bearing on a *stronger* version of the question than H-SCIENCE alone tested — see `out/PAPER_FRAMING.md`.

## FILES WRITTEN

- `out/GATE5_5_MEMO.md` — this memo
- `out/PAPER_FRAMING.md`
- `out/CHARTER_AMENDMENTS.md` (updated: post-Gate-5 scope revision)
- `out/results/gate5_5_crosshost_measurement_correlation.json`, `gate5_5_per_host_ceiling.json`, `gate5_5_sigma70_revisit.json`
- `out/results/gate5_5_fourway_comparison.{json,csv}`, `gate5_5_seqonly_verdict_table.csv`, `gate5_5_seqonly_verdict_summary.json`, `gate5_5_rs241_fourway.csv`
- `out/gate5_5_seqonly_loho_results.json`, `gate5_5_seqonly_calibration_curves.json`, `gate5_5_seqonly_n100_topconv.json`, `gate5_5_seqonly_rs241_results.json` (raw)
- `out/figures/gate5_5_crosshost_measurement_correlation.png`, `gate5_5_performance_vs_ceiling.png`, `gate5_5_fourway_comparison.png`, `gate5_5_rs241_fourway.png`
- `out/models/seqonly_loho_{EC,BS,PA}_fold{0-4}.pt` (15), `seqonly_rs241_{PRIMARY,SECONDARY}_seed{0-4}.pt` (10) — 25 new checkpoints
- `scripts/57_crosshost_measurement_correlation.py`, `58_sequence_only_model.py`, `59_sequence_only_full_pipeline.py`, `60_fourway_comparison.py` (4 new scripts)
