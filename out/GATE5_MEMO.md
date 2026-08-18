# CROSSHOST — Gate 5 Memo: THE KILL GATE

**Date:** 2026-08-05/06
**Scope:** evaluate H-MAIN (both arms), H-SCIENCE, H-DIAGNOSTIC, and RS241 transfer, exactly as pre-registered in `out/PREREGISTRATION.md` (including both dated amendments). No threshold, interval, mechanism, or metric was adjusted during this gate — see "TEMPTATIONS TO ADJUST" for the specific moments that tested this.

---

## GATE DECISION

**H-MAIN is NOT MET on *B. subtilis* (the official test) and NOT MET on *P. aeruginosa* (the charter's partial-pass check), on both readouts, under both host-feature variants, under the pre-registered primary transfer mechanism (frozen-trunk). Every one of the 8 primary-mechanism comparisons resolves as NOT MET — POINT ESTIMATE FAVORS BASELINE, not as an ambiguous overlap.** Per the charter's explicit three-way branch, this is the **negative-result outcome**: the project does not proceed to Gate 6 as originally planned, and does not narrow to a host-dependent claim (that branch requires PA to pass, which it does not). **The benchmark ships with a negative result** — *"we built the first cross-host regulatory activity benchmark, evaluated genomic and physiological host representations, and none beat simply measuring N parts in your target host."* Per the charter, this is stated as a legitimate, useful contribution, not a failure of the project: the model work becomes the baseline suite, the curated benchmark and splits remain a citable asset, and H-SCIENCE and H-DIAGNOSTIC (below) are independently reportable regardless of H-MAIN's outcome. This is not softened and not inflated: the closest any primary-mechanism comparison came was *B. subtilis*-transcription-genomic, where the model's 90% CI [0.154, 0.260] overlaps the baseline's [0.160, 0.251] by 0.091 with the model's point estimate (0.213) still fractionally below the baseline's (0.218) — a genuine near-miss, not a pass, and not a comfortable margin either way.

**One complicating, genuinely positive finding that does not change the above:** RS241 zero-shot transfer to hosts with **zero** training data (*S. enterica*, *V. natriegens*, *C. glutamicum*) shows real signal — ρ up to 0.65, mostly in the 0.4–0.6 range, decisively above zero for all but one cell. This is not eligible to satisfy H-MAIN (which is specifically a data-efficiency claim against a per-host baseline that doesn't exist for these hosts), but it is a substantive, honest, separately-reportable result about the model's general cross-host transfer capability, and belongs in the paper's discussion alongside the negative kill-gate result.

---

## TASK 1 — H-MAIN, both arms

Full data: `out/results/gate5_hmain_results.json`, `gate5_hmain_table.csv`. Figure: `out/figures/gate5_hmain_comparison.png`.

### H-MAIN-TX (transcription)

| Arm | Mean ρ | 90% CI | vs. |
|---|---:|---|---|
| Baseline @ N=100 | 0.107 | [0.066, 0.150] | — |
| **Baseline @ N=3,000** | **0.218** | **[0.160, 0.251]** | (comparison target) |
| Model, genomic, **frozen-trunk (primary)** | 0.213 | [0.154, 0.260] | **NOT MET — POINT ESTIMATE FAVORS BASELINE** (overlap 0.091) |
| Model, genomic, top-conv (supplementary) | 0.170 | [0.126, 0.211] | NOT MET — POINT ESTIMATE FAVORS BASELINE |
| Model, physiology, **frozen-trunk (primary)** | **−0.059** | [−0.122, 0.004] | **NOT MET — POINT ESTIMATE FAVORS BASELINE** (no overlap) |
| Model, physiology, top-conv (supplementary) | 0.118 | [0.066, 0.166] | NOT MET — POINT ESTIMATE FAVORS BASELINE |

**Data-efficiency multiplier (transcription only, never averaged with translation):** the baseline needs **N=3,000** to match the genomic-primary model's point estimate (0.213) — i.e. under the primary mechanism, the model at N=100 does not even achieve a finite multiplier advantage; it takes the *entire* calibration curve for the baseline to reach where the model already is. For context only (supplementary mechanism, not primary): top-conv genomic needs baseline N≈1,000; physiology top-conv needs N≈300.

### H-MAIN-TL (translation)

Baseline arm here is **not** N=3,000 — it is the fold-resolved saturation ceiling (ρ=0.253, 90% CI [0.208, 0.298], at N≈300 — re-derived with raw draws this gate for a proper bootstrap, extending Gate 4's fold-resolved point estimate). This ceiling is a value the per-host baseline **cannot cross at any available N** (flat from N=300 to the per-fold max pool, 813–943 examples) — beating it would have been a stronger claim than H-MAIN-TX's, per the pre-registration's own framing. It was not beaten.

| Arm | Mean ρ | 90% CI | Verdict |
|---|---:|---|---|
| **Baseline ceiling (N≈300, uncrossable)** | **0.253** | **[0.208, 0.298]** | — |
| Model, genomic, **frozen-trunk (primary)** | 0.205 | [0.132, 0.272] | NOT MET — POINT ESTIMATE FAVORS BASELINE |
| Model, genomic, top-conv (supplementary) | 0.224 | [0.196, 0.254] | NOT MET — POINT ESTIMATE FAVORS BASELINE |
| Model, physiology, **frozen-trunk (primary)** | 0.161 | [0.044, 0.266] | NOT MET — POINT ESTIMATE FAVORS BASELINE |
| Model, physiology, top-conv (supplementary) | 0.217 | [0.135, 0.300] | NOT MET — POINT ESTIMATE FAVORS BASELINE |

No finite data-efficiency multiplier applies (baseline never reaches a comparable value at any N by construction) — and the model did not exceed the ceiling either, so this is reported as "not exceeded," not as "multiplier undefined because the model won."

### Host-level context (charter Part V rule 5): P. aeruginosa and E. coli

Same N=100-vs-N=3,000 structure applied to the other two primary hosts (neither has a documented translation-ceiling issue like *B. subtilis*, so no special TL arm is defined for them — this is context, not a second pre-registered hypothesis). **Every primary-mechanism cell for *P. aeruginosa* and *E. coli* is also NOT MET — POINT ESTIMATE FAVORS BASELINE**, several by a wide margin (*P. aeruginosa*-genomic-transcription-frozen: 0.104 vs. baseline 0.447 — the model is not even close). Full table in `out/results/gate5_hmain_results.json` under `host_level_context_PA_EC`. **This is why the charter's PARTIAL PASS branch does not apply**: it requires MET on *P. aeruginosa*, which does not happen either.

**A pattern worth flagging honestly:** under the *supplementary* (non-primary) top-conv mechanism, several *P. aeruginosa* and *E. coli* cells have the model's point estimate slightly *exceeding* the baseline (e.g. *P. aeruginosa*-genomic-transcription: model 0.481 vs. baseline 0.447), landing as "NOT MET — INTERVALS OVERLAP" rather than a clean loss. **This does not change the gate decision** — H-MAIN is scored on the primary mechanism only, per the pre-registration's explicit instruction not to substitute the supplementary result into the headline — but it is worth the paper noting that a model with more fine-tuning capacity at N=100 comes appreciably closer to, and occasionally marginally exceeds, the baseline's point estimate, even though the pre-registered comparison does not credit this.

---

## TASK 2 — H-SCIENCE, the Bernstein test

Full grid: `out/results/gate5_hscience_results.json`, `gate5_hscience_table.csv`. Figures: `gate5_hscience_comparison.png`, `gate5_extrapolation_distance.png`.

### Overall verdict: **host-dependent, mostly statistically indistinguishable**

Genomic and physiology features **do** differ measurably in point estimates (physiology wins 4 of 6 host×readout combinations under the primary mechanism), but **only one of those six differences survives the pre-registered 90% CI as distinguishable**: *B. subtilis*-transcription, where genomic (0.213, CI [0.154,0.260]) clearly beats physiology (−0.059, CI [−0.122,0.004]) with no interval overlap. Every other cell's 90% CIs overlap — the data does not support a confident claim of which representation is better for *E. coli* or *P. aeruginosa*, on either readout.

| Host | Readout | Genomic ρ | Physiology ρ | Point winner | 90%-CI distinguishable? |
|---|---|---:|---:|---|---|
| *E. coli* | transcription | 0.602 | 0.620 | physiology | No |
| *E. coli* | translation | 0.405 | 0.493 | physiology | No |
| *B. subtilis* | transcription | **0.213** | **−0.059** | **genomic** | **Yes** |
| *B. subtilis* | translation | 0.205 | 0.161 | genomic | No |
| *P. aeruginosa* | transcription | 0.104 | 0.244 | physiology | No |
| *P. aeruginosa* | translation | 0.100 | 0.125 | physiology | No |

**Direction, as pre-registered, was not predicted, and the honest answer is not a clean win for either.** Where a difference is real and defensible (*B. subtilis*-transcription), it favors **genomic** — i.e. weakly **challenges** the Bernstein claim on that one cell, since genome-encoded features outperform the proteome-derived physiology proxy specifically on the host/readout where the comparison is least confounded by *B. subtilis*'s known translation data-floor issue. Elsewhere, the data simply doesn't resolve the question at n=5 folds.

### Extrapolation distance as a covariate (Task 2.3)

**n=3 hosts — no statistical claim is supported by this. What follows is descriptive only.**

| Host | Genomic extrap. (37-D, total / per-dim) | Physiology extrap. (6-D, total / per-dim) |
|---|---:|---:|
| *E. coli* | 11.45 / 0.309 | 1.02 / 0.170 |
| *B. subtilis* | 12.21 / 0.330 | 2.24 / 0.374 |
| *P. aeruginosa* | 29.82 / 0.806 | 3.69 / 0.615 |

**Genomic features show a perfectly monotonic (n=3) relationship between extrapolation distance and performance, on both readouts**: *E. coli* (lowest distance) → best performance; *P. aeruginosa* (highest distance) → worst performance; *B. subtilis* in between for both. **Physiology features do not show this pattern on transcription**: *B. subtilis* has the *lowest* transcription performance (−0.059) despite having *lower* extrapolation distance than *P. aeruginosa* (which scores 0.244, higher performance despite being further out-of-distribution) — see `gate5_extrapolation_distance.png`, right panel, where the *B. subtilis* point visibly breaks the pattern the *E. coli*→*P. aeruginosa* trend would predict. Physiology *does* show the monotonic pattern on translation.

**This is a substantive, if narrow (n=3), finding about representation robustness**, exactly the kind the task asked to flag: genomic features degrade *predictably* with distance from the training distribution; physiology features can fail in a way that is not predicted by distance alone, at least for transcription — consistent with Gate 4.5's Task 4 finding that physiology-conditioned optimization can land in qualitatively different (sign-inverted) optima depending on host and readout, a failure mode distance-based reasoning alone would not anticipate.

### Interpretation constraint (as instructed)

**What this result establishes:** for the specific operationalization used here — 37 genome-derived features vs. 6 proteome-derived physiology proxies, a FiLM-conditioned CNN, 3 bacterial hosts, leave-one-host-out — genomic features are not clearly worse than physiology features, and on the one host/readout where a confident distinction is possible, genomic wins. This is **not** strong support for the Bernstein claim as stated (that chassis effects require experimental physiological insight beyond genome-encoded features) — if anything, the one distinguishable result points the other way.

**What this result does not establish:** it does not resolve the Bernstein claim in general. Three hosts is far below what a claim about "genome-encoded functions" broadly would need; the physiology proxy is a coarse depth-matched proteomic summary (6 dimensions), not the kind of direct physiological measurement (growth-law parameters, ribosome content under matched conditions) the Bernstein lab's own work uses; and 5/6 cells were statistically indistinguishable, meaning "genomic isn't clearly worse" is a much weaker claim than "genomic is as good as or better than physiology." **The paper should state this result as host- and readout-specific evidence bearing on, not resolving, the Bernstein claim.**

---

## TASK 3 — H-DIAGNOSTIC

Full data: `out/results/gate5_hdiagnostic_results.json`, `gate5_hdiagnostic_table.csv`. Fallback restated: mean of training-host one-hot embeddings (Gate 3), chosen because with only 2 training hosts per LOHO fold, nearest-neighbor reduces to an arbitrary pick.

**NOT MET across all 12 cells (3 hosts × 2 readouts × 2 variants).** 10/12 are "point estimate favors baseline"; the remaining 2 (*B. subtilis*-transcription-physiology, *P. aeruginosa*-translation-physiology) are "intervals overlap," never a clean model win. The free per-host embedding (identity only, zero biology) matches or beats the host-biology-conditioned model at zero-shot (N=0) in every case tested.

**This is a real, honest, and consequential diagnostic result, consistent with H-MAIN's negative finding, not contradicting it.** The model has not demonstrated it learned host *biology* that a free identity lookup doesn't already capture — this is a genuine limitation of what the paper can claim about mechanism, independent of the H-MAIN kill-gate outcome. It does not, per the charter, kill the project on its own (H-DIAGNOSTIC is explicitly a diagnostic, not a gate) — but combined with H-MAIN's own negative result, it means the paper cannot claim the model learned a useful, biology-grounded cross-host representation; it can only report that it built one, tested it rigorously, and found it did not yet clear either bar.

---

## TASK 4 — RS241: genuinely unseen hosts

Full data: `out/results/gate5_rs241_results_summary.json`, `gate5_rs241_table.csv`. Figure: `gate5_rs241_zeroshot.png`.

**N corrected to 207 recoverable sequences (not 241), per Gate 3.5.** Per-host, per-config usable N (sequence-available, transcription/translation): PRIMARY — SE 125/121, VN 99/115, CG 118/110; SECONDARY — SE 191/196, VN 146/188, CG 176/182. **These are small samples and every number below carries wide uncertainty — stated plainly, not hedged away.**

**Zero-shot (N=0) performance is markedly better than the primary-host H-MAIN picture** — ρ mostly in the 0.4–0.65 range, decisively above zero in all but one cell (PRIMARY-physiology-*C. glutamicum*-transcription, CI crosses zero: [−0.215, 0.433]). Two consistent patterns:
1. **Genomic outperforms physiology in most cells**, sometimes by a wide margin (PRIMARY-*C. glutamicum*-transcription: genomic 0.596 vs. physiology 0.097).
2. **SECONDARY (train EC+PA, drop *B. subtilis*) outperforms PRIMARY (train EC+BS+PA) in most cells** — consistent with Gate 1.5's finding that *B. subtilis*'s weak RS241 coverage dilutes rather than helps transfer to the other RS241 hosts.

**Extrapolation distance for RS241 hosts, genomic features:** *V. natriegens* is dramatically the furthest out-of-distribution of *any* host in this entire project (28.8–35.2 total z-units, more than *P. aeruginosa*'s 29.8 among primary hosts) — yet its zero-shot performance (0.34–0.56) is not correspondingly catastrophic, unlike the clean monotonic pattern seen among the 3 primary hosts. This is a genuine departure from the Task 2.3 pattern, reported honestly rather than smoothed over: **extrapolation distance predicts primary-host genomic performance well but does not extend cleanly to RS241**, at n=3 (RS241) added to n=3 (primary) — still far too few points for a real statistical claim, but worth flagging as a place the tidy story from Task 2 does not fully hold.

**RS241 results are not eligible to satisfy H-MAIN** (no per-host baseline exists for a host with zero training data — that comparison is undefined by construction, exactly why H-MAIN is scored on primary hosts). They are reported here as the purest test of general cross-host transfer this project can run, and as material for the paper's discussion section.

---

## TEMPTATIONS TO ADJUST

Four real moments, reported as instructed — not a formality.

1. **The closest H-MAIN call.** *B. subtilis*-transcription-genomic-frozen came in at 0.213 vs. baseline 0.218 — a near-miss close enough that it was tempting to ask "would a narrower interval (the original 80% t-interval, superseded by Gate 4.5's Amendment 2) have flipped this to MET?" I did not compute or check this. Amendment 2 was committed before any H-MAIN number existed under either rule, specifically to prevent this exact kind of post-hoc rule-shopping; re-litigating it now, after seeing how close the result landed, would defeat the purpose of having pre-registered it. The 90% bootstrap stands, and the verdict is NOT MET as computed.
2. **The RS241 numbers looked good enough to wonder whether they should count toward H-MAIN.** Zero-shot ρ of 0.5–0.65 is a genuinely appealing number next to primary hosts' 0.10–0.24. It was tempting to frame RS241 as "the model does beat a data-efficient baseline, just on different hosts." I did not do this — H-MAIN is specifically a data-efficiency claim against a per-host baseline, and no such baseline exists or can exist for a host with zero training data by construction. Reporting RS241 as if it satisfied H-MAIN would misrepresent what was actually tested. It is reported separately, honestly labeled as not eligible to satisfy the hypothesis.
3. **H-SCIENCE's "4 of 6 cells favor physiology" was tempting to lead with as a summary statistic.** A bare "physiology wins 4/6" reads like a real result. I did not lead with it — only 1 of those 6 differences survives the 90% CI as distinguishable, and reporting the 4/6 count without that qualifier would overstate what the data supports. Both numbers are reported together.
4. **The genomic-vs-physiology extrapolation-distance figure was tempting to normalize by dimensionality only (per-dim) and drop the raw totals**, since per-dim numbers look more directly comparable across a 37-D and a 6-D vector. I reported both, because the FiLM generator operates on the raw (not per-dimension-averaged) vector, so the total is not a meaningless artifact either — dropping it in favor of the more "comparable-looking" number would have been a real, if minor, framing choice made after seeing what each version implied.

**No threshold, interval, mechanism, or metric was actually changed as a result of any of the above.**

---

## WHAT I COULD NOT DO

- **Did not compute a formal statistical test for the extrapolation-distance covariate** — n=3 (primary) and n=3 (RS241) hosts cannot support one; reported ordering and monotonicity only, as instructed.
- **Did not resolve why RS241's extrapolation-distance pattern (V. natriegens far out-of-distribution but not correspondingly poor-performing) departs from the primary-host pattern** — noted honestly as an open departure from Task 2.3's cleaner story, not investigated further (would require more RS241-like hosts than exist in this project's data).
- **Did not re-examine or retrain anything based on how close the H-MAIN-TX genomic result came** — per the gate's governing principle, every number here was produced before verdicts were assigned, and no result was used to trigger a second look at the pipeline that produced it.
- **Did not attempt to identify which specific genomic or physiology features drive the *B. subtilis*-transcription distinguishable result** — a natural follow-up (feature ablation) but explicitly out of scope for Gate 5, which evaluates pre-specified hypotheses, not new exploratory analysis.
- **Did not extend H-SCIENCE's formal (bootstrap-CI) comparison to RS241 hosts** — RS241's genomic-vs-physiology pattern is reported descriptively in Task 4 but was not folded into Task 2's pre-registered 3-primary-host comparison, since that comparison was specifically scoped to primary hosts in the pre-registration.

## CONTRADICTIONS WITH THE CHARTER

1. **The charter's Gate 5 framing anticipated a cleaner three-way outcome** (MET on *B. subtilis*, MET on *P. aeruginosa* only, or NOT MET on both) than what the data actually produced: a clean NOT MET under the primary mechanism, sitting alongside a genuinely encouraging supplementary-mechanism result (several PA/EC cells with the model's point estimate slightly exceeding baseline) and an even more encouraging RS241 zero-shot result. The charter's negative-result fallback text ("none beat simply measuring N parts in your target host") is accurate for the primary-mechanism, primary-host comparison — but does not capture the fuller picture this gate produced. The paper should report the full picture, not just the headline that triggers the negative-result branch.
2. **H-DIAGNOSTIC failing across the board is a real, if secondary, contradiction of what a successful project trajectory would look like** — the charter frames H-DIAGNOSTIC as evidence the model learned biology vs. identity; here it shows no clear evidence of the former. Combined with H-MAIN's result, the model built in Gate 4 has not, on the evidence gathered, demonstrated it does anything a much simpler free-embedding baseline doesn't already do at zero-shot on primary hosts — a materially different position than the charter's own framing of the model as a meaningful secondary asset ("proves the benchmark is answerable, sets a baseline others must beat").
3. **Gate 4.5's decision-rule widening (Amendment 2) was directly load-bearing** — under the original 80% t-interval, several near-miss comparisons might have resolved differently (not checked, per Temptation 1 above, but the CIs are visibly narrower under a t-interval at 80% than a bootstrap at 90%, so this is a real, not hypothetical, possibility). This is exactly the scenario Amendment 2's conservative widening was designed to guard against, and it did its job — worth stating plainly as vindication of that gate's own methodological caution, mirroring how Gate 4's FiLM sanity check vindicated itself.

## FILES WRITTEN

- `out/GATE5_MEMO.md` — this memo
- `out/results/gate5_hmain_results.json`, `gate5_hmain_table.csv`
- `out/results/gate5_hscience_results.json`, `gate5_hscience_table.csv`, `gate5_extrapolation_table.csv`
- `out/results/gate5_hdiagnostic_results.json`, `gate5_hdiagnostic_table.csv`
- `out/results/gate5_rs241_results_summary.json`, `gate5_rs241_table.csv`
- `out/gate5_bs_translation_ceiling_raw.json`, `gate5_rs241_results.json` (raw)
- `out/figures/gate5_hmain_comparison.png`, `gate5_hscience_comparison.png`, `gate5_extrapolation_distance.png`, `gate5_rs241_zeroshot.png`
- `out/models/rs241_{PRIMARY,SECONDARY}_{genomic,physiology}_seed{0-4}.pt` (20 new checkpoints)
- `scripts/50_rs241_train_and_eval.py`, `51_bootstrap_utils.py`, `52_hmain_evaluation.py`, `53_hscience_evaluation.py`, `54_hdiagnostic_evaluation.py`, `55_gate5_figures.py`, `56_rs241_evaluation_summary.py` (7 new scripts)
