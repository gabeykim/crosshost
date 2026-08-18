# CROSSHOST — Gate 7 Memo: Ceiling Correction, Calibration, and Framing Lockdown

**Date:** 2026-08-14. **Outcome: PASS.** Task 1's verdict, stated first as required: **the percent-of-ceiling metric is RETIRED.** It does not appear in the paper. The ρ=0.75 (EC–PA) vs. ρ=0.16–0.26 (BS pairs) cross-host measurement-correlation finding stands on its own and carries the full interpretive weight the retired metric was meant to carry.

---

## Session note: a mid-gate crash, and a real engineering bug found because of it

This gate ran across a session interruption. A diagnostic health-check (user-initiated, separate from Gate 7's own tasks) established: no data corruption, Task 1 fully intact and independently re-verified, Tasks 2–4 at zero on-disk output (safe to rerun, no training involved). That diagnostic also surfaced the single most valuable finding of this gate — see Task A.

---

## TASK A — The batching bug: found, fixed, verified

**The bug.** `predict_seqonly()` (`scripts/58`) and `predict_film()` (`scripts/40`) each ran a single unbatched forward pass over the entire input array. Measured directly (not estimated): 5,810 rows in 28.5s, 23,232 rows (4× the data) in **8,454s (2.35 hours)** — roughly 300× the cost for 4× the data, not linear. This is the same shape of pathology Gate 4 documented for MPS batch sizes above ~2048 (a ~30× cliff); Gate 7 confirms it also exists on CPU. Not a slow machine — a missing loop.

**The fix.** Both functions now chunk internally at `batch_size=1024` (overridable), matching this project's already-established safe size.

**Numerical identity, verified not assumed** (`scripts/audit_batching_fix.py`): chunked vs. unbatched outputs on a 2,000-row input, both models (sequence-only and FiLM-conditioned genomic), all 4 output heads each. **Max absolute difference: 7.15×10⁻⁷** (float32 precision noise, not a real difference). Verdict: NUMERICALLY IDENTICAL.

**Benchmark** (chunked path measured fresh; unbatched baseline reused from the diagnostic run that found the bug, per the calculation note in `scripts/audit_batching_fix.py` — re-measuring an already-confirmed 2.35-hour result would have cost another 2+ hours for no new information):

| | 5,810 rows | 23,232 rows | Scaling |
|---|---|---|---|
| Unbatched (measured) | 28.5s | 8,454s | ~300× for 4× data |
| Chunked, CPU | 26.25s | 114.83s | **73.6× speedup** at 23,232 rows; 4.37× time for 4× rows (linear would be 4.00×) |
| Chunked, MPS | 1.26s | 1.34s | effectively flat in this range |

**Call-site audit — every prior-gate call checked, not assumed safe.** Searched all of `scripts/` for calls to either function: every call site in Gates 1–6 passes at most **5,818 rows** (a single LOHO test fold — the largest fold is 5,818 rows — or the 207-row RS241 set, or a 500-row sanity-check slice in `scripts/41`). None ever pooled multiple folds or hosts into one array before calling. **Zero prior-gate results are affected** — they ran within the safe zone by construction, since LOHO evaluation is always one held-out fold at a time. The only call site that ever exceeded the safe zone is this gate's own new conformal-calibration script, which deliberately pools a held-out host's other 4 folds (~23,232 rows) as a calibration set — a usage pattern no prior gate needed. Prior gates ran slower than they could have (all within already-documented per-fit timescales), not incorrectly.

**Standing rule added** (`out/CHARTER_AMENDMENTS.md`, SR7): all model inference is internally chunked at ≤1024 rows going forward; this failure mode and its fix are documented with the measured numbers so future code doesn't need to rediscover it.

---

## TASK B — Provenance triage (bounded, ~30 min)

`scripts/audit_provenance.py` flags 23 orphan files with no committed producing script (24 before Task 1's crosshost-correlation figure was recommitted in Gate 7, now down to 23 after that figure was given a script — see prior gate work). All 23 are cited somewhere; none are uncited. Classified and verified in `scripts/70_provenance_triage.py` / `out/results/gate7_provenance_triage.json`:

| Category | Count | Files |
|---|---|---|
| **FIGURE, REGENERABLE** | 14 | Calibration-curve PNGs (baselines/ and figures/, 10 total), `gate5_5_fourway_comparison.png`, `gate5_5_performance_vs_ceiling.png` (moot — underlying metric retired), `gate5_5_rs241_fourway.png`, `gate5_rs241_zeroshot.png`. Underlying data all reproducible from committed scripts; only the plotting code wasn't saved. Low risk. |
| **DATA, CITED, VERIFIED** | 9 | `gate5_hmain_table.csv`, `gate5_hscience_table.csv`, `gate5_hdiagnostic_table.csv` (**pre-registered verdicts — prioritized**), `gate5_extrapolation_table.csv`, `gate5_5_sigma70_revisit.json`, `gate5_5_rs241_fourway.csv`, `cnn_architecture_validation.json`, `gate4_bs_translation_ceiling_fold_resolved.json`, `unit_test_summary.json` |
| **DATA, UNCITED** | 0 | none |

**Every DATA-CITED file verified exact-match against an independently committed source** — the three pre-registered-verdict tables diffed directly against their committed JSON (`scripts/52`–`54`'s outputs), exact to float precision; `sigma70_revisit.json` and `rs241_fourway.csv` recomputed independently and matched exactly; the three single-run artifacts (`cnn_architecture_validation.json`, the BS translation ceiling curve, the unit-test summary) cross-checked against `out/state.json`'s already-committed citations of the same numbers, consistent to reported precision. **No pre-registered verdict rests on an unverifiable orphan.** No DECISION NEEDED trigger.

---

## TASK C — Conformal calibration

`out/results/gate7_conformal_calibration.{json,csv}`, `out/figures/gate7_reliability_diagrams.png`, `out/figures/gate7_conformal_coverage.png`.

**Split-conformal regression intervals** (sequence-only LOHO zero-shot vs. B2 per-host N=3,000 baseline, both readouts, all 3 hosts, 80%/90% target coverage): empirical coverage stayed close to nominal for both systems — sequence-only 77–94%, B2 78–94% — consistent with conformal prediction's distributional guarantee, which holds regardless of how well-calibrated the underlying point predictor is.

**Reliability diagrams and ECE for the active/inactive classifier tell a sharply different, more informative story.** Sequence-only's zero-shot classifier is well-calibrated on *E. coli* (ECE 0.053 tx / 0.099 tl) but severely miscalibrated on *B. subtilis* (0.429 / 0.344) and *P. aeruginosa* (0.474 / 0.431) — roughly 4–8× worse. The simple per-host-trained B2 baseline calibrates far better throughout (ECE 0.03–0.23). **This is a real, benchmark-relevant finding stated plainly per the task's instruction:** a model that wins on rank correlation (Spearman rho) can still have badly miscalibrated absolute probabilities cross-host — the two are different properties, and this benchmark's own headline "sequence-only wins" result does not imply its probabilities are trustworthy. The reliability diagrams show the mechanism visually: for *B. subtilis*/*P. aeruginosa*, predicted probabilities in the 0.2–0.5 range correspond to empirical active rates near 0.03–0.08 (badly overconfident in the middle of the range), while the highest-confidence predictions (~0.95) correspond to only ~0.68–0.98 empirical — a genuine S-shaped miscalibration, not noise.

---

## TASK D — Host-feature ablation

`out/results/gate7_host_feature_ablation.json`, `out/figures/gate7_feature_ablation.png`. One representative configuration (genomic variant, *B. subtilis* held out, the charter's official H-MAIN stringent test host), both readouts, all 5 folds. Method: inference-time mean-ablation (zero out each feature group's z-scored columns, no retraining) on the already-trained, already-saved LOHO checkpoints — bounded to one working session as instructed.

**37 genomic columns; 32 covered by the task's 5 named groups (sigma-factor: 4, anti-SD: 1, tAI/codon: 24, RNAP: 1, chaperone/heme: 2); 5 always-retained, ungrouped:** `genome_size_bp`, `gc_content`, `n_16s_rrna_copies`, `n_trna_genes_total`, `n_ribosomal_proteins` — general genome-architecture features the task's 5-group list does not name.

**VERDICT: NO FEATURE GROUP CARRIES DETECTABLE SIGNAL.** Baseline rho = 0.131 (tx, fold-std 0.152) / 0.151 (tl, fold-std 0.056); every group's ablation delta is smaller than one fold-to-fold standard deviation (max |delta| = 0.089). One pattern worth noting without overclaiming it: dropping the sigma-factor group *increased* rho in both readouts (+0.089 tx, +0.053 tl) — the largest and most consistent delta, in the "removing it helps" direction — but it does not clear the noise threshold, so it is reported as an observation, not a finding. Per the task's own framing, this is the expected, clean supporting result given Gate 5.5's central finding (conditioning underperforms sequence-only in 17/18 distinguishable comparisons) — no full grid was run to search for something once the preliminary result showed noise-level effects throughout.

---

## TASK E — Final locked framing

`out/PAPER_FRAMING.md` rewritten in full (not incrementally patched) to reflect Gates 5, 5.5, 6, and 7 together. Contents: central claim in one sentence; the 8-item evidence hierarchy exactly as specified (the user's 5-item ordering adopted as given, with H-SCIENCE, the ECE/calibration finding, and the ablation finding appended as items 6–8 rather than dropped); unhedged scope limits with the n_hosts≤6 ceiling stated for the abstract; Evo 2 as a documented infrastructure limitation (specific technical reason, not an effort gap); the percent-of-ceiling retraction stated plainly with no trace of the old framing surviving; the Bernstein relationship including the new Gammaproteobacteria scoping observation (verified: *E. coli*, *P. aeruginosa*, and the Bernstein lab's own tested panel — Pseudomonads, Halopseudomonads, *Stutzerimonas* — are all Gammaproteobacteria; *B. subtilis* is a Firmicute; this lines up exactly with where this project's own measurement-correlation data shows chassis effects are smallest vs. largest); and the full provenance breakdown (pre-registered / post-hoc / predicted-then-corrected).

---

## The ceiling retraction — the two corrections it took to get here, tied together

**Correction 1 (Gate 5.5 → Gate 7):** the *original* file, `out/results/gate5_5_per_host_ceiling.json`, was re-verified in Gate 7 and found **correct all along** — its numbers match the genomic host-conditioned CNN's zero-shot Spearman rho (plain mean across 5 LOHO folds) to 6 decimal places against `out/gate4_loho_results.json`. Its own column header says "Actual zero-shot (genomic)."

**Correction 2 (Gate 6's own error, caught in Gate 7):** Gate 6 diffed that file against the *sequence-only* model's rho instead — a different model — found a mismatch in 5 of 6 cells, and wrongly concluded the original file was buggy, publishing `gate5_5_per_host_ceiling_CORRECTED.json` as a "fix." That Gate 6 correction is itself superseded; it compared the wrong baseline, not a real bug. Documented with a dated addendum in `out/GATE6_MEMO.md`.

**The real, non-buggy reason to retire the metric anyway:** the sequence-only model — Gate 5.5's *strongest* model, not a weak one — genuinely exceeds the correctly-computed genomic ceiling for *B. subtilis* on both readouts (102–111%, verified via `scripts/66`'s Test D: a direct three-way comparison of ceiling, genomic rho, and sequence-only rho on the same host/readout cells). A ceiling a valid model can exceed is not functioning as a ceiling. Mechanistically: raw pairwise measurement correlation is attenuated by measurement noise in *both* hosts being correlated, while a model trained on thousands of pooled examples is attenuated only by noise in the one host it's evaluated against — there is no general guarantee the latter can't exceed the former, and *B. subtilis* (smallest N, 89.9%-floor-pinned translation column, lowest activity rates of any primary host) is exactly the host with the most measurement noise to produce this effect. A supporting empirical test (`scripts/66` Test A): restricting a foundation model's evaluation to the exact same double-active subset the ceiling was computed on closes most of the apparent gap (DNABERT-2 translation: 104%→88.5% of ceiling), showing population-mismatch between the ceiling's narrow subset and the model's broader eval set is a real, additional contributing factor, on top of the attenuation argument.

Both corrections are reported together here because they were found in the same Gate 7 investigation and separating them would obscure that one was a self-caught process error (wrong baseline, Gate 6) and the other a genuine scientific finding (a real model exceeding a flawed ceiling concept, Gate 5.5's data). Full math and both dated corrections: `out/GATE5_5_MEMO.md`, `out/GATE6_MEMO.md`, `scripts/66_ceiling_metric_investigation.py`.

---

## TEMPTATIONS TO ADJUST

1. **Whether to quietly retrain sequence-only checkpoints for Task C rather than reuse the existing Gate 5.5 ones, to get a "cleaner" fresh run after the crash.** Not done — the existing checkpoints (`out/models/seqonly_loho_*.pt`) were verified present and loadable; retraining would have discarded good, already-verified work for no benefit and reintroduced hours of unnecessary compute.
2. **Whether to soften the ECE finding (Task C) since it complicates the "sequence-only wins" headline.** Not done — reported plainly, with the reliability diagrams shown, and framed explicitly as a caveat to finding 2 in `PAPER_FRAMING.md`, not buried in a subordinate clause.
3. **Whether to note the sigma-factor ablation's directional pattern (+0.089/+0.053, both "removing helps") as a real finding rather than noise, since it was the largest and most consistent delta.** Not done — it does not clear the fold-to-fold noise threshold by the pre-stated criterion, and is reported explicitly as an observation not cleared for a claim, not upgraded because it looked suggestive.
4. **Whether to keep the old `gate5_5_per_host_ceiling_CORRECTED.json` framing partially alive (e.g., "corrected for sequence-only, original for genomic") to avoid a second reversal in as many gates.** Not done — retiring the metric outright, for both models, is the more honest resolution; keeping a patched-together framing to avoid looking like the project reversed itself twice would be optimizing for appearances over correctness.

---

## WHAT I COULD NOT DO

1. **Evo 2 remains unevaluated.** No new attempt was made this gate — per the user's explicit instruction, Gate 6 is closed as final and the NVIDIA-key path is not being chased further. Documented as an infrastructure limitation in the final `PAPER_FRAMING.md`, not silently dropped.
2. **The 14 FIGURE-REGENERABLE orphans were not individually re-scripted.** Their underlying data is verified reproducible; only the plotting code is missing. Left as a Gate 8 packaging hygiene item, consistent with the ~30-minute bound on Task B.
3. **The exact root cause of the original Gate 6 "wrong baseline" error was not further forensically traced** beyond identifying what happened (compared sequence-only's rho against a file whose own column header said "genomic," without checking the label). No deeper process failure was investigated beyond that.
4. **The sigma-factor ablation's suggestive-but-not-significant directional pattern was not followed up with a larger N or a different held-out host** — out of scope for a bounded, one-session ablation task, and would require exactly the "full grid to find something" the task instructed against running given the preliminary result already showed noise-level effects.

---

## CONTRADICTIONS WITH THE CHARTER

1. **Charter Part V rule 6 ("verify before trusting a stated number") was violated twice in this project's own history before Gate 7 caught both instances** — once in the original unsaved ad-hoc computation that (as it turns out) was actually fine, and once in Gate 6's own "correction" of it, which was not. Both are now disclosed, dated, and fixed; the pattern itself (an agent's own correction being wrong) is worth flagging as a standing risk for any future gate that "fixes" a prior gate's number — the fix itself needs the same verification discipline as the original claim, not less.
2. **This gate found and fixed a real performance bug in shared infrastructure (`scripts/40`, `scripts/58`) that has been in the codebase since Gate 4/5.5** — a ~300× cliff sitting latent in code every prior gate called, undetected because no prior gate's usage pattern happened to trigger it. This is disclosed as a genuine gap in this project's own engineering review process up to this point, not just a one-off bug.

---

## FILES WRITTEN

- `out/GATE7_MEMO.md` (this file), `out/PAPER_FRAMING.md` (final, locked)
- `scripts/58_sequence_only_model.py`, `scripts/40_film_cnn_model.py` (batching fix)
- `scripts/audit_batching_fix.py`, `scripts/audit_provenance.py`, `scripts/70_provenance_triage.py`
- `scripts/68_conformal_calibration.py`, `scripts/69_host_feature_ablation.py`, `scripts/71_ablation_figure.py`
- `out/results/gate7_batching_fix_verification.json`, `gate7_provenance_triage.json`, `gate7_conformal_calibration.{json,csv}`, `gate7_host_feature_ablation.json`
- `out/figures/gate7_reliability_diagrams.png`, `gate7_conformal_coverage.png`, `gate7_feature_ablation.png`
- `out/CHARTER_AMENDMENTS.md` (SR7 appended)
- `out/GATE5_5_MEMO.md`, `out/GATE6_MEMO.md` (dated corrections, from Task 1 work carried into this gate)
- `out/state.json` (`gate_7` block)
