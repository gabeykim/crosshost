# CROSSHOST — Known Issues

For users of this benchmark, not reviewers. Read this before drawing conclusions from any number in this package. Every item below is disclosed because we found it and think you should know, not because a reviewer made us write it down.

## 1. n_hosts ≤ 6 — the generalization unit is the host, not the sequence

Three primary hosts with dense coverage (*E. coli*, *B. subtilis*, *P. aeruginosa*), three more at reduced N (~207 sequences) via RS241 (*S. enterica*, *V. natriegens*, *C. glutamicum*). **Any claim about "genome-encoded functions" or "bacterial regulatory prediction" as a general class is not supported by 6 data points.** If you bootstrap over sequences while holding the host fixed, you are measuring the wrong source of uncertainty — every interval in this project's own results bootstraps over folds/hosts, not sequences, for exactly this reason.

## 2. The translation floor-value artifact

`protein_log10` is pinned at a per-host floor/pseudo-value for a large fraction of nominally "usable" rows — 66.5% (*E. coli*), 89.9% (*B. subtilis*), 10.1% (*P. aeruginosa*). This is almost certainly the original paper's below-detection-limit convention, not real quantitative measurements. **Corrected usable-for-regression N, after excluding floor-pinned rows:** EC 9,146 / BS 1,101 / PA 17,630 — an order of magnitude smaller than the raw "usable" column implies for *B. subtilis* specifically (11,564 raw vs. 1,101 real). If you compute your own translation metrics from `data/core/three_host.parquet`, use `{host}_tl_usable` AND check for the floor value, not `{host}_tl_usable` alone.

## 3. *B. subtilis* measurement noise, and what we could establish about it

*B. subtilis* has the lowest activity rate (17.9% transcription, 10.1% translation-above-floor) and smallest cross-host-comparable N of any primary host. Gate 8 checked whether this noise explains the central finding (EC–PA cross-host correlation ρ≈0.75 vs. *B. subtilis* pairs ρ≈0.16–0.26) via measurement-error disattenuation correction, anchored by an empirically-measured *E. coli* transcription reliability of 0.91 (from 5 independent growth-condition measurements — the only host/readout combination with usable replicate-like structure in the released data). **The gap survives correction, essentially unchanged at the empirically-grounded reliability level.** Even under deliberately aggressive assumptions about *B. subtilis*'s own reliability (as low as 0.5), corrected *B. subtilis*-pair correlations never exceed ~53% of EC–PA's. Full analysis: `out/GATE8_MEMO.md` Task 1. **What we could NOT establish:** a direct reliability estimate for *B. subtilis* or *P. aeruginosa* specifically, or for translation in any host — no replicate/condition-series data exists for these in the released tables. The sensitivity analysis bounds the problem; it does not close it with certainty.

## 4. Fold variance exceeds draw variance, and grows with N — single-fold results are unreliable

Gate 3.5 found fold-to-fold standard deviation exceeds within-fold draw-to-draw standard deviation in 15 of 23 (host, readout, N) baseline combinations checked, with the ratio growing sharply as N increases (e.g. *E. coli* transcription at N=3,000: fold std 0.060 vs. draw std 0.0007, a 91× ratio). **A single fixed test fold can look highly precise (tight draw-to-draw error bars) while being an unrepresentative sample of true host-level performance.** If you evaluate a new model on this benchmark, use the full 5-fold rotation and report fold-to-fold variance, not a single fold's draw-to-draw variance — every number in `baselines/master_baselines.csv` does this; don't regress to a single-fold protocol when extending it.

## 5. The batching cliff — a real ~300× performance bug, now fixed, but know it exists

`crosshost.evaluate()`'s underlying model-inference helpers (and this project's own `predict_seqonly`/`predict_film`) are chunked internally at 1,024 rows. **This is not cosmetic.** An earlier, unbatched version of this code measured 5,810 rows in 28.5s and 23,232 rows (4× the data) in 8,454s (2.35 hours) — a ~300× cost for 4× the data, not linear, on both CPU and Apple Silicon MPS. If you write your own inference loop against this benchmark's models or a similar CNN/embedding-based architecture and pass a large pooled array (multiple folds/hosts at once) to a single unbatched forward call, you may hit the same cliff. Chunk your inference at a few thousand rows or fewer. Full writeup: `out/GATE7_MEMO.md` Task A.

## 6. Evo 2 was not evaluated — a specific technical reason, not a gap in effort

Evo 2 (up to 40B parameters, the largest and most directly relevant genomic foundation model) could not be run in this project's environment. This is a verified hard constraint: the official `evo2` package requires CUDA, Flash Attention, and (for the largest checkpoints) Transformer Engine with FP8 on a Hopper GPU — confirmed against the ArcInstitute repository's README and GitHub issue #67, which shows the library raising `ValueError: Expected a cuda device, but got: cpu` at initialization on non-CUDA hardware. No HuggingFace-compatible port exists. Two smaller genomic foundation models (DNABERT-2, 117M; PromoGen2, 148M) were evaluated instead and also did not beat the sequence-only baseline — this weakens, but does not fully close, the "the model was just too small" objection to the central finding. A free hosted-API path was identified but requires an external account this project's tooling could not create; it remains available as a future option, not chased further per an explicit project decision.

## 7. The retired ceiling metric — a reversal of a reversal, disclosed in full

An earlier "percent of cross-host measurement-correlation ceiling" metric went through two rounds of correction before being retired outright:
1. A Gate 5.5 computation (genomic-model-based) was actually correct, but
2. a Gate 6 "correction" compared it against the wrong baseline model (sequence-only instead of genomic) and published an incorrect fix, which
3. Gate 7 caught, restoring the original numbers as valid — and separately found the real, non-buggy reason to retire the metric anyway: the sequence-only model genuinely exceeds the (correct) ceiling for *B. subtilis* (102–111%), because raw pairwise measurement correlation is attenuated by noise in a way a model trained on thousands of examples is not.

**The metric does not appear anywhere in this package's baseline tables or `PAPER_FRAMING.md`.** Raw cross-host measurement correlation (item 3 above) carries the interpretive weight instead, without a derived ratio. If you see "percent of ceiling" language anywhere referencing this project, it is describing a retracted computation — check the date and gate number before trusting it. Full history: `out/GATE5_5_MEMO.md`'s dated correction, `out/GATE6_MEMO.md`'s dated addendum, `out/GATE7_MEMO.md` Task 1, `out/GATE8_MEMO.md` Task 1.

## 8. Cross-host calibration failure — the sequence-only model is confidently wrong off-distribution

Gate 7 found the sequence-only model's zero-shot active/inactive classifier is well-calibrated on *E. coli* (Expected Calibration Error 0.05–0.10) but severely miscalibrated on *B. subtilis* and *P. aeruginosa* (ECE 0.34–0.47) — a 4–9× degradation. **A model that wins on rank correlation is not automatically trustworthy in absolute probability terms cross-host.** If you deploy any model from this suite on a host meaningfully different from its training hosts, recalibrate on whatever labeled data you can get before trusting its raw predicted probabilities. See `baselines/README.md` and `baselines/gate7_conformal_and_ece.json`.

## 9. The held-out evaluation split is not cryptographically enforced

`held_out_eval/test_inputs.parquet` withholds labels for fold 4, but `data/core/three_host.parquet` (the public model-development data) still contains fold 4's rows with their real labels, since it's the same table used for development elsewhere in this project. A sufficiently motivated user could look up the answers rather than submitting predictions honestly. This benchmark has no hosted server and no technical barrier against that — like most small academic benchmarks, the withheld-eval protocol relies on the convention being respected.
