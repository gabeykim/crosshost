# CROSSHOST — Gate 3.5 Memo: Fixing the Kill-Gate Foundation

**Date:** 2026-08-04
**Scope:** (1) determine whether *B. subtilis* translation data can support H-MAIN as specified; (2) full 5-fold rotation for all Gate 3 baselines; (3) confirm the Gate 4 compute path; (4) record two Gate 3 findings as standing charter amendments.

**Bottom line up front: PASS**, with one significant premise correction (Task 1) and one very good surprise (Task 3).

---

## Task 1 — Can *B. subtilis* translation support H-MAIN? [LOAD-BEARING]

### The premise needed correcting, with numbers

The concern motivating this task was that Gate 3's floor-value fix left *B. subtilis* with only **1,101** usable translation-regression examples fold-wide, and that H-MAIN's N≥3,000 comparison arm therefore "does not exist" for this readout.

**This is not what the numbers show.** Gate 3's B2 baseline draws its nominal *N* from the full `tl_usable` pool (floor-pinned + genuinely-active rows mixed — **9,300 examples per fold** for *B. subtilis* translation, `scripts/28`), not from the 1,101-example (fold-wide) / ~943-example (per-fold train-pool) `tl_usable_for_regression` subset. The latter is a *downstream* quantity — how many of the *N* measured parts turn out to be non-floor and therefore useful for fitting the strength regressor — not a cap on how many parts can be drawn and measured in the first place. Verified directly from the archived `out/baselines/baseline2_calibration_raw.json` (not re-derived): Gate 3's N=3,000 draws for *B. subtilis* translation used **`n_actual_train=3000`, from a 9,300-example pool, in all 10 draws** (`scripts/35`, "premise correction" section). **The N=3,000 arm was computable, and Gate 3 already ran it.** This holds on every fold (classification-pool sizes: 9,169–9,300 across the 5 folds, all comfortably above 3,000).

This matters because it means the charter's own framing of B2 ("a practitioner can simply measure N parts in their host") is the correct one to apply, and it already accounts for the fact that many measured parts turn out to be floor-pinned/uninformative — that low informative-yield rate is real biology (or real assay-floor behavior), not a flaw in the baseline's construction.

### What the numbers actually establish

| Quantity (fold 0, holding fold 0 out) | Value |
|---|---:|
| `tl_usable` train pool (= B2/H-MAIN's real N-bound) | 9,300 |
| `tl_usable_for_regression` train pool (informative subset) | 943 |
| N=3,000 classification/discovery arm computable | **YES**, all 5 folds |
| N=3,000 direct-regression-only arm computable | NO, max 943 (any fold) |
| Fold-wide totals: `tl_usable` / `tl_usable_for_regression` / active / inactive | 11,564 / 1,101 / 1,101 / 10,463 |

**Reconciling Gate 3's reported numbers** (nominal N, from `tl_usable`, matches archived results): N=100 → ~11 effective regression rows, ρ=0.071±0.177 (huge run-to-run spread — the number is genuinely unstable at this N); N=1,000 → ~100 effective rows, ρ=0.277±0.032; N=3,000 → ~301 effective rows, ρ=0.289±0.010 (stable). *B. subtilis*'s translation-activity/discovery rate is only **~10.1%** per fold — an order of magnitude lower than *P. aeruginosa*'s (~90%) — which is why small nominal-N draws are so noisy for this specific host/readout.

**A second, complementary curve** (Task 1.3, as literally requested — sampling directly from the 943 regression-informative rows, isolating the regression task's own signal ceiling from the discovery-rate effect): rho rises steeply from N=10 (0.070±0.209) to N=100 (0.283±0.012), then **plateaus**: N=300 → 0.290±0.009, N=943 (the full available pool) → 0.289 (single draw, no more data to sample). Full curve, all draws: `out/gate3_5_bs_translation_feasibility.json`.

### Interpretation of the 24.5% figure: two different answers to two different questions

- **Discovery/classification sub-task (will a randomly-drawn part even be translation-active in *B. subtilis*?):** genuine headroom. The activity rate itself (~10%) is what makes small nominal-N samples so uninformative — a cross-host model that can predict *which* sequences are likely non-floor in *B. subtilis*, using patterns learned from *E. coli*/*P. aeruginosa* (much higher activity rates, much richer training signal), could plausibly do meaningfully better than a per-host model starved of examples at small N. This is exactly the kind of headroom H-MAIN was designed to test, and it is real.
- **Strength-regression sub-task (given a part IS active, how strong?):** looks like a **data floor / feature ceiling within the range we can observe** — the direct-regression curve above shows no improvement from 300 to 943 real examples (0.290 → 0.289). Sequence-only k-mer/linear features appear to have an intrinsic ceiling around ρ≈0.29 for this readout, and it is reached by ~100–300 real examples, not thousands. **Caveat, stated honestly:** 943 is still a small number to conclude a *permanent* ceiling from — we cannot rule out that ρ would climb further given, say, 5,000 real examples; we can only say it did not climb from 300→943, the only range we have data for.

**Net verdict: not simply "data floor" as originally framed.** The classification/discovery side shows genuine, exploitable headroom; the regression-precision side shows a plateau within the observed range that a same-type (sequence-only, no host feature) model is unlikely to push past with more of the same kind of data — but a cross-host model isn't limited to "more of the same kind of data," so this doesn't foreclose H-MAIN either.

### H-MAIN respecification: options with numbers, no choice made

Because the N=3,000 arm is in fact computable and already run, strict computability does **not** force a respecification. Two milder, still-real concerns remain, worth deciding on:

1. **Reliability, not existence.** At N=100 (nominal), *B. subtilis* translation's own per-host baseline is enormously noisy (ρ=0.071, std=0.177 — the std exceeds the mean). A single-draw H-MAIN comparison at Gate 5 could go either way by chance alone on this specific host/readout. **Recommend:** require the eventual cross-host model to also be evaluated over multiple seeds/draws for this comparison, not a point estimate, and report both arms' uncertainty side by side.
2. **Whether N≥3,000 is the right ceiling to quote**, given it's already close to asymptotic (ρ 0.277 at N=1,000 → 0.289 at N=3,000, a small further gain). This is not a strong objection — 3,000 remains a reasonable, already-computed number — but it means quoting a strict "30×" data-efficiency multiplier is defensible only if the comparison is fair on the classification/discovery axis, not the regression-ceiling axis.

Presented for Gabriel's decision (not chosen here), in case the reliability concern above is judged sufficient to warrant a change independent of computability:

| Option | Transcription arm | Translation arm | Multiplier vs N=100 |
|---|---:|---:|---|
| **A — transcription primary** | N=3,000 (`tx_usable`=15,697/host total, ample) | reported secondary, N-limitation stated | 30× (unchanged, transcription only) |
| **B — lower translation arm to the regression ceiling** | N=3,000 | N≈943 (or the point of plateau, N≈300) | 9.4× (at 943) or ≈3× (at 300) |
| **C — dual-arm, reported separately** | N=3,000 | N≈943 | two multipliers, not pooled |

Full numeric backing: `out/gate3_5_bs_translation_feasibility.json`. **This is a DECISION NEEDED — not resolved here per instruction.**

---

## Task 2 — Full 5-fold rotation for B1–B4

`scripts/36_fold_rotation_all_baselines.py` re-ran B1/B2/B3/B4 with every fold in turn as the held-out test set (13.8s wall-clock total — affordable specifically because the k-mer/linear substitute, not the CNN, was used throughout, exactly as Gate 3 did). Full results: `out/baselines/fold_rotation_all_baselines.json`, `out/baselines/master_baseline_results.csv` (now fold-resolved).

### Headline: Gate 3's single-fold picture holds up under rotation

| Host | Readout | B2 ρ(N=100), 5-fold mean±std | B2 ρ(N=3000), 5-fold mean±std | B3 ρ, 5-fold mean±std | B4 ρ, 5-fold mean±std |
|---|---|---:|---:|---:|---:|
| *E. coli* | transcription | 0.614±0.060 | 0.623±0.060 | 0.631±0.059 | 0.475±0.030 |
| *E. coli* | translation | 0.481±0.045 | 0.497±0.045 | 0.507±0.042 | 0.049±0.032 |
| *B. subtilis* | transcription | 0.107±0.043 | 0.218±0.060 | 0.211±0.057 | 0.223±0.037 |
| *B. subtilis* | translation | 0.107±0.061 | 0.255±0.058 | 0.258±0.054 | −0.020±0.045 |
| *P. aeruginosa* | transcription | 0.431±0.059 | 0.447±0.059 | 0.453±0.058 | 0.434±0.037 |
| *P. aeruginosa* | translation | 0.206±0.035 | 0.258±0.040 | 0.250±0.042 | −0.011±0.021 |

*E. coli*/*P. aeruginosa* still saturate close to their own N=100 ceiling; *B. subtilis* still shows the largest N=100→N=3000 gap of the three hosts on both readouts. **Gate 3's core finding was not a fold-0-specific fluke.**

### The methodologically important new finding: fold variance vs. draw variance

**Fold-to-fold standard deviation exceeds within-fold draw-to-draw standard deviation in 15 of 23 (host, readout, N) combinations checked for B2**, and the ratio grows sharply with N — precisely the opposite of what "more draws = more precision" would suggest, because averaging 10 draws shrinks draw noise while doing nothing about which 80/20 split you trained/tested on:

| Example | Fold std | Draw std | Ratio |
|---|---:|---:|---:|
| *E. coli* transcription, N=3000 | 0.0601 | 0.0007 | **91×** |
| *E. coli* transcription, N=100 | 0.0598 | 0.0060 | 10× |
| *P. aeruginosa* transcription, N=3000 | 0.0586 | 0.0014 | 42× |
| *B. subtilis* translation, N=3000 | 0.0581 | 0.0076 | 7.6× |
| *B. subtilis* translation, N=100 | 0.0611 | 0.1551 | 0.39× (draw dominates here) |

**At small N, draw variance dominates (the draw itself is the main source of noise). At large N, fold variance dominates (10-90×) — meaning Gate 3's single-fold, tightly-error-barred N=3,000 numbers were more precise than they were accurate.** This is exactly the finding the task asked to check for, and it is real: **single-fold results in this dataset should not be quoted as final without fold rotation.** Recorded as standing rule SR5 in `out/CHARTER_AMENDMENTS.md`.

Fold-resolved calibration-cost figures (error bands = std across folds): `out/baselines/calibration_curve_transcription_fold_resolved.png`, `calibration_curve_translation_fold_resolved.png`.

### B3 extended to RS241 hosts (*S. enterica*, *V. natriegens*, *C. glutamicum*)

**Data gap found and documented, not silently worked around:** RS241 lists 241 oligo IDs with expression values, but only **207/241 (85.9%)** have recoverable regulatory-sequence text in any released table (checked against both the processed library table and the raw Supplementary Data Table 1 directly) — the remaining 34 are not reconstructable from public data. This reduces every RS241 config's usable-N by ~14–15% for sequence-based baselines specifically (`scripts/37`, `out/baselines/baseline3_rs241_extension.json`).

| Config | Held out | Readout | N test (seq-available) | AUC | Spearman ρ |
|---|---|---|---:|---:|---:|
| PRIMARY (train EC+BS+PA) | *S. enterica* | transcription | 125 | 0.791 | 0.433 |
| PRIMARY | *V. natriegens* | transcription | 99 | 0.757 | 0.627 |
| PRIMARY | *C. glutamicum* | transcription | 118 | 0.820 | 0.695 |
| SECONDARY (train EC+PA) | *S. enterica* | transcription | 191 | 0.740 | 0.373 |
| SECONDARY | *V. natriegens* | transcription | 146 | 0.763 | 0.604 |
| SECONDARY | *C. glutamicum* | transcription | 176 | 0.808 | 0.663 |

Translation rho values (0.35–0.44) are similarly substantial; MCC/AUC are not reported for translation here because the RS241 usability flags don't carry a floor-corrected active/inactive label (a simplification disclosed in `scripts/37`, out of scope to redo for a diagnostic baseline given time available). **A free per-host embedding, under its documented fallback, transfers real predictive signal even to genuinely-unseen hosts with no fold structure at all** — worth keeping in mind for Gate 5's H-DIAGNOSTIC test. Full table: `out/baselines/baseline3_rs241_extension_table.csv`.

---

## Task 3 — Gate 4 compute path: a materially better answer than Gate 3 assumed

### The finding: this machine has an unused GPU

Gate 3's "no acceleration" conclusion was correct for the CPU backend (`torch.backends.mkl`/`mkldnn` both `False`, confirmed again here) but **never checked Apple's Metal backend.** `torch.backends.mps.is_available()` → **`True`**. This is a real, zero-cost, zero-setup accelerator that was sitting unused throughout Gate 3.

### Measured cost, at Gate-4-realistic scale (not Gate 3's N≤3,000 calibration values)

Gate 4's actual per-fit training set is a two-primary-host pooled LOHO training pool — **37,000–46,000 examples**, not the ≤3,000 used in Gate 3's B2 calibration grid. Benchmarked at N=45,000 (`scripts/38`, both a raw forward+backward step timing and one real `train_model()` call with real data and early stopping, for cross-validation of the estimate):

| | CPU | MPS (Apple GPU) | Speedup |
|---|---:|---:|---:|
| Raw step cost, 60-epoch fit, N=45,000, batch=1024 | 327.9 min | **11.25 min** | **29×** |
| Real end-to-end fit (early-stopped), N=2,000, real *B. subtilis* transcription data | 824.0 s (13.7 min) | **46.7 s** | **17.6×** |

Both the synthetic raw-step benchmark and the real early-stopped fit point the same direction and roughly the same magnitude — good agreement between a clean per-epoch estimate and actual convergence behavior.

### Extrapolation

| Fit count | CPU | MPS |
|---|---:|---:|
| 12 (Gate 4, LOHO × 2 host-feature variants, single readout head shared) | 2.75 days | **2.2 hr** |
| 24 (Gate 4, both readouts) | 5.5 days | **4.5 hr** |
| 60 (Gate 7, small ensembles) | 13.75 days | **11.2 hr** |
| 120 (Gate 7, full ensembles) | 27.5 days | **22.5 hr** |

**CPU is infeasible for any of these at this environment's scale — confirms and sharpens Gate 3's finding, not just for the N≤3,000 grid but for realistic full-scale fits, which are even slower than Gate 3's own worst-case estimate implied.** MPS changes the picture entirely: **Gate 4 (12–24 fits) is comfortably feasible in a single working session or overnight run, at zero additional cost.** Gate 7 (60–120 fits) is feasible but long (11–22.5 hours) — doable as a background/overnight run, though a cloud GPU would compress this to under an hour if faster iteration on ensemble configuration turns out to matter.

### Cloud GPU feasibility, checked not assumed

Network egress is reachable from this environment (`curl` to `aws.amazon.com` → HTTP 200). **Actual provisioning was not attempted** — creating a cloud account, attaching payment, and launching a billed instance are account-level, real-money actions outside what this session should do without explicit authorization, and are unnecessary now that MPS resolves the core feasibility question. If Gate 7's 11–22.5 hour MPS runtime later proves inconvenient, a small spot GPU instance (e.g. a T4/A10-class instance, roughly $0.50–1.50/hr, a few hours) would fit comfortably in the ~$150 remaining project budget.

### Recommendation (not a decision — Gabriel's call per charter Part VI)

- **Gate 4: run locally on this machine via MPS.** 12–24 fits at ~11–13 min/fit is a few hours, not days. Zero cost, zero new setup — just pass `device="mps"` (now supported in `scripts/29_model.py`, added this gate, backward-compatible default `device="cpu"`).
- **Gate 7: also feasible locally via MPS** as an overnight/background run (11–22.5 hr). A cloud GPU is optional insurance for faster turnaround, not a requirement — worth deciding once Gate 7's exact ensemble design is known (this gate does not decide it).

Full benchmark data: `out/gate3_5_cnn_compute_benchmark.json`.

---

## Task 4 — Two findings recorded as standing charter amendments

Both appended to `out/CHARTER_AMENDMENTS.md` under a new "FINDINGS — Gate 3" section, plus two new standing rules (SR5, SR6) capturing this gate's own methodological lessons (fold-vs-draw variance; nominal-vs-effective N):

1. **Translation floor-value artifact is a first-class data limitation** — floor-corrected N (9,146 / 1,101 / 17,630 for EC/BS/PA) must be used and reported in every subsequent gate's translation work, not the raw "usable" counts.
2. **Mechanistic asymmetry between readouts** — sigma-70 motif score carries real transcription signal (ρ=0.24–0.43); ΔG alone carries none for translation (|ρ|<0.02). A reportable result, and a prior for where host conditioning is likelier to show a clean win.

---

## WHAT I COULD NOT DO

- **Did not re-verify the N=3,000-arm correction by re-running Gate 3's B2 script from scratch** — relied on the already-archived `out/baselines/baseline2_calibration_raw.json` as source of truth (cross-checked its `meta.n_train_pool` and per-draw `n_actual_train` fields directly) rather than re-executing, since it is a read of already-computed, unmodified Gate 3 output. If that file were ever regenerated with different logic, this correction would need re-verification.
- **Did not extend the fold-rotation to RS241 hosts** — RS241 has no fold structure by design (it's a wholly reserved set), so "5-fold rotation" doesn't directly apply there; instead ran both frozen configs (PRIMARY/SECONDARY) once each, training on all pooled primary-host data. This is the correct analog, not a shortfall, but is worth stating explicitly since it looks different from the primary-host fold rotation.
- **Did not resolve the 34/241 RS241 sequences with no recoverable text** — confirmed absent from both the processed and raw released tables; no further retry path exists short of contacting the original authors.
- **Did not attempt to push the *B. subtilis* translation-regression ceiling further** (e.g. richer features, non-linear model) to test whether ρ≈0.29 is a true asymptote or just this feature set's ceiling — out of scope for a baseline-diagnosis gate; flagged as an open question in Task 1's writeup.
- **Did not provision cloud GPU compute** — network reachability confirmed, but actual account/billing setup was judged an authorization-requiring action outside this gate's scope, and unnecessary now that MPS resolves feasibility. See Task 3 recommendation.
- **Did not re-tune B2's k-mer/logistic-regression hyperparameters for the 5-fold rotation** — reused Gate 3's single fixed settings (k=4, LogReg C=1.0, Ridge alpha=10.0) across all folds, consistent with Gate 3's own "set once, applied uniformly" approach, not re-optimized per fold.
- **Did not run the actual CNN for the fold rotation** — Task 2 explicitly specified using the fast k-mer/linear models "since this is affordable precisely because they are ~4 orders of magnitude faster," so this is compliance with the instruction, not a gap, but it means fold-rotation results for the *eventual* CNN-based Gate 4 model remain unverified — only the k-mer/linear proxy was rotated.

## CONTRADICTIONS WITH THE CHARTER

1. **The premise carried into this gate's own Task 1 (from the prior session's framing) was itself incorrect on a load-bearing point** — "N=3,000 arm does not exist for B. subtilis translation" conflated two different quantities (measured-parts N vs. regression-informative N). Corrected here with numbers (see Task 1). This is exactly the kind of self-correction the charter's SR9 (anti-retrofit: "if evidence points somewhere other than this plan, say so") calls for, applied to this project's own prior output, not just the literature or upstream data.
2. **Gate 3's compute-infeasibility conclusion for the CNN was CPU-specific and incomplete** — it never checked for Apple's MPS backend, which is available at zero cost and delivers a 17-29× speedup, changing Gate 4 from "infeasible in this environment" to "feasible in a few hours." This doesn't retroactively invalidate Gate 3's decision to use the k-mer/linear substitute (that was still the right call given the time available in that session), but it changes what should be recommended going into Gate 4.
3. **Gate 3's realistic full-scale CNN cost was substantially worse than its own stated estimate** — Gate 3 quoted "5-7 minutes" for a 60-epoch fit, but that was for N=3,000 (the calibration-curve grid); a Gate-4-scale fit (N≈37,000-46,000, pooling two full hosts) takes ~5.5 hours on CPU, not minutes. This would have been a significant surprise if Gate 4 had proceeded on the Gate 3 estimate without this gate's re-measurement at the right scale.

## FILES WRITTEN

- `out/GATE3_5_MEMO.md` — this memo
- `out/gate3_5_bs_translation_feasibility.json`, `out/gate3_5_cnn_compute_benchmark.json`
- `out/baselines/fold_rotation_all_baselines.json`, `fold_rotation_b2_raw.json`
- `out/baselines/baseline3_rs241_extension.json`, `baseline3_rs241_extension_table.csv`
- `out/baselines/master_baseline_results.csv` (regenerated, fold-resolved)
- `out/baselines/calibration_curve_{transcription,translation}_fold_resolved.png`
- `out/CHARTER_AMENDMENTS.md` (updated: 2 findings + SR5/SR6)
- `scripts/35_bs_translation_feasibility.py`, `36_fold_rotation_all_baselines.py`, `37_baseline3_rs241_extension.py`, `38_cnn_compute_benchmark.py`, `39_compile_fold_resolved_results.py`
- `scripts/29_model.py` (updated: added optional `device` parameter to `train_model`/`predict`, backward-compatible default `cpu`)
