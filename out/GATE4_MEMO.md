# CROSSHOST — Gate 4 Memo: The Host-Conditioned Model

**Date:** 2026-08-04/05 (spanned overnight — LOHO training alone took 9.4 hours)
**Scope:** Build and validate a FiLM-conditioned CNN; train it leave-one-host-out with both host-feature variants; produce calibration curves; spot-check CNN fold variance. **Does not evaluate H-MAIN** — that is Gate 5, by design, per `out/PREREGISTRATION.md` (committed before any training code ran).

**Bottom line up front: PASS** on the FAIL conditions as stated (trained successfully on MPS within a working session; FiLM conditioning is real, not inert; leakage audit does not regress). Along the way: two real numerical bugs were found and fixed before they could corrupt the real training grid, and the results surface genuine, load-bearing fold-instability and transfer-mechanism findings that matter for Gate 5's design.

---

## Task 0 — Pre-registration

Committed to `out/PREREGISTRATION.md` before any training code executed. Freezes: the dual-arm H-MAIN respecification (H-MAIN-TX at N=3,000 vs. H-MAIN-TL at its N≈300 plateau, ρ=0.2527±0.0595 fold-resolved — re-derived fold-resolved here since Gate 3.5's ρ≈0.29 figure was fold-0-only), the decision rule (80% Student-t interval, n=5 folds, non-overlapping-CI requirement, reasoned explicitly), H-SCIENCE's exact comparison protocol, H-DIAGNOSTIC's fallback restatement, and the RS241 207/241 sequence-availability correction.

---

## Task 1 — FiLM-conditioned CNN: architecture, two real bugs, and the sanity check

### Architecture

`scripts/40_film_cnn_model.py`. 4 conv layers (128/128/64/64 channels, kernel widths 15/9/5/3 — unchanged from the Gate-1-verified 165bp part length), BatchNorm after every conv, FiLM conditioning at layers 2 and 3 (small MLP: host vector → hidden(32) → per-channel γ,β, applied as `γ⊙x+β` post-BN pre-ReLU), a residual block wrapping conv4 (skip = conv3's post-FiLM output — a disclosed interpretation of "residual on layers 3-4," see script docstring for the channel-mismatch reasoning), two jointly-trained output heads (transcription, translation) each two-stage (active-logit + strength regressor) = 4 scalar outputs from one shared trunk. **227,844 params (genomic variant) / 226,852 (physiology variant)** — within the ~300k ceiling.

**Compute-efficient grid design, disclosed:** because the architecture trains both readouts jointly from one trunk, the "12 configurations" (3 hosts × 2 readouts × 2 variants) require only **30 actual model fits** (3 hosts × 2 variants × 5 folds), each read off into 2 of the 12 reported rows — not 60. This is a direct consequence of the joint-head design the task itself specifies, not a scope cut.

### Two real bugs found and fixed before the expensive grid ran

1. **MPS boolean-mask-indexing overhead.** The first working loss implementation (`preds[...][mask]`) cost ~875ms/batch of pure gather overhead on this backend. Fixed by computing an elementwise (`reduction='none'`) loss and multiplying by a float mask instead of indexing — cut the masking overhead to ~200ms/batch.
2. **Silent target corruption via `np.nan_to_num`'s default `posinf`.** `tx_norm` contains literal `+Inf` at 398 (EC) / 613 (PA) rows — all correctly mask-excluded, but `np.nan_to_num(..., nan=0.0)` alone leaves `+Inf` mapped to its *default* replacement, ~3.4e38 (float32 max), not 0. Squaring that in MSE overflows to non-finite, and `0 × non-finite = NaN` poisons the whole batch loss — training loss was NaN from epoch 0 in the first real run, despite every tensor individually reading as "finite" under a naive check. Fixed by explicitly zeroing `nan`, `posinf`, **and** `neginf`. Both bugs are documented in `scripts/40`'s own code comments at the exact lines they were fixed, not just here.

**Also found: a sharp, non-monotonic MPS batch-size cliff for this architecture** — batch=1024 measured at 0.52ms/example; batch=2048 at 1.00ms/example (2× worse); batch=4096/8192 at 16-18ms/example (30-35× worse). Not the usual "bigger batch amortizes overhead" curve — likely a memory-pressure-triggered kernel fallback on this machine's 8GB unified memory. batch_size=1024 is used everywhere as a result.

### FiLM sanity check — **VERDICT: FiLM IS CONDITIONING**

Trained matched FiLM-enabled and FiLM-disabled (host vector structurally ignored) models on identical pooled EC+PA data, then swapped the host vector (EC's vs. PA's genomic features) for the same fixed sequences and measured the prediction change:

| Output head | FiLM-enabled \|Δ\| | FiLM-disabled \|Δ\| |
|---|---:|---:|
| tx_active_logit | 1.544 | 0.000000 |
| tx_strength | 0.552 | 0.000000 |
| tl_active_logit | 1.544 | 0.000000 |
| tl_strength | 0.298 | 0.000000 |

The disabled control is *exactly* 0 by construction (proof the test itself is valid); the enabled model shows large, consistent sensitivity. **γ/β statistics after training** (evaluated at EC/BS/PA's genomic vectors): γ2 mean=-0.006, std=0.373, range[-1.12,1.01]; β2 mean=-0.139, std=0.310; γ3 mean=-0.023, std=0.476; β3 mean=0.062, std=0.380 — real spread, not collapsed to γ≈1/β≈0. Full data: `out/gate4_film_sanity_check.json`.

### MPS timing at realistic Gate-4 scale

32.0s/epoch at N=58,084 (pooled EC+PA, both readouts). A single realistic LOHO fit (N≈46,500, 35 epochs) took **552–4,184s** in the actual grid (huge spread — see fold-variance discussion below for why runtimes varied so much; unrelated to correctness, related to which specific batches triggered more expensive paths). **Full 30-fit grid: 33,745.9s = 9.37 hours.**

---

## Task 2 — Leave-one-host-out training, both variants

Full results: `out/gate4_loho_results.json`, `out/gate4_loho_results_table.csv`. 35 epochs, no early stop triggered in any of the 30 fits (training loss kept improving through the full budget every time).

### Fold-resolved LOHO results (N=0, zero-shot on held-out host)

| Host | Variant | Readout | Spearman ρ (mean±std, 5 folds) | AUC (mean±std) |
|---|---|---|---:|---:|
| *E. coli* | genomic | transcription | 0.371±0.200 | 0.617±0.154 |
| *E. coli* | genomic | translation | 0.155±0.275 | 0.637±0.143 |
| *E. coli* | physiology | transcription | 0.292±0.162 | 0.699±0.041 |
| *E. coli* | physiology | translation | 0.255±0.264 | 0.691±0.037 |
| *B. subtilis* | genomic | transcription | 0.131±0.152 | **0.830±0.028** |
| *B. subtilis* | genomic | translation | 0.151±0.056 | 0.794±0.008 |
| *B. subtilis* | physiology | transcription | 0.239±0.027 | **0.845±0.015** |
| *B. subtilis* | physiology | translation | 0.238±0.050 | 0.811±0.039 |
| *P. aeruginosa* | genomic | transcription | 0.338±0.130 | 0.595±0.131 |
| *P. aeruginosa* | genomic | translation | 0.131±0.110 | 0.572±0.053 |
| *P. aeruginosa* | physiology | transcription | **-0.061±0.398** | 0.615±0.200 |
| *P. aeruginosa* | physiology | translation | 0.252±0.064 | 0.518±0.083 |

### Headline finding: fold stability is dramatically host-dependent

***B. subtilis* is the most stable held-out host by far** — lowest std across every variant/readout (0.027–0.152), all-positive-except-one-fold correlations, and the tightest, highest AUCs (0.79–0.88) of the three hosts. This is reassuring specifically for H-MAIN, which is defined on *B. subtilis*.

***P. aeruginosa* is the least stable, catastrophically so for the physiology variant on transcription**: ρ ranges from -0.475 to +0.500 across the 5 folds (std=0.398), and one fold's AUC (0.228) is *worse than random guessing*. This is plausibly connected to an independent, twice-established finding from Gate 1/1.5: *P. aeruginosa* has the weakest-characterized physiology data of all 6 hosts (29–57% relative disagreement between reasonable dataset choices, vs. 9–30% for other hosts) — this gate adds a *third*, model-training-based signal pointing the same direction, though it is not proof the physiology feature vector itself is the sole cause (PA-genomic is also more variable than BS-genomic, just less extremely).

### H-SCIENCE preview (numbers only — Gate 5 does the hypothesis test)

`out/figures/hscience_genomic_vs_physiology.png`. **Genuinely mixed, not a clean win for either representation:** physiology has a higher mean ρ in 4/6 host×readout combinations (both translation cases for all 3 hosts, plus *B. subtilis* transcription) and is *more stable* for *B. subtilis* specifically (std 0.027 vs. 0.152 on transcription) — but is dramatically *less* stable for *P. aeruginosa* transcription, as above. Genomic wins outright only on *E. coli* and *P. aeruginosa* transcription. **This is Gate 4's job — produce the numbers; no H-SCIENCE verdict is stated here per the pre-registration's own rule against computing hypothesis comparisons in this gate.**

---

## Task 3 — Calibration curves

Full grid: 2,100 fine-tune fits (3 hosts × 2 variants × 5 folds × 7 N/mechanism cells × 10 draws), 9,900.4s (2.75 hr). Results: `out/gate4_calibration_curves.json`, `out/gate4_calibration_results_table.csv`, figures at `out/figures/calibration_curve_{host}_{readout}.png`.

### Headline finding: the transfer mechanism matters enormously when the base model is already unstable

For *P. aeruginosa*-physiology folds with a bad zero-shot base (per Task 2), the **frozen-trunk mechanism (head_only, N≤300) inherits and is stuck with the base model's bad calibration** — e.g. fold 2's transcription ρ stays pinned around -0.32 across N=10 through N=300 (frozen). The moment **conv4 unfreezes (top_conv, same N=300)**, the same fold jumps to ρ≈+0.41 — and climbs further to +0.47 by N=3,000. This pattern repeats across the unstable PA-physiology folds (1, 2, 4): frozen-trunk fine-tuning cannot fix a badly-calibrated trunk; unfreezing even one conv layer can. This is exactly why the task specified running both mechanisms at N=300 — the crossover is real and large, not a marginal effect.

### Second finding, found while building the figures honestly: small-N fine-tuning can actively *hurt* the physiology variant relative to its own zero-shot baseline

*B. subtilis*-physiology-transcription's zero-shot (N=0) result is ρ=0.239 (good, stable — see Task 2). Fine-tuning on just 10-300 of *B. subtilis*'s own examples (frozen trunk) **drops this to slightly negative** (-0.02 to -0.04, flat across N=10/30/100/300) before recovering to ρ≈0.15-0.25 once top_conv unfreezes at N≥300 (`out/figures/calibration_curve_BS_transcription.png` — the N=0 zero-shot point sits visibly *above* the entire frozen-trunk curve). Checked for a bug, not just accepted at face value: per-draw values at N=10 are genuinely noisy in both directions (fold0 draws range -0.13 to +0.26; fold1 draws are consistently negative, -0.13 to -0.17) — this is real high-variance behavior, not a broken pipeline. Mechanistic reading: with only the held-out host's *own fixed* feature vector as conditioning input (constant across all N examples, since they're all from the same host), the FiLM generator's gradient signal at tiny N comes from a very low-diversity sample, and can drift the generator's output for *that specific host* away from the more-regularized zero-shot value learned from the general 2-host cross-training signal. The genomic variant does **not** show this dip — its N=10 point is already above N=0 and rises smoothly. Not root-caused further (would need controlled ablation, e.g. varying only the FiLM generator's learning rate at tiny N) — flagged as a real, reportable finding for Gate 5/7, not resolved here.

---

## Task 4 — CNN fold-variance spot check

Configuration requested: *B. subtilis* held out, transcription, genomic features, N=100.

| | Fold-to-fold std | Mean within-fold draw-to-draw std | Ratio (fold/draw) |
|---|---:|---:|---:|
| CNN (this config) | 0.0683 | 0.0329 | 2.08 |
| k-mer/linear (Gate 3.5, same config) | 0.0426 | 0.1142 | 0.37 |

**On this specific, narrow config, CNN fold variance (0.068) is ~1.6× the k-mer model's (0.043) — real but not dramatic**, and my pre-registered 2× threshold does not trigger. **However, this specific config is misleading in isolation**: it happens to land on *B. subtilis*, the CNN's most stable host. The broader evidence from Task 2's full grid — *P. aeruginosa*-physiology-transcription's fold std of 0.398 (9× the k-mer baseline's typical range) and one fold with AUC below random chance — shows CNN fold variance can be **far more severe** than this narrow spot-check captures. Full data: `out/gate4_cnn_fold_variance_spotcheck.json`.

**DECISION NEEDED, based on the full-grid evidence, not just the literal spot-check cell:** CNN fold-to-fold variance is comparable to the k-mer baseline for the *B. subtilis* configuration H-MAIN is actually defined on, but is substantially worse for at least one other host/variant combination (*P. aeruginosa* physiology). Since H-MAIN's decision rule (pre-registered, 5-fold, 80% CI) is specifically scoped to *B. subtilis*, this does **not** immediately invalidate the pre-registered rule — but Gate 5 should re-check this exact ratio for *B. subtilis*-translation (the other H-MAIN arm) before trusting the CI-based comparison there, since it was not separately spot-checked here.

---

## What Gate 4 deliberately does NOT do

No H-MAIN comparison is computed anywhere in this memo or its outputs — per the pre-registration, that is Gate 5's job, and mixing it in here would undermine the reason the pre-registration exists.

---

## WHAT I COULD NOT DO

- **B. subtilis's `growth_rate_mu_h` (one of 6 physiology features) is imputed**, not measured — the Zhu et al. 2025 PNAS supplementary value has been inaccessible since Gate 1 and remains so. Imputed as the mean of the other 5 hosts' raw values (1.852 h⁻¹), disclosed in `scripts/40`'s docstring and printed at every run. This affects every physiology-variant result touching *B. subtilis* (as training host or held out) — a real, carried-forward limitation, not a new one, but newly load-bearing here since Gate 4 is the first gate to actually train on this feature.
- **Did not root-cause the small-N physiology fine-tuning degradation** (Task 3, second finding) beyond a plausible mechanistic reading — would need a controlled ablation (e.g., separate learning rates for the FiLM generator vs. output heads at tiny N) that was out of scope for this gate's time budget.
- **Did not investigate why LOHO fit wall-clock time varied so much** (552s to 4,184s across nominally-identical-cost fits) — plausibly thermal throttling or memory pressure on a consumer laptop running a 9.4-hour job, not investigated further since it didn't affect correctness, only total wall-clock (which was still within a working session/overnight run).
- **Did not save the 2,100 individual fine-tuned model checkpoints from Task 3** — only their evaluation metrics were retained (saving all of them was judged unnecessary disk usage for artifacts that are cheap to regenerate from the 30 saved LOHO base checkpoints plus the fixed random seeds).
- **Did not extend the fold-variance spot check to more than the one requested configuration** before writing this memo — used the full LOHO grid's own variance patterns (already computed, no extra cost) as broader context instead, which is arguably more informative than additional narrow spot-checks would have been, but is a different kind of evidence than a second formal spot-check.
- **Did not attempt to fix or reduce PA-physiology's instability within this gate** — Task 2 explicitly asks for the model's own performance, not model improvement; flagged for Gate 5/7 to weigh when interpreting H-SCIENCE.

## CONTRADICTIONS WITH THE CHARTER

1. **Two real numerical bugs (MPS masking overhead, `nan_to_num` posinf corruption) were present in the first working implementation of code this gate itself wrote** — not inherited from an earlier gate. Caught before the expensive grid ran only because the FiLM sanity check (explicitly required by the task) was run and its output inspected carefully rather than trusted at face value (the first "successful" sanity-check run had `loss=nan` printed at every logged epoch, which would have been easy to skip past given the final gamma/beta numbers looked superficially reasonable). This is a direct vindication of the task's own instruction to "report the FiLM sanity-check result prominently" and treat non-finite loss as a stop-and-fix — worth stating because it demonstrates the validation step earned its place in the task list, not a formality.
2. **CNN fold-to-fold variance is real and, for at least one host/variant combination, severe** — worse than anything seen in Gate 3.5's k-mer baseline analysis. The pre-registered H-MAIN decision rule (5 folds, 80% CI) was designed anticipating baseline-level fold variance; Gate 5 should re-verify it holds up for the specific *B. subtilis* translation arm before trusting a MET/NOT MET call, per Task 4's DECISION NEEDED above.
3. **Small-N transfer learning is not monotonically "more data is better" for the physiology variant** — a genuinely unexpected, charter-relevant finding, since the charter's whole H-MAIN framing assumes calibration performance improves smoothly with N. It does for genomic features; it does not for physiology on at least one host/readout. This should inform Gate 5/7's interpretation of any physiology-variant calibration curve, not just this gate's own numbers.

## FILES WRITTEN

- `out/PREREGISTRATION.md`
- `out/GATE4_MEMO.md` — this memo
- `out/gate4_bs_translation_ceiling_fold_resolved.json`, `gate4_film_sanity_check.json`, `gate4_loho_results.json`, `gate4_loho_results_table.csv`, `gate4_calibration_curves.json`, `gate4_calibration_results_table.csv`, `gate4_cnn_fold_variance_spotcheck.json`
- `out/models/loho_{host}_{variant}_fold{N}.pt` (30 files), `film_sanity_check_model_on.pt`
- `out/figures/calibration_curve_{host}_{readout}.png` (6 files), `hscience_genomic_vs_physiology.png`
- `scripts/40_film_cnn_model.py` (architecture + training utilities), `41_film_sanity_check.py`, `42_loho_training.py`, `43_calibration_curves.py`, `44_cnn_fold_variance_spotcheck.py`, `45_compile_gate4_results.py`
