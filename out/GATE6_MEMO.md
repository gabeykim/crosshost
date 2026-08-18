# CROSSHOST — Gate 6 Memo: Foundation Model Head-to-Head

**Date:** 2026-08-07. **Outcome: PASS.** Two foundation models (DNABERT-2, PromoGen2) were evaluated end-to-end under the identical frozen-fold, identical-metric, 90%-bootstrap protocol used throughout this project. Evo 2 was not evaluated — see Task 1 and WHAT I COULD NOT DO.

---

## GATE DECISION, stated plainly

**The central claim survives intact.** Neither foundation model tested beats the 213,956-parameter sequence-only CNN in a way that opens a capacity-based objection to Gate 5's negative result. Of 48 total FM-vs-sequence-only comparisons across primary hosts and RS241 (24 + 24), 37 are statistically distinguishable at 90% confidence; of those, sequence-only wins 35 and a foundation model wins only 2 (both DNABERT-2/PromoGen2 at *E. coli* transcription zero-shot specifically). At the pre-committed zero-shot-only threshold from `out/GATE6_EXPECTATIONS.md` Section 2 ("the claim narrows if a foundation model distinguishably beats sequence-only on at least half of the six primary host×readout cells"), **neither model reaches even one-sixth of that bar** (1 of 6 cells each, both the same cell). This is consistent with DART-Eval and the ICLR 2025 "Specialized Foundation Models Struggle to Beat Supervised Baselines" finding, and with the honest prior stated before any evaluation ran.

A significant correction to a Gate 5.5 number was found and fixed while producing this gate's percent-of-ceiling figures — see Section 3 and the dated addendum in `out/GATE6_EXPECTATIONS.md`. It changes the *shape* of the ~50%-ceiling story but not the gate's bottom line.

---

## TASK 0 — Pre-stated expectations, outcome

Full document: `out/GATE6_EXPECTATIONS.md`, written and committed before any FM evaluation ran, per Task 1's "assess feasibility before committing" instruction (feasibility research is not evaluation).

**The ~50%-ceiling prediction's outcome, stated as required:** the premise itself needed correcting (see Section 3) — sequence-only's transcription percent-of-ceiling is not a consistent ~50% but ranges 48.7%–102.1% across hosts once computed from verified inputs. Against the corrected numbers, DNABERT-2 and PromoGen2 track a *similar host-by-host pattern* to sequence-only (both roughly track EC~60%, BS~85-95%, PA~62-65% for transcription) rather than uniformly landing above or below it — **my pre-stated prediction (FMs roughly match or underperform sequence-only, medium-high confidence) was directionally right, but the specific "near ~50%" framing it was hung on was not**, because that framing rested on a number that turned out to be wrong. Being wrong about the framing while right about the substantive prediction is reported here exactly as required, not smoothed into "I was right."

**My honest prior (Section 3 of the expectations doc) was: medium-high confidence that FMs would roughly match or underperform sequence-only.** Verdict: **correct**, and correct by a wider margin than the medium-high confidence implied (35 of 37 distinguishable comparisons favor sequence-only, not merely "mostly").

---

## TASK 1 — Evo 2

**Not evaluated. Feasibility fully assessed; execution blocked on an external account this agent cannot create.**

[VERIFIED, primary sources: `github.com/ArcInstitute/evo2` README + issue #67] Evo 2 hard-requires CUDA + Flash Attention at every released size (1b_base, 7b, 7b_base, 7b_262k, 40b, 40b_base); the 1B/20B/40B checkpoints additionally require Transformer Engine + FP8 on a Hopper GPU. It throws `ValueError: Expected a cuda device, but got: cpu` at initialization on non-CUDA hardware. This machine (8GB unified-memory Apple Silicon, MPS only, no CUDA) cannot run any Evo 2 checkpoint locally, at any size. No HuggingFace `transformers`-compatible port exists for Evo 2.

**Hosted-API path — approved by Gabriel (2026-08-07), not yet executed.** NVIDIA Build (`build.nvidia.com`) hosts a free-tier API exposing Evo2-40B's forward-pass layer embeddings via a `/biology/arc/evo2/forward` NIM endpoint — this would be the largest tractable option by a wide margin (40B vs. the 117-148M alternatives actually run), at effectively $0 against the project's <$200 budget. Account creation requires phone-verified signup; this agent has no browser automation tool and cannot complete that signup itself. Gabriel was asked to create the account and provide an API key; as of this memo, the key has not arrived. **This is not a silent gap** — see WHAT I COULD NOT DO. If a key is provided after this memo is delivered, Evo 2 can be slotted into the existing pipeline (`scripts/61_fm_embeddings.py`, `scripts/63_fm_loho_calibration_rs241.py`, `scripts/64_gate6_full_comparison.py` all already support an `evo2` tag / auto-detect its results files) without re-running DNABERT-2 or PromoGen2, and the memo/paper-framing updated as an addendum.

**Per-charter fallback question ("would a cheap cloud instance make it feasible, and at what cost") is answered by the above:** a paid cloud GPU instance is not the cheapest path — the free hosted API is strictly better if authorized, since it requires no compute provisioning at all. No cloud instance has been provisioned, per the charter's explicit "do not provision anything" instruction.

---

## TASK 2 — DNABERT-2 (primary) and PromoGen2 (supplementary)

**NTv3 was assessed and not used** — [INFERRED, medium-high confidence, from InstaDeep's own NTv3 paper/model card] its base corpus (OpenGenome2) does include bacterial sequence, so it is not fully excluded, but its post-training stage fine-tunes on only 24 animal/plant species, its architecture (U-Net, 1Mb context, single-nucleotide eukaryotic regulatory landscapes) is a poor structural match for 165bp bacterial elements, and it requires a non-`transformers`-standard loader with unconfirmed CPU/MPS compatibility. Full reasoning in `out/GATE6_EXPECTATIONS.md` Section 4.

**DNABERT-2** (`zhihan1996/DNABERT-2-117M`, [VERIFIED by direct load] 117,068,544 params, 547× the sequence-only CNN) was selected as the primary Task 2 model — named explicitly as a candidate in the charter's own text, pretrained on 135 species across 6 clades including bacteria (32.5B bases). Loading it required three documented, disclosed compatibility patches against this environment's 2026-era `transformers`/`torch` stack, none of which touch weights or computation (see `scripts/61_fm_embeddings.py` module docstring for the full mechanism of each): (1) bypassing a stale `triton` import pre-check that the model's own code already guards with a try/except; (2) fetching one file (`bert_padding.py`) the pre-check had prevented `transformers` from downloading; (3) forcing two ALiBi-tensor constructors explicitly onto `cpu` because `transformers`' newer meta-device fast-init path left them on `meta` by default, causing a device-mismatch crash unrelated to the model's actual computation.

**PromoGen2** (`jinyuan22/promogen2-base`, [VERIFIED by direct load] 147,877,120 params) was run as a supplementary check — not a charter requirement, but found during the Task 2 search and worth the low marginal cost (loaded cleanly via standard `transformers.AutoModelForCausalLM`, zero compatibility patches needed). It is pretrained *exclusively* on prokaryotic promoter/cis-regulatory sequence (17,000 bacterial genomes, 59M promoters → 1.4M curated training sequences) — the best domain match of any candidate found, directly relevant to the charter's "search rather than assuming" instruction for Task 2.

**Protocol, matching Gate 5 exactly** (`scripts/62_fm_head_model.py`, `scripts/63_fm_loho_calibration_rs241.py`): frozen embeddings (mean-pooled final hidden state, masked by the tokenizer's attention mask) cached once per unique 165bp sequence (29,042 library + 207 RS241, ~3–20 min per model depending on architecture depth), then a small MLP head (`FMHeadMLP`, embedding→64→4 output heads) trained per LOHO fold/N-cell/RS241-seed using the *identical* `compute_loss`/`evaluate` functions from `scripts/40`/`scripts/42` — same metrics (mcc, auc, spearman_rho), same fold semantics, same 90% percentile bootstrap. **One documented protocol difference, not silently matched:** the CNN's calibration curve has two transfer mechanisms (frozen-trunk `head_only` vs. unfrozen-conv4 `top_conv`) because its own trunk is adapted during fine-tuning; a frozen FM embedding has no trunk to unfreeze, so only one mechanism (`head_only`-equivalent) exists and is reported for FM systems at N=100.

**MPS memory bug found and fixed during PromoGen2 extraction:** the first extraction attempt called `model(..., output_hidden_states=True)` and read `hidden_states[-1]`, which materializes all 31 layers' tensors per forward pass; these were not released between batches on MPS's caching allocator and allocation grew unbounded across ~800 batches until it crashed (`MPS backend out of memory... allocated 8.69 GiB`). Fixed by calling `model.transformer(...)` directly (base GPT2Model, no LM head) for only the final hidden state, plus explicit `torch.mps.empty_cache()` per batch and a smaller batch size (32, down from 128). Documented in `scripts/61_fm_embeddings.py`.

---

## TASK 3 — The full comparison

Full table: `out/results/gate6_full_comparison.{json,csv}` (84 rows: 3 hosts × 2 readouts × 6 systems × {zero-shot, N100 head_only, N100 top_conv where applicable}). RS241: `out/results/gate6_rs241_comparison.{json,csv}`. Figures: `out/figures/gate6_full_comparison.png` (all 6 systems, zero-shot, all host×readout cells), `out/figures/gate6_pct_of_ceiling_corrected.png`.

### 1. Does any foundation model beat the sequence-only model? Per configuration and overall.

**No, not overall, and only in one specific cell per configuration.**

| Comparison set | Total | Distinguishable (90% CI) | FM wins | Sequence-only wins |
|---|---|---|---|---|
| Primary hosts (zero-shot + N100) | 24 | 14 | 2 | 12 |
| RS241 (zero-shot, both configs) | 24 | 23 | 0 | 23 |
| **Combined** | **48** | **37** | **2** | **35** |

Both FM wins are the same cell for both models: *E. coli* transcription, zero-shot (DNABERT-2 ρ=0.458 vs. sequence-only ρ=0.367; PromoGen2 ρ=0.495 vs. 0.367). Every other distinguishable comparison — every N=100 cell without exception, every RS241 cell but one (SE transcription PRIMARY vs. DNABERT-2, not statistically separable) — favors sequence-only, often by a wide and clearly separated margin. At N=100 specifically, sequence-only's head_only mechanism (its trunk frozen too, the fairest direct comparison to a frozen-FM-embedding head) beats both FMs in all 8 of 8 distinguishable comparisons.

### 2. Where do the FMs land relative to the ~50%-ceiling prediction?

See the correction above. Using the corrected, verified percent-of-ceiling metric (`out/results/gate5_5_per_host_ceiling_CORRECTED.json`):

| Host | Readout | Sequence-only | DNABERT-2 | PromoGen2 |
|---|---|---|---|---|
| EC | transcription | 48.7% | 60.8% | 65.7% |
| BS | transcription | 102.1% | 84.2% | 94.2% |
| PA | transcription | 63.6% | 61.6% | 65.0% |
| EC | translation | 75.8% | 67.3% | 70.3% |
| BS | translation | 110.6% | 104.0% | 101.1% |
| PA | translation | 52.5% | 39.8% | 45.1% |

FMs track the same host-by-host shape as sequence-only (BS highest, EC/PA lower) rather than uniformly exceeding or undershooting it — consistent with percent-of-ceiling being driven mostly by properties of the *data pair* (BS's low-N, low-correlation ceiling specifically) rather than by which model computes the numerator. Neither FM systematically clears 100% the way sequence-only does on BS; where FMs are below sequence-only, they are usually below on the same hosts sequence-only is also (relatively) weak on.

### 3. Does the central claim survive?

**Yes, plainly, and by a wide margin — see GATE DECISION above.** Per the pre-committed Section 2 threshold in `out/GATE6_EXPECTATIONS.md`, the claim would have narrowed only if a foundation model distinguishably beat sequence-only on ≥3 of 6 primary host×readout zero-shot cells; each model reached exactly 1 of 6, and that cell is a genomic-vs-species idiosyncrasy (*E. coli* specifically, transcription specifically) rather than a broad pattern. No conditioned-FM variant was tested here (Task 1/2 specify frozen embeddings + head, structurally parallel to the sequence-only ablation — no host information reaches either FM), so Gate 6 does not directly re-test whether conditioning helps; it tests whether capacity/pretraining was the small model's bottleneck, and the answer is no — bigger, pretrained models did worse, not merely differently. The negative H-MAIN result from Gate 5 is not attributable to the demonstration model being underpowered.

---

## TASK 4 — Framing update

`out/PAPER_FRAMING.md` updated in place (see file for full diff): a new ranked supporting-result entry for the Gate 6 foundation-model comparison, the bacterial-coverage-gap finding folded into the "what this result establishes" scope section, and the corrected percent-of-ceiling figures replacing the retracted ones everywhere they appeared.

---

## Bacterial coverage gap in genomic foundation models — a reportable finding in itself

Per Task 2's instruction ("if NTv3 is inapplicable, that is a reportable finding in itself"): the field's most prominent recent genomic foundation models split into two groups relevant here. **NTv3's *base* pretraining corpus does include bacterial sequence** (OpenGenome2, shared with Evo 2) — so the strong "excludes bacteria entirely" framing would overstate it — but its architecture and post-training are eukaryote-specialized, and no other major, easily-loadable, `transformers`-standard genomic FM was found that was built *for* bacterial regulatory sequence specifically, except PromoGen2 (a narrow, promoter-focused model, not a general genomic FM in the Evo2/NT sense). This is consistent with, and adds direct evidence to, the charter's own framing: a cross-host *bacterial* regulatory-activity benchmark fills a real gap that the largest available genomic FMs were not built to address, independent of whether those FMs are good at the task once applied to it.

---

## TEMPTATIONS TO ADJUST

1. **Whether to quietly patch `gate5_5_per_host_ceiling.json` in place rather than document the discrepancy prominently.** Not done — the correction is dated, explained, and the old file is left in place (not deleted) with the corrected file alongside it under an explicitly different name, so the record of what changed and why is preserved rather than erased.
2. **Whether to drop PromoGen2 from the reported comparison since it wasn't charter-required and complicates a clean "one alternative model" narrative.** Not done — it is reported as a supplementary check throughout, never substituted for DNABERT-2 as the primary Task 2 answer, and its numbers are shown alongside DNABERT-2's everywhere rather than cherry-picked into or out of the headline table.
3. **Whether to keep waiting indefinitely for the Evo 2 API key before delivering this memo.** Not done — the charter's FAIL condition ("no foundation model can be evaluated by any route") was already cleared by DNABERT-2 and PromoGen2; holding the entire gate hostage to an external, uncontrollable dependency would create exactly the kind of open-ended blocking the charter's gate structure is designed to avoid. Evo 2 remains addable as an addendum.
4. **Whether to reframe Task 0's prediction after the fact to match the correction more comfortably** (e.g., quietly dropping the "~50%" framing and only ever presenting the corrected numbers as if that had been the plan all along). Not done — Section 1 of `out/GATE6_EXPECTATIONS.md` is left as originally written, with the correction appended as a dated addendum, so the actual sequence of what was predicted, when, and what turned out to be wrong about the premise is auditable.

---

## WHAT I COULD NOT DO

1. **Evo 2 was not evaluated.** Local execution is hard-blocked by a CUDA requirement this environment cannot satisfy (verified via primary sources, not merely inferred). The hosted-API fallback was approved but requires an external account this agent cannot create (no browser automation tool available) and Gabriel has not yet supplied a key. This is the gate's single largest gap.
2. **The free-tier NVIDIA Build credit ceiling was never verified against NVIDIA's own documentation** (secondary sources report ~1,000 inference credits; NVIDIA's own docs do not state a number). If a key arrives, the actual limit — and whether the sequence-packing workaround discussed in `out/GATE6_EXPECTATIONS.md` is even necessary — will need to be pilot-tested empirically before committing to a full-scale extraction plan.
3. **NTv3 was not actually loaded or run** — the decision not to pursue it rests on documentation/model-card research (medium-high confidence), not a direct local-loading attempt the way DNABERT-2's and PromoGen2's feasibility was confirmed. If this project later needs a stronger "NTv3 specifically fails" claim (rather than "NTv3 was a poor prospective fit"), an actual load attempt would be needed.
4. **No other candidate genomic FMs (Carbon-500M, the original NT 2.5B/500M multi-species checkpoints, M5, ARSENAL) were attempted.** DNABERT-2 + PromoGen2 were judged sufficient to answer Task 3's questions decisively; the marginal value of a third or fourth FM given how consistent the DNABERT-2/PromoGen2 pattern already is was judged low relative to the time cost, but this is a scope choice, not a completeness claim.
5. **The `gate5_5_per_host_ceiling.json` correction was not traced to its exact root cause.** I know its `actual_zeroshot_mean` values don't match the officially bootstrapped source and that no script for it exists to inspect; I do not know what specific computation error produced those particular wrong numbers, because there is no code left to examine. This is disclosed as a real gap in the correction's forensic completeness, not just a footnote.

---

## CONTRADICTIONS WITH THE CHARTER

1. **The charter's Gate 6 framing ("here is a 250k-parameter host-conditioned model against a 40B-parameter foundation model") assumed Evo 2 would be the headline comparison.** It is not evaluated in this memo. The comparison that *is* delivered (117–148M-parameter FMs) still answers the charter's underlying question (does capacity explain the small model's ceiling) directionally and decisively, but the specific "250k vs. 40B" contrast the charter anticipated as the most quotable result is not yet available.
2. **The charter's Part V rule 6 ("verify before trusting a stated number") was violated somewhere upstream of this gate** — `gate5_5_per_host_ceiling.json` was cited in `GATE5_5_MEMO.md` and implicitly shaped `PAPER_FRAMING.md`'s "remarkably consistent ~50%" language without its own inputs having been checked against the project's other, more carefully bootstrapped results. This gate's own Task 0 prediction inherited the error before catching it. The catch happened here, not earlier, and is now fixed and disclosed — but it is a real instance of the exact failure mode the standing rule exists to prevent, not merely a hypothetical risk.
3. **This gate's own protocol deviates from the charter's literal Task 1 request** ("Extract Evo 2 embeddings... Report your figures against those [Evo 1] as reference points") by substituting DNABERT-2 and PromoGen2 where Evo 2 was specified as the primary target. This is disclosed as a scope substitution driven by hard infeasibility, not a preference, and is reversible if Evo 2 becomes available.

---

## DATED ADDENDUM — 2026-08-08, Gate 7: this memo's own ceiling "correction" was itself wrong

Gate 6 (above, Task 3 point 2 and CONTRADICTIONS #2) reported that `out/results/gate5_5_per_host_ceiling.json` was unreproducible, diffed it against the sequence-only model's officially bootstrapped zero-shot rho, found a mismatch in 5 of 6 cells, and published `gate5_5_per_host_ceiling_CORRECTED.json` as a fix.

**That diff compared the wrong baseline.** `gate5_5_per_host_ceiling.json` was never about the sequence-only model — its own column header (`out/GATE5_5_MEMO.md` Task 2.2) explicitly says "Actual zero-shot **(genomic)**". Re-verified in Gate 7 to 6 decimal places against `out/gate4_loho_results.json`: the original file's numbers are exactly the genomic host-conditioned CNN's zero-shot Spearman rho. It was correct and reproducible all along; Gate 6 introduced a new error by not reading the label and substituting a different model's numbers in its place.

**`gate5_5_per_host_ceiling_CORRECTED.json` (Gate 6's output) is therefore superseded, not `gate5_5_per_host_ceiling.json`.** The sequence-only, DNABERT-2, and PromoGen2 percent-of-ceiling numbers Gate 6 computed are not themselves wrong as *numbers* (same bootstrap methodology, correctly computed) — they were mischaracterized as "corrections to a bug" when they were actually new computations for models the original file never covered.

**This does not change Gate 6's actual verdict.** The central FM-vs-sequence-only comparison (35/37 distinguishable comparisons favor sequence-only; neither FM reaches the pre-committed zero-shot threshold) was computed via direct bootstrap CI comparison (`scripts/64`'s `build_comparison_table`/verdict logic), entirely independent of the ceiling metric — the ceiling analysis was supplementary context, not load-bearing for the gate's pass/fail decision.

**The percent-of-ceiling metric is retired in Gate 7 for an unrelated, genuine reason**: the sequence-only model (not a bug, a real result) exceeds the *correct, original* genomic-based ceiling for *B. subtilis*, both readouts. Full derivation, tests, and verdict in `out/GATE7_MEMO.md` Task 1 and `out/GATE5_5_MEMO.md`'s own dated correction. Both corrections are reported here together because they were found and resolved in the same Gate 7 investigation, and separating them would obscure that one was a self-caught process error (wrong baseline) and the other a genuine scientific finding (a real model exceeding a flawed ceiling concept).

---

## FILES WRITTEN

- `out/GATE6_EXPECTATIONS.md` (Task 0, plus dated correction addendum)
- `out/GATE6_MEMO.md` (this file)
- `data/fm_embeddings/{dnabert2,promogen2}_{library,rs241}.npz`
- `scripts/61_fm_embeddings.py`, `scripts/62_fm_head_model.py`, `scripts/63_fm_loho_calibration_rs241.py`, `scripts/64_gate6_full_comparison.py`, `scripts/65_gate6_figures.py`
- `out/gate6_{dnabert2,promogen2}_{loho_results,calibration_curves,rs241_results}.json`
- `out/results/gate6_full_comparison.{json,csv}`, `out/results/gate6_fm_verdict_{table.csv,summary.json}`
- `out/results/gate6_rs241_comparison.{json,csv}`, `out/results/gate6_rs241_fm_verdict_{table.csv,summary.json}`
- `out/results/gate5_5_per_host_ceiling_CORRECTED.json`
- `out/figures/gate6_full_comparison.png`, `out/figures/gate6_pct_of_ceiling_corrected.png`
- `out/PAPER_FRAMING.md` (updated)
- `out/state.json` (`gate_6` block)
