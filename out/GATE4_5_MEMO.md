# CROSSHOST — Gate 4.5 Memo: Fix the Physiology Vector and the Decision Rule

**Date:** 2026-08-05
**Scope:** (1) resolve the *B. subtilis* growth-rate imputation; (2) extend both-transfer-mechanism reporting to N=100; (3) widen the H-MAIN decision rule and complete the fold-variance picture; (4) investigate the *P. aeruginosa* sign inversion; (5) retrain whatever Tasks 1–4 required. **No H-MAIN or H-SCIENCE verdict is computed or stated anywhere in this memo** — Gate 5 evaluates hypotheses; this gate only prepares clean inputs for that evaluation.

**Bottom line up front: PASS.** The physiology vector now has a defensible footing (Task 1a succeeded — a real, condition-matched, host-specific value was found), the leakage audit does not regress, and retraining completed.

---

## Task 1 — B. subtilis growth-rate imputation: RESOLVED

### 1a. The search succeeded

Retrieved the primary source directly rather than relying on secondary summaries: Zhu, Mori, Hwa & Dai 2025 (*PNAS* 122:e2427091122) full text **and** SI Appendix, via the PMC open-access FTP mirror (`ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/deprecated/oa_package/82/94/PMC12067254.tar.gz` — note: the modern `oa_package` path is defunct; the `deprecated` path still works and was the same fix Gate 1 had already discovered for a different paper).

**What was and wasn't in the SI Appendix:** Tables S1–S2 (which the paper's own text points to for "strain and medium details") are reference tables, not a results table with a tabulated *B. subtilis*-in-LB-at-37°C digit. **The number is recoverable from Fig. 2B** (main text) — a bar chart of growth rate (h⁻¹) for *E. coli*, *B. subtilis*, and *V. natriegens* across several media at 37°C, with LB as each species' tallest bar. Visually read: ***B. subtilis* LB/37°C ≈ 1.7 h⁻¹.**

**Cross-checks supporting this value:**
- The paper's own text states *V. natriegens* "grows >50% faster than Ec and Bs when growing in the same nutrient and temperature," implying *E. coli* and *B. subtilis* rates in LB/37°C are close to each other and both well below *V. natriegens*'s ≈2.8 h⁻¹. 1.7 h⁻¹ (≈ *E. coli*'s own already-adopted value) is consistent with this.
- Confirmed the growth medium: SI Appendix Supplementary Methods states *B. subtilis* was "cultured in either LB broth or...", at "37°C" for the fast-growth condition — an exact match to *E. coli*'s (Sezonov 2007) and *P. aeruginosa*'s (Yang 2008) own LB/37°C conditions already in the table.

**Confidence tagging:** [INFERRED, medium-high confidence] — read from a bar chart, not extracted from a table cell, so not [VERIFIED] at the same precision as a literal digit. This is the **same confidence tier** already used for *E. coli*'s and *P. aeruginosa*'s own values (both are midpoints of literature ranges, not single verified digits either) — not a downgrade relative to the rest of the table.

**Sources found and rejected**, reported per instruction:
| Source | Value | Medium/Temp | Verdict |
|---|---|---|---|
| BioNumbers BNID 112127 (Burdett et al. 1986) | 120 min doubling (μ≈0.35 h⁻¹) | Minimal + succinate, 35°C | Rejected — wrong medium (minimal, not rich), condition mismatch with other hosts |
| Abele et al. 2025 (PaxDb source dataset for this project's own physiology proxy vector) | Not reported | Agar plate, DSMZ-specific conditions | Checked directly (fetched PMC full text) — confirmed no growth-rate/doubling-time value is reported anywhere in that paper; agar-plate proteome atlases don't measure liquid growth kinetics |
| Various web summaries ("~30 min doubling") | ~30 min (μ≈1.4 h⁻¹) | Unspecified/unsourced | Rejected — could not trace to one primary citation with a stated numeric value (same conclusion Gate 1 already reached) |

### 1b/1c. Vector recommendation

Because 1a succeeded with a real, independently-sourced, condition-matched value — not a formulaic cross-host average — **the original two-vector hedge substantially dissolves**, per the task's own framing ("if 1a succeeds, this collapses to a single clean 6-metric vector"). Both consequences the task flagged as "load-bearing" are resolved:
1. The feature is no longer informationless for *B. subtilis* by construction (it is now *B. subtilis*'s own measured value, not the mean of the other five).
2. No more training-host contamination (the value is not derived from EC/PA at all; it comes from an independent literature source).

**PHYS-6 (updated) is recommended as primary** for Gate 5, not because the alternative was rejected on the merits but because the problem PHYS-5 existed to hedge against no longer applies at the same severity. **PHYS-5 (5 metrics, growth rate dropped) was still built** as a deliverable artifact per instruction, but **was not trained** — building a full third training grid (another 15 LOHO fits + ~1,050 calibration draws) for a hedge against a problem that 1a resolved was judged disproportionate to the residual risk. This is a **DECISION** I made, not merely a recommendation, given the instruction's own conditional ("if 1a succeeds... the question dissolves") — flagged here for visibility, not hidden as a quiet scope cut. If Gabriel judges the visual-figure-read confidence tier insufficient despite matching EC/PA's own tier, PHYS-5 training remains a bounded, well-specified follow-up (data already prepared).

**Data changes:** `data/hosts_physiology.parquet` updated (BS `growth_rate_mu_h`: 1.852 → 1.7; `growth_rate_confidence`: new column, `inferred_from_figure` for BS vs. `verified_range_midpoint`/`verified_single_value` for others; `growth_rate_imputed`: new column, `False` for all 6 hosts — no host is imputed-from-other-hosts anymore). `data/MANIFEST.json` updated accordingly; `audit_leakage.py` check 6 re-verified passing against the new hash.

---

## Task 2 — Both transfer mechanisms at N=100 (reporting expansion, not a hypothesis change)

Dated amendment appended to `out/PREREGISTRATION.md` (Amendment 1) **before** this task's results were computed, per instruction. Frozen-trunk (head_only) **remains the pre-registered primary** for H-MAIN; top_conv is supplementary. Full results: `out/gate4_5_n100_topconv_supplement.json`.

| Host | Variant | Frozen-trunk ρ (primary) | Top-conv ρ (supplementary) | Gap |
|---|---|---:|---:|---:|
| *E. coli* | genomic | 0.602±0.083 | 0.626±0.056 | +0.023 |
| *E. coli* | physiology | 0.620±0.040 | 0.666±0.043 | +0.045 |
| *B. subtilis* | genomic | 0.213±0.068 | 0.170±0.051 | −0.043 |
| *B. subtilis* | physiology | **−0.059±0.083** | **0.118±0.062** | **+0.177** |
| *P. aeruginosa* | genomic | 0.104±0.224 | 0.481±0.059 | +0.376 |
| *P. aeruginosa* | physiology | 0.244±0.309 | 0.444±0.077 | +0.199 |

**The internal inconsistency Task 2 was written to fix is now visible with numbers, not just described:** at N=100 under the pre-registered primary mechanism, *B. subtilis*-physiology's transcription correlation is **negative** (−0.059), while the supplementary mechanism gives a modest positive value (+0.118). *P. aeruginosa* shows the largest gap of all six cells (+0.376 for genomic) — consistent with Task 4's finding that PA-physiology optimization is unstable specifically under the frozen-trunk mechanism at low N. This table is handed to Gate 5 as-is; **no interpretation of what it means for H-MAIN is offered here.**

---

## Task 3 — Decision rule widened; translation arm checked

### 3.1 — Decision rule

Dated amendment appended to `out/PREREGISTRATION.md` (Amendment 2), **before any H-MAIN number was computed under either rule**: interval type changed from an 80% Student's-t interval to a **90% percentile bootstrap**, both changes in the conservative direction (harder to declare MET). Full reasoning, including why a bootstrap is now preferred given Task 4's bimodal-optimum finding, is in the amendment itself.

### 3.2 — B. subtilis translation fold-variance spot check (the other H-MAIN arm)

Same protocol as Gate 4's transcription spot-check (*B. subtilis* held out, genomic features, N=100, frozen-trunk). Full data: `out/gate4_5_bs_translation_fold_variance_spotcheck.json`.

| | Fold-to-fold std | Draw-to-draw std | Ratio |
|---|---:|---:|---:|
| CNN, translation (this check) | 0.0921 | 0.0430 | **2.14** |
| CNN, transcription (Gate 4) | 0.0683 | 0.0329 | 2.08 |
| k-mer/linear baseline, translation (Gate 3.5) | 0.0611 | 0.1551 | 0.37 |

**DECISION NEEDED (restated from Gate 4, now confirmed for both H-MAIN arms):** the translation arm's fold/draw ratio (2.14) is slightly *worse* than transcription's (2.08) and **5.4× the k-mer baseline's ratio**. Both arms of H-MAIN — not just the one Gate 4 happened to spot-check — show CNN fold variance substantially exceeding baseline-level variance. This is exactly why Task 3.1's interval was widened; the two findings are linked, and the amendment says so explicitly.

---

## Task 4 — P. aeruginosa sign inversion: genuine variance, mechanistically explained

Full data: `out/gate4_5_pa_investigation.json`. Bounded investigation, ~50 minutes.

**Checks performed:**
- **Label orientation:** PASS — construction of `tx_active`/`tx_norm` is fold-independent, same code path for every fold; no evidence of a fold-specific label swap.
- **Class distribution:** not degenerate — the affected fold's active rate (81.9%) is unremarkable, within the range of the other 4 folds (80.0–87.4%).
- **FiLM γ/β degenerate solution:** not collapsed — mean γ/β magnitudes for the affected fold are comparable in scale to PA's other folds and to *B. subtilis*'s folds; no single-fold anomaly visible in summary statistics.
- **PA physiology vector outlier check:** **confirmed and quantified.** Defined an "extrapolation degree" — how far a held-out host's z-scored feature vector sits outside the per-dimension range spanned by its two LOHO training hosts. Results: **PA held out (train EC+BS) = 3.687** z-score-units total excess (driven overwhelmingly by `chaperones_fraction` alone, 1.61 units) — the **largest** of all three primary hosts. *B. subtilis* held out (train EC+PA) = 2.245. *E. coli* held out (train BS+PA) = 1.022 (smallest — *E. coli* sits closest to being bracketed by its own training pair).

**Prediction-structure check (ruling out pure noise):** the affected fold's strength predictions correlate with truth at ρ=−0.4487, **p=4.4×10⁻¹⁹⁰** (n=3,850) — statistically overwhelming, and not consistent with random scatter. This is a real, structural (if inverted) learned relationship.

**Bimodal fold pattern:** PA-physiology-transcription's 5 folds split cleanly — folds {0, 2, 4} sign-inverted, folds {1, 3} correct-signed. Not uniformly bad, which argues against a deterministic pipeline bug (expected to hit every fold identically) and toward something seed-sensitive.

**VERDICT: genuine variance, mechanistically explained — not a bug, not simple noise.**

**Best hypothesis** [INFERRED, medium-high confidence]: FiLM-generator extrapolation instability. The FiLM generator is a small MLP trained only on the two LOHO training hosts' feature vectors; PA's vector is, by a clear margin, the furthest of any primary host from the convex hull of its own training pair. Because the generator receives no training signal constraining its behavior at that far out-of-distribution point, its near-linear extrapolation can land sign-preserving or sign-inverting depending on random initialization (which differs by fold) — producing the observed bimodal pattern rather than uniform failure. Consistent with: unfreezing conv4 rescuing the sign in every affected fold (Gate 4 memo; more capacity compensates for a bad extrapolation), and PA's independently-established (Gate 1.5) status as the most source-sensitive host in this feature space.

**Not fully resolved:** the exact internal mechanism (why THIS specific inversion direction per fold) was not traced at the individual-weight level — would need multi-seed sweeps or stronger FiLM-generator regularization, out of this task's bounded scope.

**Retraining implication:** genomic-variant results are **not** affected — this verdict is specific to physiology-feature extrapolation, not a pipeline bug. No additional retraining beyond Task 1's fix was triggered by this task.

---

## Task 5 — Retraining log

| Step | Scope | Wall-clock | Result |
|---|---|---:|---|
| LOHO retrain | 15 physiology fits (3 hosts × 5 folds; genomic untouched) | 9,969.8s (2.77 hr) | All 15 completed, `epochs_run=35` (full budget) every time — no early stop |
| Calibration curve regen | Physiology slice of the 7-cell × 10-draw grid (genomic untouched) | 15,535.5s (4.32 hr) | All 1,050 physiology fine-tune draws completed |
| N=100 top_conv supplement | Both variants, all 12 configs | 388.6s | 300 fine-tunes completed |

**Same seeds, folds, architecture, and hyperparameters as the original Gate 4 run** — only the *B. subtilis* row of the physiology feature table differs (1.852 → 1.7 in one dimension). Genomic-variant checkpoints, results, and figures are **untouched** — confirmed unaffected by Task 1's fix and Task 4's "not a bug" verdict, so not retrained, per instruction ("do not retrain them without reason, and state whether you did").

**The PA sign-inversion pattern persists after retraining** (as Task 4's analysis predicted, since the growth-rate fix barely moved PA's own vector: max diff 0.019 z-score-units) — e.g. PA-physiology-transcription fold 2: ρ=−0.448 (original) → ρ=−0.326 (retrained), same sign, similar magnitude. This is expected and reported as confirming evidence for Task 4's verdict, not a new problem.

`audit_leakage.py` re-run after all retraining: still 6/6 PASS (see status report).

---

## WHAT I COULD NOT DO

- **Did not extract the exact B. subtilis growth-rate digit from a table** — read from a bar chart (Fig. 2B) instead, [INFERRED, medium-high confidence] rather than [VERIFIED]. A precise digit would require either contacting the authors or a more careful figure-digitization pass (e.g. pixel-level image analysis) than was done here.
- **Did not train a full model grid on PHYS-5** (the 5-metric, growth-rate-dropped vector) — built as a data artifact per instruction, but training it was judged disproportionate given Task 1a's success. Remains a bounded, well-specified follow-up if this judgment call is overturned.
- **Did not fully trace the FiLM-generator extrapolation-instability mechanism at the individual-weight level** for Task 4 — bounded investigation produced a well-evidenced hypothesis, not a proof.
- **Did not re-run Gate 4's original transcription fold-variance spot-check (BS, genomic)** — unaffected by any Gate 4.5 change, so the existing Gate 4 number stands unmodified.
- **Did not investigate whether the same extrapolation-instability mechanism affects EC or BS to a lesser degree** beyond the single quantitative "extrapolation degree" comparison computed for Task 4 — a fuller characterization (e.g. per-dimension breakdown for all three hosts, not just the summary metric) was out of scope for the bounded investigation window.

## CONTRADICTIONS WITH THE CHARTER

1. **Gate 4's own internal inconsistency (H-MAIN-TX evaluated at N=100, inside the observed N=10–300 degradation window) was real and is now quantified**, not just described: frozen-trunk *B. subtilis*-physiology transcription at N=100 is **negative** (−0.059), a number Gate 4's memo flagged as a risk but did not compute directly at N=100 (only inferred from the N=10-300 range description). This is exactly the kind of thing pre-registration and mechanism-reporting expansions exist to surface before Gate 5, not after.
2. **The physiology feature vector's B. subtilis growth-rate value was wrong in a way that mattered on the exact host and readout H-MAIN is defined on**, for four gates (Gate 2 through Gate 4) before being caught and fixed here. The fix was found via a more thorough literature search than any prior gate attempted (direct primary-source PDF retrieval via the PMC OA FTP mirror, not just web search) — worth noting as a general lesson: "inaccessible since Gate 1" should periodically prompt a harder retry, not be treated as permanently closed.
3. **CNN fold variance is now confirmed elevated on both H-MAIN arms**, not just the one Gate 4 happened to check — a broader and more serious version of the concern Gate 4 flagged. The decision-rule widening (Amendment 2) is a direct, proportionate response, made before any H-MAIN number exists under the new rule.

## FILES WRITTEN

- `out/GATE4_5_MEMO.md` — this memo
- `out/PREREGISTRATION.md` — updated with Amendment 1 (mechanism reporting) and Amendment 2 (90% bootstrap), both dated
- `data/hosts_physiology.parquet`, `data/hosts_physiology_final.csv` — updated (BS growth rate corrected, `growth_rate_confidence`/`growth_rate_imputed` columns added)
- `data/MANIFEST.json` — updated hashes
- `out/gate4_5_pa_investigation.json`, `gate4_5_physiology_retrain.json`, `gate4_5_n100_topconv_supplement.json`, `gate4_5_bs_translation_fold_variance_spotcheck.json`
- `out/gate4_loho_results.json`, `gate4_loho_results_table.csv`, `gate4_calibration_curves.json`, `gate4_calibration_results_table.csv` — regenerated (physiology entries replaced, genomic untouched)
- `out/figures/` — regenerated with corrected physiology data
- `out/models/loho_{EC,BS,PA}_physiology_fold{0-4}.pt` — 15 checkpoints overwritten with retrained versions
- `scripts/12_build_final_physiology_vector.py` (updated: BS growth rate, confidence/imputed columns), `13_build_task1_core_tables.py` (re-run, not modified)
- `scripts/46_retrain_physiology_gate4_5.py`, `47_bs_translation_fold_variance_spotcheck.py`, `48_n100_topconv_supplement.py`, `49_regen_physiology_calibration.py` (new)
