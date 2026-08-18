# CROSSHOST — Gate 8 Memo: Attenuation Check and Packaging

**Date:** 2026-08-14. **Outcome: PASS.** Task 1's verdict, stated first as required: **the central mechanism claim SURVIVES INTACT.** Measurement-noise attenuation does not explain a material share of the EC–PA (ρ≈0.75) vs. *B. subtilis*-pairs (ρ≈0.16–0.26) gap. At the one empirically-grounded reliability estimate available (*E. coli* transcription, 0.91, from five independent growth-condition measurements), the corrected gap ratio is **identical** to the observed one (0.341 both before and after correction). Even under deliberately aggressive assumptions favoring the alternative hypothesis (*B. subtilis* reliability as low as 0.5), corrected *B. subtilis*-pair correlations never exceed ~53% of EC–PA's. The >100%-of-ceiling anomaly Gate 7 left open **is resolved**: at every tested reliability, the corrected ceiling meets or exceeds the sequence-only model's actual performance.

---

## TASK 1 — Measurement-reliability attenuation check

### What replicate structure exists, and what does not

Gate 1 already established the manuscript-described barcode-replicate QC set is not reconstructable from released tables (`CHARTER_AMENDMENTS.md` C1) — closed for all three hosts. What Gate 8 found instead: `raw/NIHMS945382-supplement-5.xlsx` ("Supplementary Data Table 3") contains **E. coli transcription measured across five growth conditions** (LB_exp, NaCl_exp, Fe_exp, LB-stat, M9-exp) for the same 29,249 sequences — the same RNA_count/DNA_count/tx_norm structure as the main library, usability filter applied identically (`dna_count!=0 & (rna+dna)>=15`, active = usable & rna_count>0). **This is EC-transcription-only** — no equivalent series exists for *B. subtilis*, *P. aeruginosa*, or translation in any host (verified: the "Robust RSs" sheet's `Protein (a.u.)` column is a single value per sequence, not a repeated-measures series).

**Methodological caveat, stated before any number was computed and carried through the whole analysis:** growth conditions are not technical replicates. Some cross-condition variation is genuine condition-dependent regulation, not measurement noise. Treating cross-condition correlation as a reliability proxy therefore *bundles real biology into what classical test theory calls noise* — this **underestimates** true single-condition reliability, which means any disattenuation computed from it is an **upper bound** on how much of the cross-host gap measurement error alone could explain. Every number below should be read with this bias in mind — a bias that, if anything, favors finding *more* gap-closing than is really justified, making the "gap survives" result more robust, not less.

### Reliability estimate

Mean pairwise Spearman correlation across all 10 condition pairs (n=10,480–14,372 each): **0.912** (range 0.883–0.938). Under classical test theory with equal per-condition reliability, this mean *is* the reliability estimate directly. Cross-checked by a second, independent method (intraclass correlation on the "Robust RSs" 100-sequence subset's log2 mean activity / stdev): **0.929** — consistent within 0.02, despite that subset being selection-biased toward low variance (expected to run high, and does, but not wildly so).

### Disattenuation and the gap-survival check

Standard Spearman correction: `corrected = observed / sqrt(reliability_A × reliability_B)`, clipped at 1.0. Two scenarios computed (full grid: `out/results/gate8_attenuation_analysis.{json,csv}`): (1) *E. coli* uses its own empirical 0.91–0.93 estimate, the other host in each pair swept across a 0.5/0.7/0.9 sensitivity grid; (2) fully symmetric grid, both hosts in every pair at the same assumed reliability.

| | Transcription, reliability=0.9 (closest to the real EC estimate) | | |
|---|---|---|---|
| | EC–PA corrected | *B. subtilis*-pairs corrected | Ratio (BS/EC–PA) |
| Observed | 0.754 | 0.257 | 0.341 |
| Corrected | 0.838 | 0.286 | **0.341** |

The ratio is unchanged to three decimal places. Translation shows the same pattern (0.286 observed → 0.286 corrected at reliability=0.9). At the most aggressive tested assumption (reliability=0.5, applied symmetrically, favoring gap closure as much as this analysis allows), the ratio moves to at most ~0.52 — real narrowing, but *B. subtilis*-pairs still sit at roughly half of EC–PA, nowhere near parity. **A clipping artifact is worth naming explicitly**: at low assumed reliability, EC–PA's corrected value hits the mathematical ceiling of 1.0 (implying the "true" correlation would need to exceed the maximum possible value) — this is itself evidence that low reliability values (0.5, 0.7) are not self-consistent for EC/PA specifically, further supporting that the true reliability is closer to the empirically-grounded ~0.9 than to the pessimistic end of the sensitivity range.

### Does this explain the >100%-of-ceiling anomaly? Yes, at every tested reliability.

| Readout | Ceiling pair | Observed ceiling | Sequence-only model rho | Corrected ceiling @ 0.5 / 0.7 / 0.9 |
|---|---|---|---|---|
| Transcription | EC–BS | 0.258 | 0.263 (exceeds) | 0.515 / 0.368 / **0.286** (all ≥ model) |
| Translation | BS–PA | 0.263 | 0.291 (exceeds) | 0.527 / 0.376 / **0.293** (all ≥ model, narrowly at 0.9) |

At every reliability level tested, the corrected ceiling meets or exceeds the model's actual performance — the anomaly that looked like a model "beating a hard limit" is exactly what measurement-noise attenuation predicts. Figure: `out/figures/gate8_ceiling_anomaly_resolution.png`.

### Verdict: SURVIVES INTACT

Per the pre-specified three-way rubric: not "weakened" (the gap is materially unchanged at the empirically-grounded estimate, and stays substantial under every tested alternative), and certainly not "does not survive" (BS-pairs never approach EC–PA under any tested scenario). **The central mechanism claim survives intact.** `PAPER_FRAMING.md` finding 1 updated accordingly, with the full disattenuation result folded in as supporting evidence rather than a separate claim.

---

## TASK 2 — Packaging

Built in `package/`, a self-contained, installable, tested directory.

**2a. Package.** `pip install -e .` — **actually tested** in a fresh, isolated venv (`python3 -m venv`, then a separate Python process, different working directory) — not just imported in-place. 9 pytest tests, all passing: data loaders, `evaluate()` on synthetic data (random predictions correctly score near-zero rho, near-0.5 AUC), RS241 leakage check (both a real-data pass case and an injected-leak failure case, both verified), held-out-inputs label-absence check.

**2b. Held-out evaluation.** RS241's labels have already been published throughout this project (Gates 5, 5.5, 6, 7 all cite RS241 zero-shot numbers) — it cannot serve as a genuine withheld set. Fold 4 (5,813 rows) is designated instead: public inputs in `package/held_out_eval/test_inputs.parquet` (no labels), true labels in `out/gate8_held_out_eval_private_labels.parquet` (outside the distributable package, never to be published). `score_submission.py` — **tested end-to-end** with a synthetic random-prediction submission, correctly scored near-zero rho across all hosts/readouts. **Disclosed honestly, not hidden**: this is not cryptographically enforced — `data/core/three_host.parquet` (the public dev table) still contains fold 4 with real labels, since it's the same table used for development. See `KNOWN_ISSUES.md` item 9.

**2c. Baseline scores.** All 9 systems (B1 mean/majority, B2 per-host-N, B4 biophysical, free-embedding B3, genomic, physiology, sequence-only, DNABERT-2, PromoGen2), fold-resolved, identical 90% bootstrap protocol throughout, one consolidated `master_baselines.csv`/`.json` (108 rows) plus Gate 7's conformal/ECE results shipped prominently in `baselines/README.md`'s own headline finding, not an appendix.

**2d. Licensing quarantine.** Structural, not a README note: `data/core/` (MIT + attribution), `data/licensed/dnabert2_derived/` (Apache-2.0, **verified against the source GitHub repo's license API**, not assumed from the HF card which carries no license tag), `data/licensed/promogen2_derived/` (CC-BY-NC-4.0, non-commercial), `data/licensed/deepcross/` (empty, pre-provisioned per charter, disclosed as never-populated rather than silently omitted). Raw Johns et al. supplementary tables are **not redistributed** — `data/core/CITATION.md` documents provenance and points to the publisher.

**2e. Reproducibility.** Root-level `Makefile` (`audit`, `test`, `package`, `reproduce`, `reproduce-full`). **Precisely what was tested, not overclaimed:** `make audit`, `make test`, `make reproduce`, and `make package` were **all actually run this session and passed** — pasted output in this memo's status report. `make reproduce-full` (the complete from-raw-data pipeline, correct dependency graph, every target a real working script) was **not** re-executed end-to-end — doing so would repeat the cumulative ~40+ hours of CNN/FM training this project's Gates 4–8 already spent. This is disclosed in the Makefile's own header comment and in `README.md`, not buried. All 13 of Gate 7's 14 FIGURE-REGENERABLE orphans now have committed producing scripts (the 14th, the retired ceiling-metric figure, is intentionally not regenerated — documented, not silently dropped); `audit_provenance.py` now reports **0 orphans**, re-verified at the end of this gate.

**2f. Archive.** `package/ARCHIVE_METADATA.md` — draft Zenodo deposit metadata and HF dataset-card front-matter, with an explicit decision flagged for the maintainer (whether the CC-BY-NC-4.0 PromoGen2 content needs a separate Zenodo deposit, since Zenodo licenses apply per-deposit not per-file). **Nothing published, nothing submitted** — artifacts only, as instructed.

---

## TASK 3 — Known issues

`out/KNOWN_ISSUES.md` (mirrored into `package/KNOWN_ISSUES.md`), 9 items: the n_hosts≤6 ceiling, the translation floor-value artifact and corrected N, *B. subtilis* measurement noise (with Task 1's finding folded in), fold-vs-draw variance, the batching cliff, Evo 2's specific unevaluated reason, the ceiling metric's reversal-of-a-reversal history, the cross-host calibration failure, and the held-out split's non-cryptographic enforcement. Written for users, specific rather than generic, each pointing to the gate memo with the full derivation.

---

## TEMPTATIONS TO ADJUST

1. **Whether to present the disattenuation gap-survival result using only the most favorable (empirically-anchored, gap-unchanged) reliability level and omit the sensitivity range that shows some narrowing.** Not done — both the empirically-anchored result and the full 0.5/0.7/0.9 sensitivity grid are reported, including the honest observation that the gap does narrow somewhat under aggressive assumptions, just nowhere near closing.
2. **Whether to claim the held-out evaluation split is more secure than it is**, since the whole point of the exercise is credibility. Not done — `KNOWN_ISSUES.md` item 9 states plainly that fold 4's labels remain visible in the public dev table, and that there is no technical barrier against looking them up.
3. **Whether to claim `make reproduce-full` was tested end-to-end** to make the reproducibility story cleaner. Not done — the Makefile header, `README.md`, and this memo all state precisely what was and was not re-executed, and why (compute-time cost, not doubt about correctness).
4. **Whether to omit the DNABERT-2 HF-card license-tag gap** (no `license` field on the model card itself; Apache-2.0 confirmed only via the separate GitHub repo) since it complicates a clean "verified" claim. Not done — stated exactly as found in `data/licensed/dnabert2_derived/LICENSE`, including which source actually carried the confirmation.

---

## WHAT I COULD NOT DO

1. **No direct reliability estimate exists for *B. subtilis* or *P. aeruginosa* specifically, or for translation in any host** — only EC-transcription has usable replicate-like structure in the released data. The sensitivity analysis bounds this; it does not resolve it with certainty.
2. **`make reproduce-full` was not re-executed end-to-end** — the dependency graph is correct and every target maps to a real, previously-run script, but a fresh, unattended, from-raw-data-to-final-numbers run was not performed in this packaging pass (would cost tens of hours).
3. **Nothing was actually tested from a real `git clone` of a pushed remote** — no remote repository exists yet. The clean-venv `pip install -e .` test is a real but narrower verification than a true clone-and-build.
4. **The Zenodo/HF metadata YAML was not validated against either platform's live, current schema** — both platforms update periodically; the maintainer should check current docs before submitting.
5. **Held-out fold 4's non-cryptographic enforcement was not engineered around** (e.g., by actually stripping fold 4 from the public dev table) — doing so would require restructuring the fold system this project's whole history of published results depends on. Disclosed as a structural limitation instead of silently worked around.

---

## CONTRADICTIONS WITH THE CHARTER

1. **The charter's Gate 8 packaging spec assumes a benchmark being packaged around a positive or at-least-fully-resolved result.** This benchmark ships around a negative headline result (H-MAIN) with a mechanistic explanation that itself required two rounds of self-correction (Gates 6→7→8) to arrive at reliably. The packaging reflects this honestly (`KNOWN_ISSUES.md`, the retired-metric history) rather than smoothing it into a cleaner-sounding release note.
2. **DreamCROSS/DREAM-challenge-style withheld evaluation assumes a submission barrier the charter's "no hosted server yet" constraint cannot actually provide** — Task 2b's design is the honest version of this given real constraints, not the idealized version the charter's phrasing (borrowing ProteinGym/DREAM as models) might imply.

---

## FILES WRITTEN

- `out/GATE8_MEMO.md` (this file), `out/KNOWN_ISSUES.md`, `out/PAPER_FRAMING.md` (updated)
- `scripts/72_attenuation_analysis.py`, `scripts/73_attenuation_figures.py`
- `scripts/74_build_held_out_eval.py`, `scripts/75_master_baseline_export.py`
- `scripts/76_regen_calibration_figures.py`, `scripts/77_regen_remaining_figures.py`
- `out/results/gate8_attenuation_analysis.{json,csv}`, `out/gate8_held_out_eval_private_labels.parquet` (NOT for distribution)
- `out/figures/gate8_disattenuation.png`, `gate8_ceiling_anomaly_resolution.png`
- `package/` — the full distributable benchmark (pyproject.toml, LICENSE, README.md, KNOWN_ISSUES.md, ARCHIVE_METADATA.md, Makefile, environment.lock.txt, `crosshost/` package code, `data/core/` + `data/licensed/*`, `baselines/`, `held_out_eval/`, `figures/`, `tests/`)
- `Makefile` (project root)
- `out/state.json` (`gate_8` block)
