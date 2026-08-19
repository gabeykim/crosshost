# Reproducibility Statement

## What `make reproduce` regenerates (fast path — actually run this gate, exit code 0, zero errors)

`make reproduce` regenerates every figure and table in this manuscript from already-computed intermediate results (`out/results/*.json`, `data/*.parquet`) — it does not retrain any model or re-run inference. **Run in full during this gate** (2026-08-19) after extending it to cover Gates 8.5–10.5, which a prior packaging pass had left out of the fast path: it completed with exit code 0, and its own printed output reproduced every GC-confound number in Table S3/S4 of this package bit-for-bit against the numbers already written into the manuscript — an independent regeneration check, not just a read of cached files.

Target now covers, in order: both audit scripts; Gate 6 figures; the crosshost-correlation figure; the feature-ablation figure; the attenuation figures; calibration-figure regeneration; the remaining Gate-7-orphan regeneration; provenance triage; the master baseline export; then (added this gate) the two-stage-transfer table, the conditioning-mechanism comparison, the shift-prediction summary, the co-active range-restriction check, the FiLM fold-variance figure, the Gate 8.6 figures, DRAFTS parsing/joining/correlation/modality/GC-and-floor-check/feasibility, the Gate 10 figures, the in-vivo GC-confound analysis, and the Gate 10.5 figures.

**Deliberately NOT in the fast path** (need a trained-model checkpoint for live inference, or an actual training run — not just cached results): `scripts/78_conditioning_mechanisms.py` (trains concatenation and per-host-head models), `scripts/80_coactive_subset.py` Part 2 (re-scores existing checkpoints via live inference — no training, but does need `.pt` files and a forward pass, categorically different from the pure-table/JSON scripts above), `scripts/81_shift_prediction.py` (trains the shift-prediction architecture), `scripts/84_shift_regression_to_mean_control.py` (imports `scripts/81` as a module; conservatively excluded from the fast path since its import-time cost was not independently verified this gate). Their outputs (CSVs/JSONs already in `out/results/`) are shipped and used by every downstream figure/table script that needs them — you do not need to re-run these four to reproduce anything in this manuscript, only to regenerate their own underlying numbers from scratch.

## What `make reproduce-full` does NOT do (stated plainly, not implied)

`make reproduce-full` re-runs the entire pipeline from raw data, including all CNN and foundation-model training (Gates 4–8). **This was not re-executed end-to-end in this packaging pass, or in any prior one.** It would cost the cumulative 40+ hours of training this project's Gates 4–8 already spent once. The Makefile's dependency graph is accurate and every target maps to a real, working script — each individual script has been run and produced its shipped output at least once (that is how `out/` and `data/` were populated) — but running the full chain consecutively, in one sitting, from raw data alone, has not been attempted through Gate 11. If you need to verify the full pipeline yourself, budget accordingly; do not assume `make reproduce-full` has been smoke-tested end-to-end just because `make reproduce` has.

## What was tested from a genuine fresh clone (Gate 8.5, not re-tested this gate)

`git clone` of the pushed repository, `pip install -e ./package`, `pytest` (9/9 passed), `make audit`, `make reproduce` — all in a fresh clone, not just a fresh venv on the existing working directory (`out/state.json` gate_8_5.task4_packaging). **This gate re-ran `make reproduce` in place (not from a fresh clone) after extending its target list** — the fresh-clone claim above still describes Gate 8.5's narrower target list, not this gate's extended one. A fresh-clone re-verification of the extended target list has not been done and is a reasonable next check before actually posting, not assumed equivalent to this gate's in-place run.

## Audits

`scripts/audit_leakage.py` (6 checks: fold-splitting integrity, RS241 exclusion, max train-test identity, exact-duplicate co-assignment, fold-local normalization, manifest hash match) and `scripts/audit_provenance.py` (orphan-file detection, heuristic — see `FIGURE_AUDIT.md` for its disclosed limitation and what this gate found by hand-checking past it) both run as part of `make audit` and `make reproduce`, and both passed cleanly at the start and end of this gate (see the Part VI status report for exact counts).

## Environment

`environment.lock.txt` — pinned dependency versions matching the environment this project's results were produced in. Not tested on a different OS/Python version than the one used throughout (macOS, Apple Silicon MPS — see Limitations item 6 for why this excluded Evo 2).
