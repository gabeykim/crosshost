# CROSSHOST held-out evaluation

Following the ProteinGym / DREAM promoter challenge convention: fixed folds, out-of-fold predictions are the evaluation, and the test labels are not distributed publicly.

## What's withheld

Fold 4 of the frozen 5-fold split (5,813 rows, spanning all 3 primary hosts) is the official held-out test fold. `test_inputs.parquet` in this directory has the sequences and host metadata for these rows — **no activity labels.** Folds 0–3 remain fully public in `data/core/three_host.parquet` for model development.

**Honest limitation, stated plainly (see `KNOWN_ISSUES.md`):** this is not cryptographically enforced. `data/core/three_host.parquet` is the full public library and does still contain fold-4 rows with their real labels, since it is the same table used for model development elsewhere in this project. A sufficiently motivated user could look up fold 4's labels directly rather than submitting predictions honestly. This benchmark has no hosted server and no way to prevent that. The withheld-eval protocol exists to support honest, good-faith comparison — like most small academic benchmarks, it relies on the convention being respected, not on a technical barrier.

## How to submit

There is no hosted submission server yet. To get your model scored:

1. Predict on every row in `test_inputs.parquet`.
2. Write a parquet file with an `'OLIGO ID'` column matching `test_inputs.parquet`, plus, for each readout you want scored: `'{readout}_active_prob'` (predicted P(active), only needed for transcription's classifier metrics) and `'{readout}_strength_pred'` (continuous strength prediction). `readout` is `transcription` or `translation`.
3. Send the maintainer your submission file (see the repository's contact info) or, if you have access to the private labels file yourself (maintainers only), run:
   ```bash
   python held_out_eval/score_submission.py path/to/your_submission.parquet
   ```

## How scoring works

`score_submission.py` merges your predictions against the private label file (`out/gate8_held_out_eval_private_labels.parquet` in the source repository, never distributed), applies the same usable/active filtering as every other result in this project (`{host}_tx_usable`, `{host}_tl_usable`), and reports Spearman rho, MCC, and AUC per host and readout via `crosshost.evaluate()` — the identical metric suite and 90% bootstrap protocol used throughout `PAPER_FRAMING.md`. **Verified working**: tested end-to-end with a synthetic random-prediction submission during Gate 8 packaging (correctly scored near-zero rho, as expected for random noise — see `out/GATE8_MEMO.md`).

## Why fold 4, and not RS241

RS241 looks like the natural "true held-out" set, but its labels have already been used to compute and publish results throughout this project (Gates 5, 5.5, 6, 7 all report RS241 zero-shot numbers in `PAPER_FRAMING.md`). Withholding already-published labels from future submitters would not be a real evaluation barrier. Fold 4 was chosen instead — an arbitrary but genuinely-unpublished-at-the-per-row-prediction-level partition of the primary-host library.
