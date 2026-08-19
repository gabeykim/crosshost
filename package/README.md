# CROSSHOST

A benchmark for predicting bacterial regulatory DNA activity across host species. Measures whether a model — from a lightweight per-host baseline to a 148M-parameter genomic foundation model — can predict how a 165bp regulatory sequence behaves in a bacterial host it has never seen.

**Central finding:** host conditioning (genomic features, physiology features, or two genomic foundation models orders of magnitude larger than this project's own architecture) does not measurably improve cross-host prediction over a host-agnostic sequence model. The likely reason is visible directly in the raw data: phylogenetically close hosts share far more cross-host regulatory signal than distant ones. Full writeup: [`PAPER_FRAMING.md`](./PAPER_FRAMING.md). Known limitations: [`KNOWN_ISSUES.md`](./KNOWN_ISSUES.md) — **read this before drawing conclusions from any number in this package.**

**Scope ceiling, stated up front:** three primary hosts (*E. coli*, *B. subtilis*, *P. aeruginosa*) with dense coverage, plus three more (*S. enterica*, *V. natriegens*, *C. glutamicum*) at ~207 recoverable sequences via RS241. **n_hosts ≤ 6.** This is not a large-scale benchmark in the ImageNet sense; it is a small, carefully-audited one, and any claim about "genome-encoded functions" as a general class should be read with that ceiling in mind.

## Install

```bash
pip install -e .                 # core package, data loaders, evaluation API
pip install -e ".[models]"       # + torch/transformers, needed only to re-run the FM baselines
pip install -e ".[dev]"          # + pytest, for running tests/
```

Tested in an isolated venv (`python3 -m venv`, `pip install -e .`, fresh process, different working directory) — see `out/GATE8_MEMO.md` for the exact verification steps. **Tested from a fresh `git clone` of the pushed repository as of Gate 8.5** (`pip install -e ./package`, `pytest` 9/9 passed, `make audit`, `make reproduce` — see `out/state.json` gate_8_5.task4_packaging) — see Reproducibility below for exactly what was and wasn't verified.

## Quickstart

```python
import crosshost

library = crosshost.load_three_host_library()      # 29,249 rows, 165bp sequences + per-host activity
rs241 = crosshost.load_rs241()                       # held-out-host eval set, never train on this
genomic = crosshost.load_genomic_features()          # 37-D host feature vector, 6 hosts
physiology = crosshost.load_physiology_features()     # 6-D depth-matched physiology proxy
splits = crosshost.load_splits()                       # frozen fold assignment, 0..4

# score your own predictions with the project's exact metric suite + 90% bootstrap CI
result = crosshost.evaluate(predictions_per_fold, targets_per_fold)
print(result["transcription"]["spearman_rho"])   # {'mean': ..., 'lower': ..., 'upper': ..., 'n_folds': 5}
```

## Licensing quarantine

| Directory | License | Contents |
|---|---|---|
| `crosshost/` (code) | MIT | This package's own code |
| `data/core/` | MIT (data), with attribution required | Derived activity values, splits, genomic/physiology features — see `data/core/CITATION.md` |
| `data/licensed/dnabert2_derived/` | Apache-2.0 | DNABERT-2 frozen embeddings (permissive, verified against the model's source repo) |
| `data/licensed/promogen2_derived/` | **CC-BY-NC-4.0, non-commercial** | PromoGen2 frozen embeddings — do not use commercially without clearing terms with the PromoGen2 authors |
| `data/licensed/deepcross/` | N/A (empty) | Pre-provisioned quarantine per the project charter; no DeepCROSS data was ever incorporated |

**Never merge `data/licensed/*` content into `data/core/`.** The directory boundary is the enforcement mechanism, not a README note.

Raw Johns et al. 2018 supplementary tables (Springer Nature copyright) are **not redistributed** — only derived values, with attribution. See `data/core/CITATION.md`.

## Held-out evaluation

See [`held_out_eval/README.md`](./held_out_eval/README.md) for the withheld-label evaluation protocol (fixed folds, out-of-fold scoring, no hosted server yet — manual submission process documented there).

## Baseline scores

Every model from Gates 3–8, fold-resolved with 90% bootstrap intervals: [`baselines/`](./baselines/). Includes the sequence-only ablation, both foundation models, and — prominently, not in an appendix — **Gate 7's calibration failure**: the sequence-only model's zero-shot classifier is well-calibrated on *E. coli* (ECE 0.05–0.10) but severely miscalibrated on *B. subtilis* and *P. aeruginosa* (ECE 0.34–0.47), a 4–9× degradation off its training distribution. A cross-host model can be confidently wrong. See `baselines/README.md`.

## Reproducibility

- `Makefile` — `make reproduce` regenerates all figures/tables from the shipped intermediate results (fast, actually tested — see below). `make reproduce-full` re-runs the full training pipeline from raw data (slow — many hours of CNN/FM training across the original project's Gates 4–8 — **not re-executed in this packaging pass**; the Makefile's dependency graph is correct and each target maps to a real, working script, but end-to-end multi-hour execution was not repeated here. Pretrained checkpoints and result files are shipped directly so `make reproduce` doesn't require it.)
- `environment.lock.txt` — pinned dependency versions, matching the environment this project's own results were produced in.
- `scripts/audit_leakage.py` and `scripts/audit_provenance.py` — both included, both part of `make reproduce`.
- **What was actually tested, stated precisely:** `pip install -e .` in a fresh venv (works); data loaders and `evaluate()` against real and synthetic data (work, both positive and negative test cases); `scripts/audit_leakage.py` and `scripts/audit_provenance.py` (both pass, re-run at the end of every gate through Gate 11). **Tested as of Gate 8.5:** a literal `git clone` of the pushed repository, followed by `pip install -e ./package`, `pytest` (9/9 passed), `make audit`, and `make reproduce` — all in a fresh clone, not just a fresh venv on the existing working directory. **Still not tested:** the full multi-hour `make reproduce-full` training path — it would cost the cumulative 40+ hours of CNN/foundation-model training this project's Gates 4–8 already spent, and re-running it end-to-end in one sitting from raw data alone has not been attempted in any packaging pass through Gate 11. Every individual script in its dependency graph has been run and produced its output at least once (that is how `out/` and `data/` were populated); the full chain has not been re-run consecutively. Say this plainly rather than imply full-pipeline reproducibility that was never checked.

## Archive

Prepared for Zenodo (DOI metadata) and Hugging Face Datasets (dataset card) — see `ARCHIVE_METADATA.md`. **Nothing has been published.** Both are artifacts for the maintainer to review and submit.
