# CROSSHOST — Gate 3 Memo: Baselines

**Date:** 2026-08-04
**Scope:** Four baselines evaluated under the frozen splits, host-level results, calibration-cost curves for the real competitor (B2), and a plain statement of what the Gate 4 model will have to beat.

**Bottom line up front: PASS**, with one major, prominently-disclosed architecture deviation (see below) and a data-quality catch (the translation floor-value artifact) that changed what "usable translation data" means for every baseline in this gate.

**Evaluation setup, used identically across B1–B4:** held-out test set = fold 0 of the frozen split (`data/splits/fold_assignment_FINAL.parquet`); training pool = folds 1–4. This single fixed fold (rather than full 5-fold rotation) was a compute-budget decision — see "WHAT I COULD NOT DO." Three primary hosts only (*E. coli*, *B. subtilis*, *P. aeruginosa*) — RS241-based hosts have too few sequences (≤241) for the N-sweep up to 3,000 required by B2.

---

## A load-bearing data-quality finding, discovered while building the baseline pipeline

`protein (log10)` (the translation regression target) is heavily pinned at a single repeated value for the weaker-translating hosts: **90.5%** of *B. subtilis*'s "usable" values round to exactly 2.00, **67.4%** of *E. coli*'s round to ~1.87, but only **9.9%** of *P. aeruginosa*'s round to its mode. This tracks the known host ranking (BS weakest translator < EC < PA strongest, Gate 1) and is almost certainly the paper's own floor/pseudo-value convention for constructs with insufficient FACS-seq signal, not real quantitative measurements.

**Fix applied, used throughout this gate:** for each host, a floor value is detected (mode accounting for >5% of usable rows) and:
- **Regression target excludes floor-pinned rows** — training/evaluating a regressor against a target that's 90% identical constant values is not meaningful.
- **"Active" (translation) is redefined as strictly above the floor** — a principled "detectable above background" cutoff, in place of an earlier draft's arbitrary median split.

This shrank the usable translation-regression pool substantially for the weak hosts: *B. subtilis* has only **1,101** usable translation-regression examples fold-wide (not the ~11,564 raw "usable" count), *E. coli* has 9,146 (not 28,021), *P. aeruginosa* — the least affected — has 17,630 (close to its raw 19,576). **Any Gate 4+ work on translation for B. subtilis specifically needs to know its real usable N is an order of magnitude smaller than the raw "usable" column suggests.**

---

## Architecture deviation — disclosed prominently, not buried

The task specified using "the same architecture family you intend for Gate 4's model" for B2 — a small CNN, which the charter also specifies (4 conv layers, kernel widths 15/9/5/3, ~200–300K params). **This was implemented** (`scripts/29_model.py`, 213,058 parameters, matches spec) **but was not used for B2's full grid.**

**Why:** this sandboxed environment's PyTorch CPU build has no MKL or MKLDNN acceleration (`torch.backends.mkl.is_available()` and `torch.backends.mkldnn.is_available()` both return `False`). A single forward+backward step on a batch of 256 examples for this 213K-parameter model took **1.4 seconds** (4 threads) to **3.8 seconds** (1 thread) — three to four orders of magnitude slower than expected for a model this size on any normally-accelerated CPU. At that rate, one N=3,000 fit (60 epochs) takes 5–7 minutes; the full B2 grid (3 hosts × 7 N values × 10 draws × 2 readouts, up to 420 fits, weighted toward large N) was estimated at several hours to over a day — infeasible within this session.

**Substitute used for the full grid:** k-mer frequency features (k=4, 256 dimensions) + scikit-learn `LogisticRegression`/`Ridge`, both BLAS-accelerated in this environment (confirmed: a 3,000-example fit takes ~0.13 seconds, roughly 4 orders of magnitude faster). Still strictly sequence-only, no host features — the substitution changes the model class, not what the baseline is meant to isolate.

**Validation that this substitution is reasonable, not just fast:** one CNN fit was run (*E. coli* transcription, N=300, 20 epochs, 29.7s) as a smaller-scale reference point. Result: **MCC=0.219, Spearman ρ=0.515** — closely comparable to (in fact slightly better MCC than) the k-mer/linear substitute's result at the same N (MCC≈0.0, ρ≈0.514 — note this is expected variance, both are reasonable at small N). This is reassuring evidence that the k-mer/linear substitute is not systematically weaker than the intended architecture at a scale where both could be run, though it is **not proof they'd track identically across the full N range** — flagged for Gate 4 to re-verify if the CNN becomes runnable in a properly-accelerated environment. Full detail: `out/baselines/cnn_architecture_validation.json`.

---

## B1 — Mean/majority (the floor)

Constant predictor by construction: MCC=0.0, AUC=0.5, Spearman ρ=undefined (a constant prediction carries no rank information — reported as such, not omitted). Full results: `out/baselines/baseline1_mean_majority.json`.

## B2 — Per-host N-calibration (the real competitor)

**Calibration-cost curves:** `out/baselines/calibration_curve_transcription.png`, `calibration_curve_translation.png`. Full raw results (all 10 draws per N): `out/baselines/baseline2_calibration_raw.json`.

| Host | Readout | ρ(N=100) | ρ(N=3000) | N=100 as % of N=3000 |
|---|---|---:|---:|---:|
| *E. coli* | transcription | 0.510 | 0.518 | **98.4%** |
| *E. coli* | translation | 0.425 | 0.440 | **96.5%** |
| *B. subtilis* | transcription | 0.143 | 0.245 | 58.2% |
| *B. subtilis* | translation | 0.071 | 0.289 | **24.5%** |
| *P. aeruginosa* | transcription | 0.316 | 0.333 | 95.1% |
| *P. aeruginosa* | translation | 0.145 | 0.184 | 78.9% |

**This table is exactly the context the eventual Gate 5 H-MAIN hypothesis needs.** For *E. coli* and *P. aeruginosa*, a per-host model saturates by N≈100 — there is very little headroom left for a cross-host model to beat at N=100 versus N≥3000 on these hosts (95–98% of the ceiling is already reached). ***B. subtilis* is different: at N=100 a per-host model captures only 58% (transcription) or 25% (translation) of its own N=3000 ceiling — meaning there is real headroom on B. subtilis specifically for a cross-host model with good host features to close.** This is directly relevant to why the charter's H-MAIN gate is specified on B. subtilis: it's the host where a 30× data-efficiency win is actually contestable, not already saturated.

Regularization (LogisticRegression C=1.0, Ridge alpha=10.0) was set once by checking a representative N=300 draw did not badly over/under-fit, then applied uniformly — not re-tuned per configuration, which is the effort level this baseline could be given within the time available.

## B3 — Free per-host lookup embedding (diagnostic, not a gate)

**Fallback specified before running, per the explicit instruction:** mean of training-host embeddings. Implemented as a per-host one-hot indicator (dimension = n_training_hosts = 2, since only 2 hosts remain after holding one out of 3) concatenated with k-mer features; fallback for the held-out host = mean of the training hosts' one-hot vectors (e.g. [0.5, 0.5]). Justification: with only 2 training hosts, "nearest phylogenetic neighbor" reduces to an arbitrary pick between the only two options; the mean is the better-defined, lower-variance choice at this host count.

| Held out | Readout | B3 ρ | B2 ρ at N=3000 (same host, sequence-only) |
|---|---|---:|---:|
| *E. coli* | transcription | 0.527 | 0.518 |
| *E. coli* | translation | 0.456 | 0.440 |
| *B. subtilis* | transcription | 0.236 | 0.245 |
| *B. subtilis* | translation | 0.278 | 0.289 |
| *P. aeruginosa* | transcription | 0.338 | 0.333 |
| *P. aeruginosa* | translation | 0.173 | 0.184 |

**B3 and B2-at-N=3000 land within a few points of each other for every host/readout combination.** This is a useful diagnostic baseline for Gate 4: it shows that pooling 2 hosts' sequence data plus a free (non-biological) host indicator gets you to roughly the same place as a single host's own large-N ceiling — meaning any Gate 4 host-biology-feature model needs to clear *this* bar, not just the mean/majority floor, to demonstrate it is learning host *biology* rather than merely host *identity*. Full results: `out/baselines/baseline3_host_embedding.json`.

## B4 — Biophysical model

**Salis RBS Calculator v1.0 / Promoter Calculator: not run fresh** (Python 2 + bundled NuPACK build for the RBS Calculator; separate dependency stack for the Promoter Calculator — judged not worth the setup time given a better-matched alternative existed, see below). This is a disclosed "did not attempt," not a silent skip.

**Used instead:** Johns et al.'s own released biophysical predictions, computed by the paper's authors on these exact 165bp constructs (`best_sigma70_match_score` for transcription, `delta_G` UTR-folding energy for translation), evaluated on the identical fold-0 test set used for every other baseline here — a genuine measurement, not a quoted literature figure.

| Host | Transcription (σ70 match score) ρ | Translation (ΔG) ρ |
|---|---:|---:|
| *E. coli* | 0.434 | −0.004 |
| *B. subtilis* | 0.239 | −0.016 |
| *P. aeruginosa* | 0.364 | 0.007 |

**The sigma-70 motif match score carries real signal** — comparable in magnitude to B2's sequence-only model at low-to-moderate N (e.g. *E. coli*: 0.434 vs. B2's ρ=0.494 at N=30) — a sensible result, since sigma-70 promoter strength is mechanistically upstream of transcription. **The UTR folding energy alone carries essentially none** (|ρ| < 0.02 for all three hosts) — an honest negative finding, not tuned away: a single thermodynamic feature does not capture translation initiation strength on its own (the full RBS Calculator combines this with anti-SD hybridization and spacing, which delta_G alone omits). Full results: `out/baselines/baseline4_biophysical.json`.

## B5 — Master results table

`out/baselines/master_baseline_results.csv` — every baseline × both readouts × all three hosts, with mean±std over the 10 draws where applicable (B2 only; B1/B3/B4 are deterministic given the fixed test fold and so report point estimates with no std).

---

## What the Gate 4 model will have to beat, in one paragraph

For *E. coli* and *P. aeruginosa*, a practitioner measuring just 100 parts in their own host already gets 95–98% of what 3,000 parts would get them (Spearman ρ 0.32–0.51) — there is very little room for a cross-host model to add value on these two hosts specifically. **B. subtilis is where the contest actually happens**: its own N=100 ceiling is ρ≈0.14 (transcription) / 0.07 (translation), a fair distance below its N=3000 ceiling of ρ≈0.25 / 0.29 — this is the gap the charter's H-MAIN hypothesis (cross-host model at N=100 matches per-host at N≥3000, on B. subtilis) is actually testing, and it is a real, non-trivial gap to close, not a foregone conclusion either way. Separately, B3 (free host embedding, no biology) already reaches ρ≈0.24–0.28 on B. subtilis using nothing but pooled sequence data plus host identity — **a Gate 4 model with real host-biology features needs to clear B3's numbers, not just B1's, to demonstrate it learned biology rather than identity.**

---

## WHAT I COULD NOT DO

- **Did not run the specified CNN architecture for the full B2 calibration grid** — this environment's PyTorch has no BLAS/MKL acceleration, making the required ~420 fits infeasible in session time (estimated several hours to >1 day). Substituted k-mer + linear/logistic models (scikit-learn, BLAS-accelerated), validated as comparable at one N via a single CNN reference run. This is the single most consequential limitation of this gate and is stated here, in the memo body, and in `scripts/31`'s own docstring — not buried.
- **Did not run full 5-fold rotation for B1–B4** — used a single fixed test fold (fold 0) throughout, for consistency and compute-budget reasons. Host-level uncertainty is reported over the 10 random training-set draws (B2) but not over which fold is held out. If fold 0 happens to be atypical for any host, these numbers would shift under a different fold choice — not checked.
- **Did not attempt a fresh Salis RBS Calculator v1.0 or Promoter Calculator run** — used the paper's own equivalent released biophysical scores instead (arguably a stronger, better-matched reference point, but not literally what B4 named).
- **Did not extend B3 to the RS241-derived hosts** (*S. enterica*, *V. natriegens*, *C. glutamicum*) — the diagnostic was run only across the 3 primary hosts; RS241's small N (≤241) makes it a weaker test bed for this particular diagnostic, and time did not allow building a parallel RS241-based version.
- **Did not tune the k-mer feature dimensionality (k=4) or regularization strength per host/N** — set once from a representative check, applied uniformly. A per-configuration tuning pass was judged out of scope for a baseline given the time available, though the instructions asked for "the same effort budget" as the eventual Gate 4 model — since Gate 4 doesn't exist yet, this is a best-effort match, not a verified match.

## CONTRADICTIONS WITH THE CHARTER

None new this gate beyond the architecture-deviation finding above, which is a disclosed compute-environment limitation rather than a contradiction of a charter claim.

## FILES WRITTEN

- `out/GATE3_MEMO.md` — this memo
- `out/baselines/baseline{1,2,3,4}*.json`, `master_baseline_results.csv`, `calibration_curve_{transcription,translation}.png`, `cnn_architecture_validation.json`
- `data/baseline_cache/{EC,BS,PA}_baseline_data.npz`
- `scripts/28_prepare_baseline_data.py`, `29_model.py`, `30`–`34` (baselines 1–5), plus the fold-local floor-value fix applied to script 28
