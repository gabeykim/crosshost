# Baseline scores

`master_baselines.csv` / `.json` — every model, fold-resolved, 90% percentile bootstrap intervals, identical protocol throughout (`crosshost.evaluate.bootstrap_ci_90`, `N_BOOT=10000`, `SEED=0`).

| system | what it is | source gate |
|---|---|---|
| `B1_mean_majority` | constant predictor (the floor) | Gate 3 |
| `B2_per_host_N` | per-host model trained on N examples from the SAME host — "just measure N parts in your host," the project's real competitor | Gate 3/3.5 |
| `B4_biophysical` | sigma-70 motif score (transcription) / ΔG (translation), no learning | Gate 3 |
| `free_embedding_B3` | per-host lookup embedding, fallback = mean of training-host embeddings on unseen hosts | Gate 3 |
| `genomic` | FiLM-conditioned CNN, 37-D genomic host feature vector | Gate 4 |
| `physiology` | FiLM-conditioned CNN, 6-D physiology proxy vector | Gate 4 |
| `sequence_only` | same CNN trunk, NO host-conditioning pathway at all | Gate 5.5 |
| `dnabert2` | frozen DNABERT-2 (117M) embeddings + lightweight head | Gate 6 |
| `promogen2` | frozen PromoGen2 (148M) embeddings + lightweight head — see licensing note below | Gate 6 |

Columns: `host`, `readout`, `eval_point` (`single` / `N0_zeroshot` / `N100` / `N3000`), `mechanism` (`n/a` / `per_host_trained` / `head_only` / `top_conv`), `system`, `rho_mean`, `rho_lower90`, `rho_upper90`.

**Licensing note:** `promogen2` scores were computed from CC-BY-NC-4.0-licensed model embeddings (`data/licensed/promogen2_derived/`) — non-commercial use only for anything derived from them. `dnabert2` scores are Apache-2.0 (permissive). See the root `README.md` licensing table.

## Gate 7 calibration — shipped prominently, not in an appendix

`gate7_conformal_and_ece.json` — split-conformal prediction intervals (80%/90% target coverage) and Expected Calibration Error for the active/inactive classifier, `sequence_only` vs. `B2_per_host_N` (N=3,000).

**Headline finding: `sequence_only`'s zero-shot classifier is well-calibrated on *E. coli* (ECE 0.05–0.10) but severely miscalibrated on *B. subtilis* and *P. aeruginosa* (ECE 0.34–0.47) — a 4–9× degradation.** Regression-interval coverage stays close to nominal regardless (conformal prediction guarantees this independent of the underlying model's calibration) — the failure is specifically in the classifier's absolute probability estimates, not in the interval widths. See `figures/gate7_reliability_diagrams.png`.

**If you deploy `sequence_only` (or any cross-host model in this suite) on a host meaningfully different from its training hosts, do not trust its raw predicted probabilities without recalibrating on whatever labeled data you can get for that host.** A model that wins on rank correlation can still be confidently wrong in absolute terms.
