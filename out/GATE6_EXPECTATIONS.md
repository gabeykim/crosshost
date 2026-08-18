# CROSSHOST — Gate 6 Expectations (written before any foundation-model evaluation)

**Committed:** 2026-08-07, before `scripts/63_fm_loho_calibration_rs241.py` was run for any model. This is not a formal pre-registration under `out/PREREGISTRATION.md` — Gate 6 tests capacity/architecture, not a new scientific hypothesis with a pass/kill threshold — but the same discipline that caught real problems twice on this project (the *B. subtilis* growth-rate gap in Gate 4.5, the *P. aeruginosa* sign inversion investigated the same gate) applies here: state the prediction before seeing the number, so a wrong prediction is informative rather than quietly forgotten.

Feasibility assessment (below) was completed before this document was written, per Task 1's explicit instruction to assess feasibility before committing. No LOHO/calibration/RS241 evaluation of any foundation model has been run yet.

---

## 1. The ~50%-of-ceiling prediction

**Background [VERIFIED, `out/state.json` gate_5_5.task2_signal_location.per_host_ceiling]:** the sequence-only CNN (213,956 params, no host-conditioning pathway) reaches a strikingly consistent fraction of its cross-host measurement-correlation ceiling for transcription across all three primary hosts — EC 49.2%, BS 50.7%, PA 44.8% — despite the ceilings themselves varying enormously (EC–PA measurement ρ=0.75 vs. BS-pair ρ=0.16–0.26, `gate5_5_crosshost_measurement_correlation.json`). Translation does not share this consistency (EC 20.9%, BS 57.4%, PA 17.6%) and is treated separately below.

**The falsifiable prediction:**

- **If ~50% is a property of the *problem*** — i.e., a ceiling on how much of the cross-host-shared signal any reasonable sequence model can extract from 165bp regulatory DNA at this dataset's noise level — **foundation models should also land near it**, regardless of parameter count, because more capacity cannot extract signal that the data does not carry.
- **If ~50% is a property of *this project's specific 214k-parameter architecture's capacity*** — undertrained conv filters, insufficient receptive field, insufficient depth — **foundation models pretrained on billions of bases of genomic sequence should exceed it**, because they bring pretrained representational capacity this project's from-scratch CNN does not have.

This is a clean two-outcome test with no ambiguous middle reading: I commit now to reporting whichever the data shows, plainly, in Task 3.

**Translation is explicitly excluded from this specific prediction — reasoning stated before the fact:** translation's per-host ceiling pattern (57.4% BS vs. 20.9%/17.6% EC/PA) was already flagged in Gate 5.5 as "unexplained," not consistent like transcription's. I am not making a directional prediction for translation; I will report where FMs land on it descriptively, without pretending a null hypothesis existed.

---

## 2. What would change the paper's central claim, and what would not — a specific threshold

**The central claim (`out/PAPER_FRAMING.md`):** host conditioning — genomic, physiology-proxy, or free identity — does not measurably improve cross-host regulatory-activity prediction over a host-agnostic sequence model, on this benchmark, with this architecture.

**Gate 6 does not directly re-test this claim.** No foundation-model variant tested here receives host information (Task 1/2 specify frozen embeddings + a lightweight head, structurally identical in role to the sequence-only ablation — see `scripts/62_fm_head_model.py` docstring). Gate 6 tests a different axis: whether *capacity/pretraining*, not conditioning, was the bottleneck.

**Threshold for "the claim narrows":** if a foundation model **substantially and distinguishably** (90% CI non-overlapping, consistent with the rest of this project's decision rule) beats the sequence-only CNN's zero-shot LOHO Spearman ρ on **at least half of the six primary host×readout cells**, that is evidence capacity — not merely conditioning-absence — was doing real work in the small model's ceiling. It would not falsify "conditioning doesn't help" (no conditioned-FM variant is tested to compare against), but it would open, without answering, the question of whether a *large conditioned* model might show a conditioning effect the small one could not detect. The paper would need to say so explicitly rather than let the small-model result stand as if capacity were ruled out.

**Threshold for "the claim survives intact":** if foundation models match or underperform the sequence-only CNN — consistent with DART-Eval (NeurIPS 2024) and "Specialized Foundation Models Struggle to Beat Supervised Baselines" (ICLR 2025), both of which found genomic/regulatory FMs losing to lightweight supervised baselines — there is no route by which "a bigger model would show conditioning matters" remains open, because the bigger model did not even improve on the *unconditioned* baseline. This is the stronger and more useful outcome for the paper's defensibility, and Gate 5.5's `PAPER_FRAMING.md` already anticipated it as the likelier case.

---

## 3. Honest prior, with reasoning

**Prior: foundation-model embeddings will roughly MATCH or UNDERPERFORM the sequence-only CNN on both readouts, landing at or below the ~50%-of-ceiling mark for transcription. Confidence: medium-high, not high — genuinely could be wrong, see below.**

Reasoning:
1. **Field-level prior.** DART-Eval and the ICLR 2025 "Specialized Foundation Models Struggle to Beat Supervised Baselines" paper both found genomic FMs losing to small supervised models on regulatory-activity-adjacent tasks. This project's own charter cited this as the working expectation before Gate 6 was elevated in priority.
2. **Evo 1's own reference numbers are not obviously stronger than this project's sequence-only result.** Evo 1's promoter-activity Spearman was reported at ~0.43 mean, 0.61 on Kosuri promoter+RBS — both single-host, likelihood-based, not cross-host zero-shot transfer. This project's sequence-only CNN reaches zero-shot cross-host ρ=0.370 (tx) / 0.415 (tl) mean across the three primary hosts (`gate5_5_sigma70_revisit.json`) — a harder task (genuinely unseen host, not within-host) landing in a broadly similar range, not obviously behind.
3. **Task-length and objective mismatch.** DNABERT-2/NT-500M/Evo2 are pretrained on masked-token or next-token objectives over long genomic context; 165bp discrete regulatory elements with continuous FACS-seq-derived activity labels is a short, narrow, and differently-shaped target than what these models were optimized for. A frozen embedding is not guaranteed to linearly encode activity-relevant structure just because the model has seen more sequence.
4. **Counter-consideration, stated honestly:** DNABERT-2's pretraining corpus explicitly includes bacterial genomes (135 species across 6 clades, `zhihan1996/DNABERT-2-117M` model card) — closer to this task's domain than a purely eukaryotic corpus would be — so it is not a foregone conclusion. PromoGen2, if used as a supplementary check, is pretrained specifically on prokaryotic promoter/regulatory sequence — the best domain match of anything available, and a real chance to beat the prior. I am not assuming the prior is right; I am stating it and will report the actual result without adjusting the prediction after the fact.

---

## 4. Feasibility assessment (completed before this document; Task 1/Task 2 methodology)

**Evo 2 — cannot run locally at any released size.** [VERIFIED, primary sources: `github.com/ArcInstitute/evo2` README + issue #67] The official `evo2` package hard-requires CUDA + Flash Attention on Linux; 1B/20B/40B additionally require Transformer Engine + FP8 on a Hopper GPU. It throws `ValueError: Expected a cuda device, but got: cpu` at initialization on non-CUDA hardware — not a slow fallback, a hard failure. No HuggingFace `transformers`-compatible port exists for Evo 2 (unlike Evo 1). This machine has no CUDA GPU (8GB unified-memory Apple Silicon, MPS only). **Largest tractable local option: none.**

**Evo 2 — hosted-API path, in progress.** NVIDIA Build (`build.nvidia.com`) hosts a free-tier API exposing Evo2-40B's forward-pass layer embeddings via a `/biology/arc/evo2/forward` NIM endpoint — no local compute needed. This requires an account (phone-verified signup, outside this agent's tool access — no browser automation available) and has an unconfirmed credit ceiling (secondary sources report ~1,000 inference credits on the free tier, unverified against NVIDIA's own documentation). Gabriel has approved using this route (2026-08-07) and is obtaining the API key; embedding extraction for Evo 2 is designed but not yet run, and is scoped to whatever the real credit/rate limit turns out to be once the key is available and pilot-tested. **This section will be updated with the verified limit before any Evo 2 embeddings are extracted at scale.**

**NTv3 — poor fit, not pursued as the primary Task 2 model.** [INFERRED, medium-high confidence, from InstaDeep's own NTv3 paper/model card] NTv3's *base* pretraining corpus (OpenGenome2) does include bacterial/archaeal sequence, so "NTv3 excludes bacteria" would be an overstatement — but its post-training stage fine-tunes on ~16,000 functional tracks from only 24 animal/plant species, and its architecture (U-Net, 1Mb context, single-nucleotide-resolution eukaryotic regulatory landscapes) is a poor structural match for 165bp bacterial elements regardless. It also requires a custom, non-`transformers`-standard loader with unconfirmed CPU/MPS compatibility. **Not selected.**

**DNABERT-2 — selected as the Task 2 model.** [VERIFIED by direct local load, `scripts/61_fm_embeddings.py`] `zhihan1996/DNABERT-2-117M`, 117,068,544 params (547× the sequence-only CNN's 213,956), pretrained on 135 species across 6 clades including bacteria, 32.5B bases total. Named explicitly as a candidate in the charter's own Task 2 text. Loads via standard `transformers.AutoModel(..., trust_remote_code=True)` after three documented compatibility patches for this 2023-era model file against a 2026-era `transformers`/`torch` stack (bypassing a stale `triton` import pre-check that the model's own code already guards with try/except; fetching one file the pre-check had prevented `transformers` from downloading; forcing two ALiBi-tensor constructors onto `cpu` explicitly because `transformers`' newer meta-device fast-init path left them on `meta` by default) — none of which touch the model's actual weights or computation, all documented in `scripts/61_fm_embeddings.py`. Benchmark: ~7.3ms/sequence at batch=128 on MPS; full-library embedding extraction (29,042 + 207 sequences) is running now (started before this document, ~3–4 minutes total, well within "feasibility assessment" rather than "evaluation").

**PromoGen2 — considered as a secondary/supplementary check, not the headline Task 2 model.** [NOT FOUND / not attempted] 149M-param GPT-2-architecture model pretrained specifically on prokaryotic promoter/cis-regulatory sequence — the best domain match found, but newer and less independently validated than DNABERT-2. Will be attempted after DNABERT-2 and Evo 2 are complete if time and compute budget allow; not required for Task 3's comparison table to be complete, since the charter asks for one "NTv3-or-alternative" slot and DNABERT-2 fills it.

---

## ADDENDUM — dated 2026-08-07, appended after DNABERT-2/PromoGen2 evaluation, correcting Section 1's premise

**Section 1's "remarkably consistent 49.2%/50.7%/44.8%" premise is not reproducible and appears to rest on a computational error in a prior gate.** While computing percent-of-ceiling for the foundation models (Task 3, `scripts/64_gate6_full_comparison.py`), I diffed `out/results/gate5_5_per_host_ceiling.json`'s `actual_zeroshot_mean` field against the project's own officially bootstrapped sequence-only zero-shot Spearman ρ (`out/results/gate5_5_fourway_comparison.json`, the file `GATE5_5_MEMO.md`'s own headline table is built from). They match in only 1 of 6 (host, readout) cells (EC transcription, within rounding); the other 5 diverge substantially — e.g. BS transcription: 0.1305 in the ceiling file vs. the official 0.2630. No script for `gate5_5_per_host_ceiling.json` exists anywhere in `scripts/`, meaning it was produced by an unsaved ad-hoc computation in a prior session — exactly the failure mode the charter's own standing rule 6 ("verify before trusting a stated number") exists to catch. The cross-host measurement correlations themselves (the ceiling values) were independently re-run via `scripts/57` and reproduce bit-for-bit; only the "actual" side of the ratio was wrong.

**Corrected figures** (`out/results/gate5_5_per_host_ceiling_CORRECTED.json`, `scripts/64`): sequence-only transcription percent-of-ceiling is **48.7% (EC) / 102.1% (BS) / 63.6% (PA)** — not a consistent ~50%. Translation is **75.8% (EC) / 110.6% (BS) / 52.5% (PA)**. BS exceeds its own nominal ceiling on both readouts; interpreted below, not treated as paradoxical (BS's ceiling is estimated from a much smaller double-active subset — n=2099 transcription / n=314 translation, versus EC–PA's n=9741/3826 — so it is a noisier, and on this evidence probably a downward-biased, estimate for BS specifically; also a different sequence subset than the LOHO zero-shot evaluation set, so "exceeding the ceiling" is a real possibility, not a logical contradiction).

**Consequence for the Section 1 prediction:** the clean two-outcome test as originally framed ("if ~50% is a property of the problem, FMs should also land near it; if of capacity, they should exceed it") assumed a premise — tight, host-independent clustering around 50% — that does not actually hold even for the sequence-only model itself once correctly computed. The test does not fail outright, but it is weaker than stated: percent-of-ceiling now varies by host (48–110%+) enough that "near ~50%" is not a sharp target. Section 3 of `out/GATE6_MEMO.md` reports where DNABERT-2/PromoGen2 actually land on the corrected metric, host-by-host, without pretending the original clean framing survived.

**This is a correction, not a temptation acted upon** — it was found by verifying an input against a reproducible source before using it, not by adjusting a result because it was inconvenient. It is reported here per the same discipline TEMPTATIONS TO ADJUST exists to enforce: visible, dated, not smoothed over.

---

## 5. What I will NOT do in response to whatever Gate 6 shows

Per the charter's own governing principle for Gate 5 (retained as a standing requirement): if a result is inconvenient — an FM beating the sequence-only model, or landing suspiciously exactly at the predicted ~50%, or Evo 2 turning out infeasible even via the hosted API — I will report it plainly in `out/GATE6_MEMO.md`'s TEMPTATIONS TO ADJUST section rather than adjust the evaluation protocol, the comparison set, or the central-claim threshold stated in Section 2 above after the fact.
