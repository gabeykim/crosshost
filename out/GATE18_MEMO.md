**Do 93.6% and 26.0% hold, and is the regime description correct? BOTH FIGURES HOLD — 807/862 = 93.62% cell-free, 3,668/14,088 = 26.04% in vivo, recomputed from primary data and matching the committed outputs exactly. BOTH REGIME DESCRIPTIONS ARE CORRECT as stated, with two precisions worth adding. One separate finding you should know before rewriting: none of these claims are in the manuscript yet — `out/PREPRINT/MANUSCRIPT.md` was last changed in Gate 14 and contains no retention figures, no regime names, and no Gate 15 numbers at all.**

# GATE 18 MEMO — Verify two manuscript figures

Script: `scripts/100_gate18_coactive_retention.py`. Outputs: `out/results/gate18_coactive_retention.csv` (96 rows), `out/results/gate18_regime_verification.json`. No manuscript file touched, nothing trained.

---

## Premise check, before anything else

The gate describes claims "in the abstract and §3.2" and "in §3.3". **They are not there.** `out/PREPRINT/MANUSCRIPT.md` and `out/MANUSCRIPT.md` were last modified by commit `bcc660f` (Gate 14); `git log` shows no manuscript change since. Searching both files finds no `93.6`, no `26.0%` retention claim, no `0.509`, no `0.386`, no "active-fraction matched", no "uniform quantile". The only occurrence of "co-active" in the manuscript is the pre-existing phrase "n=9,741 co-active pairs" in §3.1.

Those claims currently live in **`out/GATE15_MEMO.md`** and the `gate15_*` result files. The 96 MB `out.zip` at the repo root is also Gate-14 era (its `MANUSCRIPT.md` is dated 08-19-2026) and contains no newer draft.

So either the Gate 16/17 rewrite happened somewhere not in this repository, or it has not happened yet. **Nothing here is wrong or lost from this repo's perspective** — but if you have a rewritten manuscript elsewhere, it is uncommitted and unpushed, and this verification was performed against the Gate 15 source numbers rather than against that draft's rendering of them. Worth resolving before the rewrite, since verifying numbers against a draft I cannot see is exactly the gap this gate was opened to close.

The verification itself is unaffected: the numbers checked are the ones the rewrite would be drawing on.

---

## TASK 1 — Co-active retention rates: both figures CONFIRMED

Every number recomputed from `data/three_host.parquet` and `data/drafts.parquet`, replicating the filter definitions in `scripts/80` (in vivo) and `scripts/98` (cell-free), then cross-checked against those scripts' committed outputs.

**Cross-check: 51 of 51 agree.** All 6 in-vivo pair×readout n-pairs match `gate8_5_coactive_correlation.csv` exactly; all 45 cell-free pairs match `gate15_coactive_modality.csv` exactly. The hand-derived figures were derived correctly.

| claim | recomputed | verdict |
|---|---|---|
| cell-free EC–BS co-active retention **93.6%**, 807 of 862 | 807 / 862 = **93.62%** | **CONFIRMED** |
| in-vivo EC–BS co-active retention **26.0%**, 3,668 of 14,088 | 3,668 / 14,088 = **26.04%** | **CONFIRMED** |

### All pairs and readouts

**In vivo** (full library, n=29,249):

| pair | readout | co-active / pooled | retention |
|---|---|---|---|
| EC–BS | transcription | 3,668 / 14,088 | **26.04%** |
| EC–PA | transcription | 9,741 / 19,643 | 49.59% |
| BS–PA | transcription | 2,099 / 11,969 | 17.54% |
| EC–BS | translation | 866 / 866 | **100.00%** ← see caveat |
| EC–PA | translation | 3,826 / 4,513 | 84.78% |
| BS–PA | translation | 314 / 352 | 89.20% |

**Cell-free** (45 pairs, full 1,047-sequence library): min **53.9%**, median **95.9%**, max **100.0%**. The three lowest all involve *S. enterica* (Ec–Se 53.9%, Se–Bs 54.1%, Se–Pa 54.1%) — the one species with a large silent class (45.9% inactive, Gate 15).

### Caveat that matters if a retention column goes into the paper

**Translation retention is not comparable to transcription retention, and EC–BS translation's 100% is an artifact of the definition, not a biological fact.** For transcription, "pooled" means DNA-adequate in both hosts with `rna==0` rows kept as `tx_norm = 0`. For translation, `scripts/80` builds "pooled" as *already floor-corrected* — floor-pinned rows are excluded from the denominator before anything else happens (detected floors: EC 1.87 pinning 67.1%, BS 2.00 pinning 90.2%, PA 2.23 pinning 9.9%). The co-active filter then removes only rows *below* the floor, of which EC–BS has none. Hence 866/866.

A reader seeing "EC–BS translation: 100% co-active retention" next to "EC–BS transcription: 26%" would reasonably conclude translation is far better conserved. It is not — the two denominators mean different things. **If a retention table is published, transcription and translation must not share a column without this stated.**

---

## TASK 2 — Regime descriptions: BOTH CORRECT

Verified by recomputing each regime independently of `scripts/98` and reproducing its exact n.

**"Active-fraction matched (per-host)" — ρ = 0.509, n = 218. Description CORRECT.**
Each cell-free host is restricted to the top fraction matching that host's own in-vivo active fraction, then the two are intersected. Confirmed: EC target 15,040/24,613 = 61.11% (threshold activity ≥ 0.1342), BS target 4,435/15,848 = 27.98% (threshold ≥ 0.1922), intersection n = **218**.

**"Top 26% of both (uniform quantile)" — ρ = 0.386, n = 150. Description CORRECT.**
The pooled in-vivo EC–BS co-active retention rate (3,668/14,088 = 26.04%) is applied to both hosts as a hard quantile, then intersected. Confirmed: EC threshold ≥ 0.7955, BS threshold ≥ 0.2263, intersection n = **150**.

### Two precisions to add

1. **Both regimes take their quantiles over the DNA-adequate (pooled-equivalent) population, not the co-detected population.** The description doesn't say which population the "top fraction" is a fraction *of*, and it is not the obvious one. Both use `>=` (inclusive).

2. **The "26%" regime does not retain 26% — it retains 17.4%** (150 of 862). Applying a 26% cut to each host independently and intersecting lands well below 26%. By contrast the per-host regime retains **25.3%** (218 of 862), which is almost exactly the in-vivo 26.04% it was built to mirror. This is a substantive point in the per-host regime's favour that the current description leaves out: it reproduces the in-vivo retention rate as an *emergent consequence* of matching per-host fractions, whereas the uniform-quantile regime overshoots the restriction by a third. If only one matched regime is reported, the per-host one is the better-justified choice, and the uniform one is best framed as the deliberately harsh bound it is.

### One inconsistency across the six regimes, reported because it is not visible from the table

The two matched regimes threshold on the **DNA-adequate** distribution with `>=`. The `both_above_p25` and `both_above_med` regimes threshold on the **usable** (co-detected-eligible) distribution with strict `>`. Each choice is defensible on its own, but the six rows are therefore not all on identical footing, and a reader comparing 0.647/0.592 against 0.509/0.386 is comparing thresholds struck over two different populations. Not an error; a labelling obligation.

---

## TASK 3 — Figure audit

**Main text: clean.** All eight figures are cited and all eight are present; nothing cited is missing, nothing present is uncited.

| | result |
|---|---|
| Figures cited in manuscript | 1, 2, 3, 4, 5, 6, 7, 8 (as "Figures 1-2", "Figure 3", "Figures 4-5", "Figure 6", "Figure 7", "Figure 8") |
| Files in `out/PREPRINT/figures/` | 9 files for 8 figures — Figure 2 is two panels (`Figure2a_gc_control.png`, `Figure2b_phylum_stratified.png`) |
| Cited but missing | **none** |
| Present but uncited | **none** |
| Producing scripts committed | **8 of 8** — scripts 55, 65, 68, 73, 86, 87, 91, 96, all present and git-tracked |
| Captions stating interval type | **8 of 8** |
| Captions stating N | **7 of 8** ← one gap |

**The gap: Figure 8 (shift-prediction retraction) states its interval type but not its N.** Its caption says "No interval shown — single point estimates per cell on the held-out fold" and cites "10 of 11 cells," which is a count of *cells*, not of sequences. The per-cell sequence N is available and ranges **294–9,609** depending on host pair and readout (e.g. EC→BS transcription n = 3,589; EC→PA transcription n = 9,609; BS→PA translation n = 294) — `out/results/gate8_5_shift_prediction_summary.csv`, column `n_total`. **Not fixed here**, per the instruction to touch no manuscript file; the legends file is yours to edit.

**Supplementary, for completeness:** 19 files implementing Figures S1–S10 exist and are documented in `FIGURE_LEGENDS.md`, but **none is cited from the manuscript body.** The manuscript's "Fig. S13", "Fig. S15" and "Appendix Fig. S5" references are to *Johns et al.'s* and *DRAFTS's* own supplementary figures, not to this package's. Worth a glance before submission: this package's own **Figure S5** (genomic feature-group ablation) and DRAFTS's **Appendix Fig. S5** (cited in §3.2) share a label, which is the same species of collision this project has caught twice before.

---

## Verdict

Both Task 1 figures hold. Both Task 2 descriptions are correct. Task 3 finds one incomplete caption and no missing or orphaned main figures. Nothing in this gate blocks the rewrite; three things are owed to it — the translation-retention caveat, the two regime precisions, and Figure 8's N.
