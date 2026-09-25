**VERDICT: SURVIVES — but the manuscript's framing needs correcting, and one number it has never reported must be added. Cell-free EC–BS stays above the in-vivo co-active 0.258 under every restriction regime tested (0.386–0.677, minimum 1.5× the in-vivo value). The challenge's premise — that the published comparison puts a pooled cell-free number against a co-active in-vivo number — is incorrect, but for a non-obvious reason worth stating plainly: DRAFTS and Johns use the word "usable" to mean different things, and DRAFTS-usable is already the co-active operation. However, the challenge was right that something is missing: pooled-versus-pooled shows NO modality difference (0.616 cell-free vs 0.655 in vivo), and the manuscript nowhere reports this.**

# GATE 15 MEMO — Co-active cell-free correlation

Scripts: `scripts/98_gate15_coactive_modality.py`, `scripts/99_gate15_regime_audit.py`.
Outputs: `out/results/gate15_coactive_modality.{json,csv}`, `gate15_drafts_activity_distribution.csv`, `gate15_restriction_regime_audit.csv`.
No training, no new data, no manuscript edits.

---

## TASK 1 — Is "active" definable for DRAFTS transcription?

**Yes, it is definable — and the definition was hiding in the unusable-reason column all along.**

DRAFTS marks three unusable reasons. Two (`no_DNA_counts`, `low_DNA_counts`) mean the template was absent or insufficient — a failure to *measure*, the analogue of Johns's not-usable. The third, **`no_RNA_counts`, means template present but zero RNA detected** — which is exactly Johns's `rna == 0` inactive class. The difference is what each dataset then does with those rows:

| | Johns (in vivo) | DRAFTS (cell-free) |
|---|---|---|
| template present, RNA detected | usable **and** active | **usable** |
| template present, **no RNA** | usable but **inactive** (`tx_norm` exactly 0, row kept) | **excluded** as `no_RNA_counts` |
| template absent/low | not usable | not usable |

**So DRAFTS-"usable" already applies the `rna>0` filter. It is the analogue of Johns-ACTIVE, not of Johns-usable.** This is a naming collision between two datasets from the same laboratory — the same class of trap as the *P. agglomerans* / *P. aeruginosa* "Pa" collision caught in Gate 10, and it is exactly why the `n_shared_usable` column name in `gate10_crosshost_correlations.csv` invites the reading that prompted this gate. The column name is misleading; the computation underneath it was not.

**Distribution shape.** Among RNA-detected rows the activity distribution is **continuous, not bimodal**: zero exact zeros in any species, minima of 0.0045–0.057, and no detection-floor pile-up — the most common single rounded value accounts for at most **0.60%** of a species' rows (vs. the in-vivo `protein_log10` floor, which pins 66.5%/89.9%/10.1% of rows for EC/BS/PA). There is no cell-free analogue of the in-vivo translation floor.

**The inactive class exists but is nearly vacuous — and this is the substantive finding, not a technicality:**

| | inactive fraction |
|---|---|
| DRAFTS *Cg*, *Ll* | 0.0% |
| DRAFTS *Pa*, *Vn*, *Ef*, *Pp* | 0.3–1.5% |
| DRAFTS *Bs*, *Ko*, *Ec* | 3.7–4.6% |
| DRAFTS *Se* | **45.9%** ← the one exception |
| **in vivo EC** | **38.9%** |
| **in vivo PA** | **16.7%** |
| **in vivo BS** | **72.0%** |

Cell-free lysates transcribe essentially everything with template present; living cells do not. For 9 of 10 DRAFTS species the silent class is under 5%. **The blanket statement "cell-free has no inactive class" would be wrong — *S. enterica* has 45.9%** — so the accurate statement is that it is near-vacuous for 9 of 10 species, *Se* excepted.

**Consequence for co-active restriction:** it is definable and I computed it, but in cell-free it removes ~5% of rows where in vivo it removes 74%. The two filters are the same *operation* applied to populations with radically different composition — which is why the comparison needs all thresholds shown, not one.

---

## TASK 2 — EC–BS under every restriction regime

Full shared set (1,047-sequence DRAFTS library) and the 82-sequence Johns-overlap set. In-vivo reference in the last column.

| regime | cell-free ρ (n), full set | cell-free ρ (n), overlap | in-vivo ρ (n) |
|---|---|---|---|
| **pooled-equivalent** (DNA-adequate both; inactives carried as 0) | **0.616** (862) | 0.551 (85) | **0.655** (14,088) |
| **co-detected** = the co-active operation — *what the manuscript used* | **0.677** (807) | **0.597** (82) | **0.258** (3,668) |
| both above own 25th percentile | 0.647 (499) | 0.593 (48) | — |
| both above own median | 0.592 (287) | 0.477 (28) | — |
| **matched to each host's own active fraction** (EC 61.1%, BS 28.0%, then intersect — the faithful structural analogue) | **0.509** (218) | 0.621 (22) | 0.258 (3,668) |
| **matched top-26% of both** (harshest; see caveat) | **0.386** (150) | 0.318 (11) | 0.258 (3,668) |

**The in-vivo 0.597-vs-0.258 comparison in the manuscript is regime-matched.** Script 91 filters cell-free on `Ec_usable & Bs_usable` (RNA detected in both) and in-vivo on `tx_usable & tx_active` in both hosts. Same operation, both sides.

**The caveat on `matched_26pct`, stated because it cuts against reporting the largest number:** in vivo, 26% is the *intersection* of two per-host filters that individually retain 61% (EC) and 28% (BS). Demanding the top 26% of *both* species truncates far more deeply than the in-vivo filter did, and it truncates a continuous distribution rather than removing a discrete zero-mass. It is the least faithful analogue and the least favourable number; it is reported here precisely because it is the least favourable, and it still sits 50% above 0.258.

**Range-restriction diagnostic** (mirroring Gate 8.6's in-vivo check): going pooled-equivalent → co-detected, the retained subset's IQR gets *wider*, not narrower (Ec 0.558→0.681, Bs 0.242→0.299), so mechanical range restriction cannot explain any drop there — the same direction Gate 8.6 found in vivo. For `matched_26pct` the absolute IQR balloons (7.6, 13.0) because the retained top tail is heavy-tailed, but that diagnostic is not informative for a rank statistic under a deep quantile cut, and I am not leaning on it.

### The asymmetry is the actual finding

Conditioning on co-activity does opposite things in the two modalities:

- **in vivo: 0.655 → 0.258 (−61%)** — the pooled correlation was carried mostly by agreement on *which sequences are silent* (39%/72% exact zeros aligning). Remove the zeros and little graded agreement remains.
- **cell-free: 0.616 → 0.677 (+10%)** — there is almost no silent class to carry anything, so pooled ≈ graded, and the graded agreement is genuinely high.

That divergence, not the raw level, is what the modality comparison is really measuring.

### All 45 pairs, restated per regime (full shared set)

| regime | range | mean |
|---|---|---|
| pooled-equivalent | 0.364 – 0.892 | 0.710 |
| co-detected *(the published "0.623–0.911" band)* | 0.623 – 0.911 | 0.809 |
| both above p25 | 0.467 – 0.880 | 0.724 |
| both above median | 0.334 – 0.835 | 0.626 |
| matched top-26% | 0.214 – 0.804 | 0.515 |

**Honest caveat: under the harshest cut, 4 of 45 pairs fall below the in-vivo 0.2575** — Ec–Ll (0.214), Ef–Ll (0.245), Pp–Ll (0.247), Pa–Ll (0.250), all *L. lactis*, already the lowest-correlating species in Gate 10. Under every other regime, 0 of 45 fall below it.

**EC–BS is not a favourable pick.** Among the 45 pairs it ranks **3rd lowest under co-detected** (0.677) and 15th of 45 under both `both_above_med` and `matched_26pct` — mid-pack at worst, near the bottom at the regime the manuscript actually used. The published comparison uses one of the *weakest* cell-free pairs, not a flattering one.

---

## TASK 3 — Verdict

**SURVIVES.** Cell-free EC–BS never approaches the in-vivo co-active 0.258: the minimum across all six regimes is 0.386, and at the regime-matched comparison the manuscript actually used it is 0.677 (full set) / 0.597 (overlap set) — 2.6× and 2.3× the in-vivo value. The mechanism is now better understood than when the claim was written: the difference is in *restriction sensitivity*, because cell-free has almost no silent class while in vivo most sequences are silent and the two hosts disagree about which.

**Three things the manuscript must change anyway, none of which the verdict excuses:**

1. **Report pooled-versus-pooled.** 0.616 cell-free vs 0.655 in vivo — essentially identical, no modality difference at that regime. The manuscript reports neither number and a reader who computes them will conclude, reasonably, that the finding was cherry-picked by regime. It wasn't, but the paper currently gives no way to tell.
2. **Name the regime on every correlation figure.** See Task 4 — 16 figures are co-active, 6 pooled, 6 co-detected, and the manuscript labels almost none of them.
3. **Stop calling the DRAFTS column "usable" without qualification.** Add the cross-dataset definition collision to `KNOWN_ISSUES.md`, alongside the "Pa" collision — it is the same failure mode and it nearly produced a false retraction here.

**What would change the verdict:** if the in-vivo co-active filter were replaced with a graded-activity-matched filter that also removed the bottom ~74% of *cell-free* sequences by rank, cell-free EC–BS reaches 0.386 and the gap narrows to 1.5×. That is a legitimate alternative framing, it is reported above, and it does not reach parity — but a reader who prefers it should be able to find it, which is the point of reporting all six.

---

## TASK 4 — Restriction regime audit

Full table: `out/results/gate15_restriction_regime_audit.csv` (50 figures). Summary of what each manuscript number is:

| regime | count | examples |
|---|---|---|
| **co-active** (in vivo) | 16 | all six headline cross-host correlations (EC–PA 0.754/0.742, EC–BS 0.258/0.161, BS–PA 0.257/0.263), all six GC-controlled versions, both phylum strata |
| **co-active (derived)** | 5 | all four gap ratios (0.341, 0.308, 0.286, 0.240), disattenuated gap ratio |
| **pooled** (in vivo) | 6 | the pooled counterparts (EC–BS 0.655, BS–PA 0.508, EC–PA 0.621 + translation) — reported only in §3.6, never in §3.1 |
| **co-detected** (cell-free) | 6 | 45-pair range 0.623–0.911, same/cross-phylum means 0.852/0.769, *B. subtilis* mean 0.783, headline 0.597 |
| **co-detected in both modalities** | 1 | RS234 within-species 0.689–0.901 (internally regime-matched, both arms use DRAFTS's sentinel convention) |
| **both, explicitly labelled** | 3 | §3.6's pooled→co-active pairs — the only place the manuscript does distinguish |
| **n/a — model vs truth** | 7 | H-MAIN 0.213/0.218, conditioning mechanisms, shift-prediction — not measurement correlations |
| **new this gate** | 6 | the regime sweep above, not yet in the manuscript |

The pattern: §3.1 reports co-active figures without saying so; §3.6 reports the same pairs' pooled values and *does* say so; §3.2 mixes a cell-free co-detected figure with an in-vivo co-active figure, which is correct but reads as if it weren't. Every one of these is a labelling problem, not a computation problem.
