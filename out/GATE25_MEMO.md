**Was 0.597 stale? Yes — it is the Gate 10 value, DRAFTS EC–BS on the 82-sequence Johns overlap, computed before Gate 15 defined the restriction regimes, and it appears nowhere in the manuscript. How many figures disagreed with the text? Seven of eight carried at least one real discrepancy: 23 findings, 2 blocking, 11 major, 5 minor, 5 clean. Does the modality contrast survive GC control? Yes — the co-active ratio goes from 2.63× raw to 2.69× GC-controlled, because both modalities attenuate by almost the same proportion (−17.9% cell-free, −19.8% in vivo).**

# GATE 25 MEMO

Scripts: `105` (GC control), `106` (figure regeneration), `107` (text corrections), `108` (reconciliation table). Both audits clean at start and end.

**The manuscript moved under this gate.** Ten user commits landed after Gate 24 (`9803e29` "Final manuscript" through `df7ccf4`), including a Zenodo DOI, a rewritten AI-assistance disclosure, a repair of the dangling antecedent Gate 24 reported, and a complete rewrite of the manuscript's embedded figure captions — which now split Figure 2 into **Figure 2a** and **Figure 2b** as separate captioned figures. Everything below is against that version (5,691 body words at gate open).

---

## BLOCK C — the decisive result, first

**The modality contrast survives GC partial-correlation control.** `out/results/gate25_modality_gc_control.{json,csv}`.

| | n | raw ρ | GC-controlled ρ | change |
|---|---|---|---|---|
| cell-free co-active | 807 | 0.677 | **0.556** | −17.9% |
| in-vivo co-active | 3,668 | 0.258 | **0.207** | −19.8% |
| cell-free pooled | 862 | 0.616 | 0.561 | −8.9% |
| in-vivo pooled | 14,088 | 0.655 | 0.517 | −21.1% |

**Co-active ratio: 2.63× raw → 2.69× GC-controlled.** Pooled: 0.94× → 1.09×, still no modality difference. The contrast is not a GC-composition artifact: both modalities lose a similar fraction to the control, so the ratio is preserved.

Computed by the two methods Gate 10.5 used, which agree closely (formula 0.556 / residual 0.548 cell-free; 0.207 / 0.205 in vivo). The script reproduces all four published raw values exactly before controlling for anything and refuses to run if it does not. **Not written into the Results**, per instruction.

One scope limit worth stating if this goes in the paper: the two modalities are measured on different sequence populations, so GC is controlled *within* each modality separately. This is not a paired control.

---

## BLOCK A — figures

All eight regenerated from current committed outputs by `scripts/106_gate25_regenerate_figures.py`. Full table in `out/results/gate25_figure_text_reconciliation.csv` (23 rows).

### A1 Figure 3 [BLOCKING] — confirmed and fixed

0.597 is `comparison_a_own_computation.cell_free_DRAFTS_EC_BS` in `gate10_modality_comparison.json`: n=82, ρ=0.5972649, the DRAFTS–Johns overlap. **Superseded** by `gate15_coactive_modality.csv`. Panel A did show three bars (0.597, the n=15 overlap at 0.147, full-library in-vivo 0.258) while the caption described four values under two regimes, and panel B went unmentioned.

Panel A now shows the four values the text reports, each labelled with its N. The n=15 bar is gone; it remains in Limitations item 5. Panel B keeps DRAFTS's RS234 within-species agreement, now with per-species N, and the caption describes both panels.

**The manuscript's embedded Figure 3 caption was already correct** — it described pooled and co-active with all four n values. The figure was the stale half; regenerating it brings the two into line.

### A2 Figure 8 — one report refuted, one real defect found

**"12 bar groups" is wrong.** The source CSV has exactly 11 rows and the figure draws 11 groups, matching the text. But the labels abbreviated **both** "transcription" and "translation" to `tr` via `readout[:2]`, so `PA_to_EC tr with` appeared twice for two genuinely different cells. That is almost certainly what produced the miscount. Now `tx` / `tl`.

The title/text inversion is real and fixed: the title now states both directions rather than only "1 of 11 survives".

### A3 Figure 5 — confirmed on both counts

The figure shows **six** systems (sequence-only, free host embedding, genomic + FiLM, physiology + FiLM, DNABERT-2, PromoGen2) and has **no panel (b)**; fold-to-fold variance is Figure 7. Both captions — the legends file and the manuscript's own — said otherwise. Both corrected, and all six panels now share one y-scale.

**Reported, not fixed:** §3.5 cites this figure for concatenation 0.555 and per-host-heads 0.615. Neither system is plotted here. The caption now says they are tabulated rather than plotted; whether to add them is a content decision.

### A4 Figure 2b — confirmed

Four strata, not two. All values confirmed against `gate10_5_gc_confound.json`:

| stratum | EC–PA | EC–BS | BS–PA |
|---|---|---|---|
| Proteobacteria | 0.634 (n=3,282) | 0.090 (n=321) | −0.209 (n=327) |
| Firmicutes | 0.812 (n=2,666) | 0.177 (n=2,532) | 0.271 (n=1,214) |
| Actinobacteria | 0.693 (n=824) | 0.092 (n=97) | −0.174 (n=119) |
| Bacteroidetes | 0.536 (n=842) | −0.010 (n=146) | 0.018 (n=81) |

§3.1 now reports all four with N, and **"with *B. subtilis* pairs near zero" is corrected**: the range is −0.209 to +0.271. The four EC–PA values (0.536–0.812) bracket 0.754, so there is no anomaly.

**One correction to the brief:** 0.754 is the **co-active** EC–PA value, not the pooled one (pooled EC–PA is 0.621). The strata are computed on the co-active mask, so co-active is the right comparator.

**Reported, not fixed:** the results file holds **seven** transcription strata; the figure hardcodes four. Cyanobacteria clears the threshold on all three pairs (EC–BS 0.332 n=62, EC–PA 0.543 n=455, BS–PA −0.022 n=44) and could be shown. Changing which strata appear is a content decision.

### A5 Figure 1 — confirmed, with the diagnosis corrected

The bars do sit flat at 1.0 and the true values are 1.51 / 1.08 (transcription) and 1.48 / 1.06 (translation). **But this is not axis clipping.** `scripts/72_attenuation_analysis.py` stores `float(np.clip(corrected, -1.0, 1.0))` and records `clipped_at_1`, so the cap is in the data. Extending the y-axis alone changed nothing — the first regeneration looked identical, which is how I found it.

The figure now draws a dotted line at ρ = 1 and prints the uncapped value above each genuinely clipped bar. The in-figure title is descriptive rather than assertive.

**A second defect found on the way:** `clipped_at_1` is serialised as the **string** `"True"`/`"False"` by `json.dump(default=str)` on a numpy bool, so any plain truth test matches both — my first attempt annotated every bar. Reported, not fixed: rewriting the results file would change a committed artifact.

### A6 — pipeline nomenclature stripped

"Gate 8 Task 1", "Gate 6", "Gate 7 Task 2", "Gate 8.6 Task 1" removed from Figures 1, 5, 6 and 8. No other figure carried one.

### A finding no report mentioned: the figures had no producing script

**No script produced `out/PREPRINT/figures/*.png`.** The producers write `out/figures/gate*.png`; Gate 12 renamed copies into the preprint directory by hand. `audit_provenance.py` passed only because the Gate 20 embed script mentions every preprint filename — **a false pass on all nine files, every gate since Gate 20**. `scripts/106` now regenerates and places all nine, so the provenance is real.

That script also caught itself: its first version imported each producer and called `main()`, but only one of six defines `main()`, so five silently did nothing and it copied their stale output — reporting success. It now runs each producer as a subprocess and asserts the output file was rewritten.

---

## BLOCK B — text

| | status |
|---|---|
| **B1** typo | **REFUTED — no edit made.** The manuscript correctly says "BS 15,848"; `215,848` appears nowhere. The PDF breaks the page between "BS" and "15,848" and the page number **2** falls between them. Making the "fix" would have introduced an error. |
| **B2** supercoiling | **CONFIRMED wrong — wording proposed below, not applied.** |
| **B3** "independent dataset" | applied → "the cell-free dataset" |
| **B4** "checks they did not run" | applied → "not previously applied to this comparison" |
| **B5** title | **options proposed below, none chosen** |
| **B6** "opposite directions" | applied in the abstract and §3.2 |
| **B7** population framing | applied at the §3.2 headline |
| **B8** *L. lactis* | applied to §3.3 and Limitation 2 |
| **B9** prereg timestamp | applied |

### B2 — the supercoiling claim is wrong

DRAFTS Methods: *"Plasmid DNA was extracted (Zymo Midiprep) and column-purified once more (PureLink, Invitrogen) to avoid carryover of RNase A into cell-free expression systems."* The library is supplied as **midiprepped plasmid**, which is covalently closed and supercoiled. The word "supercoil" appears nowhere in the DRAFTS paper. §3.4's "DNA supercoiling, absent in linear cell-free reactions" is factually wrong.

**Proposed wording, for your decision — not applied:**

- §3.4: *"DNA supercoiling is present on the plasmid templates DRAFTS used, but a lysate has no gyrase/topoisomerase homeostasis and no chromosomal context, so the supercoiling state is fixed at whatever the preparation produced rather than dynamically regulated. El Houdaigui et al. (2019) show that a promoter's response to supercoiling depends quantitatively on its discriminator GC content, so a fixed rather than regulated state is still a removed variable — but a weaker claim than absence."*
- §3.2: replace "native supercoiling" in the list of what lysates lack with "supercoiling homeostasis".

I did not apply either, per instruction. Note this weakens one of the three named mechanisms rather than removing it.

### B5 — title options, none chosen

1. *Cell-free systems do not reproduce the bacterial chassis effect in transcription: evidence that host specificity resides in cellular context*
2. *Cell-free transcription does not reproduce the bacterial chassis effect: evidence that host specificity resides in cellular context rather than transcription machinery*
3. *Transcriptional host specificity survives in cells but not in lysates: evidence that the chassis effect is contextual*

Option 2 is the smallest change and puts the scope limit where a reader meets it first.

---

## Verification

| check | result |
|---|---|
| resource warnings | 0 |
| images | 9 |
| caption order | Figures 1, 2a, 2b, 3, 4, 5, 6, 7, 8 — all present, ascending |
| glyphs | ρ 17, ≤ 3, γ 2, β 2, ≈ 1 — all present, non-zero |
| pages | 20 |
| §3.3 regime table | renders |
| `[email]` | 0; Zenodo DOI present |
| cross-references | all resolve |
| both copies | byte-identical |
| audits | leakage ALL 6 PASSED; provenance 0 orphans (176 files) |
| `verify-citations` | 11 VERIFIED / 0 unexplained / 3 explained / **2 skipped** — the two skips are HTTP 429 rate-limiting from the citation APIs, not a content change. Baseline is 11/0/5/0/16 and no citation was touched this gate. |

**Numbers added to the manuscript:** the A4 stratum values and their N (authorised), plus **10,420** and **55** in §3.2 — both arithmetic on numbers already present (14,088 − 3,668 and 862 − 807), required to state the asymmetry B6 asks for. Flagged because they fall outside the strict "A4 and B1 only" allowance. **No number was removed.**

---

## Verdict

0.597 confirmed stale and replaced. Seven of eight figures carried a real discrepancy; all fixed or reported. The modality contrast survives GC control at 2.69×. Two reports refuted (B1, Figure 8's bar count), one diagnosis corrected (Figure 1 is data clipping, not axis clipping), and a provenance false pass that predates this gate is now closed.
