**Do all eight figures and their captions now appear in the PDF? YES — all 8 figures (9 image files; Figure 2 is two panels) and all 8 verbatim captions are in `manuscript.pdf`. Page count 17 → 23, embedded images 0 → 9, no pandoc warnings. Two defects in the gate's prescribed method had to be corrected to get there, both found by converting rather than assuming: pandoc resolves relative image paths against the working directory, not the input file, so the given build command silently dropped every figure; and the prescribed alt-text template would have rendered each caption with its label three times over. One ordering anomaly is reported, not fixed: Figure 5 appears before Figure 4.**

# GATE 20 MEMO — Embed figures and finalize the title block

Script: `scripts/101_gate20_embed_figures.py`. Both audits clean at start and end. `verify-citations` 11/0/5/0/16, unchanged since Gate 13.

---

## TASK 1 — Inventory

**`out/PREPRINT/figures/` — 9 files, all PNG, all 150 DPI.** Nothing cited is missing; nothing present is uncited.

| file | px | in @150dpi | figure |
|---|---|---|---|
| `Figure1_disattenuation.png` | 1800×825 | 12.0×5.5 | 1 |
| `Figure2a_gc_control.png` | 1800×750 | 12.0×5.0 | 2 (panel a) |
| `Figure2b_phylum_stratified.png` | 1350×750 | 9.0×5.0 | 2 (panel b) |
| `Figure3_drafts_modality_comparison.png` | 1800×750 | 12.0×5.0 | 3 |
| `Figure4_hmain_kill_gate.png` | 2100×900 | 14.0×6.0 | 4 |
| `Figure5_conditioning_mechanisms.png` | 2400×1350 | 16.0×9.0 | 5 |
| `Figure6_calibration_collapse.png` | 2250×1350 | 15.0×9.0 | 6 |
| `Figure7_film_instability.png` | 1950×750 | 13.0×5.0 | 7 |
| `Figure8_shift_prediction_retraction.png` | 1950×900 | 13.0×6.0 | 8 |

**Nine files, eight figures.** Figure 2 is a two-panel figure sharing one caption. The gate's "eight expected" image count is therefore nine; `pdfimages -list` reports 18 rows, because every PNG carries an alpha channel and so contributes an `image` row plus a `smask` row.

**Captions:** all 8 main-text captions found in `FIGURE_LEGENDS.md`, each matched to its file, each stating N and interval type (Figure 8's N was added in Gate 19).

**Citations in `MANUSCRIPT.md`:** figure numbers appear *only* in section headings — §3.1 (Figures 1–2), §3.2 (Figure 3), §3.5 (Figures 4–5), §3.6 (Figure 6), §3.7 (Figure 7), §3.8 (Figure 8). There are no inline body-text callouts. Other "Fig." strings in the manuscript belong to *other papers* (Johns et al.'s Supplementary Figs. S13/S15; DRAFTS's Fig. 1D, Fig. 3, Appendix Fig. S5) and were left alone.

**Mismatches: none in the main text.** Two pre-existing conditions, both intentional and neither blocking:

- `FIGURE_LEGENDS.md` carries a **Figure S5 entry with no file**. That is the deliberate Gate 19 placeholder vacating the slot that collided with DRAFTS's own Appendix Fig. S5; the ablation figure is S11.
- The 19 files in `figures/supplementary/` implement Figures S1–S11 and are captioned, but **none is cited from the manuscript body** — as recorded in Gate 18. Supplementary figures were *not* embedded, per instruction.

---

## TASK 2 — Embedding

Each figure is inserted immediately after the paragraph it illustrates, located by a unique anchor string rather than a line number:

| figure | inserted after the paragraph beginning | PDF page |
|---|---|---|
| 1 | "**The gap survives measurement noise.**" | 7 |
| 2 | "**The gap survives a source-composition confound.**" | 8 |
| 3 | "**Conditioning on co-activity separates them, in opposite directions.**" | 9 |
| 5 | "**Sequence-only prediction is not distinguishably beaten…**" | 13 |
| 4 | "**The pre-registered kill gate failed on all 8 primary comparisons.**" | 14 |
| 6 | "The sequence-only model's zero-shot classifier is well calibrated…" | 15 |
| 7 | "**FiLM is measurably unstable…**" | 16 |
| 8 | "It did not survive its control." | 17 |

### Three corrections the prescribed method needed

**1. The build command drops every figure.** The gate states that pandoc resolves relative image paths against the input file. It does not — it resolves them against the **working directory**. Run verbatim from the repository root, the given command emitted nine `[WARNING] Could not fetch resource figures/… : replacing image with description` lines and produced a figure-free PDF with exit status 0. Fixed by adding `--resource-path=out/PREPRINT`. The paths themselves were left relative to the manuscript (`figures/…`), not rewritten to `out/PREPRINT/figures/…`, so that `out/PREPRINT/` remains a self-contained package for the Zenodo deposit and the figures render in any Markdown viewer.

**2. The prescribed alt-text template triples the caption label.** `![**Figure 3.** <caption verbatim>](…)` combines three labels: the gate's own `**Figure 3.**`, the legend caption's own opening `**Figure 3 — The DRAFTS modality comparison.**`, and LaTeX's automatic `Figure N:` prefix on every `\caption{}`. A test build produced `Figure 1: Figure 1 — Raw cross-host measurement correlation…`. Resolved by using the verbatim legend caption alone as the alt text — it already carries the `**Figure N — Title.**` label the template was trying to add — and suppressing LaTeX's automatic label. **No caption text was cut to achieve this**; the label the gate asked for is present, once.

**3. Figures floated away from their paragraphs.** With pandoc's default `htbp` placement, a test build put Figure 2's panel (a) on page 1 and panel (b) with its caption on page 2, after intervening body text. Fixed with `\floatplacement{figure}{H}`.

Corrections 2 and 3 need preamble packages, which the fixed build command cannot supply. They are therefore carried in a **five-line YAML `header-includes` block at the top of `MANUSCRIPT.md`**, so the documented command still works verbatim and the document stays self-contained. This is a **third change to the manuscript header beyond the two Task 3 asked for**, and is disclosed as such.

### Figure 2's two panels

Panel (a) is emitted as an uncaptioned image and panel (b) carries the single verbatim caption, so the caption sits below both panels and is not split across them. Both land on page 8. Every figure is scaled to the full 4.77-inch text width, so the uncaptioned panel and the captioned one align identically. Combining the two PNGs into one file was rejected: it would have created a figure with no producing script and an orphan in the provenance audit.

### Ordering anomaly — reported, not fixed

**Figure 5 appears on page 13 and Figure 4 on page 14.** §3.5 presents its four lines of evidence in the order host-tag, ablation, conditioning mechanisms (Figure 5), kill gate (Figure 4), so placing each figure at the paragraph it illustrates — as the gate instructs — puts them out of numeric order. The three ways out each break something the gate forbids: stacking both at the section start, renumbering the figures (altering captions), or reordering the two paragraphs (a text change this gate does not authorize). Left as-is for Gabriel to decide.

---

## TASK 3 — Title block

1. `[email]` → `gabeykim@stanford.edu`. Zero `[email]` occurrences remain; `[Zenodo DOI]` is untouched.
2. The affiliation/correspondence collapse is fixed with a trailing backslash after "USA" (pandoc hard line break). Confirmed in extracted text: `Stanford University, Stanford, CA, USA` and `Correspondence: gabeykim@stanford.edu` now occupy separate lines.

---

## TASK 4 — Build

```
pandoc out/PREPRINT/MANUSCRIPT.md -o manuscript.pdf --pdf-engine=xelatex \
    -V mainfont="Times New Roman" --resource-path=out/PREPRINT
```

Recorded as `make manuscript`. Tooling: pandoc (Homebrew), XeTeX 3.141592653-2.6-0.999998 (TeX Live 2026), `xelatex` at `/Library/TeX/texbin` — **not on the default `PATH`**, which is worth knowing for anyone reproducing this.

| | before | after |
|---|---|---|
| pages | 17 | **23** |
| embedded images (`pdfimages`, `type=image`) | 0 | **9** |
| `pdftotext \| grep -c "Figure"` | 6 | **18** |
| file size | 170 KB | 1,233 KB |

Every embedded image's pixel dimensions match its source file exactly. All 8 caption opening lines extract intact.

**Unicode glyphs — no regression. Every count is unchanged or increased, and each increase is exactly accounted for by caption text:**

| glyph | before | after | captions contain | check |
|---|---|---|---|---|
| ρ | 12 | 21 | 9 | 12 + 9 = 21 ✓ |
| ≤ | 3 | 3 | 0 | ✓ |
| γ | 2 | 2 | 0 | ✓ |
| β | 2 | 2 | 0 | ✓ |
| ≈ | 1 | 2 | 1 | 1 + 1 = 2 ✓ |

The gate's expected counts (12/3/2/2/1) are the pre-figure values; embedding captions that themselves contain ρ and ≈ necessarily raises two of them.

The §3.3 table and all 29 headings survive the rebuild unchanged.

---

## TASK 5 — Consistency

- **Cited = embedded.** {1,…,8} both ways; zero cited-but-not-embedded, zero embedded-but-not-cited.
- **Captions verbatim.** Each embedded alt text was checked, after whitespace normalisation, to be a substring of `FIGURE_LEGENDS.md`. All 8 match. The only transformation is joining each caption's soft-wrapped lines with a single space, which Markdown alt text requires; no word is added, removed, or reordered. The script asserts this and refuses to write if it fails.
- **No number changed.** Against `HEAD`, **zero numbers were removed** from the manuscript, and every number added is present in the `FIGURE_LEGENDS.md` caption text. No discrepancy to report.
- **Both copies byte-identical** — `cmp` clean, md5 `9d533d917f5f4df7ad9d88e17bf092c1`.
- **`make verify-citations`: 11 VERIFIED / 0 unexplained / 5 explained / 0 unresolved / 16 total** — no regression.

### One consequence of the byte-identical requirement

`out/MANUSCRIPT.md` carries `figures/…` paths that **do not resolve from `out/`**: `out/figures/` exists but holds the pre-Gate-12 filenames (`calibration_curve_BS_transcription.png`, …), not `Figure1_disattenuation.png`. The two requirements — paths relative to the PREPRINT copy, and both copies byte-identical — cannot both hold and also leave `out/MANUSCRIPT.md` buildable. `out/PREPRINT/MANUSCRIPT.md` is the build target; `out/MANUSCRIPT.md` is its mirror. Nothing builds from `out/`, so nothing is broken today, but it is a real asymmetry and the alternative is to stop mirroring.

---

## Verdict

Eight figures, nine files, eight verbatim captions, all present in a 23-page PDF. Three method corrections, one ordering anomaly, one mirror-path asymmetry — all reported above rather than worked around.
