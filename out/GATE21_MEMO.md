**Are all nine images present, captions in ascending order, and glyphs intact? YES on all three — 9 embedded images, captions 1 through 8 at strictly ascending positions, and every glyph count unchanged from the Gate 20 build. But the gate's diagnosis was wrong and its prescribed command would have made things worse: the floats were already pinned, the ordering came from the Markdown source, and adding `-V header-includes` on the command line *replaces* the manuscript's own YAML block rather than adding to it — which would have silently reinstated the duplicate caption labels Gate 20 removed. Fixed by moving one image line. No prose, no number, no caption changed.**

# GATE 21 MEMO — Figure ordering and the final build

Script: `scripts/102_gate21_reorder_figure4.py`. Both audits clean at start and end.

---

## TASK 1 — The diagnosis was wrong, in two ways

### It was not float placement

`\usepackage{float}` and `\floatplacement{figure}{H}` have been in the preamble **since Gate 20**, carried by the YAML `header-includes` block at the top of `MANUSCRIPT.md`. Dumping the generated LaTeX confirms both at preamble lines 71–72. With `H` placement a float cannot move at all, so float placement could not have been the cause and re-pinning it could not have been the cure.

The real cause is **Markdown source order**. Gate 20's rule was "place each figure after the paragraph it illustrates," and §3.5 presents its four lines of evidence in the order host-tag → ablation → conditioning mechanisms (Figure 5) → pre-registered kill gate (Figure 4). The generated `.tex` therefore emitted `Figure5` before `Figure4`, and the PDF faithfully reproduced that. The observation that Figure 5's caption refers back to a Figure 4 the reader has not met is correct, and is the clearest symptom of the source ordering.

### The prescribed command would have reintroduced a fixed defect

Running the gate's invocation and inspecting the preamble:

```
71:\usepackage{float}\floatplacement{figure}{H}
```

— and **nothing else**. The `\usepackage{caption}` / `\captionsetup{labelformat=empty,...}` lines are gone. A `-V header-includes:` on the command line does not append to the document's YAML `header-includes`; it **replaces** it. The caption setting Gate 20 added to suppress LaTeX's automatic `Figure N:` prefix would have been dropped, and every caption would once again have rendered as `Figure 1: **Figure 1 — Raw cross-host measurement correlation…**`. Figure order would have been unchanged, since it never depended on the preamble.

The correct command is therefore the Gate 20 one, unchanged — the manuscript's own YAML already supplies float + `H` *and* the caption fix:

```
pandoc out/PREPRINT/MANUSCRIPT.md -o manuscript.pdf --pdf-engine=xelatex \
    -V mainfont="Times New Roman" --resource-path=out/PREPRINT
```

`float` was **already installed** in this TeX Live 2026 BasicTeX (`/usr/local/texlive/2026basic/texmf-dist/tex/latex/float/float.sty`). No `tlmgr install`, no `sudo`, no substitute approach.

### The fix

The only change that restores ascending order without touching prose is to move the **Figure 4 image line** above the Figure 5 image line. It now sits immediately after §3.5's opening paragraph — the section is titled "Testing the prediction (Figures 4–5)" and Figure 4 is the pre-registered test it is named for, so leading with it reads correctly, and Figure 5's "same resampling scheme as Figure 4" now resolves backwards.

The two alternatives were both rejected as out of scope: reordering the prose paragraphs is a manuscript text change this gate forbids, and renumbering the figures would alter captions.

**The whole diff is one image line and its blank line, moved.** The script asserts it, and three independent checks confirm it against `HEAD`:

| check | result |
|---|---|
| prose lines (image lines excluded) | **identical**, 172 lines |
| every number in the document | **identical** |
| set of image lines | **identical** |
| `git diff --stat` | 4 insertions, 4 deletions across the two copies — one line plus one blank, each way |

---

## TASK 2 — Verification

| # | check | result |
|---|---|---|
| 1 | resource warnings | **0** — pandoc emitted no output at all; no figure replaced with a description |
| 2 | `pdfimages -list \| awk 'NR>2 && $3=="image"' \| wc -l` | **9** ✓ |
| 3 | caption order | **1–8 ascending** ✓ — extracted-text lines 265, 293, 325, 442, 478, 510, 550, 593 |
| 4 | glyphs | ρ 21, ≤ 3, γ 2, β 2, ≈ 2 — unchanged from Gate 20; see below |
| 5 | pages | **22** (was 23); §3.3 table still renders as an aligned table |
| 6 | placeholders | `[email]` **0 occurrences**; `gabeykim@stanford.edu` present; `[Zenodo DOI]` present ✓ |

### Check 3: the gate's exact command reports 3 of 8, and that is a grep artifact

`pdftotext -layout manuscript.pdf - | grep -n "^Figure [0-9] —"` returns only Figures 3, 7 and 8. Nothing is missing. `pdftotext` prefixes the first line of every page with a form feed (`\f`), and five of the eight captions happen to start a page — an `od -c` dump shows `2 6 5 : \f F i g u r e` for Figure 1. The `^` anchor fails on those five. Dropping the anchor returns all eight, ascending. Worth fixing in any future check script; nothing to fix in the PDF.

### Check 4: the expected glyph counts are the pre-figure values

12 / 3 / 2 / 2 / 1 were the counts **before** captions were embedded in Gate 20. The captions themselves contain 9 ρ and 1 ≈, so the document now extracts 21 and 2. Both figures are unchanged from the Gate 20 build, so **nothing regressed in this gate**; the arithmetic (12+9=21, 1+1=2) was established and reported in Gate 20.

### Check 5: the page count dropped by one, and nothing was lost

23 → 22. This is pure repagination: Figure 4 (2100×900) is shorter than Figure 5 (2400×1350), so moving it earlier packs one page tighter. Verified by extracting the text of the committed Gate 20 PDF and of this one, stripping page numbers and normalising whitespace, and sorting: **639 lines each, zero diff.** No content moved out of the document.

Figures now land on pages 7, 8, 8 (both Figure 2 panels together), 9, 12, 13, 14, 15, 16.

---

## TASK 3 — Build command recorded

`make pdf` added to the Makefile, with `make manuscript` (added in Gate 20) kept as an alias so nothing that referenced it breaks. Verified by deleting `manuscript.pdf` and rebuilding it from scratch: 22 pages, 9 images, 8 captions.

The target carries two comments that cost real time to discover:

- `--resource-path` is required, because pandoc resolves relative image paths against the **working directory**, not the input file (Gate 20).
- **Do not add `-V header-includes:`** — it replaces the manuscript's YAML block rather than appending to it, and silently drops the caption setting (Gate 21).

It also notes that `xelatex` lives in `/Library/TeX/texbin`, which is not on the default macOS `PATH`.

---

## Verdict

Nine images, eight captions in ascending order, glyphs intact, prose and numbers untouched, build reproducible from a clean tree via `make pdf`. The gate's two prescriptions — re-pin the floats, add `-V header-includes` — were respectively a no-op and a regression; both are reported above rather than applied.
