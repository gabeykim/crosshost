# CROSSHOST preprint package

Assembled Gate 11. Nothing here has been submitted or posted anywhere — see `out/SHIP_CHECKLIST.md` for what remains and who does it.

| File | Contents |
|---|---|
| `MANUSCRIPT.md` | The manuscript, with inline `(Figure N)` call-outs added to Results section headers for this submission-formatted copy — `out/MANUSCRIPT.md` (the project-root working copy) has identical content and provenance but no inline figure numbering, by design; both are kept in sync, changes go in `out/MANUSCRIPT.md` first. |
| `FIGURE_LEGENDS.md` | Captions for every figure, main text and supplementary — each states N and interval type explicitly. |
| `FIGURE_AUDIT.md` | Provenance verification for all 31 candidate figures: 28 included, 3 excluded for cause (two visualize a retired metric, one is an orphaned file superseded by a regenerated version) — read this before trusting any figure's inclusion. |
| `figures/` | 9 main-text figures (`Figure1_*.png` … `Figure7b_*.png`). |
| `figures/supplementary/` | 19 supplementary figures (`FigureS1_*.png` … `FigureS10_*.png`). |
| `SUPPLEMENTARY_TABLES.md` | Reviewer-facing renderings of the load-bearing tables (GC-confound, phylum-stratified, DRAFTS correlation matrix, RS234 modality agreement), pointers to the full machine-readable CSVs/JSONs, and the complete reference list. |
| `DATA_AND_CODE_AVAILABILITY.md` | Repo URL, archive-deposit structure and status, licensing quarantine, primary-data provenance. |
| `REPRODUCIBILITY.md` | What `make reproduce` regenerates (verified this gate, extended to cover Gates 8.5–10.5) and what `make reproduce-full` does not do (never re-executed end-to-end, stated plainly). |

All figures are 150 DPI PNGs — screen/preprint-appropriate, below the ≥300 DPI some print journals require (see `FIGURE_AUDIT.md`, "Resolution"). Regenerating at higher DPI is a one-line change per producing script, not attempted this gate to avoid touching result-generating code outside this gate's assembly-only scope.
