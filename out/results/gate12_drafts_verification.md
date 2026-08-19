# Gate 12, Task 0 — DRAFTS load-bearing claim verification

**Claim being checked:** DRAFTS (Yim, Johns et al., *Mol Syst Biol* 15:e8875, 2019) does not compare between-species (cross-host) similarity between the cell-free and in-vivo modalities, and does not claim cell-free reduces host-to-host differences. What it does report, per the second-hand framing this gate opened with: within-species in-vitro/in-vivo fidelity only (Pearson 0.71-0.90 across seven species; ~94.2% of *E. coli* activities within 1-log variation), plus predictive features (sigma70 motif strength, GC content) for cross-species activity **within** the cell-free modality alone.

**Method.** Full text of `raw/drafts/msb198875_FullTextArticle.nxml` (structured XML, tags stripped) and `raw/drafts/msb198875_Appendix.pdf` (via `pdftotext`) searched exhaustively for any statement, table, or figure that compares cross-species/cross-host correlation structure *between* the two modalities, as opposed to (a) within-species modality fidelity or (b) within-modality cross-species structure. Every surface listed in the gate prompt checked individually below, findings recorded as they are found — this file is written incrementally, not assembled at the end, per the recovery instruction.

**Status: COMPLETE.**

## VERDICT: **CONFIRMED**

DRAFTS does not, anywhere in the full text, the Appendix (Methods, all 6 tables, all 14 figures), the five source-data Excel files' sheet and column structure, or the peer review correspondence, compare cross-species (between-host) correlation structure between the cell-free and in-vivo modalities. It reports two distinct, unconnected things: (1) **within-species** in-vitro/in-vivo fidelity (RS234, seven species, Pearson r 0.71–0.90 per the paper's own text; 94.2% of *E. coli* RS29249 activities within 1-log variation), and (2) **within-cell-free** cross-species structure (the 45-pair RS1383 matrix across 10 species, sigma70/RpoD and GC content as predictive features of cross-species activity). These two findings sit next to each other in the Discussion but are never combined into the comparison this project's manuscript makes. The claim that "cell-free lysates largely abolish the cross-host differences that dominate in-vivo measurements" is **this project's own inference**, built by comparing DRAFTS's cell-free cross-species numbers against this project's own separately-computed in-vivo cross-species numbers (Johns et al. 2018 data) — DRAFTS's authors never made or implied this comparison themselves, and the manuscript's existing framing already reflects this correctly (`out/DRAFTS_SECTION.md`: "we compared our in-vivo measurements against DRAFTS... an independent dataset"). Task 0 finds no reason to revise that framing; it finds strong, exhaustive, multi-surface support for it.

**Proceeding to Task 1.**

---

## Surfaces checked

### 1. Abstract (full text, re-verified after crash)

Full text quoted:

> "Interspecies analysis of transcriptional profiles from >1,000 diverse regulatory sequences reveals functional differences in promoter activity that can be quantitatively modeled, providing a rich resource for tuning gene expression in diverse bacterial species."

This is the only sentence touching cross-species structure. It describes the **cell-free interspecies analysis on its own terms** (matches this project's own 45-pair correlation matrix work in Gate 10). **No comparison to in-vivo cross-species structure is made or implied.** No mention of "context," "compress," "reduce," or any modality-comparison framing at the abstract level.

**Finding: no cross-modality cross-species claim in the Abstract.**

### 2. Discussion (full text, re-verified after crash)

Relevant sentences, quoted verbatim and in original sequence (two adjacent but topically distinct sentences — this adjacency is exactly the kind of thing a careless read could conflate, so quoting both in full and separately):

> "In vitro transcriptional measurements faithfully recapitulated in vivo expression levels, thus avoiding laborious construction of multiple species-specific libraries and transformation into some challenging or recalcitrant species."

This is the **within-species modality-fidelity claim** (a given species' own cell-free measurement vs. its own in-vivo measurement — the RS234/Fig S5 comparison). It says nothing about cross-species structure.

> "We observed increased dissimilarity of transcriptional profiles between phylogenetically more distant species, suggesting a functional evolutionary divergence of transcriptional regulation across the natural microbial biome (Iyer et al, 2004)."

This is the **within-cell-free cross-species claim** (same-phylum vs. cross-phylum correlation pattern in the 45-pair cell-free matrix — matches this project's own Gate 10 finding of same-phylum mean 0.852 vs. cross-phylum mean 0.769). It is entirely about structure *within* the cell-free modality. No in-vivo cross-species comparison is invoked, and no claim is made that this cell-free cross-species pattern does or does not resemble an in-vivo one.

**Finding: the Discussion places these two claims (within-species fidelity; within-modality cross-species divergence) back to back but never connects them into a between-modality cross-species comparison.** No sentence anywhere in the Discussion states or implies that cell-free measurements compress, reduce, preserve, or otherwise relate to in-vivo cross-species differences.

### 3. Methods — the in-vivo measurement protocol

Searched for every "in vivo" occurrence (22 total in the full text) and the Methods description of the in-vivo protocol specifically:

> "For in vivo measurements, library cultures of all bacterial species were grown in rich media (Appendix Table S1) until mid-exponential phase (OD600 ~0.2). RNA was then extracted using RNAsnap (Stead et al, 2012), and cells were lysed using prepGEM bacteria kit (MicroGEM) for amplification of input DNA library sequences."

This describes how each species' own in-vivo sample was prepared (feeding the RS234 within-species comparison, Fig S5/Appendix Fig S5) — a protocol description, not an analysis. **No Methods subsection describes a cross-species in-vivo analysis, and no cross-species in-vivo dataset broader than RS234's seven species (used only for the seven separate within-species correlations) exists in this paper.**

### 4. Figure S5 caption (RS234, the within-species modality-fidelity comparison)

Full caption, from the Appendix PDF:

> "Appendix Figure S5. Comparison of in vitro and in vivo transcriptional measurements in 7 bacterial species. Transcriptional profile (Tx) correlations between in vivo and in vitro measurements in E. coli, S. enterica, K. oxytoca, P. putida, V. natriegens, B. subtilis, and C. glutamicum. All transcriptional measurements were made using template-switching adaptor ligation. Shaded regions represent 95% confidence interval for linear regression (dashed lines). Sample sizes (n) and Pearson correlation coefficients (r) can be found in each plot. For normalization, transcription levels in log10 scale were transformed to Z-score. All measurements are based on two biological replicates."

**Structurally seven separate per-species scatter plots (one per species, in-vitro-vs-in-vivo), not a cross-species matrix.** Confirmed independently by the underlying source-data sheet structure (Surface 8, below): the "Appendix Fig. S5" sheet in `msb198875_SourceData_Fig2.xlsx` is organized per-species (`Ec_invitro_tx`, `Ec_invivo_tx`, `Se_invitro_tx`, `Se_invivo_tx`, ...) with zero cross-species columns of any kind — the raw data released for this exact figure does not even contain the ingredients for a cross-species in-vivo correlation.

**Ancillary finding, not load-bearing for the CONFIRMED/NOT CONFIRMED verdict but worth flagging for consistency:** the paper's own text states "Pearson's r between 0.71 and 0.9" for this seven-species comparison. This project's own Gate 10 recomputation (`out/GATE10_MEMO.md` Task 4.3, `out/results/gate10_modality_comparison.json`) found the range **0.69–0.90**, with *B. subtilis* specifically at 0.69 — 0.02 below the paper's stated floor of 0.71. This was already disclosed in `out/DRAFTS_SECTION.md` ("recomputed directly from DRAFTS's own released source data... not a literal quote... should closely track, but is not a literal quote of, the paper's reported values") and in `out/GATE10_MEMO.md`'s "WHAT I COULD NOT DO." No manuscript text anywhere quotes "0.71" as if it were this project's own number — the manuscript's Section 3.4 and `out/results/` consistently use the recomputed 0.69–0.90 figure. Flagged here for completeness; does not affect Task 0's verdict, since it concerns the magnitude of within-species fidelity, not whether a cross-species cross-modality comparison exists.

### 5. Figure S7 and S8 captions (the cell-free-only cross-species figures)

> "Appendix Figure S7. Comparisons of RS1383 transcriptional profiles across 10 bacterial species using DRAFTS. (a) Box plot of transcriptional activity (Tx) distributions for each species. (b) Scatter plot matrix of in vitro transcriptional activity comparisons for 10 species. (c) Principal components analysis of transcriptional profile similarity for 10 species."

> "Appendix Figure S8. RNA polymerase subunit sigma70 as an alternative metric of evolutionary divergence. (a) Unrooted phylogenetic tree of 10 bacterial species... (b) Correlation between amino acid sequence dissimilarity of RpoD and pairwise Pearson correlation of transcriptional profiles of 421 universally-active regulatory sequences between bacterial species."

Both explicitly scoped to **"in vitro"** (S7b) and cell-free transcriptional profiles (S8b) only. Neither figure has an in-vivo counterpart or comparison panel. **This is the cell-free 45-pair cross-species structure this project already used in Gate 10 (`gate10_crosshost_correlations.csv`) — confirmed as cell-free-only at the source, not compared to any in-vivo cross-species structure anywhere in the paper.**

### 6. Appendix Tables S1–S6, in full

- **Table S1** (species/strain/medium/temperature/aeration) — pure growth-condition metadata, referenced by the in-vivo Methods paragraph (Surface 3). No comparison data.
- **Table S2** (media/buffer/cell-free-component composition) — protocol recipe table. No comparison data.
- **Table S3** (cell-free system composition per species) — protocol table. No comparison data.
- **Table S4** (plasmids/libraries used) — inventory table (RS234/RS29249/RS1383/RS7003 provenance). No comparison data.
- **Table S5** (primers used) — inventory table. No comparison data.
- **Table S6** (biological-replicate reproducibility) — checked in full (all rows, all species/libraries). Every row is a **same-species, same-modality, replicate-1-vs-replicate-2 Pearson r** (e.g., "Figure 1D / E. coli in vivo / RS29249 / 0.8932"; "Figure 1D / DRAFTS E. coli / RS29249 / 0.9333"). **Zero cross-species rows, zero cross-modality rows.** This is the table this project already used for the *E. coli* replicate-reliability figures (Gate 10) — confirmed here that it contains nothing beyond that.

**No table anywhere in the Appendix contains a cross-species correlation computed on in-vivo data, or a side-by-side comparison of cross-species correlation between the two modalities.**

### 7. Appendix Figures S9–S14, in full

- **S9** (functional category vs. activity) — cell-free only (RS1383, 10 species), box plots by antibiotic-resistance gene class. No modality comparison.
- **S10** (linear regression error distribution) — cell-free only, model trained and evaluated within the cell-free RS1383/RS7003 data. No modality comparison.
- **S11** (regression model evaluated on a separate dataset) — still cell-free-to-cell-free (RS7003-trained model evaluated on RS1383). No modality comparison.
- **S12** (dual-species hybrid lysates, transcription) — cell-free only, single-species vs. hybrid-lysate Broccoli signal. Not an in-vivo comparison; a different experimental axis (mixing lysates from two species together) than this project's chassis-effect question.
- **S13** (dual-species hybrid lysates, translation) — same as S12, GFP instead of Broccoli. Cell-free only.
- **S14** (hybrid lysates and species-selectivity, GC content) — cell-free only. Compares E. coli-selective vs. B. subtilis-selective (and vs. C. glutamicum-selective) promoters' GC content **within the cell-free hybrid-lysate experiment**. This is the closest any figure in the paper comes to a "GC content relates to species-selectivity" framing, and it is exactly the finding this project's own Gate 10 GC-confound check (`gate10_gc_and_floor_check.py`) independently re-derived from the full library — but it is still entirely a cell-free-vs-cell-free comparison (original lysate vs. hybrid lysate), not a cell-free-vs-in-vivo one.

**No figure S9–S14 makes any in-vivo comparison of any kind.**

### 8. Source-data Excel files — sheet structure, all five files

All sheet names, `openpyxl`, read-only:

| File | Sheets |
|---|---|
| `Fig1.xlsx` | RS29249, Fig. 1B–1E, Appendix Fig. S2B–S2E, README |
| `Fig2.xlsx` | Fig. 2 (one sheet per of 10 species), Appendix Fig. S4A/S4B, **Appendix Fig. S5**, README |
| `Fig3.xlsx` | RS1383, Fig. 3, Fig. 3D, Fig. 3D Sequence, README |
| `Fig4.xlsx` | RS7003, Fig. 4 (one sheet per of 10 species), README |
| `Fig5.xlsx` | Fig. 5, Appendix Fig. S12, Appendix Fig. S13, README |

No sheet name anywhere suggests a cross-species in-vivo table. The one sheet that could plausibly hold it — "Appendix Fig. S5" in `Fig2.xlsx` — was opened and its column headers read directly (108 columns): every column is `{species}_{invitro|invivo}_{rep1|rep2|combined}_{DNA|RNA|tx|zscore}`, organized strictly per-species with no cross-species pairing of any kind. **The raw data released for the one figure most likely to contain this comparison structurally cannot support it** — there is no column, anywhere, that would let a cross-species in-vivo correlation be computed even by a third party trying to construct one from the released data.

### 9. Peer review correspondence (bonus surface, not explicitly requested but checked for completeness)

`raw/drafts/msb198875_ReviewProcessFile.pdf` (13,350 words, referee reports + author responses) searched for "cross-species," "between species," "host-to-host," and "chassis effect" — **zero matches for any of the four terms.** No reviewer raised, and no author response addressed, a cross-species cross-modality comparison at any point in review.

### 10. Final broad sweep

Searched all three extracted texts (full article, appendix, review correspondence — 14,789 + ~13,350 words combined) for: "reduce host," "reduces host," "narrower," "narrows," "abolish," "compress," "host-to-host," "host to host," "context-dependent," "cellular context." **Zero matches for every term, in every document.**
