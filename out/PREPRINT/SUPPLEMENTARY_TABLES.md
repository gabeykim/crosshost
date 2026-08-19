# Supplementary Tables

Full machine-readable versions of every table below are in `out/results/` (project root) or `package/baselines/` (the packaged, pip-installable distribution). This file gives readable, reviewer-facing renderings of the tables most directly load-bearing for the manuscript's claims; it does not duplicate every CSV in the repository.

## Table S1 — Master baseline suite

Every model evaluated in this project, fold-resolved, with 90% bootstrap CI: `package/baselines/master_baselines.csv` (108 rows: hosts × readouts × N-points × mechanisms × systems, spanning B1 mean/majority through both foundation models). See `package/baselines/README.md` for column definitions. Includes the calibration-failure numbers (`package/baselines/gate7_conformal_and_ece.json`) prominently, not appended as an afterthought.

## Table S2 — DRAFTS cross-host correlation matrix (all 45 pairs)

Full table: `out/results/gate10_crosshost_correlations.csv`. Summary: range [0.623, 0.911], mean 0.81; same-phylum mean 0.852 (n=17 same-phylum pairs), cross-phylum mean 0.769 (n=28 cross-phylum pairs). Highest pair: *P. agglomerans*–*V. natriegens* (0.911). Lowest pair: *E. coli*–*L. lactis* (0.623). Per-species mean correlation with the other nine species, lowest to highest: *L. lactis* 0.729, *E. coli* 0.758, *B. subtilis* 0.783, *V. natriegens* 0.812, *C. glutamicum* 0.820, *S. enterica* 0.831, *P. putida* 0.834, *P. agglomerans* 0.835, *K. oxytoca* 0.838, *E. fergusonii* 0.853.

## Table S3 — Raw vs. GC-controlled cross-host correlation (in-vivo, this project's own data)

Full table: `out/results/gate10_5_gc_confound.csv`, `.json` (Task 2). Two independent methods (closed-form partial-correlation formula; rank-residual regression) agreed to 3 decimal places in every cell; one value shown per cell.

| Readout | Pair | N | Raw ρ | GC-controlled ρ | % change |
|---|---|---|---|---|---|
| transcription | EC–BS | 3,668 | 0.258 | 0.207 | −19.8% |
| transcription | **EC–PA** | 9,741 | **0.754** | **0.718** | **−4.7%** |
| transcription | BS–PA | 2,099 | 0.257 | 0.236 | −8.2% |
| translation | EC–BS | 866 | 0.161 | 0.125 | −22.5% |
| translation | **EC–PA** | 3,826 | **0.742** | **0.738** | **−0.6%** |
| translation | BS–PA | 314 | 0.263 | 0.230 | −12.7% |

**Gap ratio (BS-pairs mean / EC–PA), the decisive summary statistic:**

| Readout | Raw | GC-controlled | % change | Verdict |
|---|---|---|---|---|
| Transcription | 0.341 | 0.308 | −9.7% | SURVIVES INTACT |
| Translation | 0.286 | 0.240 | −15.9% | SURVIVES, ATTENUATED |

## Table S4 — Phylum-stratified cross-host correlation (assumption-free confirmation)

Full table: `out/results/gate10_5_gc_confound.json` (Task 4). `null` = stratum N too small for a stable estimate, reported as such rather than shown.

**Transcription:**

| Pair | Proteobacteria (n) | Firmicutes (n) |
|---|---|---|
| EC–BS | 0.090 (321) | 0.177 (2,532) |
| **EC–PA** | **0.634 (3,282)** | **0.812 (2,666)** |
| BS–PA | −0.209 (327) | 0.271 (1,214) |

**Translation:**

| Pair | Proteobacteria (n) | Firmicutes (n) |
|---|---|---|
| EC–BS | null (18) | 0.118 (739) |
| **EC–PA** | **0.627 (1,107)** | **0.713 (1,652)** |
| BS–PA | −0.360 (31) | 0.245 (235) |

## Table S5 — DRAFTS RS234 in-vitro-vs-in-vivo agreement, per species (recomputed from released source data, not a literal quote — see `out/GATE10_MEMO.md` Task 4.3)

| Species | N | Spearman ρ (in-vitro vs. in-vivo) |
|---|---|---|
| *E. coli* | 176 | 0.90 |
| *S. enterica* | 226 | 0.76 |
| *P. putida* | 226 | 0.74 |
| *V. natriegens* | 222 | 0.79 |
| *K. oxytoca* | 204 | 0.69 |
| *B. subtilis* | 174 | 0.69 |
| *C. glutamicum* | 230 | 0.80 |

*P. aeruginosa* is not one of RS234's seven in-vivo comparison species (and is separately confirmed absent from DRAFTS entirely — see Section 3.2/Limitations item 11).

## Table S6 — GC-vs-activity, both datasets (the confound Table S3/S4 control for)

| Host/species | Modality | Readout | Spearman ρ (GC%, activity) |
|---|---|---|---|
| *E. coli* | in vivo | transcription | −0.614 |
| *E. coli* | in vivo | translation | −0.494 |
| *B. subtilis* | in vivo | transcription | −0.200 |
| *B. subtilis* | in vivo | translation | −0.266 |
| *P. aeruginosa* | in vivo | transcription | −0.434 |
| *P. aeruginosa* | in vivo | translation | −0.230 |
| DRAFTS, 10 species | cell-free | transcription | −0.49 to −0.74 (all 10 hosts) |

Full DRAFTS per-species table: `out/results/gate10_gc_confound.json`.

## Complete reference list

See `out/PREPRINT/MANUSCRIPT.md`, "References" section — reproduced there in full, split into primary-source-verified and not-independently-re-verified citations, per Gate 9's original disclosure policy and unchanged in this gate except for the added Yim et al. (2019) DRAFTS citation.
