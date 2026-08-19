# GATE 10 MEMO — DRAFTS Acquisition and Feasibility

**Did the data acquire cleanly? Yes, via a fallback route (PMC's official bulk-download API, after both EMBO/Springer Link and PMC's own per-file download were blocked by bot-detection walls) — every file hash-verified.**
**Do the IDs join to Johns? No — the two libraries' OLIGO ID ranges are entirely disjoint. Sequence-text join works: 112 of 1047 DRAFTS sequences (10.7%) match a Johns three-host sequence, 1 matches an RS241 sequence.**
**Is a host-count sweep feasible at this N? Technically constructible (421–891 sequences depending on host count), but the sweep's relevance to the question that motivated this gate is now in serious doubt — DRAFTS's cell-free cross-host structure is measurably a different, more homogeneous phenomenon than Johns's in-vivo structure (Task 4.3). See Task 5 for the full reasoning; this is a DECISION NEEDED, not a clean go/no-go from this memo alone.**

**This gate acquired and characterized data only. No model was trained, no hypothesis was evaluated, and the frozen splits were not touched**, per explicit instruction.

---

## TASK 1 — Data acquisition

**Paper confirmed:** Yim, Johns, et al., "Multiplex transcriptional characterizations across diverse bacterial species using cell-free systems," *Molecular Systems Biology* 15:e8875 (2019). DOI 10.15252/msb.20198875. PMC6692573.

**Route 1 (EMBO/Springer Link): BLOCKED.** `doi.org/10.15252/msb.20198875` now redirects to Springer Link (EMBO journals are Springer-Nature-published); both the Springer Link article page and the direct `embopress.org` mirror serve a JavaScript "Client Challenge" bot-detection interstitial instead of content, to both `curl` (browser User-Agent) and WebFetch. The live EMBO supplement-download endpoint pattern resolves but serves the same challenge page instead of a binary file. A Wayback Machine snapshot recovered the 2019 article HTML (used to help identify supplementary file names) but never archived the binary supplementary files themselves (zero CDX hits).

**Route 2 (PMC): SUCCESS, via a fallback within the route.** `pmc.ncbi.nlm.nih.gov/articles/PMC6692573/` loads cleanly (200 OK, no CAPTCHA) and lists every supplementary file, but the direct per-file download links serve a proof-of-work JavaScript challenge (`cloudpmc-viewer`) instead of the file. **Workaround: NCBI's official PMC Open Access Web Service** (`www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=PMC6692573`) returned a bulk-package FTP link (`ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/deprecated/oa_package/53/6b/PMC6692573.tar.gz`), which downloaded cleanly and contains every file the PMC page lists.

**Route 3 (bioRxiv): BLOCKED.** Both the preprint page and its source XML returned HTTP 429 (rate-limited) on every attempt, including one retry after the stated backoff. The separate bioRxiv metadata API (not rate-limited) confirmed the preprint's own license is `cc_by_nc_nd` — more restrictive than the final published article's CC BY 4.0 — but no content was retrieved from bioRxiv.

**Route 4 (raw reads accession): FOUND, quoted exactly, not guessed.** Verified directly in the article XML's data-availability section: *"Sequencing data associated with this study are available at NCBI SRA under **PRJNA509603**."* No GEO accession is stated. Raw reads were not downloaded — the processed source-data tables (below) are sufficient for this gate's purposes, per the task's own framing ("processed matrices are likely sufficient").

**Files acquired** (all from the PMC OA bulk package; **license CC BY 4.0**, confirmed from the article XML's `<license license-type="creativeCommonsBy">` tag — I found no separate/distinct CC0 label specifically for the Source Data files, so CC BY 4.0, the article's own stated license, is what is recorded, not an assumed CC0):

| File | Size | SHA256 (first 16 hex) | Contents |
|---|---|---|---|
| `msb198875_SourceData_Fig3.xlsx` | 1,083,660 B | `a5f9ab026a931eee` | **Key file**: 10-species x 1047-sequence cell-free matrix ("Fig. 3" sheet) |
| `msb198875_SourceData_Fig2.xlsx` | 328,242 B | `e002ced3f158eb50` | RS234 in-vitro/in-vivo comparison, 7 species ("Appendix Fig. S5" sheet) |
| `msb198875_Appendix.pdf` | 4,520,346 B | `8bb7e233846b04e0` | Methods, Tables S1-S6, Figure legends |
| `msb198875_FullTextArticle.{pdf,nxml}` | 2.9MB / 178KB | — | Full text (license + data-availability verification) |
| `msb198875_SourceData_Fig{1,4,5}.xlsx`, `_Appendix.zip`, `_ReviewProcessFile.pdf`, `PMC6692573.tar.gz` | — | — | Acquired, not used this gate (RS29249/RS7003/hybrid-lysate data — out of scope for Tasks 2-5) |

Full URL/size/SHA256/license record for every file: `data/MANIFEST.json` (26 total entries, up from 12 before this gate — 3 new `data/*.parquet` outputs plus 10 `raw/drafts/*` source files added). All hashes independently re-verified against the manifest as part of writing it (not just trusted from the acquisition pass).

---

## TASK 2 — Inventory

Script: `scripts/88_parse_drafts.py`. Outputs: `data/drafts.parquet` (1047 rows x 90 cols), `data/drafts_rs234_invivo.parquet` (234 rows x 72 cols).

**Ten species, CONFIRMED — with one load-bearing correction to a natural assumption.** The abstract's own species list (verified against Appendix Table S1 and the Results text, both independently, exact string matches): seven Proteobacteria (*E. coli, E. fergusonii, S. enterica, **P. agglomerans**, K. oxytoca, P. putida, V. natriegens*), two Firmicutes (*B. subtilis, L. lactis*), one Actinobacterium (*C. glutamicum*). **"Pa" in every DRAFTS table is Pantoea agglomerans, NOT Pseudomonas aeruginosa** — confirmed three independent ways: Appendix Table S1's species/growth-condition table lists "P. agglomerans" at the exact column position "Pa" occupies; the Results section names all ten species by full binomial including "P. agglomerans"; a full-text-and-appendix grep for "aeruginosa" returns **zero** matches anywhere in either document. **This project's *P. aeruginosa* — one of its three primary hosts — is entirely absent from DRAFTS.** Everything downstream in this memo that discusses "shared hosts" reflects this correction; it changes the shape of Tasks 3 and 4 substantially from what a naive reading of the abstract would suggest.

**The 421-sequence matrix: exists, and 421 is independently recomputed, not taken on faith.** Every one of the 10 `{species}_zscore_union` columns has exactly 421 non-null values, and the 10-way intersection of "non-null in all 10" is also exactly 421 — the paper's own headline figure, reproduced exactly by two different counting methods.

**A real data-quality trap, caught before it could propagate:** the `{sp}_tx`/`{sp}_zscore` columns are NOT null for unusable rows — they contain the literal strings `'no_DNA_counts'`, `'no_RNA_counts'`, or `'low_DNA_counts'` in place of a number. A naive `.notna()` check would have reported 1047/1047 "usable" for every species (100%, obviously wrong) rather than the true per-species usable counts below. Fixed via explicit `pd.to_numeric(errors='coerce')` before any usability logic.

**Per-species usable N** (numeric value present, not a sentinel string), out of 1047:

| Species | Usable N | Unusable reasons (count) |
|---|---|---|
| Ec (E. coli) | 897 | no_DNA:83, no_RNA:43, low_DNA:24 |
| Ef (E. fergusonii) | 828 | no_DNA:156, low_DNA:54, no_RNA:9 |
| Se (S. enterica) | 475 | no_RNA:403, no_DNA:142, low_DNA:27 |
| Ko (K. oxytoca) | 752 | no_DNA:236, no_RNA:32, low_DNA:27 |
| Pa (P. agglomerans) | 876 | no_DNA:125, low_DNA:43, no_RNA:3 |
| Pp (P. putida) | 802 | no_DNA:192, low_DNA:41, no_RNA:12 |
| Vn (V. natriegens) | 822 | no_DNA:178, low_DNA:40, no_RNA:7 |
| Bs (B. subtilis) | 830 | no_DNA:141, low_DNA:44, no_RNA:32 |
| Cg (C. glutamicum) | 908 | no_DNA:102, low_DNA:37 |
| Ll (L. lactis) | 932 | no_DNA:95, low_DNA:20 |

*S. enterica* is the clear outlier (475/1047, dominated by `no_RNA_counts`) — worth flagging for Gate 11 if *S. enterica* is included in any subset.

**The RS234 in-vivo arm: located, N confirmed.** `msb198875_SourceData_Fig2.xlsx`, sheet "Appendix Fig. S5" — 234 sequences x 7 species (Ec, Se, Pp, Vn, Ko, Bs, Cg — **not** Pa/Ef/Ll, matching the abstract's "seven species for in-vivo comparison"), with per-species in-vitro AND in-vivo columns (2 replicates each, DNA/RNA/tx/zscore). Per-species N usable in both modalities: Ec 176, Se 226, Pp 226, Vn 222, Ko 204, Bs 174, Cg 230 (all /234).

**Units/normalization:** `_tx` = a linear-scale ratio (right-skewed, e.g. *E. coli* median 0.203, max 101.5 — consistent with an RNA/DNA-count-derived ratio, not yet log-transformed). `_zscore` = z-score of log10(tx), per-species, matching the paper's stated methodology ("transcription levels in log10 scale were transformed to Z-score"). `_zscore_union` = the same, restricted to the 421-sequence 10-way-usable set.

**Replicate structure: recoverable, and directly useful — a second, independent reliability estimate.** `msb198875_Appendix.pdf`, Appendix Table S6 ("Reproducibility between biological replicates"), Figure 1D row, states directly (**transcribed, not computed**): Pearson r between biological replicate 1 and 2 = **0.8932** for *E. coli* in vivo (RS29249) and **0.9333** for DRAFTS *E. coli* cell-free (RS29249). This is a genuinely new, independent *E. coli* reliability estimate — via true biological replicates rather than Gate 8's five-growth-condition proxy (0.912/0.929) — landing in the same low-to-mid-0.9s range. Table S6 also gives per-species replicate-r for the RS1383/RS7003 libraries (all 10 species, all in the high 0.8s-0.98 range) — a further, broad confirmation that DRAFTS's own measurements are internally reproducible, independent of the cross-host or cross-modality questions below.

---

## TASK 3 — Join to Johns tables [LOAD-BEARING]

Script: `scripts/89_join_drafts_johns.py`. Full report: `out/results/gate10_join_report.json`.

**ID join: does not work.** DRAFTS's Oligo ID range is 13097-14477; this project's `three_host.parquet` range is 14478-43726; RS241's range is 14772-48469. **Zero overlap with either**, and the ranges are contiguous (DRAFTS ends exactly one below where three_host begins) — consistent with the same overall synthesis/ID-numbering pipeline having allocated non-overlapping ID batches to the two studies (both from the Wang lab, both built on the same 184-genome mining effort), not two unrelated numbering schemes.

**Sequence-text join: works, and is the join key.** 112 of DRAFTS's 1047 sequences (10.7%) exact-match a `three_host.parquet` sequence; 1 of 1047 matches an RS241 sequence (via a library-ID lookup, since RS241's own table has no sequence-text column — 207/241 RS241 sequences are recoverable this way, matching Gate 3.5's prior finding exactly).

**Fold compatibility: PARTIALLY COVERED, 111/112 (99.1%), reasonably distributed.** Of the 112 matched sequences, 111 already sit inside the frozen 5-fold assignment (fold distribution: 24/23/6/31/27 — fold 2 is thin at 6, worth noting for anyone reusing this specific 112-row overlap directly). The 1 uncovered row is confirmed to be exactly the 1 sequence that is also an RS241 member (excluded from the main fold assignment for the same reason every other RS241 ID was). **Caught and fixed a real bug while computing this:** `fold_assignment_FINAL.parquet`'s `OLIGO ID` column is stored as string dtype, not int; a first pass using an int-keyed reindex silently returned 0/112 matches (a dtype mismatch, not a real non-coverage finding) before being caught and fixed by matching dtypes explicitly.

**This is a small overlap (112 of 1047) and a small fraction of the frozen splits (111 of 29,042 rows, 0.4%)** — worth stating plainly: the 935 DRAFTS sequences with no Johns match at all have no existing fold assignment and would need entirely new genome-blocked splits (with the same near-duplicate safety-net treatment Gate 2 built) for any DRAFTS-native modeling. **This is exactly the "new split work needed" case the task asked to distinguish** — not fully covered, not fully uncovered, but a small, usable partial overlap plus a much larger unsplit remainder. **No frozen split was modified this gate.**

**Primary-host overlap, corrected for the P. agglomerans/P. aeruginosa finding:** of this project's three primary hosts (EC, BS, PA), only **EC and BS** appear in DRAFTS. Of the 112 matched rows: 94 usable in Ec, 83 usable in Bs, 82 usable in both.

**Source-genome metadata: carried natively by DRAFTS**, does not need to be joined from Johns's tables — `source_species`/`source_genus`/`source_family`/`source_order`/`source_class`/`source_phylum`/`source_genome` are present directly in the Fig. 3 sheet for every one of the 1047 sequences.

---

## TASK 4 — Characterization on DRAFTS's own terms

### 4.1/4.2 — Cross-host correlation, all 45 pairs, and the phylogenetic pattern

Script: `scripts/90_drafts_crosshost_correlations.py`. Full table: `out/results/gate10_crosshost_correlations.csv`. Figure: `out/figures/gate10_crosshost_correlation_matrix.png`.

**Every one of the 45 pairwise cell-free correlations is high**: range **[0.623, 0.911]**, mean 0.81. Same-phylum pairs average 0.852; cross-phylum pairs average 0.769 — **the phylogenetic pattern holds in direction** (same-phylum higher), but the effect size is small compared to Johns's in-vivo pattern. The weakest single pair (*E. coli*-*L. lactis*, 0.623) is still more than **double** the strongest *B. subtilis* in-vivo pair from Johns (≈0.26).

**Does *B. subtilis* behave as an outlier in cell-free the way it does in vivo? No.** Per-species mean correlation with the other 9 species: *L. lactis* is actually the lowest (0.729), not *B. subtilis* (0.783, roughly mid-pack, close to *E. coli*'s 0.758). *B. subtilis*, *L. lactis*, and *C. glutamicum* (the two Firmicutes and the one Actinobacterium — the non-Proteobacteria) do NOT collectively stand out as low-correlating outliers; *C. glutamicum* actually has one of the higher mean correlations (0.820), above four of the seven Proteobacteria.

### 4.3 — The modality comparison [LOAD-BEARING]

Script: `scripts/91_drafts_modality_comparison.py`. Full output: `out/results/gate10_modality_comparison.json`. Figure: `out/figures/gate10_modality_comparison.png`.

**P. aeruginosa is confirmed absent from DRAFTS** (Task 2) — the only primary-host pair testable both in-vivo and cell-free is **EC-BS**.

**Comparison A (this project's own computation, EC-BS pair, same underlying biology, two modalities):**

| | n | Spearman rho |
|---|---|---|
| Cell-free (DRAFTS), EC-BS | 82 | **0.597** |
| In-vivo (Johns), EC-BS, on the same 112-seq overlap | 15 | 0.147 *(n too small to trust)* |
| In-vivo (Johns), EC-BS, full library (reference) | 3,668 | **0.258** |

**The cell-free EC-BS correlation (0.597) is more than double the in-vivo EC-BS correlation (0.258, the trustworthy reference figure).** This is the single most important number in this gate.

**Comparison B (DRAFTS's own paired RS234 in-vitro-vs-in-vivo design, SAME species, different modalities):** recomputed directly from DRAFTS's own released source data (`data/drafts_rs234_invivo.parquet`) — **not** a value transcribed from the paper's text, since no clean extractable per-species table exists in the Appendix beyond the replicate-reproducibility table (Table S6, which is a different quantity — see Task 2). Per species, in-vitro vs in-vivo Spearman rho: Ec 0.90, Se 0.76, Pp 0.74, Vn 0.79, Ko 0.69, **Bs 0.69**, Cg 0.80 — all in a fairly narrow, fairly high [0.69, 0.90] band.

**WHAT I COULD NOT DO:** the exact per-species in-vitro-vs-in-vivo r values as the paper itself would report them (likely annotated directly on the Appendix Figure S5 scatter-plot images) were not extractable via `pdftotext`, which cannot read text baked into a plot image. The recomputed values above use exactly the data DRAFTS released for this comparison and should track closely, but are not a literal quote — labeled as such throughout.

**Putting A and B together, the answer to "is the cell-free host effect the same object as the in-vivo host effect" is: NO, not for the cross-host comparison that matters to this project.** Within a single species, DRAFTS's cell-free measurement agrees reasonably well with an in-vivo measurement of that SAME species (Comparison B: 0.69-0.90) — DRAFTS is a reasonably faithful proxy for a given host's own regulatory behavior. But the CROSS-HOST DIFFERENCES — the actual chassis effect, which is what this whole project studies — are not preserved: cell-free shows every pair correlating strongly and *B. subtilis* behaving unremarkably (4.1/4.2), while in-vivo shows a dramatic split with *B. subtilis* as a severe outlier. **A conditioning result trained and tested on DRAFTS's cell-free data would be answering a real but different question — does host-count matter for cross-species cell-free TXTL transcription — not the in-vivo chassis-effect question this project's central claim is about.** This is stated as an observation from directly computed correlations, not a model result, per the gate's own scope limit.

### 4.4 — The GC/phylum confound

Script: `scripts/92_drafts_gc_and_floor_check.py`. Figure: `out/figures/gate10_gc_vs_activity.png`.

**Confirmed, strongly, in every one of the 10 recipient hosts.** Spearman(source-genome GC%, log10 activity) ranges from **-0.49 to -0.74** across all 10 species — lower-GC source sequences are more active, in every recipient, with no exception. Separately: Firmicutes-derived sequences show higher mean log-activity than Proteobacteria-derived sequences in **every one of the 10 recipient hosts** — directly reproducing the paper's own stated finding ("sequences derived from Firmicutes were more active overall than those from Proteobacteria, likely reflecting genomic GC content") in this project's own re-derived numbers, not just quoted from the text. **This is a real, universal confound for any future host-count sweep or cross-host model on DRAFTS**: source-sequence GC/phylum composition must be controlled for (e.g., stratified sampling across source phyla, or including GC% as a covariate) or apparent host effects may partly be source-composition artifacts, exactly as the task warned.

### 4.5 — Floor-value artifact check

Script: `scripts/92_drafts_gc_and_floor_check.py`. No Johns-style floor artifact found: the most common single rounded `tx` value in any species accounts for at most **0.60%** of that species' usable rows (*B. subtilis*, the highest) — nowhere close to Johns's 89.9% *B. subtilis*-translation pile-up. **DRAFTS's transcription-only cell-free readout does not show an analogous detection-floor convention.** Unusable fractions (sentinel-string rows) range 11.0%-54.6% per species (*S. enterica* highest, consistent with its low usable-N flagged in Task 2) but do not concentrate at a single pinned value the way Johns's translation floor did.

---

## TASK 5 — Feasibility verdict for a host-count sweep

Script: `scripts/93_drafts_feasibility.py`. Full output: `out/results/gate10_feasibility.json`.

**1. Maximum training hosts: 10** (all of DRAFTS's cell-free panel). N usable in ALL k hosts simultaneously, by k (best achievable subset at each k, exhaustively checked through k=6, 300-subset sample for k=7-9):

| k | N range across subsets | Best subset | Best N |
|---|---|---|---|
| 2 | 430-891 | Cg, Ll | 891 |
| 3 | 424-862 | Pa, Cg, Ll | 862 |
| 4 | 422-828 | Ec, Pa, Cg, Ll | 828 |
| 5 | 421-798 | Ef, Pa, Vn, Cg, Ll | 798 |
| 6 | 421-778 | Ef, Pa, Pp, Vn, Cg, Ll | 778 |
| 7 | 421-762 | Ec, Ef, Pa, Vn, Bs, Cg, Ll | 762 |
| 8 | 421-750 | Ec, Ef, Pa, Pp, Vn, Bs, Cg, Ll | 750 |
| 9 | 421-717 | (all except Se) | 717 |
| 10 | 421 (fixed) | all | 421 |

**2. A volume-controlled sweep is constructible**, but the total N it can hold constant across the FULL k=2..10 range is modest: **1,782 total rows**, limited by k=2's ceiling (891/host x 2 = 1,782 is the smallest of the per-k maximum totals — every larger k could in principle support a bigger total, but the design needs ONE total-N that works everywhere). Concretely: k=2 → 891/host, k=5 → 356/host, k=10 → 178/host (subsampled from the 421 available), all summing to ≈1,780.

**3. Is N enough for fold-resolved 90% bootstrap intervals? Marginal, not clearly too small — but this project's own established variance behavior (SR5: fold-to-fold variance often exceeds draw-to-draw, growing with N) means intervals at the higher host counts (k=7-10, per-host-N 178-254) should be expected to be WIDE.** This is comparable in scale to configurations this project has already run and found noisy-but-informative (e.g. Gate 3.5's N≈300 direct-regression arms) — not the clear "too small to distinguish anything" verdict that would save a gate outright, but not a generously-powered design either. **Compounding factor, not accounted for in the N figures above: DRAFTS has no existing fold structure of its own** (Task 3) — 935 of 1047 sequences are entirely unsplit, and Gate 11 would need to build new genome-blocked splits (with the same near-duplicate safety net Gate 2 built) before any of the N figures above could be used in a real 5-fold rotation. That is real, non-trivial infrastructure work, not a data-availability problem, but it belongs in the feasibility accounting.

**4. Host-subset sampling:** distinct k-host subsets available: k=2: 45, k=3: 120, k=4: 210, k=5: 252, k=6: 210, k=7: 120, k=8: 45, k=9: 10, k=10: 1. Ample subset diversity exists at every k from 2 through 8 to separate "how many hosts" from "which hosts" (e.g., at k=5, 252 distinct 5-host combinations are available to sample from) — this part of the design is not constrained by DRAFTS's size.

**5. Recommendation: NOT a clean go, NOT a clean no-go — DECISION NEEDED.** Two independent concerns, of different character:
   - **Sample size** is marginal but workable with disclosed caveats (wide intervals likely at high k; new splits required).
   - **Relevance is the dominant concern, and it is more serious than sample size.** Task 4.3 found DRAFTS's cross-host structure is measurably NOT the same phenomenon as the in-vivo chassis effect this project's central claim is about — *B. subtilis* is not an outlier in cell-free, all cross-host correlations are compressed into a much higher, narrower band, and the EC-BS pair specifically transfers more than twice as well cell-free as in-vivo. **A host-count sweep on DRAFTS would be well-powered to answer "does training-host count matter for cell-free cross-species transcription transfer" — a real, publishable question in its own right — but would NOT, by itself, resolve whether more training hosts would have fixed FiLM's underdetermination on the ORIGINAL in-vivo problem**, which is the objection that motivated this gate in the first place. **If Gate 11 proceeds, its pre-registration should state this scope limit explicitly up front, not discover it after results are in** — exactly the failure mode Gate 8.5/8.6 already had to correct for once (the shift-prediction retraction) and should not repeat by implication here.

---

## A self-caused audit regression, found and fixed before delivery

Adding the `raw/drafts/*` source files to `data/MANIFEST.json` (as instructed) broke `audit_leakage.py`'s check 6, which had always resolved every manifest key relative to `data/` — it reported all 10 new raw-file entries as "missing" even though they were present on disk at their correct `raw/drafts/` path. **This was a real regression, introduced by this gate's own work, not a pre-existing bug.** Fixed in `scripts/audit_leakage.py` by resolving manifest keys prefixed `raw/` or `out/` relative to the project root instead of `data/`, backward-compatible with every existing bare-filename entry. Rerun after the fix: 26/26 files checked, 0 mismatches, all 6 checks pass. Disclosed here rather than silently fixed and left unmentioned — this is exactly the kind of thing the standing audit practice exists to catch.

## WHAT I COULD NOT DO

1. Routes 1 (EMBO/Springer) and 3 (bioRxiv) for data acquisition — both blocked (bot-detection wall; rate-limiting), worked around via Route 2's PMC OA bulk-package fallback, which delivered every needed file.
2. Raw sequencing reads (SRA PRJNA509603) — not downloaded; the processed source-data tables were sufficient and the task itself flagged raw reads as optional.
3. The exact per-species RS234 in-vitro-vs-in-vivo Pearson r values as the paper's own Figure S5 panels would report them — not extractable from the PDF (baked into plot images); recomputed instead from DRAFTS's own released source data for the same comparison, clearly labeled as recomputed, not quoted.
4. A rigorous phylogenetic tree ordering for the correlation-matrix figure — used a coarse manually-assigned order (Proteobacteria, then Firmicutes, then Actinobacteria) rather than an actual 16S/rpoD-distance-based ordering (DRAFTS's own Fig. 3D sheet has pairwise 16S-identity/rpoD-distance data that could build a proper tree — not used here, flagged as available for Gate 11 if needed).
5. New genome-blocked splits for the 935 DRAFTS-only sequences with no existing fold assignment — explicitly out of scope this gate ("do not modify the frozen splits"); would be required infrastructure work before any DRAFTS-native LOHO evaluation.

## TEMPTATIONS TO ADJUST

1. **Reporting "10 shared species" or "3 shared primary hosts" without catching the P. agglomerans/P. aeruginosa distinction** — the natural, easy misreading given the abbreviation "Pa" and the charter's own established use of "PA" for P. aeruginosa throughout this project. Caught via direct primary-source verification (three independent checks) before it could propagate into Task 3/4/5's framing.
2. **Presenting the modality-comparison finding (Task 4.3) more tentatively than the numbers support**, to avoid appearing to pre-judge Gate 11's design. Not done — the numbers are what they are (0.597 vs 0.258, more than 2x) and are reported plainly as the gate's central finding, with the appropriately bounded claim ("would not by itself resolve..." rather than "DRAFTS is useless").
3. **Treating Task 5's feasibility question as purely a sample-size question**, since that is the easier, more mechanical thing to answer. Not done — the relevance concern from Task 4.3 is reported as the dominant issue, ahead of (not instead of) the sample-size analysis.
4. **Silently substituting a text-quoted number for the RS234 in-vitro-vs-in-vivo comparison** when a clean one could not be extracted from the PDF, to avoid the extra step of computing it from the source data directly. Not done — computed it properly instead, which turned out to be strictly better (exact N, verifiable) than a transcribed figure-legend number would have been.
