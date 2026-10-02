**Did Tasks 1 and 2 apply, and what does the DRAFTS source say about the clustering? Both applied. On Task 3 the external read is correct and the manuscript was wrong: Yim et al. describe the clustering as "distinct gram-negative and gram-positive groups", not phylum. Their term is now adopted. A second error surfaced in the same check — our sentence attributed both their clustering and their phylogenetic gradient to a bare "Fig. 3", when the clustering is Fig. 3C and the gradient is Fig. 3D. Both panels are now cited correctly. No number changed.**

# GATE 29 MEMO

Script: `scripts/113_gate29_corrections.py`. Both audits clean at start and end. Abstract stays at 250 words; body 6,650 → 6,717.

---

## TASK 1 — Attenuation claim scoped [APPLIED]

Arithmetic recomputed from `gate25_modality_gc_control.json` rather than taken on report:

| regime | cell-free | in vivo | apart by |
|---|---|---|---|
| co-active | 0.677 → 0.556 = **17.9%** | 0.258 → 0.207 = **19.8%** | 1.11× |
| pooled | 0.616 → 0.561 = **8.9%** | 0.655 → 0.517 = **21.1%** | **2.36×** |

Confirmed: the generalisation is false for the pooled half.

**Before:**
> Both modalities attenuate by nearly the same proportion, so the contrast is not an artifact of source-genome GC composition.

**After:**
> Under co-active restriction the two modalities attenuate by nearly the same proportion, so the co-active contrast — the one this paper's claim rests on — is not an artifact of source-genome GC composition. The pooled attenuations differ considerably, but the pooled comparison shows no modality difference before or after the control, so nothing in the argument turns on them.

The pooled figures and the 0.94× → 1.09× ratio are untouched, as instructed. No percentage is restated, so no number is duplicated into the paragraph.

---

## TASK 2 — Three checks, two confounds [APPLIED]

**Before:**
> What can be said is that the two confounds a careful reader proposes first have each been tested and neither survives contact with the data.

**After:**
> What can be said is that the two confounds a careful reader proposes first have each been tested — **GC twice, on the in-vivo pairs and on the modality contrast** — and neither survives contact with the data.

The "surviving three does not rule out a fourth" caveat in the preceding sentence is untouched. The rejected "tested three ways" phrasing was not used, for the reason given: the three checks are not three methods, they are two methods with one applied to two different comparisons.

---

## TASK 3 — Clustering terminology [SOURCE CHECKED; the read was right]

### What DRAFTS actually says

Main text, verbatim from `raw/drafts/msb198875_FullTextArticle.pdf`:

> "Pearson correlations of these pairwise comparisons showed varying levels of transcriptional concordance, which when clustered further revealed **distinct gram-negative and gram-positive groups (Fig 3C)**. Principal component analysis of the RS1383 transcription profiles also revealed **distinct groups for gram-negative and gram-positive species** (Appendix Fig S7C)."

And, separately:

> "Pairwise comparisons of 16S rRNA and sigma70 protein sequence similarity showed that **more phylogenetically related species tend to share a more similar transcription profile (Fig 3D**, Appendix Fig S8)."

**They say gram-stain for the clustering.** Their term is adopted.

### A second error the same check turned up

Our sentence attributed *both* results to a bare "Fig. 3". The caption separates them:

- **Fig. 3C** — "Pairwise comparison of transcriptional profiles… between bacterial species", the heat map that clusters.
- **Fig. 3D** — "Correlation between evolutionary divergences (16S rRNA percent identity) and pairwise Pearson correlation of transcriptional profiles. Dashed line represents linear regression." This is where r = 0.73 lives.
- **Fig. 3E** — "Activity profiles (Tx) of regulatory sequences from **donor phyla** Proteobacteria and Firmicutes", which is about the phyla the *sequences were mined from*, not how the recipient species group.

Gate 28's memo put r = 0.73 in panel E. That was wrong, and panel E is where the word "phyla" appears in their caption — which is plausibly where the manuscript's "clusters by phylum" came from in the first place. The number and its direction were right; the panel was not. Panel D is the only panel with a 16S axis and a linear regression, so r = 0.73 is unambiguously D's.

### Before and after

**Before:**
> Yim et al.'s own Fig. 3 reports that cell-free cross-species correlation **clusters by phylum** and tracks phylogeny, at r = 0.73 against 16S rRNA similarity across the species pairs.

**After:**
> Yim et al.'s own Fig. 3 reports that cell-free cross-species correlation **clusters into distinct gram-negative and gram-positive groups (their Fig. 3C)**, and that more phylogenetically related species share more similar transcription profiles, at r = 0.73 against 16S rRNA similarity **(their Fig. 3D)**.

**The substance of the paragraph is unaffected.** Our argument is that their result describes how agreement *varies* with phylogenetic similarity while this paper's claim is about its *level* under restriction. That holds whether the clustering is labelled by gram stain or by phylum — the correction is terminological and attributional, exactly as anticipated.

r = 0.73 and the "similarity" direction were verified in Gate 28 and were not rechecked, per instruction.

---

## Numbers

**None changed.** The guard compares numeric-token counts in both directions and refuses to write on any difference. It fired once, on `3` going 9 → 11: that is `Fig. 3C` and `Fig. 3D` replacing one bare `Fig. 3`. Panel letters, not quantities — allowed through a documented exception rather than by weakening the check. Every measured value in the manuscript is bit-for-bit as it was.

---

## TASK 4 — Verification

| check | result |
|---|---|
| resource warnings | **0** |
| images | **9** |
| captions | Figures 1, 2a, 2b, 3, 4, 5, 6, 7, 8 — ascending |
| glyphs | ρ 17, ≤ 3, γ 2, β 2, ≈ 1 |
| pages | **22** |
| regime table | renders |
| "BS 15,848" | intact |
| `[email]` | 0 |
| cross-references | all resolve |
| Limitations | 12 items |
| both copies | byte-identical |
| abstract | **250 words** (unchanged) |
| body | 6,650 → **6,717** |
| audits | leakage ALL 6 PASSED; provenance 0 orphans (177 files) |
| `verify-citations` | 11 VERIFIED / 0 unexplained / 4 explained / **1 skipped** — `dnabert2: network/lookup error: The read operation timed out`. Transient, same class as the Gate 26–27 HTTP 429s; VERIFIED and unexplained-mismatch counts are unchanged. |

---

## Verdict

Two corrections applied as specified. The source check confirmed the external read and turned up a panel misattribution of my own from Gate 28; both are fixed. No number changed.
