**Did all four apply, and what is the abstract's final word count? All four applied. The abstract is exactly 250 words, down from 252 — the 12-word GC clause in, the two proposed cuts out. One number was added to the manuscript: r = 0.73, verified in the DRAFTS source before use. Nothing was removed.**

# GATE 28 MEMO — Four approved edits

Script: `scripts/112_gate28_approved_edits.py`. Both audits clean at start and end. Body 6,215 → 6,650 words; 20 → 22 pages.

---

## TASK 1 — Abstract [APPLIED, 252 → 250]

**Added (12 words):** "…a 2.6-fold gap **that survives partial-correlation control for source-genome GC composition (2.63× raw, 2.69× controlled)**."

**Cut (14 words):** "a barrier practitioners describe as unresolved for engineering non-model hosts" (10) and "over all measured sequences" (4).

Both ratios already appear in §3.3, so this adds no unverifiable figure to the abstract. The three evidence lines, the silent-class sentence and the single-species-pair limit are untouched, as recommended.

### Final abstract

> A regulatory DNA sequence characterized in one bacterial species often behaves differently in another — the chassis effect. Comparing cell-free (DRAFTS, Yim et al. 2019) against in-vivo (Johns et al. 2018) measurements of the same 165 bp library, we find that cross-host agreement responds to conditioning on co-activity very differently in the two modalities. Pooled, *E. coli*–*B. subtilis* transcription correlates at **0.616 cell-free and 0.655 in vivo** — no modality difference. Restricted to sequences detectably active in both hosts, the in-vivo correlation falls to **0.258** while the cell-free correlation rises to **0.677**, a 2.6-fold gap that survives partial-correlation control for source-genome GC composition (2.63× raw, 2.69× controlled). In living cells, cross-host agreement is carried substantially by agreement on which sequences are silent; lysates silence almost nothing. **This comparison rests on a single species pair**, the only one present in both datasets. We interpret it as evidence, not proof, that host specificity resides substantially in cellular context rather than in the transcriptional machinery lysates retain. That interpretation predicts host-descriptor features encoding machinery composition should carry no signal, and we tested that prediction directly: a pre-registered hypothesis failed on all 8 primary comparisons; an arbitrary host-identity tag matched or beat both a 37-feature genomic and a 6-feature physiology vector in all 12 cells; and sequence-only prediction was not distinguishably beaten in 16 of 18 tests across three conditioning mechanisms. **n_hosts ≤ 6** throughout. We release the benchmark, splits, evaluation code, and baseline suite, including the negative result and a retracted rescue attempt.

---

## TASK 2 — Yim et al. Fig. 3 as an objection [APPLIED to §3.3]

One paragraph, placed at the end of §3.3. It states the objection at full strength before answering it, distinguishes their *gradient* from this paper's *level*, quotes the 45-pair band against the in-vivo 0.258, and names the recomputation as the open check — including that the gradient could strengthen rather than weaken under co-active restriction. The recomputation was **not** run.

**One precision against the source.** The brief describes r = 0.73 as a gradient "against 16S distance". DRAFTS Fig. 3E plots it against **16S rRNA similarity (%)**, so the correlation is positive — closer species agree more. The paragraph says similarity. Verified by extracting the figure panel text from `raw/drafts/msb198875_FullTextArticle.pdf`; the same panel also carries r = 0.91 and r = 0.61 for its two sub-panels, and 0.73 is the one for the whole set.

**r = 0.73 is the only number this gate added to the manuscript.** It is quoted from DRAFTS, not computed here, and the text attributes it to them.

---

## TASK 3 — Fitness decoupling [APPLIED to §3.4 and Limitations item 2]

One paragraph at the end of §3.4, directly after the "hypothesis consistent with the data, not a demonstrated causal claim" paragraph, so the hedge lands on a named competitor rather than on nothing.

**A precision that matters, and I have flagged it in the text rather than smoothing it over.** Yim et al.'s sentence is:

> "the distribution of **DNA libraries** was much more uniform in vitro than in vivo, mostly due to the decoupling of gene expression from cell fitness in vitro, which improves accuracy by minimizing noisy measurements from low-abundance library members"

Their claim is about **library representation at the DNA level**, not about silent sequences at the RNA level. The extension to activity is the natural one and it is the one that makes the account a competitor — but it is **our** extension, not their claim, and the paragraph says so explicitly. Writing it as though they had made the stronger claim would have misattributed it.

With that stated, the paragraph does what was asked: it gives their account plainly, says both accounts predict the same thing for every quantity available here, and says that separating them needs per-sequence fitness or burden data neither dataset provides.

**Limitations item 2** gained a matching clause, so the limitation is findable from the Limitations section and not only from §3.4. Item count unchanged at 12.

---

## TASK 4 — Provenance disclosure [APPLIED verbatim]

Added to Data and Code Availability, after the reproducibility paragraph where the audits are described, exactly as specified — "during development" rather than gate numbers.

---

## Numbers

| | |
|---|---|
| added | **0.73** — Yim et al.'s reported gradient, quoted from DRAFTS Fig. 3E and verified in the source PDF before use |
| removed | **none** |
| 2.63× / 2.69× in the abstract | already present in §3.3; not new to the manuscript |
| 45-pair band 0.623–0.911, in-vivo 0.258 in the new §3.3 paragraph | already present in §3.3; restated, not introduced |

The script refuses to write if any numeric token disappears.

---

## TASK 5 — Verification

| check | result |
|---|---|
| resource warnings | **0** |
| images | **9** |
| captions | Figures 1, 2a, 2b, 3, 4, 5, 6, 7, 8 — present, ascending |
| glyphs | ρ 17, ≤ 3, γ 2, β 2, ≈ 1 — all present |
| pages | **22** (from 20; the two new paragraphs) |
| regime table | renders |
| "BS 15,848" | intact on one line |
| `[email]` | 0 |
| cross-references | all resolve |
| Limitations | 12 items |
| both copies | byte-identical |
| abstract | **250 words** |
| body | 6,215 → **6,650** |
| audits | leakage ALL 6 PASSED; provenance 0 orphans (177 files) |
| `verify-citations` | **11 / 0 / 5 / 0 / 16 — back to baseline.** The two HTTP 429 skips of Gates 26–27 have cleared on their own, which confirms they were rate-limiting and not a content problem. |

---

## Verdict

Four edits applied, abstract at exactly 250, one number added and verified at source, none removed. Two places where the brief's paraphrase of DRAFTS was looser than DRAFTS — 16S similarity rather than distance, and library representation rather than silent sequences — are corrected in the text rather than reproduced.
