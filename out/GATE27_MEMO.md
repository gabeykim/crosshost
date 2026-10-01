**Are the figures current or stale, and which account was wrong? The figures are current. All nine match their embedded captions, verified by extracting the images from `manuscript.pdf` itself. The review was wrong on all three counts — Figure 2b shows five strata with N on every bar, Figure 5 shows nine systems on a shared y-axis, and Figure 3 shows the four-bar panel with no 0.597 and panel B present. Gate 26's verification was correct. The reviewer read a stale PDF.**

# GATE 27 MEMO

Scripts: `110` (text corrections), `111` (figure-state table). Both audits clean at start and end.

---

## TASK 1 — Figure ground truth

### Why the two accounts could both look plausible

The working tree was **clean** at gate open, so the figures on disk are byte-for-byte what commit `3cfbca6` holds. Timestamps:

- all nine figures written **2026-09-30 21:07**
- `manuscript.pdf` built **21:08:34** — after every figure
- commit `3cfbca6` at **21:11:46** — after both

So the committed PDF was built from the committed figures. There is no window in which a stale figure could have reached it.

### What the PDF actually contains

I did not re-read the source PNGs — that would not settle the question, since the dispute is about what the *PDF* shows. I extracted the embedded images with `pdfimages -png manuscript.pdf` and inspected them.

| figure | caption promises | PDF shows | |
|---|---|---|---|
| **2b** | five strata, per-stratum N on every bar | **five strata** (Proteobacteria, Firmicutes, Actinobacteria, Bacteroidetes, Cyanobacteria), **N on all 15 bars** | **MATCH** |
| **5** | nine systems incl. three conditioning mechanisms, one shared y-scale | **nine systems**, six panels, **y-axis shared** (identical −0.2 to 0.7 range, ticks on leftmost only); EC transcription shows concatenation **0.555** and per-host-heads-avg **0.615** | **MATCH** |
| **3** | panel A four bars (0.616/862, 0.655/14,088, 0.677/807, 0.258/3,668); panel B RS234, seven species | **four bars with exactly those values and N**; **no 0.597**; **no n=15 bar**; **panel B present**, seven species with N | **MATCH** |
| **8** | 11 bar groups, tx/tl distinguished, title stating both directions | title reads "Reference-value alone matches or beats the original model in 10 of 11; 1 of 11 survives"; **11 bar groups**; labels **tx/tl distinct** | **MATCH** |

The remaining five (1, 2a, 4, 6, 7) are **byte-identical** to the Gate 25 commit — `md5` equal across `ef23f32`, `3cfbca6` and disk — so regeneration is deterministic and they need no re-inspection; each was read directly in Gate 20 or Gate 25.

Full table in `out/results/gate27_figure_state.csv`. **Nine of nine MATCH.**

### Which account was wrong

**The review.** Every specific claim it makes is contradicted by the PDF in the repository:

- "Figure 2b does not match its caption" — it does; five strata, N on every bar.
- "Figure 5 does not match its caption" — it does; nine systems, shared y-axis.
- "Figure 3 is unchanged" — it was changed in Gate 25; 0.597 is gone and panel A has four bars.

The most likely explanation is a PDF downloaded before Gate 25 (`ef23f32`, 2026-09-30) or built from an older checkout. **Gate 26's verification stands.** Nothing was regenerated this gate, because nothing was stale.

One thing worth saying plainly: Gate 25 *did* contain exactly the failure this gate was opened to catch — a regeneration script that reported success while five of six producers never ran. That is presumably why the review was credible. But it was caught and fixed inside Gate 25, and the fix held.

---

## TASK 2 — Three text corrections [APPLIED]

**2a.** The pooled ratio is now as checkable as everything around it.

> before: "Pooled, it moves from 0.94× to 1.09×, still no modality difference."
> after: "Pooled, the cell-free correlation falls from 0.616 to 0.561 (8.9%) and the in-vivo from 0.655 to 0.517 (21.1%), moving the ratio from 0.94× to 1.09× — still no modality difference."

**2b.** Two checks → three, all named.

> before: "The co-active contrast survived two checks — disattenuation for measurement error, and partial-correlation control for source-genome GC composition, confirmed independently by phylum stratification (Section 3.1). Neither was pre-registered; both were motivated after the fact… and surviving both does not rule out a third."
> after: "The co-active contrast survived three checks — disattenuation for measurement error; partial-correlation control for source-genome GC composition on the in-vivo pairs, confirmed independently by phylum stratification (Section 3.1); and that same GC control applied to the modality contrast itself (Section 3.3). None was pre-registered; all three were motivated after the fact, by a subsequent review, by a finding in the cell-free dataset, and by a selection-artifact objection respectively, and surviving three does not rule out a fourth."

**2c.** Drafting residue removed.

> before: "The GC partial-correlation control applied to the in-vivo pairs in Section 3.1 **was not originally applied to this comparison**, which is the check that speaks most directly to a selection artifact."
> after: "The GC partial-correlation control applied to the in-vivo pairs in Section 3.1 **applies equally to this comparison, where it is the check that speaks most directly to a selection artifact, so we apply it here.**"

**Numbers changed: four added — 0.561, 0.517, 8.9, 21.1 — all from Task 2a, all read from `gate25_modality_gc_control.json` rather than typed. None removed.**

---

## TASK 3 — The GC control in the abstract [ASSESSED, NOT APPLIED]

**Proposed clause (12 words),** attached to the existing gap sentence rather than added as a new one:

> …while the cell-free correlation rises to **0.677**, a 2.6-fold gap **that survives partial-correlation control for source-genome GC composition (2.63× raw, 2.69× controlled)**.

**The abstract is already at 252 words**, so it is over the 250 ceiling before anything is added. Adding 12 takes it to 264; 14 words must come out.

**What I would cut, in order of preference:**

1. **"a barrier practitioners describe as unresolved for engineering non-model hosts" (10 words).** Pure framing. The chassis effect is named in the same sentence and the Introduction makes this point properly. Cheapest 10 words in the abstract.
2. **"over all measured sequences" (4 words)** from "Pooled over all measured sequences". "Pooled" is defined in the next clause by contrast with "restricted to sequences detectably active in both hosts", so the gloss is redundant.

Together that is exactly 14 words and lands the abstract back at **250**.

**What I would not cut, and why:** the three evidence lines (8 primary comparisons, 12 cells, 16 of 18) are the negative result, which is the paper's main contribution; the silent-class sentence is the mechanism; and "**This comparison rests on a single species pair**" is the scope limit that makes the rest honest. Cutting any of those to make room for a robustness check would be a bad trade.

---

## TASK 4 — Two carried-over items [ASSESSED, NOT APPLIED]

### 4a. Yim et al. Fig. 3 as an objection

**One paragraph if cited as-is; roughly half a gate if engaged properly.**

Their result is that cell-free cross-species correlation structure clusters by gram stain and tracks phylogeny (r = 0.73 against 16S distance). Read naively it says cell-free *does* reproduce host-to-host structure, which is the title's opposite.

The reconciliation is available and does not need new data: both can be true because they measure different things. Their gradient is about **how cross-host correlation varies with phylogenetic distance**; our claim is about **its absolute level**. Cell-free can preserve the ordering — closer species agree more — while compressing the whole band upward, and our own numbers say it does: the 45-pair cell-free band is **0.623–0.911 co-active** against an in-vivo co-active **0.258** for the one pair measured both ways. A phylogenetic gradient within a band that never drops near the in-vivo value is not evidence that lysates reproduce the chassis effect.

- **One paragraph in §3.3**, citing their Fig. 3, stating the distinction and quoting our 45-pair band: this is writing, not analysis, and nothing in it is unverified.
- **Properly engaged** would mean recomputing the phylogenetic gradient *under our restriction regimes* — we have all 45 pairs in `gate10_crosshost_correlations.csv` and DRAFTS released 16S distances — and reporting whether r = 0.73 holds co-active. That is a real analysis with a real risk: if the gradient strengthens under co-active restriction, it is a harder objection than it looks now.

**My read: the paragraph is v1, the recomputation is v2** — unless you want the stronger claim, in which case it should be done before submission rather than after a reviewer asks.

### 4b. The fitness-decoupling account

**One paragraph plus one Limitations clause. No new analysis, because our data cannot settle it.**

Yim et al. attribute greater in-vitro uniformity to expression being decoupled from cell fitness. That is the simplest competing explanation for the missing silent class: in a living cell a strongly expressed sequence can be costly and is selected against or shut down, while a lysate has no fitness to couple to. It predicts the same observation §3.4 explains by context — and it is a **better-motivated** competitor than the confounds already tested, because it comes from the dataset's own authors.

**We cannot distinguish the two.** Separating them needs a fitness or burden measurement per sequence per host, which neither dataset has. What we have — GC, activity in the other host, the silent fraction — is consistent with both accounts.

So the honest treatment is to name it in §3.4 as a competing explanation of equal standing, and add a Limitations clause saying we cannot separate them and what would. That is a strengthening, not a concession: §3.4 already frames itself as a hypothesis rather than a demonstrated cause, and naming the strongest competitor is what makes that framing credible.

**My read: both are v1.** 4b especially — a reviewer who knows the DRAFTS paper will raise it, and the paper is better for having named it first.

---

## Verification

| check | result |
|---|---|
| resource warnings | 0 |
| images | 9 |
| captions | Figures 1, 2a, 2b, 3, 4, 5, 6, 7, 8 — present, ascending |
| glyphs | ρ 17, ≤ 3, γ 2, β 2, ≈ 1 |
| pages | 20 |
| regime table | renders |
| `[email]` | 0 |
| "BS 15,848" | still intact on one line |
| cross-references | all resolve |
| both copies | byte-identical |
| audits | leakage ALL 6 PASSED; provenance 0 orphans (177 files) |

Body 6,172 → 6,215 words.

---

## Verdict

Figures current, review wrong, Gate 26 correct. Three text corrections applied, four numbers added, none removed. Tasks 3 and 4 assessed and not applied.
