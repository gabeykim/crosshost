**Did all five decisions apply cleanly? Yes. Does the PDF still break "BS" from "15,848"? No — fixed. A non-breaking space between each host label and its count prevents the line break, so "BS 15,848, PA 21,473" now sits intact at the top of page 3 and the misreading is no longer possible. The number itself was not touched.**

# GATE 26 MEMO — Applying the Gate 25 decisions

Script: `scripts/109_gate26_apply_decisions.py`. Both audits clean at start and end. Body 5,938 → 6,172 words.

---

## TASK 1 — Supercoiling [APPLIED to §3.2 and §3.4]

**§3.2, before:**
> …which retain RNA polymerase, sigma factors, and ribonucleotides but lack an intact membrane, **native supercoiling**, macromolecular resource competition, and growth-phase physiology.

**§3.2, after:**
> …which retain RNA polymerase, sigma factors, and ribonucleotides but lack an intact membrane, **supercoiling homeostasis**, macromolecular resource competition, and growth-phase physiology.

**§3.4, before:**
> **DNA supercoiling, absent in linear cell-free reactions**, acts as a global transcriptional regulator whose effect on a promoter depends quantitatively on its discriminator GC content (El Houdaigui et al., 2019).

**§3.4, after:**
> **DNA supercoiling** acts as a global transcriptional regulator whose effect on a promoter depends quantitatively on its discriminator GC content (El Houdaigui et al., 2019). **DRAFTS supplied its library as midiprepped plasmid, so those templates are supercoiled; what a lysate lacks is the gyrase/topoisomerase homeostasis and chromosomal context that regulate supercoiling in a living cell, leaving it fixed at whatever the preparation produced rather than dynamically maintained. That is a weaker claim than absence, and it is the accurate one.**

### A third mention exists, and I did not change it

**§2.3** also names supercoiling: "…closer in kind to the genomic vector's machinery annotations than to the membrane potential, **supercoiling**, resource competition, or growth-phase transitions implicated by the mechanistic hypothesis in Section 3.4."

This is a list of what the hypothesis *concerns*, not a claim that lysates lack supercoiling, so it stays true under the correction and needs no edit. Reporting it because the instruction was to stop if the fix required a third passage — it does not, but you should know the mention is there.

---

## TASK 2 — Title [APPLIED, option 2]

> **Cell-free transcription does not reproduce the bacterial chassis effect: evidence that host specificity resides in cellular context rather than transcription machinery**

Confirmed in the rendered PDF.

---

## TASK 3 — GC control into §3.3 [APPLIED]

One paragraph, placed after the regime-dependence material and before the Pandi comparison, so it sits with the other scope checks rather than in the headline. Carries all five required elements and the scope limit: cell-free 0.677 → 0.556 (17.9%), in vivo 0.258 → 0.207 (19.8%), ratio 2.63× → **2.69×**, pooled 0.94× → 1.09×, both methods agreeing (0.556/0.548; 0.207/0.205), and that GC is controlled within each modality separately because the two populations differ — not a paired control.

---

## TASK 4 — Figure 5 [SYSTEMS ADDED; the panel stays readable]

Added concatenation, per-host-heads-averaged and per-host-heads-nearest. **Nine systems per panel, and it is still legible** at 30° rotated labels, so no text change was needed. EC transcription now shows concatenation **0.555** and per-host-heads-averaged **0.615** against sequence-only 0.367 and FiLM 0.371 — exactly the values §3.5 cites. The citation resolves.

The merge is legitimate rather than convenient: `gate8_5_conditioning_mechanisms.csv` and `gate6_full_comparison.csv` are on identical footing — same five folds, same 90% bootstrap scheme, and `sequence_only` 0.367130 and `film_genomic`/`genomic` 0.371166 agree to six decimal places across the two files. The caption records this.

---

## TASK 5 — Figure 2b [ADDED Cyanobacteria, rather than explaining the gap]

**Chose to add it.** Stating the threshold in the caption would have required explaining why a stratum that *clears* the stated threshold is absent, and there is no principled reason — the four-phylum list was hardcoded. Adding it makes figure, caption, text and results file all say the same thing.

Cyanobacteria: EC–PA 0.543 (n = 455), EC–BS 0.332 (n = 62), BS–PA −0.022 (n = 44). §3.1 now reports five strata. Two strata (Euryarchaeota, Planctomycetes) remain unplotted because at least one pair falls below threshold; the caption says so.

**This slightly weakens the stratified claim, and the text now says so.** Cyanobacteria has the highest EC–BS of any stratum at 0.332, so the EC–PA margin over the *B. subtilis* pairs is narrowest there (0.543 against 0.332) on the smallest N. §3.1 states this rather than leaving a reader to notice it. The ordering still holds in every stratum.

---

## TASK 6a — Provenance disclosure [WORDING PROPOSED, NOT APPLIED]

**Proposed, for Data and Code Availability** (one sentence, appended to the paragraph describing the audits):

> The provenance audit matches output filenames against script text, which is a substring check rather than a dependency check; from Gate 20 to Gate 24 it returned no orphans for the nine preprint figures although no committed script produced them, because the figure-embedding script mentions every filename. Producing scripts were added in Gate 25 and the audit now covers them genuinely.

Alternative, shorter, **for Limitations** (as a clause on the existing audit item) if you would rather it not sit in Data and Code Availability:

> The provenance audit's filename matching gave a false pass on the nine preprint figures from Gate 20 to Gate 24 — the embedding script mentions every filename — until producing scripts were added in Gate 25.

I lean to Data and Code Availability: that is where the reproducibility claim is made, so it is where the caveat belongs. Not applied either way, per instruction.

---

## TASK 6b — The page break [FIXED]

**Before:** page 2 ended "…Per-host usable transcription counts: EC 24,613, BS" and page 3 opened "15,848, PA 21,473", with the page number **2** printed between them. Three readers read that as "BS 215,848".

**After:** a non-breaking space (U+00A0) between each host label and its count. TeX cannot break there, so the whole group moves together: page 2 now ends "…EC 24,613," and page 3 opens "**BS 15,848, PA 21,473.**"

Chosen over `\needspace` because it needs no preamble package and cannot push a figure around; chosen over reflowing the sentence because that would have meant editing text to fix typesetting. Applied to all three labels so the list behaves consistently. **No number changed** — only the character class of three spaces.

---

## Numbers that changed, and why

| number | why |
|---|---|
| 0.205, 0.548, 0.556, 1.09, 19.8, 2.63, 2.69 | Task 3, the GC-control paragraph |
| 0.022, 0.332, 0.543, 44, 62, 455 | Task 5, Cyanobacteria's three values and their N |
| 0.33 | Task 5, the upper end of the *B. subtilis*-pair range, recomputed over five strata |
| **0.27 (removed)** | the previous upper end of that range over four strata. It is a recomputed summary, and every underlying per-stratum value is still printed in the same sentence. This is the only token to leave the manuscript, and the script required an explicit documented exception to allow it. |

Nothing else changed. The guard refuses to write if any other numeric token disappears.

---

## Verification

| check | result |
|---|---|
| resource warnings | 0 |
| images | 9 |
| captions | Figures 1, 2a, 2b, 3, 4, 5, 6, 7, 8 — present, ascending |
| glyphs | ρ 17, ≤ 3, γ 2, β 2, ≈ 1 — all non-zero |
| pages | 20 |
| regime table | renders |
| `[email]` / Zenodo DOI | 0 / present |
| cross-references | 9 refs, all resolve |
| Limitations | 12 items |
| both copies | byte-identical |
| audits | leakage ALL 6 PASSED; provenance 0 orphans (176 files) |
| `verify-citations` | 11 / 0 / 3 / **2 skipped** — still HTTP 429 rate-limiting, as predicted. Noted, not chased. |
| "BS 15,848" | intact on one line, top of page 3 |

---

## Verdict

All five decisions applied. One third supercoiling mention found and left in place with reasons. Task 5 resolved by adding the stratum rather than explaining its absence, which narrows one margin — stated in the text. Task 6a proposed in two forms and not applied. Task 6b fixed with a non-breaking space and no change to the number.
