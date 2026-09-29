# GATE 23 TASK 2 — Cut candidates

**Nothing in this file has been applied.** Task 2 is a report; the manuscript on disk carries Task 1 only.

Body word count as it stands (excluding references and figure captions): **6240**.

- Apply the 10 **CUT** recommendations: **5915** words (−325).
- Apply CUT + the 6 **BORDERLINE** ones as well: **5756** words (−484).

Every quoted string below was verified to appear exactly once in `out/PREPRINT/MANUSCRIPT.md`; word counts are of the quoted span, not estimates.

## Recommended cuts

Each of these fails at least one **Cut if** test and passes none of the **Keep if** tests.

| § | words | why | text |
|---|---|---|---|
| 1. Introduction | 38 | restates | and that "it is impractical to postulate that the observable chassis-effect between a given set of hosts can be explained by a single or even a set … |
| 2.3 | 23 | restates | That distinction bears on how these results relate to prior physiology-based claims (Discussion) and on why this vector's failure does not … |
| 2.6 | 36 | process | **Amendment 1** added the supplementary unfrozen-conv4 transfer mechanism alongside the pre-registered frozen-trunk one at N=100, after validation … |
| 3.1 | 35 | restates | The two regimes answer different questions. Co-active asks how well two hosts agree on the relative strength of sequences that work in both; pooled … |
| 3.5 | 26 | restates | The hypothesis was committed before any host-conditioned model was trained; both amendments are dated, preceded the results they touch, and make the … |
| 3.6 | 12 | presentation | We quantify the magnitude rather than claiming it as a novel warning. |
| 3.7 | 18 | restates | , which matters for any model conditioned on a handful of domains independent of this paper's negative result |
| 4. Discussion | 35 | restates | The co-active EC–PA-versus-*B. subtilis* contrast has survived two attacks aimed at explaining it away — disattenuation for measurement error, and … |
| 4. Discussion | 42 | restates | The modality contrast survived a third check, the one that came closest to overturning it: pooled against pooled, the two modalities are … |
| 4. Discussion | 60 | restates | **On the form of the negative result.** The evidence in Section 3.5 was run because Sections 3.2–3.4 predicted, before any model was trained, that … |

**§1. Introduction (38 words, restates).** Second Chan & Bernstein quote. The clause immediately after it ('that host physiology explains chassis effects and genome sequence alone should not be expected to') is a complete paraphrase, and para 1 already cites the same paper for the same point. One quote carries the hypothesis; two restate it.

**§2.3 (23 words, restates).** A forward pointer to an argument the Discussion makes in full ('What this says about the physiology hypothesis'). The preceding sentence already establishes the distinction; this one only announces that it matters later.

**§2.6 (36 words, process).** Process detail about a mechanism that by its own admission never entered the headline result. The operative facts -- both amendments predate evaluation, both make MET harder -- are stated in the surrounding sentences and survive. Nothing downstream cites unfrozen-conv4.

**§3.1 (35 words, restates).** Section 2.2 defines both regimes, including that pooled includes sequences measured as inactive. Keep the operative sentence that follows ('We use co-active as primary because it isolates strength from silence').

**§3.5 (26 words, restates).** Every clause is already in section 2.6: the commit date, that both amendments preceded evaluation, and that both make MET harder to declare.

**§3.6 (12 words, presentation).** Explains a framing choice rather than making one. The LaFleur comparison that follows does the work of showing this is a known effect being quantified.

**§3.7 (18 words, restates).** The Discussion says this in its own paragraph ('worth carrying into any host-conditioned model built on a handful of domains'). The generalisation belongs there, not twice.

**§4. Discussion (35 words, restates).** A recapitulation of section 3.1's two robustness checks. What is NOT in 3.1 -- that neither correction was pre-registered, that both were motivated after the fact, and that surviving two checks does not rule out a third -- is the rest of the paragraph and must stay.

**§4. Discussion (42 words, restates).** Restates section 3.2's pooled result, and its second sentence duplicates 3.2's own 'We report this first because it is the comparison a reader is most likely to compute independently.' Cut here; keep the statement where the choice is actually made.

**§4. Discussion (60 words, restates).** Near-verbatim restatement of the Introduction, which already says the prediction was pre-registered, that four independent attempts followed, that two looked like counterexamples, that one survives narrowly and one failed its control.

## Borderline

I would keep some of these. Flagged because you asked for every candidate, with my own lean stated.

| § | words | why | text |
|---|---|---|---|
| 1. Introduction | 12 | presentation | This is a bounded claim about five named suites, each read directly. |
| 2.2 | 6 | presentation | Section 3.1's primary figures are co-active. |
| 2.5 | 32 | presentation | We report this as a measured number rather than as a claim that clustering was applied, since approximate clustering carries no completeness … |
| 3.2 | 12 | double hedge | The claim should be read as one well-characterized instance, not a survey. |
| 3.3 | 65 | process | **A host-count sweep on DRAFTS was considered and not run.** DRAFTS's ten-host panel could in principle test whether the conditioning failure in … |
| 3.3 | 32 | double hedge | That concerns a different quantity — absolute yield within a single host — but supports the general caution this section sharpens by identifying … |

**§1. Introduction (12 words, presentation).** Explains the presentation rather than making a claim. The preceding sentence already names all five suites, which is the bound. Borderline because the explicitness is a defence against an over-reading a reviewer might make.

**§2.2 (6 words, presentation).** A pointer, not a finding. Cheap to keep at 7 words and it pre-empts a real ambiguity, but it is a presentation statement.

**§2.5 (32 words, presentation).** Explains a reporting choice. Borderline, and I would keep it: it tells a reader why to trust 0.8485 over 'we clustered', which is the difference between a measured guarantee and an asserted one.

**§3.2 (12 words, double hedge).** Hedges a claim already hedged twice: the paragraph's own bold lead ('rests on a single species pair') and Limitations item 2. Borderline because it is the sentence a reader quotes back at you.

**§3.3 (65 words, process).** Describes work not done. The reasoning is sound and pre-empts a reviewer, but it answers a question no reader has asked yet and changes nothing they should believe or do. If kept, it compresses to one clause.

**§3.3 (32 words, double hedge).** Concedes the citation does not bear on the claim, then asserts that it does. Keep Pandi's 0.41 and the quotation; this sentence is the hedge on the hedge.

## Scanned and kept

Sections you flagged where I did not find a cut worth making.

| § | words | why | text |
|---|---|---|---|
| 1. Introduction | 8 | scope | Part of the field routes around the question. |

**§1. Introduction (8 words, scope).** Scanned and kept. The orthogonal-machinery paragraph (78 words) is not a restatement and it scopes who the paper is for, which changes what a reader should do with it.

## Where the words are, and where they are not

The Discussion is the densest source: **three of the sixteen candidates are there and they account for 178 of the 325 recommended words** — more than half. All three are recapitulations of Results or the Introduction. The Discussion's own arguments (the physiology-hypothesis paragraph, what would settle the open questions, the epistemic frame on the robustness checks) are not candidates and should not be touched.

Two places you flagged yielded less than expected:

- **§2.6** gives up only the Amendment 1 sentence (42 words). Amendment 2's content — the 90% percentile bootstrap, 10,000 resamples, hierarchical resampling — is cited by §3.5 and by the captions of Figures 4 and 5, so it cannot be compressed without breaking them.

- **§3.3** yields two borderline candidates and no clean cut. Its bulk is the six-row regime table and the paragraph that makes the table checkable; the enumeration of what was searched in DRAFTS is what makes the novelty claim verifiable, so it stays by the **Keep if** test on checkability.

Nothing in the Abstract, §2.1, §2.4, §3.4, §3.8 or §3.9 reached candidate status. §3.8 is almost entirely controls that killed our own results, which the **Keep if** tests protect explicitly.

