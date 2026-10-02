**Did Tasks 1–4 apply, and what is the abstract's final word count? All four applied. The abstract is exactly 250 words — unchanged, with the 93.6% / 26.0% contrast and the practitioner takeaway added and 30 words cut from elsewhere to pay for them. No measured value changed; the only numbers added anywhere are 93.6 and 26.0, both already in §3.2.**

# GATE 30 MEMO — Framing edits from a cold read

Script: `scripts/114_gate30_framing.py`. Both audits clean at start and end. Body 6,717 → 6,794 words.

---

## TASK 1 — §3.4: the competing account falls on the same side of the dichotomy [APPLIED]

**Before** (paragraph ended):
> …We name it as a competing explanation of equal standing rather than one this paper rules out.

**After** (two sentences added):
> …We name it as a competing explanation of equal standing rather than one this paper rules out. **Both accounts nonetheless place the difference in cellular context rather than in transcriptional machinery — expression coupled to growth is a property of a living cell, not of its polymerase. What the data available here cannot settle is which contextual mechanism dominates, not whether the difference is contextual.**

This is the highest-value edit in the gate: the paragraph previously conceded a competitor without noting that the competitor is on the paper's own side of the title's claim. A reader reaching it had no way to see that.

---

## TASK 2 — Abstract [APPLIED, 250 → 250]

**Added:**
- 2a: "…on which sequences are silent: **26.0% of these sequences are co-active in vivo against 93.6% in lysate**." The sentence previously ended "lysates silence almost nothing" with no number.
- 2b: "Screening constructs for one fixed chassis via lysate is therefore reasonable; ranking a part across candidate hosts is not."

**Cut, to pay for them — 30 words, none from the protected set:**

| cut | from | to | saved |
|---|---|---|---|
| opening framing clause: "we find that cross-host agreement responds to conditioning on co-activity very differently in the two modalities" | 36 | 21 | **15** |
| pooled sentence, compressed to "Pooled, the two modalities are indistinguishable (0.616 cell-free, 0.655 in vivo)" | 17 | 11 | **6** |
| prediction lead-in: "should carry no signal, and we tested that prediction directly" → "carry no signal, and we tested it" | 18 | 16 | **2** |
| availability: "baseline suite… a retracted rescue attempt" → "baselines… a retraction" | 19 | 16 | **3** |
| (net of the two additions, +27 and +19 against the 19-word sentence replaced) | | | |

**Nothing protected was touched.** The three §3.5 evidence lines, the single-species-pair limit, the n_hosts ≤ 6 bound and the GC-control clause are all verbatim as they were — which is why this was applied directly rather than reported first. The pooled figures 0.616 and 0.655 survive the compression; only the words around them went.

---

## TASK 3 — §3.3 renamed [APPLIED]

Options considered:

1. **"Magnitude, robustness, and the published objections"** ← chosen
2. "How large the contrast is, and what it survives"
3. "The contrast under scrutiny: magnitude, controls, and prior art"

Chosen because it names the three things actually in the section — the regime table, the GC control, and the Yim/Pandi engagement — in the order they appear, and because "published objections" is the phrase that tells a skimmer the strongest prior-art challenge is answered here. Option 2 is more readable but vaguer about prior art; option 3 names it but reads like a subtitle.

Content was not moved into §3.2, as instructed. All `Section 3.3` cross-references still resolve, since the number is unchanged.

---

## TASK 4 — The other open check is now cross-referenced [APPLIED]

Appended to the Yim paragraph's flag:

> **It is one of two checks this paper leaves open; the other, cross-modality correlations for the three RS241 hosts DRAFTS covers, is stated in Limitations item 2.**

Verified: Limitations item 2 is the modality-contrast item and does contain the RS241 clause. There are 12 items, and this is now the manuscript's only `Limitations item N` reference — it resolves, and the check for that is in the verification table below.

---

## TASK 5 — Decoupling §3.5 from the modality claim [ASSESSED, NOT APPLIED]

### 1. Every place §3.5 is framed as testing the §3.4 prediction

| # | location | text |
|---|---|---|
| 1 | **Abstract** | "That interpretation predicts host-descriptor features encoding machinery composition carry no signal, and we tested it:" followed by the three evidence lines |
| 2 | **Introduction** | "Comparing the two (Section 3.2) motivates a testable prediction: if host specificity is substantially contextual rather than mechanistic, host-descriptor features… should carry no signal… We pre-registered a test of this before training any model. It held on every primary comparison." |
| 3 | **§3.4 closing** | "the link from them to the modeling failure in Section 3.5 is an inference. If host specificity resides substantially in context-dependent silencing… then host-descriptor features… should carry no signal… however well constructed." |
| 4 | **§3.5 heading** | "Testing the prediction (Figures 4–5)" |
| 5 | **§3.5 opening** | "Four independent lines of evidence test the prediction above, and each comes out as Section 3.4 says it should." |

**The Discussion no longer frames it.** The "On the form of the negative result" paragraph — which said the §3.5 evidence "was run because Sections 3.2–3.4 predicted… that it should come out this way" — was cut in Gate 24 as a restatement of the Introduction. So the Discussion is already decoupled, and four of the five remaining sites are short.

### 2. The minimal decoupling

Change sites 4 and 5 and soften 1 and 3; leave 2 intact.

- §3.5 heading → "Host-descriptor features carry no signal (Figures 4–5)"
- §3.5 opening → "Four independent lines of evidence test whether host-descriptor features carry signal for cross-host prediction. None does."
- Abstract → "…carry no signal. We tested that separately:" (the conditional stated, the confirmation not claimed)
- §3.4 → keep "the link from them to the modeling failure in Section 3.5 is an inference", which already concedes the point, and drop "however well constructed"

The sections stay adjacent and the conditional stays stated. What goes is the claim that the second result *confirms* the first.

### 3. What would be lost

**Less than it appears, and the paper survives.** §3.5 is already the stronger half on its own terms: a pre-registered kill gate that failed on all 8 primary comparisons, an arbitrary host tag matching or beating both feature vectors in 12 of 12 cells, and sequence-only not distinguishably beaten in 16 of 18 tests. None of that depends on §3.4 being right. As an independent negative result it reads: *host-descriptor features built from machinery annotations do not help at this n_hosts, by four independent tests, one of them pre-registered.* That is a complete, publishable claim.

What would be lost is the paper's narrative spine — the move from an observation about two datasets to a prediction to a pre-registered test. The title's "evidence that host specificity resides in cellular context" would then rest on §3.2–§3.4 alone, which is one species pair. That is a real cost: §3.5 currently supplies the paper's only *pre-registered* evidence for the title, and decoupling removes it from that role.

**My read: the readers are right about the inference and wrong about the remedy.** The two-host limitation does weaken the inference — §3.7 says so in the paper's own voice ("a γ/β generator fit from only two training-host vectors is a poor default in this few-domain regime"), and Limitations item 1 bounds everything at n_hosts ≤ 6. But full decoupling throws away a correct conditional to avoid overstating it.

### 4. The one-clause alternative

Add to §3.5's opening, or to §3.4 where the inference is already conceded:

> With two training hosts per fold, these features would also fail if host specificity were mechanistic but simply not learnable at this n_hosts; the prediction is consistent with the result, not confirmed by it.

**This is the option I would take.** It costs one sentence, it states the readers' objection in the paper's own voice before a reviewer makes it, it keeps the narrative, and it is more honest than either the current framing or full decoupling — because the conditional *is* consistent with the result, and that is all anyone can say.

---

## TASK 6 — Verification

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
| cross-references | all `Section X.Y` resolve; `Limitations item 2` resolves (12 items) |
| both copies | byte-identical |
| abstract | **250 words** |
| body | 6,717 → **6,794** |
| audits | leakage ALL 6 PASSED; provenance 0 orphans (177 files) |
| `verify-citations` | 11 / 0 / 4 / **1 skipped** — `dnabert2`, the same transient network timeout as Gate 29; VERIFIED and unexplained-mismatch counts unchanged |

**Numbers: 93.6 and 26.0 each +1 occurrence** (moved into the abstract from §3.2, authorised). Task 4's cross-reference adds one `2` and one `241`, both identifiers rather than measured values. **Nothing else changed**, and the guard refuses to write on any other delta.

---

## Verdict

Four framing edits applied, abstract held at exactly 250, no measured value altered. Task 5 assessed with a recommendation: the one-clause acknowledgement rather than full decoupling.
