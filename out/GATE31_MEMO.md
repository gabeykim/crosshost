**Where did the clause land, and did any number change? §3.5's opening, immediately after the sentence that makes the confirmation claim. No number changed — the clause contains no digits, and the guard confirmed numeric-token counts are identical in both directions.**

# GATE 31 MEMO — The §3.5 inference clause

Script: `scripts/115_gate31_inference_clause.py`. Both audits clean at start and end. Abstract unchanged at 250 words; body 6,794 → 6,828.

---

## TASK 1 — The acknowledgement [APPLIED at §3.5's opening]

**Before:**
> Four independent lines of evidence test the prediction above, and each comes out as Section 3.4 says it should.

**After:**
> Four independent lines of evidence test the prediction above, and each comes out as Section 3.4 says it should. **With two training hosts per fold, these features would also fail if host specificity were mechanistic but simply not learnable at this n_hosts; the prediction is consistent with the result, not confirmed by it.**

### Why §3.5's opening rather than §3.4's closing

Two reasons.

1. **It is where the claim being qualified is made.** "Each comes out as Section 3.4 says it should" is the confirmation claim; the qualification belongs against it, not two paragraphs earlier. A reader now meets the concession *before* the four evidence lines, so it frames them rather than arriving as an afterthought.
2. **§3.4's closing already concedes the same point from the other direction** — "the link from them to the modeling failure in Section 3.5 is an inference". Adding the clause there would state one concession twice in a single paragraph, which reads as anxiety rather than precision.

The four other framing sites — abstract, Introduction, §3.4 closing, §3.5 heading — are untouched, per instruction.

---

## TASK 2 — Redundancy with §3.7

**They do not read as redundant, and I propose no change.**

The two sentences share a premise and make different claims:

| | §3.5 (new) | §3.7 (existing) |
|---|---|---|
| claim | the negative result cannot discriminate between a contextual account and a mechanistic one that is merely unlearnable at this n_hosts | FiLM's fold-to-fold variance is 4.0–7.8× worse than the simpler alternatives, and a γ/β generator fit from two vectors is a poor default |
| type | a scope limit on an **inference** | a **measured** instability, with numbers |
| about | all four lines of evidence | FiLM specifically |

They are also two sections apart with §3.6 between them, so a reader does not encounter them in sequence.

**One optional addition I did not make**, since Task 2 scopes a cross-reference to the redundant case: appending "(Section 3.7 measures what that costs)" to the new clause would turn an asserted limitation into a checkable one. Two words, no number, and it would make the §3.5 concession point at the paper's own evidence for it. Say the word and it is a one-line gate.

---

## TASK 3 — Verification

### A build failure to disclose first

**The first `make pdf` of this gate failed with exit 43 — `fclose: No space left on device` — and I nearly reported a stale PDF as verified.** The volume was at 98% capacity. `make` returned non-zero, but every downstream check still produced plausible output, because they were reading the **Gate 30** `manuscript.pdf` from 17:35:58 against a manuscript last written at 20:07:15.

I caught it on the exit code, discarded that entire round of checks, and rebuilt. The second build succeeded (exit 0) with 5.8 GB free, so the condition was transient. Every number in the table below comes from the rebuilt PDF, confirmed newer than the manuscript (20:09:32 against 20:07:15) and confirmed to contain the new clause.

This is the same failure mode Gate 25 produced and Gate 27 was opened to adjudicate: a step that does nothing, reports success downstream, and leaves a stale artifact behind. The lesson it teaches is the one that caught it here — check the exit code, and check that the output is newer than its input.

Also noted: **`out.zip` no longer exists.** It has been in every status report since Gate 19 as untracked and stale; it was deleted outside this session.

| check | result |
|---|---|
| new clause present in the PDF | **yes** |
| PDF newer than the manuscript | **yes** (20:09:32 > 20:07:15) |
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
| abstract | **250 words**, unchanged |
| body | 6,794 → **6,828** |
| audits | leakage ALL 6 PASSED; provenance 0 orphans (177 files) |
| `verify-citations` | **11 / 0 / 5 / 0 / 16 — full baseline, zero skips.** The `dnabert2` timeout of Gates 29–30 has cleared. |

---

## Numbers

**None changed.** The clause contains no digits. The guard compares numeric-token counts in both directions and reported them identical, and it also asserts the abstract word count is untouched.

---

## Verdict

One clause, at the point of inference, with the other four framing sites left alone. Not redundant with §3.7. A transient disk-full build failure was caught on its exit code before its stale output could be reported as verification.
