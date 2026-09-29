**Body word count: 5,870, down from 6,240 (−370). Were any numbers lost? No. Zero numeric tokens disappeared from the manuscript; the script refuses to write if one does. The 0.912 reliability, the 0.8485 identity bound and all four disattenuation figures are present and confirmed in the rendered PDF.**

# GATE 24 MEMO — Apply the approved cuts

Script: `scripts/104_gate24_apply_cuts.py`. 13 cuts applied, 6 borderline candidates untouched. Both audits clean at start and end. `verify-citations` 11/0/5/0/16.

---

## What was applied

**C1–C10, verbatim from `out/GATE23_CUT_CANDIDATES.md`.** The candidates file carries no C-labels, so the mapping is by the order the file lists them; the ten spans sum to 38+23+36+35+26+12+18+35+42+60 = **325 words**, matching the file's stated total, which confirms the mapping.

| | § | words | what went |
|---|---|---|---|
| C1 | 1 | 38 | second Chan & Bernstein quote |
| C2 | 2.3 | 23 | forward pointer to the Discussion |
| C3 | 2.6 | 36 | the Amendment 1 sentence |
| C4 | 3.1 | 35 | re-explanation of the two regimes |
| C5 | 3.5 | 26 | restatement of the pre-registration |
| C6 | 3.6 | 12 | "We quantify the magnitude rather than…" |
| C7 | 3.7 | 18 | FiLM generalisation tail |
| C8 | 4 | 35 | recapitulation of §3.1's robustness checks |
| C9 | 4 | 42 | recapitulation of §3.2's pooled result |
| C10 | 4 | 60 | "On the form of the negative result" |
| C11 | 1 | 6 | "This paper addresses the second group." |
| C12 | 3.2 | 18 | "We report this first because…" |
| C17 | 2.7 | 24 | middle clause of the AI-assistance statement |

**Not touched:** all six borderline candidates (§1 bounded-claim, §2.2 pointer, §2.5 measured-number note, §3.2 "one well-characterized instance", §3.3 host-count sweep, §3.3 Pandi hedge tail).

**C11, C12 and C17 were not in the candidates file** — they were specified directly in this gate. C12 reverses the Gate 23 recommendation, which was to keep §3.2's sentence and cut the Discussion's duplicate (C9). Applying both removes the reporting-choice narration entirely, which is consistent with the principle stated in the gate, so both were applied as instructed.

**C13–C16 were "do not cut" instructions for spans that were never candidates** in `GATE23_CUT_CANDIDATES.md`, so they are no-ops. For the record: C13's "shift-prediction warning in §3.7" does not exist — the retraction is §3.8 and the actionable warning is Limitations item 8; C14 and C16 were argued for keeping in the candidates file's own prose and never listed; C15's listed candidate was the *hedge sentence* after Pandi, not the citation, and it is borderline and therefore untouched. Nothing was cut that the gate protected.

---

## No number lost — the explicit confirmation

**Zero numeric tokens disappeared.** Four now appear fewer times, each because a duplicate went:

| token | was | now | why, and where it survives |
|---|---|---|---|
| `2024` | 12 | 11 | C1's "(Chan & Bernstein, 2024)"; the same citation is in Introduction ¶1, the Discussion, and References |
| `100` | 3 | 2 | C3's "N=100"; survives in §2.6's H-MAIN statement and the Figure 4 caption |
| `1` | 10 | 9 | the "**Amendment 1**" label itself |
| `4` | 11 | 10 | "unfrozen-conv4" |

**The three figures you asked me to confirm explicitly, counted in the rendered PDF text, not just the Markdown:**

- **0.912** reliability — present, **2 occurrences** (§3.1, and Limitations item 3)
- **0.8485** identity bound — present, **2 occurrences** (§2.5 twice: the bound and the per-fold list)
- **disattenuation figures** — **0.341** (5), **0.515** (1), **0.308** (1), **0.240** (1), **0.286** (1), and the cross-check **0.929** (1), all present

The script asserts all eight before writing and fails closed if any is missing.

### One tokenizer defect found and fixed

The first run failed with `numeric tokens would disappear entirely: ['100,']`. That was a false alarm from my own regex: a greedy `[\d,]*` swallowed the trailing comma, so `N=100,` in the cut sentence and `N=100 ` elsewhere read as different figures. The tokenizer now requires a token to end in a digit. Worth recording because the guard would have blocked a correct cut, and the tempting fix was to delete the guard.

---

## Two seams left dangling, applied as instructed and not repaired

Both are consequences of cutting a sentence that introduced something the next sentence refers to. The gate says "Apply only these," and repairing either needs words added, so neither was touched.

**1. Discussion, after C8.** The paragraph now opens:

> **Robustness, not proof.** Neither correction was pre-registered; both were motivated after the fact…

**"Neither correction" has no antecedent.** The two corrections — disattenuation and the GC partial correlation — were introduced by the sentence C8 removed. A reader arrives at "Neither correction" cold. Minimal repair, if you want it: *"Neither of the two corrections in Section 3.1 was pre-registered."*

**2. §2.6, after C3.** The paragraph now reads:

> Two amendments, both dated 2026-08-05 and both made before H-MAIN was evaluated. **Amendment 2** widened the decision interval…

**It announces two amendments and describes one.** Amendment 1 is now never named. Minimal repair: drop the "**Amendment 2**" label to "The second widened…", or say "two amendments, of which only the second bears on any reported number."

---

## One place the stated rationale does not hold

C17's reasoning is that §2.5 already establishes the clause being dropped. It does not. The dropped clause was:

> Every reported number is produced by a committed script in the public repository and verified by the audit suite described in Section 2.5.

§2.5 is *Splitting*; it describes the leakage defences and the standing audit, but it makes no claim about the provenance of every reported number. After the cut, **the manuscript no longer states anywhere that every number it reports comes from a committed script.** The nearest survivor is Data and Code Availability, which lists "the leakage, provenance, and restriction-regime audits" as available code — naming the audit without making the claim.

This was applied as specified, because it is a content judgement and yours to make. Flagging it because it is the paper's only statement of script-level provenance, and because the audit it referred to is the one this project has run at every gate.

**Side effect:** that clause held the manuscript's only reference to Section 2.5, so **§2.5 is now cross-referenced from nowhere.** The heading remains and nothing is unresolved.

---

## Verification

| check | result |
|---|---|
| body word count | **5,870** (from 6,240, −370) |
| resource warnings | **0** |
| images | **9** |
| caption order | **1–8 ascending** (255, 283, 314, 432, 468, 500, 539, 582) |
| glyphs | ρ 21, ≤ 3, γ 2, β 2, ≈ 2 — all present, all unchanged |
| pages | **21** (unchanged; the cuts reclaimed less than a page) |
| regime table | renders as an aligned 6-row table |
| both copies | byte-identical, md5 `a6f96c999971c040632c9f9ab0f720b9` |
| cross-references | 9 `Section` refs, **all resolve**; `Section 2.5` no longer referenced (see above) |
| Limitations | 12 items, sequential |
| figure embeds | 9 |
| `verify-citations` | 11 / 0 / 5 / 0 / 16 — no regression |
| whitespace | 0 triple newlines, 0 space-before-punctuation; the only double spaces are the YAML preamble's own indentation |
| §2.7 | matches the gate's specified replacement text exactly |

---

## Verdict

13 cuts applied, 370 words out, no number lost, build clean. Two dangling antecedents and one over-broad rationale are reported above rather than patched, since each would need words added that this gate did not authorise.
