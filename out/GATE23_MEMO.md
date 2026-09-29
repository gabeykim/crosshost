**Body word count after Task 1: 6,240. Task 2 identifies 16 candidates totalling 484 words — 325 I recommend cutting, 159 borderline.**

**Two premise corrections first. The manuscript is not a hybrid: it is byte-identical to the Gate 21 commit, with the *uncompressed* Results and Methods (6,349 body words, not ~4,650), and there is no `## 6. What we could not do` section to merge — Gate 22 verified its absence, and that was one of the three markers that matched. The compressed draft has still not reached disk. Task 1 was applied anyway, because the 19-item Limitations is exactly as this gate describes it and the 12-item target list is fully specified; Task 2 is reported against the 6,240-word file that exists.**

# GATE 23 MEMO — Limitations to 12 items, and a full cut scan

Script: `scripts/103_gate23_limitations.py`. Both audits clean at start and end. `verify-citations` 11/0/5/0/16, unchanged since Gate 13.

---

## Premise check

| this gate assumes | actual state |
|---|---|
| compressed Results and Methods | **uncompressed.** Results 3.1–3.9, body 6,349 words before this gate's edit |
| a `## 6. What we could not do` section to merge in | **does not exist.** Gate 22 verified this; nothing to merge |
| §2.1 lacks the denominator note (Gate 22 Task 2 never ran) | **§2.1 already has it**, including 52.0% / 18.9% / 83.8%, the 17.9% recomputation, the 1.0-point difference, and "the percentages are not comparable across denominators" |
| the ceiling metric is absorbed into §3.7 | absorbed into **§3.8** under this numbering; §3.8 already ends "It appears nowhere in this paper's results" |
| Data Availability has an existing `reproduce-full` clause | **it does not.** One was written, carrying the numbers (see below) |

The 19-item Limitations, on the other hand, is exactly as described, and the 12-item target list specifies every number. So Task 1 went ahead. The two things that depended on Gate 22's unapplied edits both resolved in the gate's favour: §2.1 already explains the 17.9/18.9 difference, so cutting that clause loses nothing, and the Introduction already records the four attempts, so cutting the fifth-attempt item loses nothing.

---

## TASK 1 — Limitations 19 → 12

**940 words → 712 in the section. Zero numbers left the manuscript.** The script asserts this before writing: it strips list markers (so renumbering 19 items down to 12 does not read as a loss) and refuses to write if any number in the old section is absent from the new document.

| old item | disposition |
|---|---|
| 1 | → new 1, unchanged |
| 2 | → new 2, unchanged |
| 3 | → new 3, plus the Pearson clause absorbed from old 17 |
| 4, 5, 6, 7, 8 | → new 4, 5, 6, 7, 8, unchanged |
| 9 + 10 | **merged** → new 9, "Three dataset collisions" |
| 11 (ceiling metric) | **cut.** No numbers. §3.8 already states it appears nowhere in the results |
| 12 | → new 12, unchanged |
| 13 (~300× inference cliff) | **relocated to §2.4**, carrying ~300× and the 1,024-row chunk size |
| 14 + 15 | **merged** → new 10, Evo 2 and the non-native protocols |
| 16 | → new 11, unchanged |
| 17 | **split:** Pearson clause → item 3; 17.9%/18.9%/1.0-point cut, since all three are in §2.1 |
| 18 (`make reproduce-full`) | **relocated to Data and Code Availability**, carrying Python 3.14.2, the six defects, four iterations, and 40+ hours |
| 19 (no fifth attempt) | **cut.** The Introduction and Discussion both record the four attempts |

### Every number dropped, in full

**None.** Four numbers now appear fewer times because a duplicate was removed, and each is still in the manuscript:

| number | occurrences | still at |
|---|---|---|
| 17.9% | 2 → 1 | §2.1 |
| 18.9% | 2 → 1 | §2.1 |
| 1.0-point | 2 → 1 | §2.1 |
| "Section 3.8" | 4 → 3 | Introduction, Limitations item 8 |

Nothing is entirely absent and nothing is entirely new. Verified by a counted diff of every numeric token against `HEAD`.

### Two relocations, and why

The gate says both "reduce to 12 items" and "Preserve every number," and adds: if absorbing an item would lose a number that appears nowhere else, say so and keep it. Two items were designated for cutting but carried numbers found nowhere else in the manuscript — `~300×`/`1,024`, and `Python 3.14.2`/`40+ hours`/`six defects`/`four iterations`. Cutting them outright would have broken "preserve every number"; keeping them as items would have broken "12 items". Moving them to the sections where they belong satisfies both: an implementation detail into §2.4 Architecture, and a reproducibility statement into Data and Code Availability.

This is why the net body saving is only **109 words** (6,349 → 6,240) against 228 removed from Limitations: about 119 words came back as relocations. If you would rather those two clauses simply go, say so and the body drops to roughly 6,120.

---

## TASK 2 — Cut scan

In `out/GATE23_CUT_CANDIDATES.md`. **Nothing applied.** 16 candidates, every quoted span verified to appear exactly once, word counts measured rather than estimated.

- **325 words recommended** → body **5,915**
- **+159 borderline** → body **5,756**

The Discussion is the densest source: three candidates, 178 of the 325 words, all three recapitulations of Results or the Introduction. Two of the places you flagged yielded less than expected — §2.6 gives up only the Amendment 1 sentence (Amendment 2's bootstrap parameters are cited by §3.5 and by the Figure 4 and 5 captions), and §3.3 yields two borderline candidates and no clean cut, because the enumeration of what was searched in DRAFTS is what makes the novelty claim checkable.

---

## TASK 3 — Build

Built with `make pdf`. **The gate's invocation was not used**: it adds `-V "header-includes:..."`, which replaces the manuscript's own YAML `header-includes` rather than appending to it, dropping `\captionsetup{labelformat=empty}` and restoring LaTeX's automatic `Figure N:` on top of each caption's own label. Gate 21 tested and documented this, and the Makefile carries a comment saying not to add it. Float pinning is already in the manuscript's YAML, so the flag is redundant as well as harmful.

| check | result |
|---|---|
| resource warnings | **0** |
| images | **9** |
| caption order | **1–8 ascending** (267, 295, 327, 444, 480, 512, 552, 595); Figure 4 precedes Figure 5 |
| glyphs | ρ 21, ≤ 3, γ 2, β 2, ≈ 2 — all present, non-zero, unchanged |
| pages | **21** (was 22) |
| §3.3 regime table | renders as a table |
| `[email]` | 0; `gabeykim@stanford.edu` present; `[Zenodo DOI]` present |
| Limitations | 12 items, sequential 1–12 |
| cross-references | 10 `Section 3.x`/`2.x` refs, **all resolve**; no `Limitations item N` references exist |
| both copies | byte-identical, md5 `75b2c27830b57df95b3c79c4483dc994` |
| `verify-citations` | 11 / 0 / 5 / 0 / 16 — no regression |

---

## Verdict

Limitations consolidated to 12 items with no number lost from the manuscript, two clauses relocated rather than deleted, and the build clean. The cut scan is a report and changes nothing. The compressed draft is still not on disk, and when it lands this gate's Limitations block can be dropped into it unchanged — but the Task 2 scan will need redoing against it, since every word count here is measured against the uncompressed text.
