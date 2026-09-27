**Was the manuscript current, and did all six tasks complete? It was NOT current when the gate opened — both tracked copies were still the Gate 14 version, and the draft had been saved to a new path (`out/PREPRINT/MANUSCRIPT_v8.md`) rather than over them. That was reported and the gate stopped before any edit, per Task 0. After confirmation, the draft was synced to both canonical paths as a separate commit, and all six tasks then completed. No number changed except the five Task 1 additions and Task 4's N range; `verify-citations` did not regress (11/0/5/0/16, identical to Gates 13–14).**

# GATE 19 MEMO — Manuscript sync and pre-submission corrections

Two commits, as specified: `d8fb600` (sync as received, no edits) and this gate's corrections commit. Both audits clean at start and end.

---

## TASK 0 — Manuscript currency

**Initial state: stale.** Both `out/MANUSCRIPT.md` and `out/PREPRINT/MANUSCRIPT.md` were byte-identical to commit `bcc660f` (Gate 14), `git diff HEAD` empty, on-disk mtime Aug 19 16:19 — five weeks old. None of the Task 0 markers were present: old title, no §2.2, no §3.3, 8 results subsections rather than 9, 12 Limitations items rather than 19, and 0.616 / 0.677 / 93.6 all missing from §3.2.

**The draft was located at `out/PREPRINT/MANUSCRIPT_v8.md`** (Sep 26 21:59, 46,449 bytes, untracked), and passed every marker: "evidence that" title, `### 2.2 Two dataset collisions, and the restriction regimes`, 9 results subsections with 3.3 = "Scope and supporting checks", `## 5. Limitations and unresolved gaps` with 19 items already sequential, and all six §3.2 figures. It also already carried the inline (Figure N) callouts, so one file serves both copies.

**Resolution:** copied to both canonical paths, `_v8` deleted, committed as `d8fb600` with no edits, so the corrections land as a visible diff. The repository now holds one manuscript, not three.

---

## TASKS 1–5 — Corrections applied

**Task 1 — regime description.** Replaced only the first sentence, per the corrected reading (the original first sentence was a single semicolon-joined sentence covering both variants; treating it as "two sentences" would have deleted the "intermediate quantile cuts" sentence the instruction says to keep). New text adds 61.11%, 27.98%, 25.3%, 26.04%, 17.4% and the DNA-adequate population note. Everything from "The intermediate quantile cuts…" onward is unchanged.

**Task 2 — retention scoped to transcription.** §3.2 now reads "…against 26.0% (3,668 of 14,088) in vivo **(transcription)**". No retention column or cross-readout table was added anywhere, per instruction — and for the reason Gate 18 established: `scripts/80` builds translation's pooled set already floor-corrected, so EC–BS translation reads 866/866 = 100% as a definitional artifact.

**Task 3 — Figure S5 collision resolved.** Old label **Figure S5** → new label **Figure S11** (genomic feature-group ablation). File `FigureS5_feature_ablation.png` → `FigureS11_feature_ablation.png` via `git mv`. The S5 slot is now an explicit placeholder in `FIGURE_LEGENDS.md` explaining that it is left empty deliberately because §3.2 cites DRAFTS's own *Appendix Fig. S5*, with a pointer to S11. Renumbering within S1–S10 could not have fixed this — any figure of ours occupying the S5 slot collides equally — so vacating the slot is the minimal unambiguous fix. **References checked: exactly one existed** (the legends entry) plus the filename; the manuscript's only "Appendix Fig. S5" is correctly attributed to DRAFTS ("their Appendix Fig. S5") and needed no change. No dangling reference remains.

**Task 4 — Figure 8 caption.** Added "N per cell ranges 294–9,609 depending on host pair and readout (e.g. EC→BS transcription n=3,589; EC→PA transcription n=9,609; BS→PA translation n=294…)" alongside the existing interval statement. All four figures verified against `gate8_5_shift_prediction_summary.csv`. All eight main captions now state both N and interval type.

**Task 5 — disclosures.** §2.7 AI assistance added verbatim at the end of Section 2 (after §2.6, before §3). Acknowledgements placeholder replaced verbatim. Funding, Competing Interests and Author Contributions untouched.

---

## TASK 6 — Consistency pass

**6.1 No number changed that should not have.** Diffing the corrections against the as-received sync commit, the set of numbers added but not re-added-after-removal is exactly: **17.4%, 25.3%, 26.04%, 27.98%, 61.11%** (Task 1) plus **2.5** and **2.7** (the new §2.7 heading and its Section 2.5 cross-reference). Task 4's N range lives in `FIGURE_LEGENDS.md` and verified exactly against source. Total manuscript diff: 7 insertions, 3 deletions.

**6.1 Reconciliation against `out/results/`: 48 numbers checked, 0 contradictions.** Four values flagged as "not found" were each confirmed to be **absences, not errors** — numbers the draft chooses not to cite:

| flagged | status |
|---|---|
| n = 3,826 / 866 / 314 (translation co-active) | the draft cites translation ρ (0.161, 0.742, 0.263 — all three verified correct) but not their n |
| 0.597 (cell-free EC–BS, 82-seq overlap) | the draft uses the larger full-shared-set 0.677 as its headline instead, and reports the 15-sequence overlap (ρ = 0.147) in Limitations item 5 — a defensible strengthening, not an error |

Also noted, not changed: GC-controlled **translation** values (0.125 / 0.738 / 0.230) appear nowhere; the draft reports only transcription GC-controlled figures plus both gap ratios.

**6.2 Cross-references: all resolve.** Sections referenced — 2.2, 2.5, 3.1, 3.2, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9 — all exist among the 21 headings present. Zero unresolved. No "Limitations item N" references exist in the text, so none can dangle. (The new §2.7's "audit suite described in Section 2.5" was checked specifically: §2.5 Splitting does describe the standing leakage audit, so it resolves.)

**6.3 Limitations: 19 items, sequential 1–19, no duplicates, no gaps.**

**6.4 Both copies byte-identical** — `cmp` clean, md5 `7f57d6260fec26cfd8e531f2a4178e53` on both.

**6.5 `make verify-citations`: 11 VERIFIED, 0 unexplained mismatches, 5 manually-reviewed-and-explained, 0 unresolved, 16 total** — identical to Gates 13 and 14. No regression.

**Additional check (stray drafts).** Scanned every `.md` under `out/` for a full manuscript (Abstract + Results + References together), excluding the two canonical copies: **none found.** The repository will archive to Zenodo with exactly one manuscript in two identical locations.

---

## Verdict

Manuscript synced and all six tasks complete. One instruction correction applied with prior confirmation (Task 1's sentence count). Zero unintended number changes, zero unresolved cross-references, zero citation regressions, zero stray drafts.
