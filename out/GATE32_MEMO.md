**Is the repository submission-ready, and did anything fail verification? Nothing failed verification — every check in Task 3 passes, including the two Gate 31 added (exit status 0, PDF newer than both manuscript copies). The repository is submission-ready with one real gap: there is no LICENSE file at the repository root, and the manuscript states no license for the benchmark. Everything else is present and current.**

# GATE 32 MEMO — Cross-reference, verification-failure note, and submission readiness

Script: `scripts/116_gate32_crossref.py`. Both audits clean at start and end. Abstract unchanged at 250; body 6,828 → 6,834.

---

## TASK 1 — The §3.7 cross-reference [APPLIED]

**Before:**
> …the prediction is consistent with the result, not confirmed by it.

**After:**
> …the prediction is consistent with the result, not confirmed by it **(Section 3.7 measures what that costs)**.

Confirmed present in the rendered PDF. One note on how that was checked: a single-line `grep` of the extracted text returned **zero matches**, because pdftotext wraps the phrase across lines 429–430. Unwrapping the text first finds it. The naive check would have reported the edit missing from a PDF that contains it — the mirror image of the Gate 31 failure, and worth recording for the same reason.

---

## TASK 2 — Verification-failure note [WORDING PROPOSED, NOT APPLIED]

Per instruction, proposed here rather than applied.

**Preferred, one sentence, for Data and Code Availability immediately after the provenance-audit paragraph:**

> Build and regeneration steps should be checked on exit status and on whether the output is newer than its input, rather than on reported success, because a step that silently does nothing leaves its previous output in place and every downstream check will then validate a stale artifact without error.

**Alternative, two sentences, if you want the instances named** — more convincing to a sceptical reader, and more self-incriminating:

> Build and regeneration steps should be checked on exit status and on whether the output is newer than its input, rather than on reported success, because a step that silently does nothing leaves its previous output in place and every downstream check will then validate a stale artifact without error. This occurred three times during development: a figure-regeneration pass in which five of six producers never ran, the provenance false pass described above, and a document build that failed on a full disk while the checks that followed validated the previous build.

I would take the **two-sentence version**. The one-sentence version states a rule; the two-sentence version gives the evidence that the rule was earned, which is the standard the rest of this paper's disclosures are held to. It is also the only place a reader learns the figure-regeneration failure happened at all.

---

## TASK 3 — Final verification

| check | result |
|---|---|
| `make pdf` exit status | **0** |
| **PDF newer than `out/PREPRINT/MANUSCRIPT.md`** | **YES** (20:25:51 > 20:25:38) |
| **PDF newer than `out/MANUSCRIPT.md`** | **YES** |
| new cross-reference present in the PDF | **yes** (across a line break) |
| `audit_leakage.py` | ALL 6 CHECKS PASSED |
| `audit_provenance.py` | 0 orphans, 177 files |
| `make verify-citations` | **11 / 0 / 5 / 0 / 16 — full baseline, zero skips** |
| `make test` | **9 passed** |
| resource warnings | **0** |
| images | **9** |
| captions | Figures 1, 2a, 2b, 3, 4, 5, 6, 7, 8 — ascending |
| glyphs | ρ 17, ≤ 3, γ 2, β 2, ≈ 1 — all present |
| pages | **22** |
| §3.3 regime table | renders |
| "BS 15,848" | intact |
| `[email]` | 0 |
| Zenodo DOI in PDF | present, and `https://doi.org/10.5281/zenodo.23051834` resolves (302 → zenodo.org) |
| GitHub URL | resolves (200) |
| `Section X.Y` references | 10, **all resolve** |
| `Limitations item N` references | 1, resolves (12 items) |
| both copies byte-identical | **yes**, md5 `a4c591a2e0c560c3bc3efd2b3a69e751` |
| abstract | **250 words** |
| body | **6,834 words** |

**Numbers: none changed.** The cross-reference adds one `3.7`, a section pointer rather than a measured value; the guard refuses to write on any other delta.

---

## TASK 4 — Submission readiness

### Present and current

| item | where | state |
|---|---|---|
| Manuscript | `out/PREPRINT/MANUSCRIPT.md`, `out/MANUSCRIPT.md` | byte-identical, 6,834 body words, 250-word abstract |
| PDF | `manuscript.pdf` | 22 pages, tracked, built from the current manuscript with exit 0 |
| Main figures | `out/PREPRINT/figures/` | 9 files for 8 figures (Figure 2 is two panels), all with producing scripts since Gate 25 |
| Supplementary figures | `out/PREPRINT/figures/supplementary/` | 19 files, Figures S1–S11 |
| Figure legends | `out/PREPRINT/FIGURE_LEGENDS.md` | every caption states N and interval type |
| Frozen splits | `data/splits/` | `fold_assignment_FINAL.parquet` plus clusters at 30/50/70% identity |
| Evaluation code | `package/crosshost/evaluate.py` | covered by the test suite |
| Baselines | `package/baselines/` | `master_baselines.{csv,json}`, nine systems |
| Audits | `scripts/audit_leakage.py`, `audit_provenance.py`, `99_gate15_regime_audit.py`, `97_verify_citations.py` | all clean |
| Data manifest | `data/MANIFEST.json` | 26 files hash-verified, 0 mismatches |
| Reproduction | `Makefile` | `help`, `audit`, `verify-citations`, `pdf`/`manuscript`, `test`, `reproduce`, `reproduce-full`, `package`, `clean` |
| Repository README | `README.md` | added Gate 14 after a clean-clone test found the gap |
| DOI | `10.5281/zenodo.23051834` | in the manuscript; resolves |
| Tests | `package/tests/test_package.py` | 9 pass, including two leakage checks |

### Flagged

1. **No LICENSE at the repository root.** `package/LICENSE` is MIT and covers the curated distribution, and `package/data/licensed/` carries per-dataset terms for the DNABERT-2 and PromoGen2 derivatives. But a reader cloning `github.com/gabeykim/crosshost` sees no top-level licence, which defaults to all rights reserved — the opposite of what the paper's availability statement implies. **This is the one thing I would fix before submission**, and it is a one-file change, but choosing the licence is yours.

2. **The manuscript states no licence for the benchmark.** Data and Code Availability says the data are public and that derived activity values are redistributed with attribution, but names no licence for the code or the benchmark itself. A reader has to infer it from `package/pyproject.toml`.

3. **Task 2's sentence is not applied**, pending your choice of wording. One line.

4. **Disk at 99%, 3.9 GB free** — down from 5.8 GB at Gate 31, where a build failed on exactly this. The next build is one bad moment from the same failure, and it will not announce itself.

5. **One open analysis, deliberately**: recomputing Yim et al.'s phylogenetic gradient under this paper's restriction regimes. §3.3 states it as the open check and says it could strengthen the objection. Flagged v2.

6. `out.zip`, carried as untracked and stale since Gate 19, no longer exists.

### Not flagged, but worth knowing

The three self-caught verification failures are each disclosed somewhere: the provenance false pass in Data and Code Availability, the figure regeneration in `out/GATE25_MEMO.md` and `out/state.json`, and the stale-PDF build in `out/GATE31_MEMO.md`. Only the first is in the manuscript. Task 2's two-sentence variant would put all three in front of a reader who reruns the pipeline.

---

## Verdict

Nothing failed verification, including the two checks Gate 31 added. One real gap — no root LICENSE — plus one unapplied sentence awaiting your wording choice. Everything else is present, current, and verified.
