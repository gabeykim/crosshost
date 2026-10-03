**Is the repository submission-ready with no remaining gaps? The licence gap is closed and everything in Task 3 passes. One gap remains, and it is not the one this gate was opened to fix: the repository tracks five Johns et al. supplementary tables under `raw/`, three of them large, contrary to the project's own `.gitignore` policy, which says only "the two smallest xlsx supplements" are kept. That is Springer Nature copyrighted material being re-hosted, and removing it from a public repository is a decision I am not making for you.**

# GATE 33 MEMO — License and the verification-failure note

Script: `scripts/117_gate33_license_and_note.py`. Both audits clean at start and end. Abstract unchanged at 250; body 6,834 → 6,950.

---

## TASK 1 — Root LICENSE [APPLIED] and the contradiction check

### The stop condition was checked and not met

The instruction was to stop if anything in `package/data/licensed/` is contradicted by a root MIT licence. It is not, and the reason is worth stating precisely: **the restricted payloads in that directory are gitignored.** `.gitignore` excludes `package/data/licensed/*/*.npz`, so the PromoGen2 embeddings (CC-BY-NC-4.0) and DNABERT-2 embeddings (Apache-2.0) are not distributed by the repository at all. Only their `LICENSE` and `README` files are tracked. Nothing restricted ships from that directory, so nothing there is contradicted. Proceeded.

### The LICENSE as written

Canonical MIT, `Copyright (c) 2026 Gabriel Kim`, plus a short scope note. Verified:

- filename `LICENSE`, no extension
- **the MIT body is byte-identical to the canonical template** (checked by string comparison, not by eye)
- the MIT body is 71% of the file, so GitHub's detector should classify it as MIT rather than "Other"
- the MIT text matches `package/LICENSE` exactly apart from the copyright line

**One intentional difference from `package/LICENSE`:** it says `Copyright (c) 2026 CROSSHOST contributors`; the root says `Copyright (c) 2026 Gabriel Kim`, as instructed. Flagging it because the gate also asked to match `package/LICENSE`, and these two now differ. Either is defensible; they should probably agree before submission.

**The scope note is not decoration.** A bare MIT at the root would claim rights over material the project does not hold — the `raw/` publisher material below, and the quarantined model-derived data. The note names what is and is not covered and points at `package/LICENSE` and `package/data/licensed/` for the per-component breakdown.

### Manuscript clause [APPLIED]

Appended to the availability paragraph:

> …PromoGen2 and DNABERT-2 embeddings were computed from publicly released model weights and are not redistributed. **Code and benchmark artifacts are released under the MIT license; source datasets retain their original terms, documented per dataset in the repository.**

---

## TASK 2 — Verification-failure note [APPLIED, two-sentence version]

Appended to the provenance-disclosure paragraph in Data and Code Availability:

> …Producing scripts were subsequently added and the audit now covers them genuinely. **More generally, build and regeneration steps should be checked on exit status and on whether the output is newer than its input, rather than on reported success, because a step that silently does nothing leaves its previous output in place and every downstream check will then validate a stale artifact without error. This occurred three times during development: a figure-regeneration pass in which five of six producers never ran, the provenance false pass described above, and a document build that failed on a full disk while the checks that followed validated the previous build.**

No gate numbers, matching the existing disclosure's register. All three instances confirmed present in the rendered PDF.

---

## The gap this gate found: `raw/` re-hosts publisher material against the project's own policy

`.gitignore` states the policy in its own words:

> "Johns et al. 2018 supplementary materials — Springer Nature / Nature Methods copyright (see crosshost-PROJECT-CHARTER.md Part V, 'do not re-host raw tables without checking'). **The smaller derived/verification artifacts (bodytext.txt, main-text figure thumbnails, the two smallest xlsx supplements) are kept** for provenance-audit purposes; the largest and most directly re-hostable items … are excluded."

**Five supplements are tracked, not two:**

| file | size | within policy? |
|---|---|---|
| `raw/NIHMS945382-supplement-3.xlsx` | 8.4 MB | **no** |
| `raw/NIHMS945382-supplement-4.xlsx` | 5.4 MB | **no** |
| `raw/NIHMS945382-supplement-5.xlsx` | 7.6 MB | **no** |
| `raw/NIHMS945382-supplement-6.xlsx` | 46 KB | yes — one of the two smallest |
| `raw/NIHMS945382-supplement-7.xlsx` | 14 KB | yes — one of the two smallest |

So **21.4 MB of Springer Nature copyrighted supplementary tables are being re-hosted**, where the project's stated intent was to keep only the two smallest (60 KB combined) for provenance checks. `.gitignore` does exclude supplement-2 and the PMC tarball, so the policy was applied — just not completely.

**I did not remove them.** Deleting tracked files from a public repository is a judgement with consequences — the bytes stay in git history unless the history is rewritten, and a rewrite breaks the Zenodo archive's correspondence to the repository. The options, briefly:

1. **Leave them and change the policy text** — simplest, but the charter's "do not re-host raw tables without checking" then needs to have been checked.
2. **Remove from the working tree and gitignore them** — stops further distribution; the bytes remain in history.
3. **Remove from history** — complete, and invalidates the DOI's snapshot correspondence.

This is pre-existing, not caused by this gate. The new root LICENSE does not make it worse — its scope note explicitly disclaims `raw/` — but it is the reason the scope note exists.

---

## TASK 3 — Final verification

| check | result |
|---|---|
| `make pdf` exit status | **0** |
| **PDF newer than `out/PREPRINT/MANUSCRIPT.md`** | **YES** (16:39:52 > 16:39:38) |
| **PDF newer than `out/MANUSCRIPT.md`** | **YES** |
| all four new clauses present in the PDF | **yes**, confirmed by unwrapped grep |
| root `LICENSE` exists, correct filename, no extension | **yes**, 1,513 B |
| MIT body byte-identical to canonical template | **yes** |
| MIT body share of file (detectability) | **71%** |
| `audit_leakage.py` | ALL 6 CHECKS PASSED |
| `audit_provenance.py` | 0 orphans, 177 files |
| `make verify-citations` | **11 / 0 / 5 / 0 / 16 — full baseline, zero skips** |
| `make test` | **9 passed** |
| resource warnings / images | 0 / **9** |
| captions | Figures 1, 2a, 2b, 3, 4, 5, 6, 7, 8 — ascending |
| glyphs | ρ 17, ≤ 3, γ 2, β 2, ≈ 1 |
| pages | **22** |
| regime table / "BS 15,848" | renders / intact |
| `[email]` | 0 |
| Zenodo DOI | present; resolves (302 → zenodo.org) |
| GitHub URL | resolves (200) |
| `Section X.Y` references | 10, all resolve |
| `Limitations item N` references | 1, resolves (12 items) |
| both copies byte-identical | **yes**, md5 `6e5bce0c3ea0dd27b99c0a68e1b4aac6` |
| abstract | **250 words** |
| body | **6,950 words** |

**No measured value changed.** Both new clauses are digit-free; the guard compared numeric-token counts in both directions and found them identical, and separately asserts the abstract is untouched.

Disk recovered to 5.1 GB free (98%), up from 3.9 GB at Gate 32.

---

## Verdict

Licence gap closed: root `LICENSE` present and detectable, named in the manuscript, scoped so it does not overclaim. Verification-failure note applied in the two-sentence form. One pre-existing gap surfaced and left for your decision — five tracked Johns supplements against a stated policy of two.
