**Are the supplements untracked, and does the Zenodo archive contain them? Yes to both. Supplements 3, 4 and 5 are untracked, explicitly gitignored, and still present on disk; 6 and 7 remain tracked. The Zenodo archive does contain all five — record 23051835 is a single 134.59 MB GitHub snapshot of commit `2079b29`, and its central directory lists `gabeykim-crosshost-2079b29/raw/NIHMS945382-supplement-3/-4/-5.xlsx`. No takedown attempted. Two further findings below, one of which contradicts the licence this project just adopted.**

# GATE 34 MEMO — Untrack the publisher supplements

Both audits clean at start and end. **No manuscript change**: md5 `6e5bce0c3ea0dd27b99c0a68e1b4aac6`, identical to Gate 33. No number changed.

---

## TASK 1 — Untracking

Step 1 was already staged, so this gate began at step 2.

**State now:**

| file | tracked | gitignored | on disk |
|---|---|---|---|
| `supplement-3.xlsx` | no | **yes** | yes, 8.4 MB |
| `supplement-4.xlsx` | no | **yes** | yes, 5.4 MB |
| `supplement-5.xlsx` | no | **yes** | yes, 7.6 MB |
| `supplement-6.xlsx` | **yes** | no | yes, 46 KB |
| `supplement-7.xlsx` | **yes** | no | yes, 14 KB |

**`.gitignore`** — explicit entries added for 3, 4 and 5, and the comment rewritten so the stated policy matches what is tracked. It previously said the "two smallest xlsx supplements" are kept while five were tracked; it now says exactly that and names the three exclusions with their combined size. It also records that the three were tracked until this commit, that they remain in earlier commits and in the v1.0.0 Zenodo archive, and that history was deliberately not rewritten because that would break the DOI correspondence.

**`README.md`** — new section, *Publisher material not redistributed*, before Reproducibility. It names the three files with sizes, gives both retrieval routes (PMC6065261 and doi:10.1038/nmeth.4633), states that only `make reproduce-full` needs them, and names the three scripts that read them: `scripts/01_inspect_supplements.py`, `scripts/02_build_core_tables.py` (3 and 4) and `scripts/72_attenuation_analysis.py` (5).

**Verified locally:** all five are on disk and readable — `pd.read_excel` on supplement-5's "Robust RSs" sheet succeeds, which is what `scripts/72` does.

**`make reproduce` still runs: exit 0.** Zero references to the supplements anywhere in its log, and no tracked output changed. That was predictable from the dependency graph — none of the 23 scripts in the fast path is among the four that read these files — but it was run rather than argued.

One correction to the `.gitignore`'s old wording: it said to obtain the files "before running `make reproduce-full`'s Gate-1/2 stage". That is right, and `scripts/02_build_core_tables.py` — the script that most needs them — **is not invoked by any Makefile target at all**, including `reproduce-full`. It is a from-raw bootstrap run by hand. The README says which scripts read the files rather than implying a target covers them.

---

## The Zenodo archive does contain them

- Concept DOI `10.5281/zenodo.23051834` resolves to **record 23051835**, "gabeykim/crosshost: CROSSHOST benchmark v1.0.0", published **2026-09-30**.
- One file: `gabeykim/crosshost-v1.0.0.zip`, **134.59 MB**, a GitHub release snapshot of commit **`2079b29`**.
- Confirmed by an HTTP range request for the last 3 MB of the zip and reading its central directory — no 134 MB download. All five entries are present:
  `gabeykim-crosshost-2079b29/raw/NIHMS945382-supplement-3.xlsx`, `-4`, `-5`, `-6`, `-7`.

`2079b29` predates this gate, so the archive is a faithful snapshot of a repository state in which the three were tracked. **No takedown attempted**, as instructed.

### Finding: the Zenodo record's licence is CC-BY-4.0

The deposit's metadata declares **`cc-by-4.0`**. The repository's root `LICENSE` is MIT, and the manuscript now says "Code and benchmark artifacts are released under the MIT license". **Three statements, two licences.** Zenodo defaults to CC-BY-4.0 when a GitHub release has no detected licence — and at `2079b29` there was no root LICENSE, which is exactly the gap Gate 33 closed. A new release would pick up MIT.

This needs a decision, and it is adjacent to the takedown question: any new deposit version would both carry the correct licence and omit the three supplements.

### Finding: the full article PDF is also tracked

The `.gitignore` policy names three kept items — "bodytext.txt, main-text figure thumbnails, the two smallest xlsx supplements". Two tracked files are not in that list:

| file | size | in the stated policy |
|---|---|---|
| `raw/nihms945382.pdf` | 1.15 MB | **no** — the full PMC author manuscript |
| `raw/nihms945382.nxml` | 108 KB | **no** — its full text in XML |

This is the same class as the supplements and the full article PDF is the most directly re-hostable item of all. I did not untrack them: the instruction named three files, and whether the PMC author manuscript is redistributable is a different question from the Springer Nature supplementary tables — NIHMS deposits are public under the NIH Public Access Policy, which is not the same as a redistribution licence.

**For contrast, the DRAFTS files are fine.** `raw/drafts/msb198875_*` is *Molecular Systems Biology*, fully open access, and the manuscript already cites it as "the PMC6692573 open-access package". Nothing to do there.

---

## TASK 2 — Copyright line [APPLIED]

`package/LICENSE`: `Copyright (c) 2026 CROSSHOST contributors` → `Copyright (c) 2026 Gabriel Kim`. Both LICENSE files now carry the identical line, verified by string comparison.

---

## TASK 3 — Verification

| check | result |
|---|---|
| supplements 3/4/5 untracked, gitignored, on disk | **yes** |
| supplements 6/7 still tracked | **yes** |
| `.gitignore` comment matches tracked contents | **yes** |
| README documents retrieval | **yes**, both PMC and DOI routes |
| both LICENSE files name the same holder | **yes**, byte-identical copyright lines |
| `make reproduce` | **exit 0**, no supplement references, no tracked output changed |
| `make test` | **9 passed** |
| `audit_leakage.py` / `audit_provenance.py` | ALL 6 PASSED / 0 orphans, 177 files |
| `make verify-citations` | **11 / 0 / 5 / 0 / 16 — full baseline.** An earlier run in this gate showed 2 skips; a clean re-run returned zero, confirming transient network as in Gates 29–32. |
| PDF newer than both manuscript copies | **yes** (16:39:52 > 16:39:38, unchanged from Gate 33) |
| both copies byte-identical | **yes**, md5 `6e5bce0c3ea0dd27b99c0a68e1b4aac6` — unchanged |
| pages / images / captions | 22 / 9 / ascending |
| glyphs | ρ 17, ≤ 3, γ 2, β 2, ≈ 1 |
| regime table / "BS 15,848" / `[email]` | renders / intact / 0 |
| `Section X.Y` references | 10, all resolve |
| abstract / body | 250 / 6,950 words — both unchanged |

---

## Verdict

Ongoing distribution of the three supplements is stopped, the stated policy now matches the tracked contents, and a fresh clone is told how to obtain them. The Zenodo archive still contains them and still declares CC-BY-4.0 against the repository's MIT — both are reported, neither acted on.
