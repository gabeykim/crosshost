# Gate 14 — Clean-Clone Verification

**Clone location:** `/tmp/crosshost-verify` (provided pre-cloned; confirmed a genuine `git clone` of `origin` = `https://github.com/gabeykim/crosshost.git`, at commit `299bfdf` — Gate 13's own final commit — `git status` clean, "up to date with origin/main").

**Method:** worked entirely inside the clone, following only what is written in files inside the clone itself — no step taken on the basis of prior knowledge of the project's internal structure. Every command and its result recorded below, in the order run, written incrementally as the verification proceeded (not reconstructed afterward). Nothing in the clone was fixed; every defect found here is fixed only in `~/Downloads/crosshosts`, per instruction.

**Status: IN PROGRESS.**

---

## Task 1 — Environment setup from scratch

### 1.1 Confirming no venv inheritance

```
python3 -m venv .venv && source .venv/bin/activate && which python3
-> /private/tmp/crosshost-verify/.venv/bin/python3
VIRTUAL_ENV=/private/tmp/crosshost-verify/.venv
```

Confirmed: the interpreter resolves inside the clone's own `.venv`, not any pre-existing environment. (Separately noted: this tool's Bash sessions do not persist shell state — including venv activation — between separate tool invocations. Every subsequent command in this report was run with `cd /tmp/crosshost-verify && source .venv/bin/activate && ...` chained explicitly in the same call, to guarantee the venv was actually active for that command, not assumed from a prior one.)

### 1.2 What the repository documents, read verbatim

**No `README.md` exists at the repository root.** `find . -maxdepth 1 -type f` at the root lists only: `Makefile`, `crosshost-PROJECT-CHARTER.md`, `.gitignore`. The only `README.md` in the entire repository is `package/README.md`. **This is the first undocumented-step finding**: the gate's own instruction ("Read the README and Makefile") assumes a root README exists; a genuine stranger cloning this repository has no root-level entry point document at all and must discover `package/README.md` by exploring the tree — nothing at the root points to it directly except the Makefile's own header comment.

**The root `Makefile`'s header comment**, verbatim:
> "CROSSHOST reproduction Makefile. Lives at the project root because scripts/ and data/ (the things being regenerated) live here -- package/ is a curated EXPORT of selected outputs for distribution, built BY this Makefile's `package` target, not a self-contained pipeline in its own right. See package/README.md for the distinction."

This is the only pointer from the root to any setup documentation. Following it to `package/README.md`'s "Install" section, verbatim:
```bash
pip install -e .                 # core package, data loaders, evaluation API
pip install -e ".[models]"       # + torch/transformers, needed only to re-run the FM baselines
pip install -e ".[dev]"          # + pytest, for running tests/
```

**This is itself ambiguous about working directory.** Taken literally from the repository root (where a stranger following the Makefile's pointer would still be sitting, having only *read* `package/README.md`, not `cd`'d into it), `pip install -e .` targets the root — which has no `pyproject.toml` — and fails. The same README's own later sentence ("Tested from a fresh `git clone`... `pip install -e ./package`") uses the correct relative path, but the "Install" code block itself does not, and nothing in the surrounding text instructs `cd package` first. **Second undocumented-step finding.**

**No environment-variable, Python-version, or system-dependency documentation exists anywhere in the root tree.** `package/environment.lock.txt` and `package/pyproject.toml` exist, but the root Makefile never references either, and nothing at the root states a minimum/tested Python version. `python3 --version` in this environment: 3.14.2 (see 1.4 for whether this matters).

### 1.3 Following the documented path exactly

| # | Command | Result |
|---|---|---|
| 1 | `python3 -m venv .venv && source .venv/bin/activate` | OK — interpreter confirmed inside clone (1.1) |
| 2 | `make audit` (no install step attempted first — nothing in the Makefile or its `help` output states one is required before `make audit`) | **FAILED**: `ModuleNotFoundError: No module named 'numpy'`. `make audit`, the Makefile's own first-listed, simplest command, cannot run in a truly fresh clone with no undocumented setup step. |
| 3 | `pip install -e ./package` (inferred from `package/README.md`'s *second* mention of the command, since the "Install" section's own literal `pip install -e .` fails at the repo root — no `pyproject.toml` there) | OK — installed numpy 2.4.6, pandas 3.0.3, scipy 1.18.0, scikit-learn 1.9.0, pyarrow 24.0.0, matplotlib 3.11.0, `crosshost-benchmark` 1.0.0 (editable). No version conflicts against Python 3.14.2 — all pinned versions had prebuilt `cp314-macosx-arm64` wheels. |
| 4 | `make audit` (retry) | **FAILED** — see below, a real deviation, not a setup gap |

**`make audit` result on retry:**
```
[PASS] 1. No RS appears in more than one fold
       29042 rows, 0 duplicate OLIGO IDs
[PASS] 2. No RS241 sequence appears in the main-library fold assignment
       241 RS241 ids checked, 0 found in main-library folds
[PASS] 3. Max train-test identity per fold < 0.9
       per-fold max identities: fold0=0.8485, fold1=0.8395, fold2=0.8485, fold3=0.8485, fold4=0.8364; worst=0.8485
[PASS] 4. No exact-duplicate sequence pair is split across folds
       181 exact-duplicate groups checked, 0 split across folds
[PASS] 5. Fold-local normalization: fitted statistics use training-pool data only
       EC: cached=1.87, train_pool_recompute=1.87, all_data_recompute=1.87; BS: cached=2.00, ...; PA: cached=2.23, ... -- all match training-pool-only fit
[FAIL] 6. Table hashes match data/MANIFEST.json
       24 files checked against manifest, 2 mismatches: ['raw/drafts/msb198875_SourceData_Appendix.zip: file missing', 'raw/drafts/PMC6692573.tar.gz: file missing']
RESULT: AT LEAST ONE CHECK FAILED
```

**Checks 1–5 match the real repo's numbers exactly** (29,042 rows, worst identity 0.8485 — no deviation). **Check 6 does not**: the real repo reports 26/26 files, 0 mismatches; the clone reports 24/26, 2 missing. Root cause identified, not guessed: `raw/drafts/PMC6692573.tar.gz` and `raw/drafts/msb198875_SourceData_Appendix.zip` were deliberately excluded from git via `.gitignore` in Gate 10 (large, redundant bundles — their contents are also available in the other, git-tracked `raw/drafts/*.xlsx`/`.pdf` files) — **but their hashes remain in `data/MANIFEST.json`, and `audit_leakage.py` check 6 has no concept of "expected-absent in a clean checkout."** A stranger cloning the repository and running the audit exactly as documented gets a FAIL, not a PASS, on this check — every time, unconditionally, regardless of anything they do right. This is a real defect in the audit script's design, not a fresh-clone artifact to be worked around.

### 1.4 Undocumented steps required

1. **No environment/dependency setup is documented at the repository root at all.** `make audit` (the Makefile's first-listed command) fails immediately with `ModuleNotFoundError` in a truly fresh clone.
2. **No root-level `README.md` exists.** The only path to setup instructions is the Makefile's own header comment pointing to `package/README.md`.
3. **`package/README.md`'s own "Install" code block (`pip install -e .`) is wrong if followed from the repository root** — it has no `pyproject.toml` there. The correct command, `pip install -e ./package`, appears later in the same file's prose but not in the "Install" section's own code block.

---

## Task 2 — Running the checks (Iteration 1, before any fix)

### `make verify-citations`

**FAILED, all 16 citations.** Every single one: `SSL: CERTIFICATE_VERIFY_FAILED: unable to get local issuer certificate`. Root cause identified directly: `scripts/97_verify_citations.py` imports `certifi` (with a try/except fallback to the system default SSL context) to build its HTTPS request context, exactly because this same certificate problem was hit and fixed during the script's own development in Gate 13 — **but `certifi` was never added as a declared dependency anywhere** (not in `package/pyproject.toml`'s base `dependencies`, not in `dev`, not in `models`). Confirmed directly: `python3 -c "import certifi"` fails in the fresh venv; `pip list | grep certifi` shows nothing. On the original development machine, `certifi` was present incidentally (from some other, unrelated prior install), masking this gap completely — the exact "works on my machine" failure mode this gate exists to catch. Expected result (11 verified, 5 explained mismatches, 0 real errors) — **0 of 16 could even be checked.**

### `make reproduce`

**FAILED at the first step.** `reproduce: audit` in the Makefile — `make reproduce` depends on `make audit` and `make` stops at the first failing prerequisite by default. Since `make audit` fails on check 6 (above), **none of `reproduce`'s 23 figure/table-regeneration scripts ran at all.** The manuscript's reproducibility claim ("`make reproduce`... re-run in place, successfully" — Gate 11) is not reproducible by a stranger from a clean clone today: the chain never gets past the audit gate.

**Summary of Iteration 1: 2 real defects found, both blocking.** (1) `audit_leakage.py` check 6 has no concept of deliberately-`.gitignore`d manifest entries, so it always fails in a clean checkout. (2) `scripts/97_verify_citations.py` depends on an undeclared `certifi` dependency. Both must be fixed before any further comparison against expected values is even possible — `make reproduce` cannot run until (1) is fixed, since it depends on `audit`.

---

## Task 3 — Fixing defects in the real repo

**Three fixes made in `~/Downloads/crosshosts`, committed as `9d6bda5`:**

1. **Added `README.md` at the repository root** — states the actual setup path (`pip install -e ./package`, not `pip install -e .`), lists the `make` targets, points to `package/README.md` for the full project description and to `out/MANUSCRIPT.md`/`out/PREPRINT/MANUSCRIPT.md` for the writeup. This is the fix for undocumented-step findings #1–#3 (Task 1.4) — a documentation fix, no code changed, per the gate's own stated preference.
2. **Fixed `package/README.md`'s "Install" section** to state explicitly which directory its code block assumes, and cross-reference the new root README for the repo-root path. Documentation-only.
3. **Fixed `scripts/audit_leakage.py` check 6** (code change — documentation alone cannot fix this, since the check's own logic, not a missing instruction, was wrong) to read `.gitignore` and treat a gitignored manifest entry's absence as an expected skip rather than a failure, while still hash-verifying it if present. Also added `certifi` as a new `citations` extra in `package/pyproject.toml` (`scripts/97_verify_citations.py`'s own dependency, previously undeclared and masked by an incidental prior install on the original dev machine).

**Deleted `/tmp/crosshost-verify`, cloned fresh again** (`git clone https://github.com/gabeykim/crosshost.git /tmp/crosshost-verify`, confirmed at commit `9d6bda5`) — **Iteration 2:**

```
python3 -m venv .venv && source .venv/bin/activate    # OK, confirmed inside the new clone
pip install -e ./package                                # OK -- as now documented in the new root README
pip install -e "./package[citations]"                    # OK -- installs certifi
make audit                                                # PASSED, all 6 checks
make verify-citations                                     # PASSED -- 11 verified, 5 explained, 0 real errors, 0 unresolved
make reproduce                                            # (see below)
```

**`make audit`, full output, Iteration 2:**
```
[PASS] 1. No RS appears in more than one fold — 29042 rows, 0 duplicate OLIGO IDs
[PASS] 2. No RS241 sequence appears in the main-library fold assignment — 241 RS241 ids checked, 0 found
[PASS] 3. Max train-test identity per fold < 0.9 — worst=0.8485 (fold-by-fold: 0.8485/0.8395/0.8485/0.8485/0.8364)
[PASS] 4. No exact-duplicate sequence pair is split across folds — 181 groups checked, 0 split
[PASS] 5. Fold-local normalization — EC/BS/PA cached values all match training-pool-only fit
[PASS] 6. Table hashes match data/MANIFEST.json — 24 files checked against manifest, 0 mismatches;
       2 gitignored and skipped: [msb198875_SourceData_Appendix.zip, PMC6692573.tar.gz]
RESULT: ALL CHECKS PASSED
```

**Deviation from the real repo's own numbers, reported rather than smoothed over, per instruction:** the real repo reports "26 files checked against manifest, 0 mismatches" (all 26 physically present there); the clone reports "24 files checked... 0 mismatches; 2 gitignored and skipped." **This is an expected, explained difference, not an unexplained one** — 24 + 2 = 26, the same total set of manifest entries, and the 2 that differ are exactly the two files Gate 10 deliberately excluded from git. Every number that is expected to be identical between environments (rows, identities, checked-file count among files that ARE present, 0 mismatches) *is* identical. `audit_provenance.py`: 167 files, 0 orphans — matches exactly.

**`make verify-citations`, Iteration 2: matches Gate 13's own reported values exactly** — 11 verified, 0 unexplained mismatches, 5 manually-reviewed-and-explained mismatches (same 5, same explanations), 0 unresolved, 16 total.

**`make reproduce`, Iteration 2: progressed much further than Iteration 1 (which never got past `audit`), then failed on a fourth defect.** Every script through `scripts/87_gate8_6_figures.py` ran cleanly, and every number printed matches the real repo's established values exactly (spot-checked: `EC transcription concat: mean=0.555 vs seqonly=0.367`, `film_std_over_seqonly_std: 7.800`, the two-stage-transfer table, the shift-prediction summary table, the co-active range-restriction IQR ratios, the fold-variance table — all identical to the manuscript's own stated figures). Then:

```
python3 scripts/88_parse_drafts.py
...
ImportError: `Import openpyxl` failed. Use pip or conda to install the openpyxl package.
make: *** [reproduce] Error 1
```

**Root cause identified, not guessed:** `scripts/88_parse_drafts.py` calls `pd.read_excel(...)` to parse the DRAFTS source-data spreadsheets — part of the core `make reproduce` path since Gate 11 added the Gate-10/10.5 scripts to it. Pandas imports its Excel engine (`openpyxl`) lazily, inside `read_excel()` itself, so it never appears as a top-level `import openpyxl` statement anywhere in the script — a plain grep for import statements across all 23 `reproduce`-path scripts (run before this fix, to check for further hidden gaps in one pass) found nothing, confirming this is a real blind spot in that kind of check, not a scan that was run carelessly. `openpyxl` was never declared in `package/pyproject.toml`. On the original dev machine it was present incidentally (used directly, via `import openpyxl`, during Gates 12–13's own ad-hoc Excel-sheet-structure checks) — masked exactly the same way `certifi` was.

**Checked specifically for further hidden optional-dependency gaps of the same shape** (`grep -l "read_excel\|to_excel\|openpyxl"` across all 23 `reproduce`-target scripts) before fixing: only `scripts/88` uses this pattern. No other script in the fast-reproduce path calls `pd.read_excel` or otherwise pulls in an undeclared optional pandas engine.

**Fixed** (`package/pyproject.toml`: added `openpyxl==3.1.5` to base `dependencies`, not an extra, since `scripts/88` runs unconditionally as part of `make reproduce`), **committed as `87e6b61`, pushed.** With `openpyxl` installed manually as a spot-check (before the commit, to confirm the diagnosis), `make reproduce` completed end to end, exit code 0, and every number in the DRAFTS/GC-confound task sequence (scripts 88–96) matched the manuscript's own established values exactly — spot-checked in full: per-host GC-vs-activity correlations (EC −0.614/−0.494, BS −0.200/−0.266, PA −0.434/−0.230), all 6 raw-vs-GC-controlled pair correlations, the decisive gap-ratio verdicts (transcription 0.341→0.308 SURVIVES INTACT, translation 0.286→0.240 SURVIVES ATTENUATED), and the full phylum-stratified table — all bit-for-bit identical to Gate 10.5's and Gate 11's own reported figures.

**Deleted `/tmp/crosshost-verify`, cloned fresh a third time** (commit `87e6b61`) — **Iteration 3, the clean pass:**

```
python3 -m venv .venv && source .venv/bin/activate           # OK
pip install -e ./package                                      # OK -- openpyxl now installs automatically, no manual step
pip install -e "./package[citations]"                          # OK -- certifi installs automatically
make audit               # PASSED, all 6 checks, identical output to Iteration 2
make verify-citations     # PASSED, 11/0/5/0/16, identical to Iteration 2 and to Gate 13's own report
make reproduce            # PASSED -- exit code 0, 25/25 script invocations completed, zero errors/tracebacks anywhere in the full log
```

**`make reproduce`, Iteration 3: clean, end to end, exit code 0.** Every number checked matches the established values exactly, including the full DRAFTS/GC-confound sequence: `EC-PA raw=0.754 -> GC-controlled=0.718`, `film_std_over_seqonly_std: 7.800181878927151`, `EC transcription concat: mean=0.555 vs seqonly=0.367`, and the overall verdict `{'SURVIVES INTACT', 'SURVIVES, ATTENUATED'}` — all bit-for-bit identical to Gates 10.5/11/13's own reported figures. `grep -iE "error|traceback|failed"` across the entire log: zero matches.

**Iteration 3 is the clean pass. It took 3 iterations to get there:**
- **Iteration 1** (the clone as originally provided, commit `299bfdf`): 4 real defects found — no root README / ambiguous install path, `audit_leakage.py` check 6 unconditionally failing, `certifi` undeclared. Fixed, committed as `9d6bda5`.
- **Iteration 2** (fresh clone at `9d6bda5`): the 3 fixes confirmed working; `make audit` and `make verify-citations` passed; `make reproduce` progressed far further than Iteration 1 (23 scripts, not 0) before failing on a 4th defect, `openpyxl` undeclared. Fixed, committed as `87e6b61`.
- **Iteration 3** (fresh clone at `87e6b61`): `make audit`, `make verify-citations`, and `make reproduce` all passed with zero undocumented steps and zero manual intervention beyond the four commands the new root `README.md` actually documents.

---

## Task 4 — Reproducibility statement, updated

(see `out/PREPRINT/REPRODUCIBILITY.md` and `out/MANUSCRIPT.md` Section 6 — updated to state the Gate 14 clean-clone verification precisely: date, exact targets verified, environment, and the one still-open limitation, `make reproduce-full`, which remains untested end-to-end and is disclosed as such, not silently implied to be covered by this gate's work.)

**Environment used for this verification:** macOS, Apple Silicon (arm64), Python 3.14.2, `pip` 25.3. No other OS or Python version has been tested. All pinned dependency versions in `package/pyproject.toml` (including the two added this gate, `certifi` and `openpyxl==3.1.5`) had prebuilt `cp314-macosx-arm64` wheels available — no compilation-from-source step was needed or tested.

**Remaining limitation, disclosed rather than implied closed:** `make reproduce-full` (the complete from-raw-data pipeline, including all CNN/foundation-model training) was explicitly not attempted, per instruction — it would cost the cumulative ~40 hours of training this project's Gates 4–8 already spent. This gate closes the "was a genuine clean clone ever tested" gap for `make audit`, `make verify-citations`, and `make reproduce` (the fast path) specifically. It does not close, and was never asked to close, the separate `make reproduce-full` gap.
