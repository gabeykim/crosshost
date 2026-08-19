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

(populated below as each fix is made and re-verified)
