# CROSSHOST

A benchmark for predicting bacterial regulatory DNA activity across host species. See [`package/README.md`](./package/README.md) for the full project description, central finding, licensing quarantine, and quickstart — this file exists only to get a fresh clone of the *repository* (not the curated `package/` distribution) working, since nothing at the repository root pointed to that file before Gate 14 found the gap.

**Repository layout, in brief:** `scripts/` and `data/` are where the actual pipeline lives and runs; `package/` is a curated *export* of selected outputs, built by `make package`, not a self-contained pipeline in its own right (see the `Makefile`'s own header comment). Most people cloning this repo want to run `make audit` / `make reproduce`, both of which run from the repository root, not from inside `package/`.

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e ./package                # NOT `pip install -e .` -- there is no pyproject.toml at the repo root
pip install -e "./package[citations]"    # + certifi, needed only for `make verify-citations`
```

`requires-python = ">=3.10"` (declared in `package/pyproject.toml`). Verified working on Python 3.14.2, macOS/Apple Silicon, as of Gate 14 (`out/results/gate14_clean_clone_verification.md`) — no other Python version or OS has been tested.

## Checks

```bash
make audit             # leakage + provenance audits -- fast, offline
make verify-citations   # mechanically verify every manuscript citation -- live network (CrossRef/arXiv)
make test               # package test suite
make reproduce          # regenerate all figures/tables from shipped intermediate results -- fast
make reproduce-full      # full pipeline from raw/ -- MANY HOURS, not re-tested end-to-end, see Makefile header
```

Run `make help` for the same list with one-line descriptions.

## Publisher material not redistributed

Five Johns et al. (2018) files are **not** included in this repository, because they are Springer Nature / Nature Methods copyright and this project does not hold redistribution rights. Availability on PMC under the NIH Public Access Policy is not a redistribution licence.

| file | size | what it is |
|---|---|---|
| `raw/nihms945382.pdf` | 1.2 MB | the full author manuscript |
| `raw/nihms945382.nxml` | 108 KB | its full text in XML |
| `raw/NIHMS945382-supplement-3.xlsx` | 8.4 MB | supplementary table |
| `raw/NIHMS945382-supplement-4.xlsx` | 5.4 MB | supplementary table |
| `raw/NIHMS945382-supplement-5.xlsx` | 7.6 MB | supplementary table |

Download them from the publisher and place them in `raw/` under exactly those filenames:

- PMC open-access package: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6065261/> (the supplementary files are listed on that page)
- or Nature Methods: <https://doi.org/10.1038/nmeth.4633>

**You only need the three supplements for `make reproduce-full`.** `make audit`, `make test`, `make verify-citations` and `make reproduce` (the fast path) all run from shipped intermediate results and read none of these files. The scripts that do are `scripts/01_inspect_supplements.py`, `scripts/02_build_core_tables.py` (supplements 3 and 4) and `scripts/72_attenuation_analysis.py` (supplement 5). **No script reads the pdf or the nxml** — they were present for reference only.

Kept for provenance-audit purposes: `nihms945382_bodytext.txt` (the extracted body text, used as a unit-test fixture by `scripts/03` and `scripts/15`), the five main-text figure thumbnails, and the two smallest supplements (`-6`, 46 KB and `-7`, 14 KB).

## Reproducibility

Verified from a genuine fresh `git clone` (not an existing working directory, not an isolated venv reusing a checkout already on disk) as of Gate 14: `make audit`, `make verify-citations`, and `make reproduce` all pass with no undocumented step beyond what is written in this file. Full verification record, including every gap found and fixed to get there: `out/results/gate14_clean_clone_verification.md`. `make reproduce-full` was not attempted (would cost the ~40 cumulative hours of training this project's Gates 4–8 already spent) — see `out/PREPRINT/REPRODUCIBILITY.md` for what that does and does not mean.

## The manuscript

`out/MANUSCRIPT.md` (working copy) and `out/PREPRINT/MANUSCRIPT.md` (submission-formatted, with figure numbering) are the same content. Start there for the actual scientific writeup; this README is setup instructions only.
