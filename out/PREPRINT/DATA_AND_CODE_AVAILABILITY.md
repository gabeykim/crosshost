# Data and Code Availability

## Repository

Code, frozen splits, evaluation API, and derived data tables: `https://github.com/gabeykim/crosshost`. **Currently private.** Making it public is a decision for the maintainer (Gabriel), not made by this session — see `out/SHIP_CHECKLIST.md`. Once public, this line becomes the canonical citation target; until then, this statement describes the plan, not a live link a reader can follow.

## Archive deposits (prepared, not yet created)

Two Zenodo deposits, kept structurally separate so the core benchmark stays cleanly, commercially reusable:

1. **Core benchmark** (MIT license) — `data/core/`, the `crosshost` package, baseline suite, held-out evaluation split. Covers everything needed to reproduce every result in this manuscript except the two foundation-model comparisons (Appendix, Section 3.9).
2. **PromoGen2-derived embeddings** (CC-BY-NC-4.0, non-commercial, matching PromoGen2's own license) — kept out of the core deposit specifically so the core deposit's MIT license applies without exception.

Full metadata drafts (title, description, keywords, related-identifiers) for both: `package/ARCHIVE_METADATA.md`. **Neither has been created; no DOI has been minted.** This statement will need updating with the actual DOIs once that happens — do not cite a DOI from this document as if it exists yet.

## Licensing quarantine (already enforced in the repository structure, not just stated)

| Directory | License | Contents |
|---|---|---|
| `crosshost/` (code) | MIT | Package code |
| `data/core/` | MIT (data), attribution required | Derived activity values, splits, genomic/physiology features |
| `data/licensed/dnabert2_derived/` | Apache-2.0 | DNABERT-2 frozen embeddings |
| `data/licensed/promogen2_derived/` | **CC-BY-NC-4.0, non-commercial** | PromoGen2 frozen embeddings |
| `data/licensed/deepcross/` | N/A (empty) | Pre-provisioned quarantine; never populated |

`data/licensed/*` content is never merged into `data/core/` — the directory boundary is the enforcement mechanism.

## Primary data provenance

- **Johns et al. (2018)**, *Nature Methods* 15:323–329 — raw supplementary tables (Springer Nature copyright) are **not redistributed**; only derived values, with attribution (`data/core/CITATION.md`). BioProject PRJNA431139.
- **Yim, Johns et al. (2019)**, *Molecular Systems Biology* 15:e8875 (DRAFTS) — acquired via NCBI's PMC Open Access bulk-download service (PMC6692573), license **CC BY 4.0** per the article's own stated license (confirmed from the article XML's license tag, not assumed). Every acquired file's URL, size, SHA256, and license is recorded in `data/MANIFEST.json`, independently re-verified against the manifest, not just trusted from the acquiring script's own report. Raw sequencing reads (SRA PRJNA509603) were not downloaded — the processed source-data tables were sufficient for every analysis in this manuscript.

## What is and is not included in the released package

**Included:** frozen genome-blocked splits (5-fold, plus the held-out evaluation fold), derived activity tables, genomic and physiology host-feature vectors, the full baseline suite (9 systems, fold-resolved, 90% bootstrap CI), DRAFTS-derived intermediate tables (`data/drafts.parquet`, `data/drafts_rs234_invivo.parquet`), all figures and their producing scripts, both audit scripts (`audit_leakage.py`, `audit_provenance.py`).

**Not included:** raw Johns et al. supplementary files (copyright; use `data/core/` derived values instead); model checkpoints beyond what `out/models/` already contains locally (not yet decided whether these ship in the Zenodo deposit or are regenerable-only — see `out/SHIP_CHECKLIST.md`); DRAFTS's own raw supplementary files beyond what this project's `data/` pipeline outputs need (the acquired `raw/drafts/*.xlsx`/`.pdf` files are tracked in `data/MANIFEST.json` with hashes but excluded from the git repository itself via `.gitignore` for two large redundant bundles — see Gate 10's standing note).
