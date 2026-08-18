# Archive metadata — prepared, nothing published

Per Gate 8 Task 2f: artifacts prepared for review. **Nothing has been submitted or published to Zenodo or Hugging Face.** Everything below is a draft for the maintainer to check, edit, and actually submit.

**Gate 8.5 update:** the source repository now exists and is pushed —
`https://github.com/gabeykim/crosshost` (**currently private**; the
maintainer has not authorized making it public). The Zenodo
single-vs-split-deposit question flagged below as "DECISION NEEDED FROM
MAINTAINER" in Gate 8 is now **settled**: **two separate Zenodo deposits**
— one for the MIT-licensed core benchmark (`data/core/`, `crosshost/`,
baselines, held-out eval), one for the CC-BY-NC-4.0 PromoGen2-derived
embeddings (`data/licensed/promogen2_derived/`) — keeping the core
benchmark cleanly licensed for downstream commercial use. Neither has
actually been created yet; this remains prepared-not-submitted.

## Zenodo

Suggested metadata for a Zenodo deposit (upload type: dataset + software, since this package bundles both):

```yaml
title: "CROSSHOST: A Benchmark for Predicting Bacterial Regulatory DNA Activity Across Host Species"
upload_type: dataset
description: >
  A benchmark measuring whether models can predict bacterial regulatory DNA
  activity in host species unseen during training. Includes a frozen
  165bp-sequence library (E. coli, B. subtilis, P. aeruginosa; 3 additional
  hosts at reduced N via RS241), frozen splits, an evaluation API, and
  baseline scores for 9 systems from a mean/majority floor through two
  genomic foundation models. Central finding: host conditioning does not
  measurably improve cross-host prediction over a sequence-only model,
  explained by cross-host measurement correlation that varies sharply with
  host phylogenetic distance. n_hosts <= 6.
creators:
  - name: "[MAINTAINER NAME -- fill in]"
    affiliation: "[fill in]"
license: MIT
  # NOTE (SETTLED, Gate 8.5): this deposit covers data/core/ and crosshost/
  # ONLY. data/licensed/promogen2_derived/ (CC-BY-NC-4.0) ships as a SEPARATE
  # Zenodo deposit -- see below -- so this deposit's MIT license applies
  # cleanly to everything in it, with no per-file exception needed.
keywords:
  - genomics
  - synthetic biology
  - regulatory DNA
  - bacteria
  - cross-host prediction
  - benchmark
related_identifiers:
  - relation: isSupplementTo
    identifier: "10.1038/s41592-018-0021-8"  # Johns et al. 2018, Nat Methods -- VERIFY this DOI before submitting
communities: []  # e.g. add a "bioinformatics" or "synthetic-biology" community if desired
```

### Second deposit — CC-BY-NC-4.0 PromoGen2-derived embeddings (SETTLED, Gate 8.5)

```yaml
title: "CROSSHOST: PromoGen2 zero-shot embeddings (CC-BY-NC-4.0, non-commercial)"
upload_type: dataset
description: >
  Frozen PromoGen2 (jinyuan22/promogen2-base) embeddings for the CROSSHOST
  three-host library and RS241, used as one of two genomic foundation-model
  baselines in the CROSSHOST benchmark (see the core deposit / repository
  for the benchmark itself). Distributed separately under CC-BY-NC-4.0
  (non-commercial), matching PromoGen2's own license, so the core CROSSHOST
  deposit can remain MIT-licensed and commercially reusable without
  entanglement. Requires the core CROSSHOST deposit's data/splits to be
  useful on its own.
creators:
  - name: "[MAINTAINER NAME -- fill in, same as core deposit]"
    affiliation: "[fill in]"
license: "CC-BY-NC-4.0"
related_identifiers:
  - relation: isPartOf
    identifier: "[core CROSSHOST deposit DOI -- fill in once minted]"
communities: []
```

**What remains for the maintainer:**
1. Fill in creator name/affiliation/ORCID (both deposits).
2. ~~Decide the single-deposit-vs-split-deposit question~~ — **done**: two deposits, per above.
3. Verify the Johns et al. 2018 DOI is exactly right before citing it as `related_identifiers`.
4. Actually create the Zenodo deposit, upload `package/` (or a zip of it), and mint the DOI.
5. Once minted, add the DOI badge to `README.md` and cite it in any future paper.

## Hugging Face Datasets

Suggested dataset card front-matter (`package/README.md` already has the human-readable content; HF wants YAML front-matter at the top of a dataset card, which would need to be added to a copy intended specifically for HF upload):

```yaml
---
license: mit
license_details: >
  Core benchmark (data/core/) is MIT with attribution to Johns et al. 2018.
  data/licensed/promogen2_derived/ is CC-BY-NC-4.0 (non-commercial) --
  see that directory's own LICENSE file. data/licensed/dnabert2_derived/
  is Apache-2.0.
task_categories:
  - other
tags:
  - genomics
  - bacteria
  - regulatory-dna
  - cross-species
  - benchmark
size_categories:
  - 10K<n<100K
---
```

**What remains for the maintainer:**
1. Create a Hugging Face account/organization if one doesn't exist.
2. Decide whether to upload `data/core/` as one HF dataset and `data/licensed/*` as separate, clearly-labeled datasets (recommended, given the license split) or find another way to represent the quarantine on HF's platform.
3. Actually run `huggingface-cli upload` or use the web UI — not done here.

## What was verified vs. not, for both

- Verified: the package's own internal structure, licensing quarantine, and data integrity (via `make audit`, `make test`).
- Not verified: whether the exact YAML schemas above validate against Zenodo's/HF's current (2026) metadata requirements — these platforms update their schemas periodically; check the current docs before submitting rather than trusting this file blindly.
