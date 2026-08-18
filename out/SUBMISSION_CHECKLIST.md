# CROSSHOST — Venue Assessment and Submission Checklist

**Nothing has been submitted, deposited, or made public. This is
preparation and a realistic assessment, not an aspirational plan.**

## Realistic assessment

This is a three-primary-host (plus three at N≈207) benchmark and a
negative result, built from one 2018 dataset, by a solo, unaffiliated
author, using free/local compute. That is a real, useful contribution — a
carefully-controlled negative result with a public benchmark attached is
exactly the kind of thing NeurIPS Datasets & Benchmarks workshops and
*Scientific Data* exist for — but it is not a main-track NeurIPS
Datasets & Benchmarks paper on its own scale or novelty, and treating it
as one would waste a submission cycle chasing a venue this work is not
positioned to win.

## Venue plan

### Primary: MLCB (Machine Learning in Computational Biology) or a NeurIPS/ICLR workshop

**Reachable, appropriate, and fast.** MLCB and workshop tracks (NeurIPS
D&B workshop, ICLR workshop tracks with a computational-biology or
genomics fit) have review cycles measured in weeks, not months, review
solo/small-team work fairly, and are a strong fit for a well-controlled
negative result plus a reusable benchmark artifact — exactly the kind of
paper these venues exist to surface that a main track would deprioritize
for lack of a positive headline result.

**What remains:** confirm current-year submission deadlines and CFP scope
for MLCB and the relevant NeurIPS/ICLR workshops (call-for-papers pages
change year to year and were not checked this session — the manuscript is
ready, the venue-specific formatting and exact deadline are not).

### Secondary: *Scientific Data* (curated-resource descriptor)

Appropriate for the benchmark-as-resource framing specifically — a data
descriptor emphasizing the dataset, splits, evaluation API, and baseline
suite rather than the negative-result narrative as the primary
contribution. Complements, does not compete with, the MLCB/workshop
submission (different framing of the same underlying work, standard
practice for a benchmark paper). Longer review cycle than a workshop;
submit after the workshop version is in, not simultaneously to the same
audience with the same framing.

### Preprint: bioRxiv, immediately, regardless of workshop/journal timing

**This is what makes the outreach in `out/OUTREACH_DRAFTS.md` concrete** —
none of those four drafts should be sent without a citable public URL to
point to. bioRxiv has no meaningful barrier to entry for this kind of work
and no reason to wait for peer review before posting. **Recommended first
action once the manuscript is finalized: post to bioRxiv, then begin
outreach.**

### NeurIPS Datasets & Benchmarks main track: a reach, not the plan

Stated plainly rather than left ambiguous, per instruction: this is
possible but should not be treated as the primary target. Main-track D&B
review expects either substantially larger scale (more hosts, more
sequences, more tasks) or a headline positive result neither of which this
submission has — n_hosts ≤ 6 and an honest negative are real, defensible
contributions at workshop/MLCB/*Scientific Data* scale, not at main-track
scale. If workshop reviewers respond unusually well and specifically
encourage a main-track resubmission, revisit; do not target it as the
first attempt.

## Preparation status

### Preprint-ready manuscript

`out/MANUSCRIPT.md` — complete draft, all sections written, all provenance
(pre-registered / post-hoc / predicted-then-corrected / attempted-and-retracted)
labeled throughout. **Not yet formatted for bioRxiv submission** (no
LaTeX/Word conversion, no author/affiliation block beyond a placeholder,
reference list has 3 gaps disclosed in Section "References" — see
`out/MANUSCRIPT.md` for exactly which citations still need verification
before this can go out). Figures referenced in the Results text
(calibration curves, mechanism comparison, disattenuation, shift-control)
already exist as PNGs in `out/figures/` from Gates 4–8.6 but have not been
assembled into a single figure set with manuscript-style numbering/captions
— that assembly is outstanding.

### Zenodo deposit structure — settled, not executed

Per Gate 8.5/8.6: **two separate deposits.**
1. **Core benchmark** (MIT license) — `data/core/`, `crosshost/` package,
   baselines, held-out evaluation split. Metadata draft:
   `package/ARCHIVE_METADATA.md`.
2. **PromoGen2-derived embeddings** (CC-BY-NC-4.0, non-commercial,
   matching PromoGen2's own license) — kept structurally separate so the
   core deposit stays cleanly, commercially reusable. Metadata draft: same
   file, second deposit block.

Neither deposit has been created. `package/ARCHIVE_METADATA.md` documents
what remains: creator name/affiliation/ORCID, the Johns et al. DOI
double-checked before citing as a `related_identifiers` field, and actually
running the Zenodo upload once an account exists.

### Hugging Face Datasets mirror — drafted, not executed

Front-matter and structure drafted in `package/ARCHIVE_METADATA.md`.
Requires a Hugging Face account/organization (not yet created) and a
decision on whether `data/core/` and `data/licensed/*` upload as one
dataset with clear license disclosure or as separate HF datasets
(recommended, matching the Zenodo split) — not yet executed either way.

### Repository

`https://github.com/gabeykim/crosshost` — created, pushed, **currently
private**. `pip install`, the full test suite, `make audit`, and
`make reproduce` are all verified working from a genuine clean clone (Gate
8.5 Task 4). **Not made public.**

## What remains for Gabriel — nothing below has been done by this session

1. **Zenodo account** (if none exists) and actually creating both deposits,
   minting two DOIs.
2. **Hugging Face account/organization** (if none exists) and actually
   uploading the dataset mirror(s).
3. **Making the GitHub repository public** — currently private; this
   session has not done this and will not without explicit instruction.
4. **Filling in the manuscript's author/affiliation block**, converting to
   bioRxiv's submission format, and resolving the 3 disclosed reference
   gaps (DNABERT-2's exact citation, the DART-Eval and ICLR 2025 baseline
   citations, and re-verifying the BEND/Genomic Benchmarks/DART-Eval/
   DNALONGBENCH/NT-suite "human/animal/plant" claim — see
   `out/MANUSCRIPT.md` Section 6 and its References section).
5. **Assembling the manuscript figure set** with consistent numbering and
   captions from the existing `out/figures/*.png` files.
6. **Confirming current CFP deadlines** for MLCB and relevant NeurIPS/ICLR
   workshops.
7. **Actually posting to bioRxiv.**
8. **Actually sending any of the four outreach drafts** in
   `out/OUTREACH_DRAFTS.md` — every recipient contact needs independent
   verification first (each draft has a bracketed reminder), and none
   should go out without a live preprint URL to reference.
9. **Deciding submission order and timing** across MLCB/workshop vs.
   *Scientific Data* vs. the outreach emails — this document recommends
   preprint first, then outreach, then workshop submission, then
   *Scientific Data*, but the actual calendar is Gabriel's call.
