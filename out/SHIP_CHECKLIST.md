# CROSSHOST — Ship Checklist

**Nothing has been submitted, deposited, made public, or sent. This is a checklist, not a log of actions taken.**

## One-line summary of what the preprint claims (sanity-check this before posting)

A pre-registered test found host-conditioned models do not beat a sequence-only baseline at predicting bacterial regulatory DNA activity across hosts (n_hosts ≤ 6), and the sequence-only model is also badly miscalibrated off its training host; what bounds performance instead is raw cross-host measurement agreement itself (a Johns et al. 2018 pattern we quantify and robustness-check two independent ways — measurement noise, and separately, source-genome GC composition), with an independent cell-free dataset (DRAFTS) suggesting the reason is that host specificity lives in cellular context, not in the transcriptional machinery our genomic features described.

If that doesn't match your own understanding of what this project found, say so before anything below happens — the checklist assumes it's right.

---

## Ordered checklist

| # | Action | Owner | Status |
|---|---|---|---|
| 1 | Sanity-check the one-line summary above against your own understanding | **Gabriel** | Not done |
| 2 | Review `out/PREPRINT/MANUSCRIPT.md` in full, including the "What I could not do" section and the disclosed reference gaps | **Gabriel** | Not done |
| 3 | Fill in author name, affiliation, and (if desired) ORCID — currently a placeholder everywhere it appears (manuscript title block, `package/ARCHIVE_METADATA.md`) | **Gabriel** | Not done |
| 4 | Resolve the 2 remaining disclosed reference gaps before posting if you want them closed first (DNABERT-2's exact citation; the DART-Eval/ICLR-2025-baseline citations and the BEND/Genomic Benchmarks/DART-Eval/DNALONGBENCH/NT-suite task-composition claim — see `out/PREPRINT/MANUSCRIPT.md` References and Section 6) — **or post with them disclosed as-is; both are legitimate options, this is a judgment call, not a blocker** | **Gabriel** | Not done |
| 5 | Convert `out/PREPRINT/MANUSCRIPT.md` to bioRxiv's submission format (their own template/portal handles most of this; the content itself is complete) | **Gabriel** | Not done |
| 6 | Create a bioRxiv account (if none exists) and submit | **Gabriel** | Not done |
| 7 | Create the two Zenodo deposits (MIT core benchmark; CC-BY-NC-4.0 PromoGen2-derived embeddings) and mint DOIs — metadata drafted in `package/ARCHIVE_METADATA.md`, not yet executed | **Gabriel** | Not done |
| 8 | Once the bioRxiv preprint and Zenodo DOIs exist, update `out/PREPRINT/DATA_AND_CODE_AVAILABILITY.md` and `package/README.md`'s Archive section with the real URLs/DOIs (currently placeholders) | **Either — mechanical once the URLs exist** | Not done |
| 9 | Decide whether to make `github.com/gabeykim/crosshost` public — **confirm with yourself explicitly before doing this; the repo is currently private and this session will not change that visibility without your direct instruction** | **Gabriel** | Not done |
| 10 | Verify each outreach draft's recipient contact info in `out/OUTREACH_DRAFTS.md` (all 5 have bracketed verification reminders; one pair — recipients 1 and 5 — share a surname across two unrelated labs, flagged explicitly in that file) | **Gabriel** | Not done |
| 11 | Send outreach, once the preprint has a live URL — recommended order in `out/OUTREACH_DRAFTS.md`'s own sending checklist (Harris Wang lab first, since their answers could still affect the manuscript) | **Gabriel** | Not done |
| 12 | Confirm current-year CFP deadlines and scope for MLCB and relevant NeurIPS/ICLR workshops (not checked this session — see `out/SUBMISSION_CHECKLIST.md`) | **Gabriel** | Not done |
| 13 | Submit to the chosen workshop/MLCB venue, then *Scientific Data* as a secondary submission afterward (not simultaneously — see `out/SUBMISSION_CHECKLIST.md` for the reasoning) | **Gabriel** | Not done |

## What was done on this session's side, this gate (context, not a to-do)

- Manuscript integrated: DRAFTS section, GC-control robustness check, updated Limitations (12 items, synced with `out/KNOWN_ISSUES.md`), final abstract (245 words).
- Every load-bearing number in the manuscript reconciled against `out/results/*.json` and `out/state.json` — zero discrepancies found (see the Part VI status report for the full account).
- `out/PREPRINT/` assembled: figure set (9 main + 19 supplementary, audited for provenance, 3 stale/retired-metric figures found and excluded), supplementary tables, data/code availability statement, reproducibility statement.
- `make reproduce` extended to cover Gates 8.5–10.5 (previously missing) and verified end-to-end this gate — exit code 0, numbers independently reproduced bit-for-bit against the manuscript.
- `package/README.md`'s reproducibility claims corrected (was stale, predating Gate 8.5's clean-clone verification).
- Outreach drafts rewritten in full (5 recipients, reordered by value, DRAFTS/GC-control findings incorporated, a real naming-collision risk between two "Wang lab" recipients flagged explicitly).

## What remains unresolved and is not this session's call

Nothing above is blocked on further analysis — everything remaining is either a account-creation/submission action, a judgment call on timing/venue, or a decision (repo visibility, whether to resolve the 2 disclosed reference gaps first) that belongs to Gabriel per the charter's own standing rule that decisions with consequences outside this session are not this session's to make unilaterally.
