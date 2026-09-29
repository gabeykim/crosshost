#!/usr/bin/env python3
"""Gate 23 Task 1: consolidate Limitations from 19 items to 12.

Two clauses are relocated rather than deleted, so that no number leaves the
manuscript:

  * the ~300x unbatched-inference cost cliff and the 1,024-row chunk size move to
    Section 2.4 (Architecture), where an implementation detail belongs;
  * the `make reproduce-full` / clean-clone / platform clause moves to Data and
    Code Availability.

Everything else is a merge of existing items.  The script asserts that every
number present in the old Limitations section is still somewhere in the new
document before it writes, and prints the 12-item roll-up.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PREPRINT = ROOT / "out" / "PREPRINT" / "MANUSCRIPT.md"
MIRROR = ROOT / "out" / "MANUSCRIPT.md"

LIM_HEAD = "## 5. Limitations and unresolved gaps"
AVAIL_HEAD = "## Data and Code Availability"

NEW_LIMITATIONS = """## 5. Limitations and unresolved gaps

1. **n_hosts ≤ 6.** The generalization unit is the host, not the sequence: three primary hosts with dense coverage, three more at reduced coverage. Claims about "genome-encoded functions" or "bacterial regulatory prediction" as a general class are not supported by six data points. Every interval bootstraps over folds and hosts, not sequences.

2. **The modality contrast rests on one species pair and is regime-dependent.** *E. coli*–*B. subtilis* is the only pair present in both datasets. Cell-free EC–BS ranges 0.386–0.677 across six restriction levels against an in-vivo co-active 0.258 — a ratio between 1.5× and 2.6×. Under the harshest cut, 4 of 45 DRAFTS pairs fall below 0.258, all involving *L. lactis*. The pooled cell-free figures require imputing activity 0 for `no_RNA_counts` rows, a definitional assignment rather than a measurement. We did not compute cross-modality correlations for the three RS241 hosts DRAFTS covers, which is the most direct route to widening the comparison.

3. **No direct reliability estimate exists for *B. subtilis* or *P. aeruginosa*, or for translation in any host.** No replicate or condition-series data exists for them in the released tables. The one measured reliability (0.912, EC transcription) is for the host and readout the central claim depends on least; *B. subtilis*'s own reliability is bounded by a sensitivity grid, which is a different strength of evidence than a measurement. Relatedly, we did not re-derive Johns et al.'s per-cell Pearson values from their Supplementary Figs. S13/S15 — only the captions and main-text summary were available in the source materials, and the qualitative claim is confirmed from quoted text.

4. **The translation floor artifact.** `protein_log10` is pinned at a per-host floor for 66.5% / 89.9% / 10.1% of nominally usable rows (EC/BS/PA). Corrected usable-for-regression N: 9,146 / 1,101 / 17,630.

5. **No within-sequence paired in-vivo/cell-free comparison at adequate N.** The DRAFTS–Johns overlap usable in both *E. coli* and *B. subtilis* simultaneously is 15 sequences (ρ = 0.147, reported for completeness and not used as evidence). Section 3.2 compares two populations measured on the same library by the same laboratory. DRAFTS measures transcription only, so those conclusions are scoped to transcriptional host-specificity.

6. **Fold variance exceeds draw variance and grows with N** — in 15 of 23 baseline combinations checked, reaching a 91× ratio at *E. coli* transcription N = 3,000. A single fixed test fold can look precise while being unrepresentative; use the full five-fold rotation when extending this benchmark.

7. **Cross-host calibration failure.** ECE 3.5–9× worse off the training host (Section 3.6). Recalibrate before trusting raw predicted probabilities cross-host.

8. **Regression to the mean is a serious confound for any shift-prediction extension of this dataset.** Run a reference-value-only baseline and a shuffled-sequence control before trusting such a result (Section 3.8).

9. **Three dataset collisions, each of which silently corrupts a naive join (Section 2.2).** "Usable" differs between the two datasets, so joining on `n_shared_usable` compares mismatched restriction regimes. The ID spaces are disjoint, so any join must use sequence text. And DRAFTS's `Pa` is *Pantoea agglomerans*, not *P. aeruginosa* — which is absent from DRAFTS entirely — so a script matching on the bare two-letter code merges two unrelated organisms.

10. **Evo 2 was not evaluated, and the two foundation models that were did not run under their native protocols.** The official package requires CUDA, Flash Attention, and Transformer Engine on a Hopper GPU (verified against the ArcInstitute repository and issue #67); this work ran on Apple Silicon, MPS-only. For the two models evaluated instead, a uniform frozen-embedding head-to-head was chosen for comparability between architecturally different models, but it measurably cost PromoGen2 relative to its published native numbers (Section 3.9), and the comparison was not re-run under that protocol. The capacity objection is weakened, not closed.

11. **We could not separate two explanations for the transcription/translation asymmetry** — that translation is fundamentally less cross-host-conserved, versus that the FACS-seq readout is too noisy to support this analysis. Our working view is a mix weighted toward noise for *B. subtilis* specifically, stated as inference rather than a measured result.

12. **The held-out evaluation split is not cryptographically enforced.** The public model-development table contains the withheld fold's labels, since it is the same table used for development; the protocol relies on convention, not a technical barrier.
"""

ARCH_ANCHOR = ("Those three across six (host, readout) cells give the 18 comparisons "
               "in Section 3.5; FiLM is the comparator, not one of the three.")
ARCH_ADD = (" Unbatched inference over a large pooled array showed a ~300× non-linear "
            "cost cliff on both CPU and Apple Silicon MPS; inference is chunked at 1,024 "
            "rows by default, numerically verified identical to the unbatched path.")

AVAIL_ANCHOR = ("Derived activity values from Johns et al. are redistributed with "
                "attribution; PromoGen2-derived content is held in a separate deposit "
                "under CC BY-NC-4.0.")
AVAIL_ADD = ("\n\n`make reproduce` (the fast path), `make audit`, and `make verify-citations` "
             "were verified from a fresh clone on macOS, Apple Silicon, Python 3.14.2 — a "
             "test that surfaced six undocumented defects, each fixed and the test re-run "
             "until it passed with zero undocumented steps (four iterations). "
             "`make reproduce-full`, the complete from-raw-data pipeline including all model "
             "training, was not re-executed end to end, as it would cost the 40+ cumulative "
             "hours already expended; every script in the dependency graph has run and "
             "produced its output at least once.")

NUM = re.compile(r"\d[\d,]*\.?\d*")


def main() -> int:
    doc = PREPRINT.read_text(encoding="utf-8")

    if LIM_HEAD not in doc or AVAIL_HEAD not in doc:
        print("FAIL: could not find the Limitations or Availability heading")
        return 1
    head, rest = doc.split(LIM_HEAD, 1)
    old_lim_body, tail = rest.split("\n---\n\n" + AVAIL_HEAD, 1)
    old_lim = LIM_HEAD + old_lim_body

    n_old = len(re.findall(r"^\d+\. \*\*", old_lim, flags=re.M))
    if n_old != 19:
        print(f"FAIL: expected 19 items, found {n_old}")
        return 1

    for anchor in (ARCH_ANCHOR, AVAIL_ANCHOR):
        if doc.count(anchor) != 1:
            print(f"FAIL: relocation anchor not unique: {anchor[:50]}")
            return 1

    new = head + NEW_LIMITATIONS + "\n---\n\n" + AVAIL_HEAD + tail
    new = new.replace(ARCH_ANCHOR, ARCH_ANCHOR + ARCH_ADD)
    new = new.replace(AVAIL_ANCHOR, AVAIL_ANCHOR + AVAIL_ADD)

    n_new = len(re.findall(r"^\d+\. \*\*", NEW_LIMITATIONS, flags=re.M))
    if n_new != 12:
        print(f"FAIL: new section has {n_new} items, expected 12")
        return 1

    # Every number that was in the old Limitations must still be somewhere in the
    # document.  The list markers ("13. ") are numbering, not content, so they are
    # stripped before the comparison -- otherwise renumbering 19 items down to 12
    # always reads as a loss.
    def content_numbers(text: str) -> set[str]:
        stripped = re.sub(r"^\d+\. ", "", text, flags=re.M)
        return set(NUM.findall(stripped))

    lost = sorted(content_numbers(old_lim) - content_numbers(new))
    if lost:
        print(f"FAIL: numbers would be lost from the manuscript: {lost}")
        return 1

    PREPRINT.write_text(new, encoding="utf-8")
    MIRROR.write_text(new, encoding="utf-8")
    print(f"Limitations {n_old} -> {n_new} items; no number lost from the manuscript")
    print(f"old Limitations words: {len(old_lim.split())}  new: {len(NEW_LIMITATIONS.split())}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
