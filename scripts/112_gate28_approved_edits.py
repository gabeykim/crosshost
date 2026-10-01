#!/usr/bin/env python3
"""Gate 28: four approved edits.

T1  abstract -- add the GC-control clause, remove two spans to stay at 250 words
T2  section 3.3 -- engage Yim et al. Fig. 3 as the strongest published objection
T3  section 3.4 -- the fitness-decoupling competing account, plus a Limitations clause
T4  Data and Code Availability -- the provenance-audit disclosure, verbatim

Both DRAFTS claims were verified against raw/drafts/msb198875_FullTextArticle.pdf
before being written into the manuscript:
  * Fig. 3E reports r = 0.73 against "16S rRNA similarity, %" -- similarity, not
    distance, so the correlation is positive.
  * "the distribution of DNA libraries was much more uniform in vitro than in
    vivo, mostly due to the decoupling of gene expression from cell fitness in
    vitro". Their claim is about LIBRARY REPRESENTATION, not about silent
    sequences; the extension to activity is ours and is labelled as such.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PREPRINT = ROOT / "out" / "PREPRINT" / "MANUSCRIPT.md"
MIRROR = ROOT / "out" / "MANUSCRIPT.md"

T2 = ("\n\n**The strongest published argument the other way.** Yim et al.'s own Fig. 3 reports that "
 "cell-free cross-species correlation clusters by phylum and tracks phylogeny, at r = 0.73 against "
 "16S rRNA similarity across the species pairs. Read as a claim about host specificity, that is "
 "the most direct published challenge to the framing here. The two results concern different "
 "quantities: their gradient describes how cross-host agreement *varies* with phylogenetic "
 "similarity within the cell-free system, while the claim here is about its *level* under matched "
 "restriction. Both can hold at once, and here they do — the 45-pair cell-free band runs "
 "0.623–0.911 co-active and never approaches the in-vivo co-active 0.258 for the one pair "
 "measured in both systems, so a phylogenetic gradient inside that band does not account for the "
 "modality difference. The direct test is to recompute their gradient under the restriction "
 "regimes used here; the released data supports it, and **we have not run it**. We flag it as the "
 "open check rather than a settled point, because the gradient could as easily strengthen under "
 "co-active restriction as weaken.")

T3 = ("\n\n**The strongest competing account comes from the dataset's own authors.** Yim et al. "
 "observe that their DNA library distribution is far more uniform in vitro than in vivo, and "
 "attribute that to the decoupling of gene expression from cell fitness in a lysate. Their "
 "statement is about library representation rather than about silent sequences, so extending it to "
 "activity is our reading, not their claim — but extended, it predicts what we observe: in a "
 "living cell, expression that is costly can be selected against or shut down, producing a large "
 "silent class, while a lysate has no fitness for expression to couple to and silences almost "
 "nothing. That account and the contextual one above make the same prediction for every quantity "
 "available here. Separating them requires a per-sequence fitness or burden measurement in each "
 "host, which neither dataset provides; GC, cross-host activity and the silent fraction are "
 "consistent with both. We name it as a competing explanation of equal standing rather than one "
 "this paper rules out.")

T4 = ("\n\nThe provenance audit matches output filenames against script text, which is a substring "
 "check rather than a dependency check. During development it reported no orphans for the nine "
 "preprint figures although no committed script produced them, because the figure-embedding script "
 "names every file. Producing scripts were subsequently added and the audit now covers them "
 "genuinely.")

LIM2 = (" The mechanism behind the contrast is also unidentified: Yim et al.'s fitness-decoupling "
 "account predicts the same observation, and separating it from the cellular-context account "
 "requires per-sequence fitness or burden data that neither dataset provides.")

EDITS = [
 # --- T1 abstract ---
 ("T1-add", "a 2.6-fold gap.",
  "a 2.6-fold gap that survives partial-correlation control for source-genome GC composition "
  "(2.63× raw, 2.69× controlled)."),
 ("T1-cut1", "the chassis effect, a barrier practitioners describe as unresolved for engineering "
  "non-model hosts.", "the chassis effect."),
 ("T1-cut2", "Pooled over all measured sequences, *E. coli*", "Pooled, *E. coli*"),
 # --- T2 section 3.3, appended after the Pandi paragraph ---
 ("T2", "supports the general caution this section sharpens by identifying which axis of fidelity "
  "holds and which does not.",
  "supports the general caution this section sharpens by identifying which axis of fidelity holds "
  "and which does not." + T2),
 # --- T3 section 3.4, appended after the hypothesis paragraph ---
 ("T3", "or measuring the same library in vivo across a matched growth-condition series in "
  "multiple hosts.",
  "or measuring the same library in vivo across a matched growth-condition series in multiple "
  "hosts." + T3),
 # --- T3 Limitations clause, appended to item 2 ---
 ("T3-lim", "We did not compute cross-modality correlations for the three RS241 hosts DRAFTS "
  "covers, which is the most direct route to widening the comparison.",
  "We did not compute cross-modality correlations for the three RS241 hosts DRAFTS covers, which "
  "is the most direct route to widening the comparison." + LIM2),
 # --- T4 Data and Code Availability ---
 ("T4", "as it would cost the 40+ cumulative hours already expended; every script in the "
  "dependency graph has run and produced its output at least once.",
  "as it would cost the 40+ cumulative hours already expended; every script in the dependency "
  "graph has run and produced its output at least once." + T4),
]

NUM = re.compile(r"\d(?:[\d,]*\d)?(?:\.\d+)?")
strip = lambda t: re.sub(r"^\d+\. ", "", t, flags=re.M)
abstract = lambda t: t.split("## Abstract")[1].split("\n---")[0].strip()
body = lambda t: len(re.sub(r"^!\[.*$", "", t.split("## References")[0].split("# Cell-free")[-1],
                            flags=re.M).split())


def main() -> int:
    doc = PREPRINT.read_text(encoding="utf-8")
    before = doc
    for label, old, new in EDITS:
        n = doc.count(old)
        if n != 1:
            print(f"FAIL: {label} matched {n} times: {old[:60]!r}")
            return 1
        doc = doc.replace(old, new, 1)

    a, b = set(NUM.findall(strip(before))), set(NUM.findall(strip(doc)))
    gone = sorted(a - b)
    if gone:
        print(f"FAIL: numeric tokens would disappear: {gone}")
        return 1

    aw = len(abstract(doc).split())
    PREPRINT.write_text(doc, encoding="utf-8")
    MIRROR.write_text(doc, encoding="utf-8")
    print("applied T1 (add + 2 cuts), T2, T3, T3-limitations, T4")
    print(f"abstract words {len(abstract(before).split())} -> {aw}")
    print(f"body words {body(before)} -> {body(doc)}")
    print(f"numbers added: {sorted(b - a) or 'none'}   removed: {gone or 'none'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
