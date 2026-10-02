#!/usr/bin/env python3
"""Gate 30: four framing edits.

T1  section 3.4 -- make explicit that the competing fitness-decoupling account
    falls on the SAME side of the title's dichotomy, so naming it does not
    weaken the title.
T2  abstract -- put the 93.6%/26.0% co-activity contrast in, add the
    practitioner takeaway, and stay at 250 words.
T3  section 3.3 -- rename to name its contents.
T4  section 3.3 -- cross-reference the other open check.

Task 5 is an assessment and applies nothing.

The only measured values added anywhere are 93.6 and 26.0, both already in
section 3.2; the guard asserts it.
"""
from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PREPRINT = ROOT / "out" / "PREPRINT" / "MANUSCRIPT.md"
MIRROR = ROOT / "out" / "MANUSCRIPT.md"

EDITS = [
 # ---------- T1: section 3.4 ----------
 ("T1",
  "We name it as a competing explanation of equal standing rather than one this paper rules out.",
  "We name it as a competing explanation of equal standing rather than one this paper rules out. "
  "Both accounts nonetheless place the difference in cellular context rather than in transcriptional "
  "machinery — expression coupled to growth is a property of a living cell, not of its "
  "polymerase. What the data available here cannot settle is which contextual mechanism dominates, "
  "not whether the difference is contextual."),
 # ---------- T2: abstract ----------
 ("T2-open",
  "Comparing cell-free (DRAFTS, Yim et al. 2019) against in-vivo (Johns et al. 2018) measurements "
  "of the same 165 bp library, we find that cross-host agreement responds to conditioning on "
  "co-activity very differently in the two modalities.",
  "We compare cell-free (DRAFTS, Yim et al. 2019) against in-vivo (Johns et al. 2018) measurements "
  "of the same 165 bp library."),
 ("T2-pooled",
  "Pooled, *E. coli*–*B. subtilis* transcription correlates at **0.616 cell-free and 0.655 in "
  "vivo** — no modality difference.",
  "Pooled, the two modalities are indistinguishable (**0.616 cell-free, 0.655 in vivo**)."),
 ("T2-silent",
  "In living cells, cross-host agreement is carried substantially by agreement on which sequences "
  "are silent; lysates silence almost nothing.",
  "In living cells, cross-host agreement is carried substantially by agreement on which sequences "
  "are silent: **26.0% of these sequences are co-active in vivo against 93.6% in lysate**. "
  "Screening constructs for one fixed chassis via lysate is therefore reasonable; ranking a part "
  "across candidate hosts is not."),
 ("T2-leadin",
  "That interpretation predicts host-descriptor features encoding machinery composition should "
  "carry no signal, and we tested that prediction directly:",
  "That interpretation predicts host-descriptor features encoding machinery composition carry no "
  "signal, and we tested it:"),
 ("T2-avail",
  "We release the benchmark, splits, evaluation code, and baseline suite, including the negative "
  "result and a retracted rescue attempt.",
  "We release the benchmark, splits, evaluation code, and baselines, including the negative result "
  "and a retraction."),
 # ---------- T3: section 3.3 title ----------
 ("T3", "### 3.3 Scope and supporting checks",
  "### 3.3 Magnitude, robustness, and the published objections"),
 # ---------- T4: cross-reference the other open check ----------
 ("T4",
  "We flag it as the open check rather than a settled point, because the gradient could as easily "
  "strengthen under co-active restriction as weaken.",
  "We flag it as the open check rather than a settled point, because the gradient could as easily "
  "strengthen under co-active restriction as weaken. It is one of two checks this paper leaves "
  "open; the other, cross-modality correlations for the three RS241 hosts DRAFTS covers, is stated "
  "in Limitations item 2."),
]

NUM = re.compile(r"\d(?:[\d,]*\d)?(?:\.\d+)?")
strip = lambda t: re.sub(r"^\d+\. ", "", t, flags=re.M)
abstract = lambda t: t.split("## Abstract")[1].split("\n---")[0].strip()
body = lambda t: len(re.sub(r"^!\[.*$", "", t.split("## References")[0].split("# Cell-free")[-1],
                            flags=re.M).split())

# T2 moves 26.0 and 93.6 into the abstract (both already in 3.2, so each count
# rises by one). T4's "Limitations item 2" and "three RS241" add one "2" and one
# "241"; those are a cross-reference and a panel-style identifier, not measured
# values. Everything else must be unchanged.
EXPECTED = {"26.0": +1, "93.6": +1, "2": +1, "241": +1}


def main() -> int:
    doc = PREPRINT.read_text(encoding="utf-8")
    before = doc
    for label, old, new in EDITS:
        n = doc.count(old)
        if n != 1:
            print(f"FAIL: {label} matched {n} times: {old[:60]!r}")
            return 1
        doc = doc.replace(old, new, 1)

    a, b = Counter(NUM.findall(strip(before))), Counter(NUM.findall(strip(doc)))
    changed = {k: (a[k], b[k]) for k in set(a) | set(b) if b[k] - a[k] != EXPECTED.get(k, 0)}
    if changed:
        print(f"FAIL: unexpected numeric token changes: {changed}")
        return 1

    aw = len(abstract(doc).split())
    if aw > 250:
        print(f"FAIL: abstract is {aw} words, over the 250 ceiling")
        return 1

    PREPRINT.write_text(doc, encoding="utf-8")
    MIRROR.write_text(doc, encoding="utf-8")
    print("applied T1, T2 (5 edits), T3, T4")
    print(f"abstract {len(abstract(before).split())} -> {aw} words")
    print(f"body {body(before)} -> {body(doc)} words")
    print(f"numeric token deltas, all expected: {EXPECTED}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
