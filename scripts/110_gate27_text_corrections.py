#!/usr/bin/env python3
"""Gate 27 Task 2: three text corrections.

2a  report the pooled GC-controlled correlations, so the pooled ratio is as
    checkable as every other number in the paragraph
2b  the Discussion says two robustness checks; there are three
2c  remove drafting residue ("was not originally applied to this comparison")

The pooled values are read from the results file, never typed.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PREPRINT = ROOT / "out" / "PREPRINT" / "MANUSCRIPT.md"
MIRROR = ROOT / "out" / "MANUSCRIPT.md"
RESULTS = ROOT / "out" / "results"


def pooled_clause() -> str:
    p = json.load(open(RESULTS / "gate25_modality_gc_control.json"))["pooled"]
    cf, iv = p["cell_free"], p["in_vivo"]
    return (f"Pooled, the cell-free correlation falls from {cf['rho_raw']:.3f} to "
            f"{cf['rho_gc_controlled_formula']:.3f} ({abs(cf['pct_change'])}%) and the in-vivo "
            f"from {iv['rho_raw']:.3f} to {iv['rho_gc_controlled_formula']:.3f} "
            f"({abs(iv['pct_change'])}%), moving the ratio from {p['ratio_raw']}× to "
            f"{p['ratio_gc_controlled']}× — still no modality difference.")


EDITS = [
 ("2a", "Pooled, it moves from 0.94× to 1.09×, still no modality difference.", None),
 ("2b",
  "The co-active contrast survived two checks — disattenuation for measurement error, and "
  "partial-correlation control for source-genome GC composition, confirmed independently by "
  "phylum stratification (Section 3.1). Neither was pre-registered; both were motivated after "
  "the fact, one by a subsequent review and one by a finding in the cell-free dataset, and "
  "surviving both does not rule out a third.",
  "The co-active contrast survived three checks — disattenuation for measurement error; "
  "partial-correlation control for source-genome GC composition on the in-vivo pairs, confirmed "
  "independently by phylum stratification (Section 3.1); and that same GC control applied to the "
  "modality contrast itself (Section 3.3). None was pre-registered; all three were motivated "
  "after the fact, by a subsequent review, by a finding in the cell-free dataset, and by a "
  "selection-artifact objection respectively, and surviving three does not rule out a fourth."),
 ("2c",
  "The GC partial-correlation control applied to the in-vivo pairs in Section 3.1 was not "
  "originally applied to this comparison, which is the check that speaks most directly to a "
  "selection artifact.",
  "The GC partial-correlation control applied to the in-vivo pairs in Section 3.1 applies equally "
  "to this comparison, where it is the check that speaks most directly to a selection artifact, "
  "so we apply it here."),
]

NUM = re.compile(r"\d(?:[\d,]*\d)?(?:\.\d+)?")
strip = lambda t: re.sub(r"^\d+\. ", "", t, flags=re.M)


def main() -> int:
    doc = PREPRINT.read_text(encoding="utf-8")
    before = doc
    for label, old, new in EDITS:
        if new is None:
            new = pooled_clause()
        n = doc.count(old)
        if n != 1:
            print(f"FAIL: {label} matched {n} times")
            return 1
        doc = doc.replace(old, new, 1)

    a, b = set(NUM.findall(strip(before))), set(NUM.findall(strip(doc)))
    gone = sorted(a - b)
    if gone:
        print(f"FAIL: numeric tokens would disappear: {gone}")
        return 1

    PREPRINT.write_text(doc, encoding="utf-8")
    MIRROR.write_text(doc, encoding="utf-8")
    body = lambda t: len(re.sub(r"^!\[.*$", "", t.split("## References")[0].split("# Cell-free")[-1],
                                flags=re.M).split())
    print("applied 2a, 2b, 2c")
    print(f"body words {body(before)} -> {body(doc)}")
    print(f"numbers added: {sorted(b - a) or 'none'}   removed: {gone or 'none'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
