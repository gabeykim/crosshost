#!/usr/bin/env python3
"""Gate 32 Task 1: point the Section 3.5 inference clause at the measurement
that backs it, converting an asserted limitation into a checkable one.

Task 2's wording is proposed in the memo and NOT applied, per instruction.
"""
from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PREPRINT = ROOT / "out" / "PREPRINT" / "MANUSCRIPT.md"
MIRROR = ROOT / "out" / "MANUSCRIPT.md"

OLD = ("the prediction is consistent with the result, not confirmed by it.")
NEW = ("the prediction is consistent with the result, not confirmed by it "
       "(Section 3.7 measures what that costs).")

NUM = re.compile(r"\d(?:[\d,]*\d)?(?:\.\d+)?")
strip = lambda t: re.sub(r"^\d+\. ", "", t, flags=re.M)
abstract = lambda t: t.split("## Abstract")[1].split("\n---")[0].strip()
body = lambda t: len(re.sub(r"^!\[.*$", "", t.split("## References")[0].split("# Cell-free")[-1],
                            flags=re.M).split())

# The cross-reference adds one "3.7". That is a section pointer, not a measured
# value; nothing else may change.
EXPECTED = {"3.7": +1}


def main() -> int:
    doc = PREPRINT.read_text(encoding="utf-8")
    before = doc
    n = doc.count(OLD)
    if n != 1:
        print(f"FAIL: anchor matched {n} times")
        return 1
    doc = doc.replace(OLD, NEW, 1)

    a, b = Counter(NUM.findall(strip(before))), Counter(NUM.findall(strip(doc)))
    changed = {k: (a[k], b[k]) for k in set(a) | set(b) if b[k] - a[k] != EXPECTED.get(k, 0)}
    if changed:
        print(f"FAIL: unexpected numeric token changes: {changed}")
        return 1
    if len(abstract(before).split()) != len(abstract(doc).split()):
        print("FAIL: abstract changed; it should not be touched")
        return 1

    PREPRINT.write_text(doc, encoding="utf-8")
    MIRROR.write_text(doc, encoding="utf-8")
    print("applied: Section 3.7 cross-reference on the 3.5 inference clause")
    print(f"numeric token delta, expected only: {EXPECTED}")
    print(f"abstract {len(abstract(doc).split())} words (unchanged)")
    print(f"body {body(before)} -> {body(doc)} words")
    return 0


if __name__ == "__main__":
    sys.exit(main())
