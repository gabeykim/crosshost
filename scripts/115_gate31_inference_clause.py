#!/usr/bin/env python3
"""Gate 31: the one-clause acknowledgement that Section 3.5 cannot discriminate
between a contextual and a mechanistic account at this n_hosts.

Placed at Section 3.5's opening, immediately after the sentence that makes the
confirmation claim, so a reader meets the qualification before the four evidence
lines rather than after them. Section 3.4's closing already concedes the same
point from the other direction ("the link from them to the modeling failure in
Section 3.5 is an inference"); putting it there too would state one concession
twice in one paragraph.

The four other framing sites are untouched, per instruction.
No numeric token changes -- the clause contains no digits.
"""
from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PREPRINT = ROOT / "out" / "PREPRINT" / "MANUSCRIPT.md"
MIRROR = ROOT / "out" / "MANUSCRIPT.md"

OLD = ("Four independent lines of evidence test the prediction above, and each comes out as "
       "Section 3.4 says it should.")
NEW = ("Four independent lines of evidence test the prediction above, and each comes out as "
       "Section 3.4 says it should. With two training hosts per fold, these features would also "
       "fail if host specificity were mechanistic but simply not learnable at this n_hosts; the "
       "prediction is consistent with the result, not confirmed by it.")

NUM = re.compile(r"\d(?:[\d,]*\d)?(?:\.\d+)?")
strip = lambda t: re.sub(r"^\d+\. ", "", t, flags=re.M)
abstract = lambda t: t.split("## Abstract")[1].split("\n---")[0].strip()
body = lambda t: len(re.sub(r"^!\[.*$", "", t.split("## References")[0].split("# Cell-free")[-1],
                            flags=re.M).split())


def main() -> int:
    doc = PREPRINT.read_text(encoding="utf-8")
    before = doc
    n = doc.count(OLD)
    if n != 1:
        print(f"FAIL: anchor matched {n} times")
        return 1
    doc = doc.replace(OLD, NEW, 1)

    a, b = Counter(NUM.findall(strip(before))), Counter(NUM.findall(strip(doc)))
    changed = {k: (a[k], b[k]) for k in set(a) | set(b) if a[k] != b[k]}
    if changed:
        print(f"FAIL: numeric tokens changed: {changed}")
        return 1

    aw_before, aw_after = len(abstract(before).split()), len(abstract(doc).split())
    if aw_before != aw_after:
        print(f"FAIL: abstract changed {aw_before} -> {aw_after}; it should not be touched")
        return 1

    PREPRINT.write_text(doc, encoding="utf-8")
    MIRROR.write_text(doc, encoding="utf-8")
    print("applied: inference clause at Section 3.5's opening")
    print("numeric tokens changed: NONE (counts identical in both directions)")
    print(f"abstract {aw_after} words (unchanged)")
    print(f"body {body(before)} -> {body(doc)} words (+{body(doc) - body(before)})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
