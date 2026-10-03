#!/usr/bin/env python3
"""Gate 33: name the licence in Data and Code Availability, and add the
verification-failure note adjacent to the provenance disclosure.

Both clauses are digit-free, so no numeric token may change; the guard asserts
counts are identical in both directions.
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
 # Task 1b -- licence named where availability is stated
 ("T1-licence",
  "PromoGen2 and DNABERT-2 embeddings were computed from publicly released model "
  "weights and are not redistributed.",
  "PromoGen2 and DNABERT-2 embeddings were computed from publicly released model "
  "weights and are not redistributed. Code and benchmark artifacts are released under "
  "the MIT license; source datasets retain their original terms, documented per dataset "
  "in the repository."),
 # Task 2 -- verification-failure note, immediately after the provenance disclosure
 ("T2-note",
  "Producing scripts were subsequently added and the audit now covers them genuinely.",
  "Producing scripts were subsequently added and the audit now covers them genuinely. "
  "More generally, build and regeneration steps should be checked on exit status and on "
  "whether the output is newer than its input, rather than on reported success, because "
  "a step that silently does nothing leaves its previous output in place and every "
  "downstream check will then validate a stale artifact without error. This occurred "
  "three times during development: a figure-regeneration pass in which five of six "
  "producers never ran, the provenance false pass described above, and a document build "
  "that failed on a full disk while the checks that followed validated the previous "
  "build."),
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

    a, b = Counter(NUM.findall(strip(before))), Counter(NUM.findall(strip(doc)))
    changed = {k: (a[k], b[k]) for k in set(a) | set(b) if a[k] != b[k]}
    if changed:
        print(f"FAIL: numeric tokens changed: {changed}")
        return 1
    if len(abstract(before).split()) != len(abstract(doc).split()):
        print("FAIL: abstract changed; it should not be touched")
        return 1

    PREPRINT.write_text(doc, encoding="utf-8")
    MIRROR.write_text(doc, encoding="utf-8")
    print("applied: licence clause + verification-failure note")
    print("numeric tokens changed: NONE (counts identical in both directions)")
    print(f"abstract {len(abstract(doc).split())} words (unchanged)")
    print(f"body {body(before)} -> {body(doc)} words")
    return 0


if __name__ == "__main__":
    sys.exit(main())
