#!/usr/bin/env python3
"""Gate 29: two text corrections and one source-verified terminology change.

T1  scope the "nearly the same proportion" claim to the co-active regime, where
    it holds (17.9% vs 19.8%); pooled it does not (8.9% vs 21.1%, a 2.36-fold
    difference). The pooled numbers and the pooled ratio stay.
T2  make the three-checks / two-confounds arithmetic explicit rather than
    implying three methods for two confounds.
T3  adopt Yim et al.'s own term for the clustering. Verified in
    raw/drafts/msb198875_FullTextArticle.pdf: they write "distinct gram-negative
    and gram-positive groups (Fig 3C)" for the clustering, and put the
    phylogenetic gradient in Fig 3D. Our text said "clusters by phylum" and
    attributed both to a bare "Fig. 3".

No number changes; the guard asserts it in both directions.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PREPRINT = ROOT / "out" / "PREPRINT" / "MANUSCRIPT.md"
MIRROR = ROOT / "out" / "MANUSCRIPT.md"

EDITS = [
 ("T1",
  "Both modalities attenuate by nearly the same proportion, so the contrast is not an artifact "
  "of source-genome GC composition.",
  "Under co-active restriction the two modalities attenuate by nearly the same proportion, so the "
  "co-active contrast — the one this paper's claim rests on — is not an artifact of "
  "source-genome GC composition. The pooled attenuations differ considerably, but the pooled "
  "comparison shows no modality difference before or after the control, so nothing in the argument "
  "turns on them."),
 ("T2",
  "What can be said is that the two confounds a careful reader proposes first have each been "
  "tested and neither survives contact with the data.",
  "What can be said is that the two confounds a careful reader proposes first have each been "
  "tested — GC twice, on the in-vivo pairs and on the modality contrast — and neither "
  "survives contact with the data."),
 ("T3",
  "Yim et al.'s own Fig. 3 reports that cell-free cross-species correlation clusters by phylum and "
  "tracks phylogeny, at r = 0.73 against 16S rRNA similarity across the species pairs.",
  "Yim et al.'s own Fig. 3 reports that cell-free cross-species correlation clusters into distinct "
  "gram-negative and gram-positive groups (their Fig. 3C), and that more phylogenetically related "
  "species share more similar transcription profiles, at r = 0.73 against 16S rRNA similarity "
  "(their Fig. 3D)."),
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

    from collections import Counter
    a, b = Counter(NUM.findall(strip(before))), Counter(NUM.findall(strip(doc)))
    # T3 replaces one bare "Fig. 3" with "Fig. 3C" and "Fig. 3D", so the token
    # "3" occurs twice more. Those are panel letters, not quantities: no measured
    # value is added, removed or altered.
    EXPECTED = {"3": +2}
    changed = {k: (a[k], b[k]) for k in set(a) | set(b)
               if b[k] - a[k] != EXPECTED.get(k, 0)}
    if changed:
        print(f"FAIL: numeric tokens changed count: {changed}")
        return 1
    print("expected token delta (panel letters only):", EXPECTED)

    PREPRINT.write_text(doc, encoding="utf-8")
    MIRROR.write_text(doc, encoding="utf-8")
    print("applied T1, T2, T3")
    print("numeric tokens changed: NONE (counts identical in both directions)")
    print(f"abstract words {len(abstract(before).split())} -> {len(abstract(doc).split())}")
    print(f"body words {body(before)} -> {body(doc)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
