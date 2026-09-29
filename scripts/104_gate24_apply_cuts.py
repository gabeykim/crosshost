#!/usr/bin/env python3
"""Gate 24: apply the approved cuts from out/GATE23_CUT_CANDIDATES.md.

C1-C10 are the ten non-borderline candidates, in the order the candidates file
lists them.  C11, C12 and C17 were specified directly in the gate and are not in
that file.  The six borderline candidates are NOT touched.

Two removals extend slightly beyond the quoted span to keep the sentence
well-formed -- C1 takes the ", and" that joined the dropped quote to the one
before it, and the two whole-paragraph cuts take their trailing blank line.  No
word is added anywhere.

The script refuses to write if any numeric token disappears from the document.
Section cross-references are numeric tokens too, so they are tracked separately
and reported rather than blocking.
"""
from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PREPRINT = ROOT / "out" / "PREPRINT" / "MANUSCRIPT.md"
MIRROR = ROOT / "out" / "MANUSCRIPT.md"

# (label, exact text to delete).  Order does not matter; each is located afresh.
CUTS = [
 ("C1", ', and that "it is impractical to postulate that the observable chassis-effect '
        'between a given set of hosts can be explained by a single or even a set of '
        'predictable genome-encoded functions without experimental insight" '
        '(Chan & Bernstein, 2024)'),
 ("C2", " That distinction bears on how these results relate to prior physiology-based "
        "claims (Discussion) and on why this vector's failure does not contradict them."),
 ("C3", " **Amendment 1** added the supplementary unfrozen-conv4 transfer mechanism "
        "alongside the pre-registered frozen-trunk one at N=100, after validation found "
        "frozen-trunk interacts with a small-N degradation; frozen-trunk remained primary "
        "and the supplementary mechanism never entered the headline result."),
 ("C4", "The two regimes answer different questions. Co-active asks how well two hosts "
        "agree on the relative strength of sequences that work in both; pooled asks how "
        "well they agree on what works *and* how strongly. "),
 ("C5", " The hypothesis was committed before any host-conditioned model was trained; "
        "both amendments are dated, preceded the results they touch, and make the bar "
        "harder to clear."),
 ("C6", "We quantify the magnitude rather than claiming it as a novel warning. "),
 ("C7", ", which matters for any model conditioned on a handful of domains independent "
        "of this paper's negative result"),
 ("C8", " The co-active EC–PA-versus-*B. subtilis* contrast has survived two attacks "
        "aimed at explaining it away — disattenuation for measurement error, and "
        "partial-correlation control for source-genome GC composition, the latter "
        "confirmed independently by phylum stratification (Section 3.1)."),
 ("C9", "The modality contrast survived a third check, the one that came closest to "
        "overturning it: pooled against pooled, the two modalities are indistinguishable, "
        "and only the response to co-activity conditioning separates them. We report that "
        "first in Section 3.2 for that reason.\n\n"),
 ("C10", "**On the form of the negative result.** The evidence in Section 3.5 was run "
         "because Sections 3.2–3.4 predicted, before any model was trained, that it "
         "should come out this way. Two of the four independent attempts to find an "
         "exception initially looked like counterexamples; one survives only in a single "
         "explainable cell and the other did not survive at all.\n\n"),
 ("C11", " This paper addresses the second group."),
 ("C12", " We report this first because it is the comparison a reader is most likely to "
         "compute independently."),
 ("C17", "Every reported number is produced by a committed script in the public repository "
         "and verified by the audit suite described in Section 2.5. "),
]

# Trailing punctuation must not join the token: "N=100," and "N=100" are the
# same figure, and a greedy [\d,]* makes them look like different ones.
NUM = re.compile(r"\d(?:[\d,]*\d)?(?:\.\d+)?")
SECREF = re.compile(r"(?:Section|Sections)\s+\d+\.\d+(?:–\d+\.\d+)?")
MUST_SURVIVE = ["0.912", "0.929", "0.8485", "0.341", "0.515", "0.308", "0.240", "0.286"]


def tokens(text: str) -> Counter:
    """Numeric tokens, with list markers and section cross-references removed."""
    t = re.sub(r"^\d+\. ", "", text, flags=re.M)
    t = SECREF.sub(" ", t)
    return Counter(NUM.findall(t))


def main() -> int:
    doc = PREPRINT.read_text(encoding="utf-8")
    before = doc

    for label, text in CUTS:
        n = doc.count(text)
        if n != 1:
            print(f"FAIL: {label} found {n} times, expected 1: {text[:60]!r}")
            return 1
        doc = doc.replace(text, "", 1)

    if "  " in "".join(l for l in doc.split("\n") if not l.startswith("|")):
        doc = re.sub(r"(?<=\S)  +(?=\S)", " ", doc)

    a, b = tokens(before), tokens(doc)
    gone = sorted(k for k in a if b[k] == 0)
    if gone:
        print(f"FAIL: numeric tokens would disappear entirely: {gone}")
        return 1
    fewer = {k: (a[k], b[k]) for k in a if 0 < b[k] < a[k]}

    missing = [v for v in MUST_SURVIVE if v not in doc]
    if missing:
        print(f"FAIL: required figures missing after cuts: {missing}")
        return 1

    PREPRINT.write_text(doc, encoding="utf-8")
    MIRROR.write_text(doc, encoding="utf-8")

    body = doc.split("## References")[0].split("# Cell-free")[-1]
    body = re.sub(r"^!\[.*$", "", body, flags=re.M)
    wb = before.split("## References")[0].split("# Cell-free")[-1]
    wb = re.sub(r"^!\[.*$", "", wb, flags=re.M)
    print(f"applied {len(CUTS)} cuts")
    print(f"body words {len(wb.split())} -> {len(body.split())} "
          f"(-{len(wb.split()) - len(body.split())})")
    print(f"numeric tokens lost entirely: none")
    print(f"numeric tokens now appearing fewer times: {fewer or 'none'}")
    print(f"all required figures present: {MUST_SURVIVE}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
