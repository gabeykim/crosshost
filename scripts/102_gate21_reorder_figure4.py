#!/usr/bin/env python3
"""Gate 21: move the Figure 4 embed above the Figure 5 embed in
out/PREPRINT/MANUSCRIPT.md, then mirror the file to out/MANUSCRIPT.md.

Gate 20 placed each figure after the paragraph it illustrates.  Section 3.5
discusses the three conditioning mechanisms (Figure 5) before the pre-registered
kill gate (Figure 4), so that rule put Figure 5 in front of Figure 4 in the PDF --
and Figure 5's caption refers back to Figure 4 ("same resampling scheme as
Figure 4"), which a reader then meets first.

This is source order, not LaTeX float placement: \\floatplacement{figure}{H} has
been in the preamble since Gate 20, so floats cannot move at all.  The only fix
that touches no prose is to move the Figure 4 image block earlier.  It goes
immediately after Section 3.5's opening paragraph -- the section is titled
"Testing the prediction (Figures 4-5)" and Figure 4 is the pre-registered test
that titles it, so leading with it reads correctly and restores ascending order.

No manuscript text is changed: the script asserts that the document with all
image lines stripped is byte-identical before and after.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PREPRINT = ROOT / "out" / "PREPRINT" / "MANUSCRIPT.md"
MIRROR = ROOT / "out" / "MANUSCRIPT.md"

FIG4 = "](figures/Figure4_hmain_kill_gate.png)"
FIG5 = "](figures/Figure5_conditioning_mechanisms.png)"
OPENING = ("Four independent lines of evidence test the prediction above, "
           "and each comes out as Section 3.4 says it should.")


def prose(text: str) -> str:
    """The document with image lines removed and blank-line runs collapsed.

    Collapsing matters: stripping an image line leaves the blank lines that
    surrounded it behind, so moving a figure necessarily moves a blank line with
    it.  Collapsing runs compares the prose itself rather than its spacing.
    """
    out: list[str] = []
    for ln in text.split("\n"):
        if ln.startswith("!["):
            continue
        if ln.strip() == "" and out and out[-1] == "":
            continue
        out.append("" if ln.strip() == "" else ln)
    return "\n".join(out)


def main() -> int:
    doc = PREPRINT.read_text(encoding="utf-8")
    before_prose = prose(doc)

    lines = doc.split("\n")
    idx = {}
    for i, ln in enumerate(lines):
        for key, needle in (("fig4", FIG4), ("fig5", FIG5), ("open", OPENING)):
            if needle in ln:
                if key in idx:
                    print(f"FAIL: {key} anchor is not unique")
                    return 1
                idx[key] = i
    missing = [k for k in ("fig4", "fig5", "open") if k not in idx]
    if missing:
        print(f"FAIL: could not locate {missing}")
        return 1

    if idx["fig4"] < idx["fig5"]:
        print("FAIL: Figure 4 already precedes Figure 5; nothing to do")
        return 1
    if not (idx["open"] < idx["fig5"] < idx["fig4"]):
        print("FAIL: unexpected layout of Section 3.5")
        return 1

    fig4_line = lines.pop(idx["fig4"])
    # drop the now-doubled blank line left behind
    if lines[idx["fig4"] - 1] == "" and lines[idx["fig4"]] == "":
        lines.pop(idx["fig4"])
    # insert after the opening paragraph (which is a single line + its blank)
    lines[idx["open"] + 1:idx["open"] + 1] = ["", fig4_line]

    out = "\n".join(lines)

    if prose(out) != before_prose:
        print("FAIL: prose changed; refusing to write")
        return 1
    order = [ln.split("](figures/")[1].rstrip(")") for ln in out.split("\n")
             if ln.startswith("![") and "](figures/" in ln]
    print("figure order now:", " ".join(o.split("_")[0] for o in order))

    PREPRINT.write_text(out, encoding="utf-8")
    MIRROR.write_text(out, encoding="utf-8")
    print("wrote both copies; prose byte-identical")
    return 0


if __name__ == "__main__":
    sys.exit(main())
