#!/usr/bin/env python3
"""Gate 20: embed the eight main-text figures, with their verbatim captions, into
out/PREPRINT/MANUSCRIPT.md.

Captions are lifted verbatim from out/PREPRINT/FIGURE_LEGENDS.md (main-text section
only).  The only transformation applied is joining a caption's soft-wrapped lines
with a single space, which a Markdown image alt-text must be.  No word is added,
removed, or reordered; the script asserts this by comparing whitespace-normalised
token streams before writing.

Figure 2 is two files (panels a and b).  Panel (a) is emitted as a plain image and
panel (b) carries the single Figure 2 caption, so the caption sits below both panels
and is not split.

Each figure is inserted immediately after the paragraph it illustrates, identified by
a unique anchor string rather than a line number.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "out" / "PREPRINT" / "MANUSCRIPT.md"
LEGENDS = ROOT / "out" / "PREPRINT" / "FIGURE_LEGENDS.md"
FIGDIR = ROOT / "out" / "PREPRINT" / "figures"

# figure number -> image files, in the order they should appear.  The caption is
# attached to the LAST file listed, so a multi-panel figure gets one caption below
# all of its panels.
FILES = {
    1: ["Figure1_disattenuation.png"],
    2: ["Figure2a_gc_control.png", "Figure2b_phylum_stratified.png"],
    3: ["Figure3_drafts_modality_comparison.png"],
    4: ["Figure4_hmain_kill_gate.png"],
    5: ["Figure5_conditioning_mechanisms.png"],
    6: ["Figure6_calibration_collapse.png"],
    7: ["Figure7_film_instability.png"],
    8: ["Figure8_shift_prediction_retraction.png"],
}

# figure number -> unique opening substring of the paragraph it illustrates.
ANCHORS = {
    1: "**The gap survives measurement noise.**",
    2: "**The gap survives a source-composition confound.**",
    3: "**Conditioning on co-activity separates them, in opposite directions.**",
    5: "**Sequence-only prediction is not distinguishably beaten",
    4: "**The pre-registered kill gate failed on all 8 primary comparisons.**",
    6: "The sequence-only model's zero-shot classifier is well calibrated",
    7: "**FiLM is measurably unstable, independent of whether conditioning helps.**",
    8: "It did not survive its control.",
}

# The YAML preamble pandoc needs so that (i) figures stay at the point of insertion
# instead of floating away from their paragraph, and (ii) LaTeX does not prepend its
# own "Figure N:" to a caption that already begins "**Figure N — ...**".
PREAMBLE = """---
header-includes: |
  \\usepackage{float}
  \\floatplacement{figure}{H}
  \\usepackage{caption}
  \\captionsetup{labelformat=empty,font=small,justification=raggedright,singlelinecheck=false}
---

"""


def read_captions() -> dict[int, str]:
    """Return {figure number: caption, joined to one line} for the main-text figures."""
    text = LEGENDS.read_text(encoding="utf-8")
    main = text.split("## Main text", 1)[1].split("## Supplementary", 1)[0]
    captions: dict[int, str] = {}
    for block in main.split("\n\n"):
        block = block.strip()
        m = re.match(r"\*\*Figure (\d+) —", block)
        if not m:
            continue
        captions[int(m.group(1))] = " ".join(block.split("\n"))
    return captions


def verbatim(caption_joined: str, legends_text: str) -> bool:
    """True if every token of the joined caption appears, in order, in the legends."""
    return " ".join(caption_joined.split()) in " ".join(legends_text.split())


def main() -> int:
    captions = read_captions()
    legends_text = LEGENDS.read_text(encoding="utf-8")

    missing = sorted(set(FILES) - set(captions))
    if missing:
        print(f"FAIL: no caption found for figure(s) {missing}")
        return 1
    for n, names in FILES.items():
        for name in names:
            if not (FIGDIR / name).is_file():
                print(f"FAIL: Figure {n} file not found: {name}")
                return 1
        if not verbatim(captions[n], legends_text):
            print(f"FAIL: Figure {n} caption is not verbatim after line-joining")
            return 1
        if "[" in captions[n] or "]" in captions[n]:
            print(f"FAIL: Figure {n} caption contains a square bracket")
            return 1

    doc = MANUSCRIPT.read_text(encoding="utf-8")
    if doc.lstrip().startswith("---\nheader-includes"):
        print("FAIL: manuscript already carries the pandoc preamble; nothing done")
        return 1
    if "](figures/" in doc:
        print("FAIL: manuscript already embeds figures; nothing done")
        return 1

    for n in sorted(FILES, key=lambda k: list(ANCHORS).index(k)):
        anchor = ANCHORS[n]
        if doc.count(anchor) != 1:
            print(f"FAIL: anchor for Figure {n} occurs {doc.count(anchor)} times")
            return 1
        start = doc.index(anchor)
        end = doc.index("\n\n", start)  # end of the anchor paragraph
        names = FILES[n]
        blocks = [f"![](figures/{name})" for name in names[:-1]]
        blocks.append(f"![{captions[n]}](figures/{names[-1]})")
        doc = doc[:end] + "\n\n" + "\n\n".join(blocks) + doc[end:]

    doc = PREAMBLE + doc

    # title block
    if doc.count("Stanford University, Stanford, CA, USA\nCorrespondence: [email]") != 1:
        print("FAIL: title block not in the expected form")
        return 1
    doc = doc.replace(
        "Stanford University, Stanford, CA, USA\nCorrespondence: [email]",
        "Stanford University, Stanford, CA, USA\\\nCorrespondence: gabeykim@stanford.edu",
    )

    MANUSCRIPT.write_text(doc, encoding="utf-8")
    print(f"embedded {sum(len(v) for v in FILES.values())} image files "
          f"for {len(FILES)} figures; title block updated")
    return 0


if __name__ == "__main__":
    sys.exit(main())
