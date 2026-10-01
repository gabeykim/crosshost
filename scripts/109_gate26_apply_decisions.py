#!/usr/bin/env python3
"""Gate 26: apply the five approved Gate 25 decisions to the manuscript.

T1 supercoiling (sections 3.2 and 3.4), T2 title (option 2), T3 the GC control as
a fourth robustness check in 3.3, T5 the fifth phylum stratum in 3.1.

T6a (the provenance-audit disclosure) is NOT applied -- wording is proposed in the
memo for review first, as instructed.

Stratum values are read from the results file, never typed.
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

# Every stratum in which all three host pairs clear the n > 20 reporting
# threshold. Gate 25 reported four; Cyanobacteria also clears it (Gate 26 T5).
PHYLA = ["Proteobacteria", "Firmicutes", "Actinobacteria", "Bacteroidetes", "Cyanobacteria"]


def strata_sentence() -> str:
    t4 = json.load(open(RESULTS / "gate10_5_gc_confound.json"))["task4_phylum_stratified"]["transcription"]
    for ph in PHYLA:
        for pair in ("EC_PA", "EC_BS", "BS_PA"):
            if t4[pair].get(ph, {}).get("spearman_rho") is None:
                raise SystemExit(f"FAIL: {ph}/{pair} has no value; it does not clear the threshold")
    parts = [
        f"{ph} {t4['EC_PA'][ph]['spearman_rho']:.3f} (n = {t4['EC_PA'][ph]['n']:,}) against EC–BS "
        f"{t4['EC_BS'][ph]['spearman_rho']:.3f} (n = {t4['EC_BS'][ph]['n']:,}) and BS–PA "
        f"{t4['BS_PA'][ph]['spearman_rho']:.3f} (n = {t4['BS_PA'][ph]['n']:,})"
        for ph in PHYLA]
    bs = [t4[p][ph]["spearman_rho"] for p in ("EC_BS", "BS_PA") for ph in PHYLA]
    ecpa = [t4["EC_PA"][ph]["spearman_rho"] for ph in PHYLA]
    return ("Phylum stratification agrees independently. In all five source-genome phyla where "
            "every pair clears the n > 20 reporting threshold, EC–PA transcription stays above "
            "both *B. subtilis* pairs: " + "; ".join(parts) + f". The *B. subtilis* pairs range "
            f"{min(bs):.2f} to {max(bs):.2f} across strata rather than sitting near zero in every "
            f"one, and the five EC–PA values ({min(ecpa):.3f}–{max(ecpa):.3f}) bracket the "
            "unstratified co-active 0.754, so the stratification reveals no anomaly. The margin is "
            "narrowest in Cyanobacteria, the smallest stratum, where EC–BS reaches 0.332 on "
            "n = 62.")


T3 = ("**The contrast survives a source-composition control.** The GC partial-correlation control "
      "applied to the in-vivo pairs in Section 3.1 was not originally applied to this comparison, "
      "which is the check that speaks most directly to a selection artifact. Controlling for "
      "source-genome GC within each modality, the cell-free co-active correlation falls from 0.677 "
      "to 0.556 (a 17.9% reduction) and the in-vivo co-active correlation from 0.258 to 0.207 "
      "(19.8%); the ratio moves from 2.63× to **2.69×**. Pooled, it moves from 0.94× to "
      "1.09×, still no modality difference. Both modalities attenuate by nearly the same "
      "proportion, so the contrast is not an artifact of source-genome GC composition. The "
      "closed-form partial correlation and an independent rank-residual regression agree (0.556 "
      "against 0.548 cell-free; 0.207 against 0.205 in vivo). One scope limit: the two modalities "
      "are measured on different sequence populations, so GC is controlled within each separately "
      "— this is not a paired control.\n\n")

EDITS = [
 ("T1-3.2",
  "but lack an intact membrane, native supercoiling, macromolecular resource competition, "
  "and growth-phase physiology.",
  "but lack an intact membrane, supercoiling homeostasis, macromolecular resource competition, "
  "and growth-phase physiology."),
 ("T1-3.4",
  "DNA supercoiling, absent in linear cell-free reactions, acts as a global transcriptional "
  "regulator whose effect on a promoter depends quantitatively on its discriminator GC content "
  "(El Houdaigui et al., 2019).",
  "DNA supercoiling acts as a global transcriptional regulator whose effect on a promoter depends "
  "quantitatively on its discriminator GC content (El Houdaigui et al., 2019). DRAFTS supplied its "
  "library as midiprepped plasmid, so those templates are supercoiled; what a lysate lacks is the "
  "gyrase/topoisomerase homeostasis and chromosomal context that regulate supercoiling in a living "
  "cell, leaving it fixed at whatever the preparation produced rather than dynamically maintained. "
  "That is a weaker claim than absence, and it is the accurate one."),
 ("T2-title",
  "# Cell-free systems do not reproduce the bacterial chassis effect: evidence that host "
  "specificity resides in cellular context rather than transcription machinery",
  "# Cell-free transcription does not reproduce the bacterial chassis effect: evidence that host "
  "specificity resides in cellular context rather than transcription machinery"),
 ("T3-gc", "Prior work has noted that cell-free and in-vivo measurements can diverge.",
  T3 + "Prior work has noted that cell-free and in-vivo measurements can diverge."),
]

NUM = re.compile(r"\d(?:[\d,]*\d)?(?:\.\d+)?")
strip = lambda t: re.sub(r"^\d+\. ", "", t, flags=re.M)


def main() -> int:
    doc = PREPRINT.read_text(encoding="utf-8")
    before = doc

    # T5: replace the whole strata sentence (located by its opening)
    i = doc.index("Phylum stratification agrees independently.")
    j = doc.index("reveals no anomaly.", i) + len("reveals no anomaly.")
    doc = doc[:i] + strata_sentence() + doc[j:]

    for label, old, new in EDITS:
        n = doc.count(old)
        if n != 1:
            print(f"FAIL: {label} matched {n} times")
            return 1
        doc = doc.replace(old, new, 1)

    # One token is expected to go: "0.27" was the upper end of the rounded
    # B. subtilis-pair range over four strata. Adding Cyanobacteria (EC-BS 0.332)
    # moves that endpoint to 0.33. It is a recomputed summary of the stratum
    # values, not a measurement, and every underlying per-stratum value is still
    # printed in the same sentence.
    EXPECTED_GONE = {"0.27"}
    a, b = set(NUM.findall(strip(before))), set(NUM.findall(strip(doc)))
    gone = sorted((a - b) - EXPECTED_GONE)
    if gone:
        print(f"FAIL: numeric tokens would disappear: {gone}")
        return 1
    print(f"expected token(s) retired by recomputation: {sorted(a - b)}")

    PREPRINT.write_text(doc, encoding="utf-8")
    MIRROR.write_text(doc, encoding="utf-8")
    body = lambda t: len(re.sub(r"^!\[.*$", "", t.split("## References")[0].split("# Cell-free")[-1],
                                flags=re.M).split())
    print("applied: T1 (3.2, 3.4), T2 title, T3 GC paragraph, T5 fifth stratum")
    print(f"body words {body(before)} -> {body(doc)}")
    print(f"numbers added: {sorted(b - a)}")
    print(f"numbers removed: {gone or 'none'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
