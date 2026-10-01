#!/usr/bin/env python3
"""Gate 25 Block A4 + Block B: manuscript text corrections.

Applies A4 (all four phylum strata with per-stratum N, and the "near zero"
correction), B3, B4, B6, B7, B8 and B9.

NOT applied, by instruction: B1 (refuted -- the manuscript is already correct),
B2 (supercoiling -- wording proposed in the memo, not applied), B5 (title options
proposed, not chosen), and Block C (reported, not written into the Results).

Stratum values are read from out/results/gate10_5_gc_confound.json rather than
typed, so the figure, the text and the results file cannot drift apart.
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

PHYLA = ["Proteobacteria", "Firmicutes", "Actinobacteria", "Bacteroidetes"]


def strata_sentence() -> str:
    t4 = json.load(open(RESULTS / "gate10_5_gc_confound.json"))["task4_phylum_stratified"]["transcription"]
    parts = []
    for ph in PHYLA:
        ecpa, ecbs, bspa = (t4["EC_PA"][ph], t4["EC_BS"][ph], t4["BS_PA"][ph])
        parts.append(f"{ph} {ecpa['spearman_rho']:.3f} (n = {ecpa['n']:,}) against EC–BS "
                     f"{ecbs['spearman_rho']:.3f} (n = {ecbs['n']:,}) and BS–PA "
                     f"{bspa['spearman_rho']:.3f} (n = {bspa['n']:,})")
    bs_vals = [t4[p][ph]["spearman_rho"] for p in ("EC_BS", "BS_PA") for ph in PHYLA]
    ecpa_vals = [t4["EC_PA"][ph]["spearman_rho"] for ph in PHYLA]
    return ("Phylum stratification agrees independently. In all four source-genome phyla where "
            "every pair clears the reporting threshold, EC–PA transcription stays well above "
            "both *B. subtilis* pairs: " + "; ".join(parts) + ". The *B. subtilis* pairs range "
            f"{min(bs_vals):.2f} to {max(bs_vals):.2f} across strata rather than sitting near zero "
            f"in every one, and the four EC–PA values ({min(ecpa_vals):.3f}–"
            f"{max(ecpa_vals):.3f}) bracket the unstratified co-active 0.754, so the stratification "
            "reveals no anomaly.")


EDITS = [
 # --- A4 ---
 ("A4", "Phylum stratification agrees independently: within Proteobacteria, EC–PA "
        "transcription ρ = 0.634 with *B. subtilis* pairs near zero; within Firmicutes, "
        "EC–PA ρ = 0.812.", None),
 # --- B4 ---
 ("B4", "two robustness checks they did not run",
        "two robustness checks not previously applied to this comparison"),
 # --- B3 ---
 ("B3", "one by a finding in an independent dataset",
        "one by a finding in the cell-free dataset"),
 # --- B6 abstract ---
 ("B6a", "cross-host agreement responds to conditioning on co-activity in opposite directions "
         "in the two modalities",
         "cross-host agreement responds to conditioning on co-activity very differently in the "
         "two modalities"),
 # --- B6 + B7 section 3.2 ---
 ("B6b", "**Conditioning on co-activity separates them, in opposite directions.** Restricted to "
         "sequences active in both hosts, the in-vivo correlation falls to **0.258** (n = 3,668), "
         "a 61% drop, while the cell-free correlation *rises* to **0.677** (n = 807), a 10% "
         "increase. That is a **2.6-fold** modality gap under matched restriction, against no gap "
         "when pooled.",
         "**Conditioning on co-activity separates them.** Restricted to sequences active in both "
         "hosts, the in-vivo correlation falls to **0.258** (n = 3,668), a 61% drop, while the "
         "cell-free correlation *rises* to **0.677** (n = 807), a 10% increase. That is a "
         "**2.6-fold** modality gap under matched restriction, against no gap when pooled. The two "
         "movements are not symmetric and should not be read as such: the in-vivo collapse "
         "discards 10,420 of 14,088 sequences and is the robust half of the contrast, whereas the "
         "cell-free rise follows from dropping 55 of 862 and is small enough to be fragile. All "
         "four figures are population Spearman correlations over the full stated N, not resampled "
         "estimates, so no interval is attached to them."),
 # --- B8 section 3.3 ---
 ("B8a", "and 4 of the 45 DRAFTS species pairs fall below 0.258 — all four involving "
         "*L. lactis*.",
         "and 4 of the 45 DRAFTS species pairs fall below 0.258 — all four involving "
         "*L. lactis*. Those four should be read with caution for a second reason: Yim et al. "
         "excluded ScrFI/CCNGG-containing sequences from their *L. lactis* analysis as "
         "restriction artifacts, so a *L. lactis* pair may carry residual artifact that has "
         "nothing to do with the modality contrast."),
 # --- B8 Limitation 2 ---
 ("B8b", "Under the harshest cut, 4 of 45 DRAFTS pairs fall below 0.258, all involving "
         "*L. lactis*.",
         "Under the harshest cut, 4 of 45 DRAFTS pairs fall below 0.258, all involving "
         "*L. lactis* — a species for which Yim et al. excluded ScrFI/CCNGG-containing "
         "sequences as restriction artifacts, so those four pairs may be artifact-affected."),
 # --- B9 ---
 ("B9", "`PREREGISTRATION.md`, committed 2026-08-04, before any host-conditioned training.",
        "`PREREGISTRATION.md`, committed 2026-08-04, before any host-conditioned training. The "
        "commit timestamp is independently verifiable in the public repository's git history."),
]


def main() -> int:
    doc = PREPRINT.read_text(encoding="utf-8")
    before = doc
    applied = []
    for label, old, new in EDITS:
        if new is None:
            new = strata_sentence()
        n = doc.count(old)
        if n != 1:
            print(f"FAIL: {label} matched {n} times: {old[:70]!r}")
            return 1
        doc = doc.replace(old, new, 1)
        applied.append((label, len(new.split()) - len(old.split())))

    # no number may vanish
    NUM = re.compile(r"\d(?:[\d,]*\d)?(?:\.\d+)?")
    strip = lambda t: re.sub(r"^\d+\. ", "", t, flags=re.M)
    a, b = set(NUM.findall(strip(before))), set(NUM.findall(strip(doc)))
    gone = sorted(a - b)
    if gone:
        print(f"FAIL: numeric tokens would disappear: {gone}")
        return 1

    PREPRINT.write_text(doc, encoding="utf-8")
    MIRROR.write_text(doc, encoding="utf-8")
    body = lambda t: len(re.sub(r"^!\[.*$", "", t.split("## References")[0].split("# Cell-free")[-1],
                                flags=re.M).split())
    for label, d in applied:
        print(f"  {label}: {d:+d} words")
    print(f"body words {body(before)} -> {body(doc)}")
    print(f"new numeric tokens added: {sorted(b - a) or 'none'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
