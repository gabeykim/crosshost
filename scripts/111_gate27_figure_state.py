#!/usr/bin/env python3
"""Gate 27 Task 1: figure ground truth.

Establishes, for each of the nine preprint figures, whether the file on disk is
current and whether the rendered PDF agrees with the embedded caption. Content
rows for Figures 2b, 3, 5 and 8 were established by extracting the images from
manuscript.pdf with `pdfimages -png` and inspecting them directly -- not by
re-reading the source PNGs, which would not have settled the question.
"""
from __future__ import annotations

import csv
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "out" / "results" / "gate27_figure_state.csv"
FIGS = ROOT / "out" / "PREPRINT" / "figures"

# figure, producing script, caption promises, figure shows, verdict, how verified
ROWS = [
 ("Figure 1", "scripts/73_attenuation_figures.py (via 106)",
  "three reliabilities, EC-PA bars clipped at rho=1 with uncapped values 1.51/1.08 and 1.48/1.06 printed",
  "as promised",
  "MATCH", "byte-identical to the Gate 25 commit; inspected in Gate 25"),
 ("Figure 2a", "scripts/96_gate10_5_figures.py (via 106)",
  "raw vs GC-controlled, three pairs, both readouts",
  "as promised",
  "MATCH", "byte-identical to the Gate 25 commit; inspected in Gate 25"),
 ("Figure 2b", "scripts/96_gate10_5_figures.py (via 106)",
  "five strata (Proteobacteria, Firmicutes, Actinobacteria, Bacteroidetes, Cyanobacteria), per-stratum N on every bar",
  "five strata, N printed on all 15 bars",
  "MATCH", "extracted from manuscript.pdf (img-004) and inspected this gate"),
 ("Figure 3", "scripts/106_gate25_regenerate_figures.py",
  "panel A four bars (0.616/862, 0.655/14,088, 0.677/807, 0.258/3,668); panel B RS234, seven species",
  "panel A four bars with exactly those values and N; no 0.597; no n=15 bar; panel B present, seven species with N",
  "MATCH", "extracted from manuscript.pdf (img-006) and inspected this gate"),
 ("Figure 4", "scripts/55_gate5_figures.py (via 106)",
  "eight primary H-MAIN comparisons, 90% bootstrap CI",
  "as promised",
  "MATCH", "byte-identical to the Gate 25 commit; inspected in Gate 20"),
 ("Figure 5", "scripts/65_gate6_figures.py (via 106)",
  "nine systems including the three conditioning mechanisms, six cells, one shared y-scale",
  "nine systems, six panels, y-axis shared (identical -0.2 to 0.7 range, ticks on leftmost only); EC transcription shows concatenation 0.555 and per-host-heads-avg 0.615",
  "MATCH", "extracted from manuscript.pdf (img-010) and inspected this gate"),
 ("Figure 6", "scripts/68_conformal_calibration.py (via 106)",
  "reliability diagrams, three hosts, both readouts",
  "as promised",
  "MATCH", "byte-identical to the Gate 25 commit; inspected in Gate 20"),
 ("Figure 7", "scripts/87_gate8_6_figures.py (via 106)",
  "fold-to-fold std, four systems, six cells",
  "as promised",
  "MATCH", "byte-identical to the Gate 25 commit; inspected in Gate 20"),
 ("Figure 8", "scripts/87_gate8_6_figures.py (via 106)",
  "11 bar groups, tx/tl distinguished, title stating both directions",
  "title: 'Reference-value alone matches or beats the original model in 10 of 11; 1 of 11 survives'; 11 bar groups; labels tx/tl distinct",
  "MATCH", "extracted from manuscript.pdf (img-016) and inspected this gate"),
]

HEADER = ["figure", "producing_script", "caption_promises", "figure_shows", "verdict",
          "how_verified", "file_mtime", "mtime_vs_commit_3cfbca6", "pdf_newer_than_figure",
          "script_run_since_gate26"]

FILES = {"Figure 1": "Figure1_disattenuation.png", "Figure 2a": "Figure2a_gc_control.png",
         "Figure 2b": "Figure2b_phylum_stratified.png",
         "Figure 3": "Figure3_drafts_modality_comparison.png",
         "Figure 4": "Figure4_hmain_kill_gate.png",
         "Figure 5": "Figure5_conditioning_mechanisms.png",
         "Figure 6": "Figure6_calibration_collapse.png",
         "Figure 7": "Figure7_film_instability.png",
         "Figure 8": "Figure8_shift_prediction_retraction.png"}


def main() -> int:
    commit_t = int(subprocess.run(["git", "show", "-s", "--format=%ct", "3cfbca6"],
                                  cwd=ROOT, capture_output=True, text=True).stdout.strip())
    pdf_t = (ROOT / "manuscript.pdf").stat().st_mtime
    import datetime
    rows = []
    for fig, script, promises, shows, verdict, how in ROWS:
        st = (FIGS / FILES[fig]).stat().st_mtime
        rows.append([fig, script, promises, shows, verdict, how,
                     datetime.datetime.fromtimestamp(st).strftime("%Y-%m-%d %H:%M:%S"),
                     "predates" if st < commit_t else "postdates",
                     "yes" if pdf_t > st else "NO",
                     "yes, during Gate 26 via scripts/106; not since (nothing changed)"])
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(HEADER)
        w.writerows(rows)
    print(f"Wrote {OUT}: {len(rows)} rows")
    print("verdicts:", {r[4] for r in rows})
    print("all figures predate the commit:", all(r[7] == "predates" for r in rows))
    print("pdf newer than every figure:", all(r[8] == "yes" for r in rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
