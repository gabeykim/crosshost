"""
GATE 8 - Task 2c: consolidated master baseline table for shipping. Extends
scripts/64's six-system comparison (sequence_only, free_embedding_B3,
genomic, physiology, dnabert2, promogen2) with B1 (mean/majority) and B4
(biophysical) at their natural single evaluation point, and B2 (per-host-N)
at N=100 and N=3000, all via the identical bootstrap_ci_90 protocol used
throughout the project -- so every number in this table is directly
comparable to every other number and to every number in PAPER_FRAMING.md.
"""
import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
boot = import_module("51_bootstrap_utils")
comparison64 = import_module("64_gate6_full_comparison")

OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"
PKG = Path(__file__).resolve().parent.parent / "package"
BASELINES_DIR = PKG / "baselines"
BASELINES_DIR.mkdir(parents=True, exist_ok=True)

HOSTS = ["EC", "BS", "PA"]
READOUTS = ["transcription", "translation"]


def b1_draws(host, readout):
    d = json.load(open(OUT / "baselines" / "fold_rotation_all_baselines.json"))
    pf = d["b1_mean_majority"]["summary"][host][readout]["per_fold"]
    return [np.array([x["spearman_rho"]]) for x in pf if x["spearman_rho"] is not None]


def b4_draws(host, readout):
    d = json.load(open(OUT / "baselines" / "fold_rotation_all_baselines.json"))
    pf = d["b4_biophysical"]["summary"][host][readout]["per_fold"]
    return [np.array([x["spearman_rho"]]) for x in pf if x["spearman_rho"] is not None]


def b2_draws(host, readout, n):
    d = json.load(open(OUT / "baselines" / "fold_rotation_b2_raw.json"))
    out = []
    for f in range(5):
        draws = d[host][readout][str(n)][str(f)]["draws"]
        out.append(np.array([x["spearman_rho"] for x in draws if x["spearman_rho"] is not None]))
    return out


def main():
    rows = []
    for host in HOSTS:
        for readout in READOUTS:
            ci = boot.bootstrap_ci_90(b1_draws(host, readout))
            rows.append({"host": host, "readout": readout, "eval_point": "single", "mechanism": "n/a",
                         "system": "B1_mean_majority", "rho_mean": ci["mean"], "rho_lower90": ci["lower"], "rho_upper90": ci["upper"]})

            ci = boot.bootstrap_ci_90(b4_draws(host, readout))
            rows.append({"host": host, "readout": readout, "eval_point": "single", "mechanism": "n/a",
                         "system": "B4_biophysical", "rho_mean": ci["mean"], "rho_lower90": ci["lower"], "rho_upper90": ci["upper"]})

            for n in [100, 3000]:
                ci = boot.bootstrap_ci_90(b2_draws(host, readout, n))
                rows.append({"host": host, "readout": readout, "eval_point": f"N{n}", "mechanism": "per_host_trained",
                             "system": "B2_per_host_N", "rho_mean": ci["mean"], "rho_lower90": ci["lower"], "rho_upper90": ci["upper"]})

    b1b2b4_df = pd.DataFrame(rows)

    six_way = pd.read_csv(RESULTS / "gate6_full_comparison.csv")

    master = pd.concat([b1b2b4_df, six_way], ignore_index=True)
    master.to_csv(BASELINES_DIR / "master_baselines.csv", index=False)
    with open(BASELINES_DIR / "master_baselines.json", "w") as f:
        json.dump(master.to_dict(orient="records"), f, indent=2, default=str)
    print(f"Wrote {BASELINES_DIR / 'master_baselines.csv'} ({len(master)} rows, "
          f"{master.system.nunique()} systems: {sorted(master.system.unique())})")

    # Gate 7 conformal/ECE, shipped alongside, prominently (per Task 2c instruction)
    conformal = json.load(open(RESULTS / "gate7_conformal_calibration.json"))
    with open(BASELINES_DIR / "gate7_conformal_and_ece.json", "w") as f:
        json.dump(conformal, f, indent=2, default=str)
    print(f"Wrote {BASELINES_DIR / 'gate7_conformal_and_ece.json'}")


if __name__ == "__main__":
    main()
