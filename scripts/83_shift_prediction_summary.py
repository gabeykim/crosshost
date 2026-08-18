"""
GATE 8.5 - Task 1B summary: bootstrap-CI'd (pair, readout, variant) table
from scripts/81's raw per-fold results. Separated into its own committed
script (rather than left as an ad-hoc analysis) so
out/results/gate8_5_shift_prediction_summary.csv has a real producing
script, per scripts/audit_provenance.py.
"""
import json
import numpy as np
import pandas as pd
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
boot = import_module("51_bootstrap_utils")

OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"


def main():
    d = json.load(open(RESULTS / "gate8_5_shift_prediction.json"))
    rows = []
    for pair in d:
        for ro in d[pair]:
            node = d[pair][ro]
            if "per_fold" not in node or len(node["per_fold"]) == 0:
                continue
            for variant in ["no_condition", "with_condition"]:
                rho_draws = [np.array([node["per_fold"][f][variant]["spearman_rho"]])
                             for f in node["per_fold"] if node["per_fold"][f][variant]["spearman_rho"] is not None]
                fve_zero_draws = [np.array([node["per_fold"][f][variant]["frac_var_explained_vs_zero"]])
                                   for f in node["per_fold"]]
                fve_mean_draws = [np.array([node["per_fold"][f][variant]["frac_var_explained_vs_mean"]])
                                   for f in node["per_fold"]]
                ci_rho = boot.bootstrap_ci_90(rho_draws)
                ci_fvz = boot.bootstrap_ci_90(fve_zero_draws)
                ci_fvm = boot.bootstrap_ci_90(fve_mean_draws)
                rows.append({
                    "pair": pair, "readout": ro, "variant": variant, "n_total": node["n_total"],
                    "rho_mean": ci_rho["mean"], "rho_lower90": ci_rho["lower"], "rho_upper90": ci_rho["upper"],
                    "fve_zero_mean": ci_fvz["mean"], "fve_zero_lower90": ci_fvz["lower"], "fve_zero_upper90": ci_fvz["upper"],
                    "fve_mean_mean": ci_fvm["mean"], "fve_mean_lower90": ci_fvm["lower"], "fve_mean_upper90": ci_fvm["upper"],
                })
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "gate8_5_shift_prediction_summary.csv", index=False)
    print(df.round(3).to_string(index=False))
    print(f"\nWrote {RESULTS / 'gate8_5_shift_prediction_summary.csv'} ({len(df)} rows)")


if __name__ == "__main__":
    main()
