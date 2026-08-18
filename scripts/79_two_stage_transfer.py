"""
GATE 8.5 - Task 1D: separate the two stages for transfer.

The 4-head joint architecture (scripts/40, scripts/58) always trained an
active/inactive CLASSIFIER and a strength REGRESSOR jointly, and every prior
gate's headline numbers (H-MAIN, the four-way comparison, Gate 6/7) reported
them together or led with the regression Spearman rho. This script reuses
the EXISTING zero-shot LOHO results already on disk -- no new training --
and simply separates the two stages that scripts/42's evaluate() function
already computes independently (mcc/auc for the classifier, spearman_rho for
the regressor on actives) into their own bootstrap-CI'd comparison, for
sequence-only vs. both FiLM variants (genomic, physiology).

Sources (all already-computed, pre-existing files, zero re-training):
  - out/gate5_5_seqonly_loho_results.json  (sequence-only, 3 hosts x 5 folds)
  - out/gate4_loho_results.json            (FiLM genomic + physiology, 3 hosts x 5 folds)
"""
import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
boot = import_module("51_bootstrap_utils")

OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

HOSTS = ["EC", "BS", "PA"]
READOUTS = ["transcription", "translation"]
CLASSIFICATION_METRICS = ["mcc", "auc"]
REGRESSION_METRIC = "spearman_rho"


def seqonly_draws(host, readout, metric):
    with open(OUT / "gate5_5_seqonly_loho_results.json") as f:
        d = json.load(f)
    out = []
    for f_ in range(5):
        val = d[host]["per_fold"][str(f_)]["eval"][readout][metric]
        out.append(np.array([val]) if val is not None else np.array([]))
    return out


def film_draws(host, variant, readout, metric):
    with open(OUT / "gate4_loho_results.json") as f:
        d = json.load(f)
    out = []
    for f_ in range(5):
        val = d[host][variant]["per_fold"][str(f_)]["eval"][readout][metric]
        out.append(np.array([val]) if val is not None else np.array([]))
    return out


def main():
    rows = []
    for host in HOSTS:
        for readout in READOUTS:
            for metric in CLASSIFICATION_METRICS + [REGRESSION_METRIC]:
                stage = "classification" if metric in CLASSIFICATION_METRICS else "regression"
                systems = {
                    "sequence_only": seqonly_draws(host, readout, metric),
                    "film_genomic": film_draws(host, "genomic", readout, metric),
                    "film_physiology": film_draws(host, "physiology", readout, metric),
                }
                for system, draws in systems.items():
                    ci = boot.bootstrap_ci_90(draws)
                    rows.append({
                        "host": host, "readout": readout, "stage": stage, "metric": metric,
                        "system": system, "mean": ci["mean"], "lower90": ci["lower"], "upper90": ci["upper"],
                        "n_folds": ci["n_folds"],
                    })

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "gate8_5_two_stage_transfer.csv", index=False)

    # verdict: for each (host, readout, stage), is FiLM (either variant) distinguishably
    # better than sequence-only? Use auc for classification (mcc is noisier / more often
    # degenerate at small n_active_eval), spearman_rho for regression.
    verdict_rows = []
    for host in HOSTS:
        for readout in READOUTS:
            for stage, metric in [("classification", "auc"), ("regression", "spearman_rho")]:
                sub = df[(df.host == host) & (df.readout == readout) & (df.stage == stage) & (df.metric == metric)]
                seq_row = sub[sub.system == "sequence_only"]
                if len(seq_row) == 0 or seq_row["mean"].isna().all():
                    continue
                seq_ci = {"mean": seq_row["mean"].values[0], "lower": seq_row["lower90"].values[0], "upper": seq_row["upper90"].values[0]}
                for variant in ["film_genomic", "film_physiology"]:
                    other_row = sub[sub.system == variant]
                    if len(other_row) == 0 or other_row["mean"].isna().all():
                        continue
                    other_ci = {"mean": other_row["mean"].values[0], "lower": other_row["lower90"].values[0], "upper": other_row["upper90"].values[0]}
                    no_overlap = (seq_ci["lower"] > other_ci["upper"]) or (other_ci["lower"] > seq_ci["upper"])
                    winner = "sequence_only" if seq_ci["mean"] > other_ci["mean"] else variant
                    verdict_rows.append({
                        "host": host, "readout": readout, "stage": stage, "metric": metric, "vs_variant": variant,
                        "seq_only_mean": seq_ci["mean"], "film_mean": other_ci["mean"],
                        "distinguishable": no_overlap, "point_winner": winner,
                        "film_wins_distinguishably": bool(no_overlap and winner == variant),
                    })
    vdf = pd.DataFrame(verdict_rows)
    vdf.to_csv(RESULTS / "gate8_5_two_stage_transfer_verdict.csv", index=False)

    n_class_film_wins = int(vdf[(vdf.stage == "classification")]["film_wins_distinguishably"].sum())
    n_reg_film_wins = int(vdf[(vdf.stage == "regression")]["film_wins_distinguishably"].sum())
    n_class_total = int((vdf.stage == "classification").sum())
    n_reg_total = int((vdf.stage == "regression").sum())

    summary = {
        "n_classification_comparisons": n_class_total,
        "n_classification_film_wins_distinguishably": n_class_film_wins,
        "n_regression_comparisons": n_reg_total,
        "n_regression_film_wins_distinguishably": n_reg_film_wins,
        "stage_specific_finding": (n_class_film_wins != n_reg_film_wins),
    }
    with open(RESULTS / "gate8_5_two_stage_transfer_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("=" * 80)
    print("TASK 1D: TWO-STAGE TRANSFER -- classification vs regression, seq-only vs FiLM")
    print("=" * 80)
    print(vdf.to_string(index=False))
    print()
    print(json.dumps(summary, indent=2))
    print(f"\nWrote {RESULTS / 'gate8_5_two_stage_transfer.csv'}, "
          f"{RESULTS / 'gate8_5_two_stage_transfer_verdict.csv'}, "
          f"{RESULTS / 'gate8_5_two_stage_transfer_summary.json'}")


if __name__ == "__main__":
    main()
