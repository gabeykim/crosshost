"""
GATE 5.5 - Task 1: four-way comparison -- sequence-only / free-embedding /
genomic / physiology -- across every configuration, fold-resolved with 90%
percentile bootstrap intervals, matching Gate 5's protocol exactly.

Evaluation points compared (where each system has a defined value):
  - N=0 zero-shot (LOHO base): sequence-only, genomic, physiology, free-embedding (B3)
  - N=100 frozen-trunk (primary mechanism): sequence-only, genomic, physiology
  - N=100 top-conv (supplementary): sequence-only, genomic, physiology
  - RS241 zero-shot (both configs): sequence-only, genomic, physiology, free-embedding (B3, primary hosts only)
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


def seqonly_zeroshot_draws(host, readout, metric):
    with open(OUT / "gate5_5_seqonly_loho_results.json") as f:
        d = json.load(f)
    out = []
    for f_ in range(5):
        val = d[host]["per_fold"][str(f_)]["eval"][readout][metric]
        out.append(np.array([val]) if val is not None else np.array([]))
    return out


def seqonly_n100_draws(host, readout, mechanism, metric):
    if mechanism == "head_only":
        with open(OUT / "gate5_5_seqonly_calibration_curves.json") as f:
            d = json.load(f)
        out = []
        for f_ in range(5):
            draws = d[host][str(f_)]["N100_head_only"]
            out.append(np.array([x["eval"][readout][metric] for x in draws if x["eval"][readout][metric] is not None]))
        return out
    else:
        with open(OUT / "gate5_5_seqonly_n100_topconv.json") as f:
            d = json.load(f)
        out = []
        for f_ in range(5):
            draws = d[host][str(f_)]
            out.append(np.array([x["eval"][readout][metric] for x in draws if x["eval"][readout][metric] is not None]))
        return out


def genomic_physiology_zeroshot_draws(host, variant, readout, metric):
    with open(OUT / "gate4_loho_results.json") as f:
        d = json.load(f)
    out = []
    for f_ in range(5):
        val = d[host][variant]["per_fold"][str(f_)]["eval"][readout][metric]
        out.append(np.array([val]) if val is not None else np.array([]))
    return out


def genomic_physiology_n100_draws(host, variant, readout, mechanism, metric):
    if mechanism == "head_only":
        with open(OUT / "gate4_calibration_curves.json") as f:
            d = json.load(f)
        out = []
        for f_ in range(5):
            draws = d[host][variant][str(f_)]["N100_head_only"]
            out.append(np.array([x["eval"][readout][metric] for x in draws if x["eval"][readout][metric] is not None]))
        return out
    else:
        with open(OUT / "gate4_5_n100_topconv_supplement.json") as f:
            d = json.load(f)
        out = []
        for f_ in range(5):
            draws = d[variant][host][str(f_)]
            out.append(np.array([x["eval"][readout][metric] for x in draws if x["eval"][readout][metric] is not None]))
        return out


def b3_draws(host, readout, metric):
    with open(OUT / "baselines" / "fold_rotation_all_baselines.json") as f:
        d = json.load(f)
    per_fold = d["b3_free_host_embedding"]["summary"][host][readout]["per_fold"]
    return [np.array([pf[metric]]) for pf in per_fold]


def main():
    rows = []
    print("=" * 80)
    print("FOUR-WAY COMPARISON: sequence-only / free-embedding / genomic / physiology")
    print("=" * 80)
    for host in HOSTS:
        for readout in READOUTS:
            # N=0 zero-shot, all 4 systems
            for system, draws_fn in [
                ("sequence_only", lambda: seqonly_zeroshot_draws(host, readout, "spearman_rho")),
                ("free_embedding_B3", lambda: b3_draws(host, readout, "spearman_rho")),
                ("genomic", lambda: genomic_physiology_zeroshot_draws(host, "genomic", readout, "spearman_rho")),
                ("physiology", lambda: genomic_physiology_zeroshot_draws(host, "physiology", readout, "spearman_rho")),
            ]:
                ci = boot.bootstrap_ci_90(draws_fn())
                rows.append({"host": host, "readout": readout, "eval_point": "N0_zeroshot", "mechanism": "n/a",
                             "system": system, "rho_mean": ci["mean"], "rho_lower90": ci["lower"], "rho_upper90": ci["upper"]})
                print(f"  {host} {readout} N=0 {system}: rho={ci['mean']} [{ci['lower']},{ci['upper']}]" if ci["mean"] is not None else f"  {host} {readout} N=0 {system}: n/a")

            # N=100, both mechanisms, 3 systems (no B3 -- not defined at N=100)
            for mechanism in ["head_only", "top_conv"]:
                for system, draws_fn in [
                    ("sequence_only", lambda mech=mechanism: seqonly_n100_draws(host, readout, mech, "spearman_rho")),
                    ("genomic", lambda mech=mechanism: genomic_physiology_n100_draws(host, "genomic", readout, mech, "spearman_rho")),
                    ("physiology", lambda mech=mechanism: genomic_physiology_n100_draws(host, "physiology", readout, mech, "spearman_rho")),
                ]:
                    ci = boot.bootstrap_ci_90(draws_fn())
                    rows.append({"host": host, "readout": readout, "eval_point": "N100", "mechanism": mechanism,
                                 "system": system, "rho_mean": ci["mean"], "rho_lower90": ci["lower"], "rho_upper90": ci["upper"]})
                    print(f"  {host} {readout} N=100 {mechanism} {system}: rho={ci['mean']} [{ci['lower']},{ci['upper']}]" if ci["mean"] is not None else f"  {host} {readout} N=100 {mechanism} {system}: n/a")

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "gate5_5_fourway_comparison.csv", index=False)
    with open(RESULTS / "gate5_5_fourway_comparison.json", "w") as f:
        json.dump(rows, f, indent=2, default=str)
    print(f"\nWrote {RESULTS / 'gate5_5_fourway_comparison.csv'} ({len(df)} rows)")

    # verdict computation: for each (host, readout, eval_point, mechanism), is sequence-only
    # distinguishable from genomic/physiology (90% CI)?
    print("\n" + "=" * 80)
    print("VERDICT: is sequence-only distinguishable from genomic/physiology?")
    print("=" * 80)
    verdict_rows = []
    for host in HOSTS:
        for readout in READOUTS:
            for eval_point, mech_list in [("N0_zeroshot", ["n/a"]), ("N100", ["head_only", "top_conv"])]:
                for mechanism in mech_list:
                    sub = df[(df.host == host) & (df.readout == readout) & (df.eval_point == eval_point) & (df.mechanism == mechanism)]
                    seq_row = sub[sub.system == "sequence_only"]
                    if len(seq_row) == 0 or seq_row.rho_mean.isna().all():
                        continue
                    seq_ci = {"mean": seq_row.rho_mean.values[0], "lower": seq_row.rho_lower90.values[0], "upper": seq_row.rho_upper90.values[0]}
                    for other_system in ["genomic", "physiology"]:
                        other_row = sub[sub.system == other_system]
                        if len(other_row) == 0 or other_row.rho_mean.isna().all():
                            continue
                        other_ci = {"mean": other_row.rho_mean.values[0], "lower": other_row.rho_lower90.values[0], "upper": other_row.rho_upper90.values[0]}
                        no_overlap = (seq_ci["lower"] > other_ci["upper"]) or (other_ci["lower"] > seq_ci["upper"])
                        winner = "sequence_only" if seq_ci["mean"] > other_ci["mean"] else other_system
                        verdict_rows.append({"host": host, "readout": readout, "eval_point": eval_point, "mechanism": mechanism,
                                              "vs_system": other_system, "seq_only_mean": seq_ci["mean"], "other_mean": other_ci["mean"],
                                              "distinguishable": no_overlap, "point_winner": winner})
                        print(f"  {host} {readout} {eval_point} {mechanism} seqonly({seq_ci['mean']:.3f}) vs {other_system}({other_ci['mean']:.3f}): "
                              f"distinguishable={no_overlap}, winner={winner}")
    vdf = pd.DataFrame(verdict_rows)
    vdf.to_csv(RESULTS / "gate5_5_seqonly_verdict_table.csv", index=False)
    n_distinguishable = vdf["distinguishable"].sum()
    n_total = len(vdf)
    n_conditioning_wins = ((vdf["distinguishable"]) & (vdf["point_winner"] != "sequence_only")).sum()
    n_seqonly_wins = ((vdf["distinguishable"]) & (vdf["point_winner"] == "sequence_only")).sum()
    print(f"\nOf {n_total} comparisons: {n_distinguishable} distinguishable ({n_conditioning_wins} conditioning wins, {n_seqonly_wins} sequence-only wins)")
    with open(RESULTS / "gate5_5_seqonly_verdict_summary.json", "w") as f:
        json.dump({"n_total": int(n_total), "n_distinguishable": int(n_distinguishable),
                    "n_conditioning_wins": int(n_conditioning_wins), "n_seqonly_wins": int(n_seqonly_wins)}, f, indent=2)
    print(f"Wrote {RESULTS / 'gate5_5_seqonly_verdict_table.csv'}")


if __name__ == "__main__":
    main()
