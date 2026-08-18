"""
GATE 8.5 - Task 1A: formal comparison of the two alternative conditioning
mechanisms (concatenation, per-host-heads) against sequence-only and FiLM,
with 90% percentile bootstrap intervals, matching scripts/60's protocol
exactly (same bootstrap utility, same fold-draws convention: each fold's
zero-shot rho is a length-1 draw array).
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


def draws_from(path, host, readout, metric, key_path):
    with open(path) as f:
        d = json.load(f)
    out = []
    for f_ in range(5):
        node = d[host]["per_fold"][str(f_)]
        for k in key_path:
            node = node[k]
        val = node[readout][metric]
        out.append(np.array([val]) if val is not None else np.array([]))
    return out


def main():
    with open(OUT / "gate4_loho_results.json") as f:
        film_d = json.load(f)

    def film_genomic_draws(host, readout):
        out = []
        for f_ in range(5):
            val = film_d[host]["genomic"]["per_fold"][str(f_)]["eval"][readout]["spearman_rho"]
            out.append(np.array([val]) if val is not None else np.array([]))
        return out

    rows = []
    for host in HOSTS:
        for readout in READOUTS:
            systems = {
                "sequence_only": lambda h=host, r=readout: draws_from(
                    OUT / "gate5_5_seqonly_loho_results.json", h, r, "spearman_rho", ["eval"]),
                "film_genomic": lambda h=host, r=readout: film_genomic_draws(h, r),
                "concat": lambda h=host, r=readout: draws_from(
                    OUT / "gate8_5_concat_loho_results.json", h, r, "spearman_rho", ["eval"]),
                "perhost_heads_avg": lambda h=host, r=readout: draws_from(
                    OUT / "gate8_5_perhost_heads_loho_results.json", h, r, "spearman_rho", ["eval_avg"]),
                "perhost_heads_nearest": lambda h=host, r=readout: draws_from(
                    OUT / "gate8_5_perhost_heads_loho_results.json", h, r, "spearman_rho", ["eval_nearest"]),
            }

            for system, draws_fn in systems.items():
                draws = draws_fn()
                ci = boot.bootstrap_ci_90(draws)
                rows.append({"host": host, "readout": readout, "system": system,
                             "rho_mean": ci["mean"], "rho_lower90": ci["lower"], "rho_upper90": ci["upper"],
                             "n_folds": ci["n_folds"]})
                print(f"  {host} {readout} {system}: rho={ci['mean']} [{ci['lower']},{ci['upper']}]"
                      if ci["mean"] is not None else f"  {host} {readout} {system}: n/a")

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "gate8_5_conditioning_mechanisms.csv", index=False)
    with open(RESULTS / "gate8_5_conditioning_mechanisms.json", "w") as f:
        json.dump(rows, f, indent=2, default=str)
    print(f"\nWrote {RESULTS / 'gate8_5_conditioning_mechanisms.csv'} ({len(df)} rows)")

    print("\n" + "=" * 80)
    print("VERDICT: does any alternative mechanism beat sequence-only AND FiLM distinguishably?")
    print("=" * 80)
    verdict_rows = []
    for host in HOSTS:
        for readout in READOUTS:
            sub = df[(df.host == host) & (df.readout == readout)]
            seq_row = sub[sub.system == "sequence_only"]
            film_row = sub[sub.system == "film_genomic"]
            if len(seq_row) == 0 or seq_row.rho_mean.isna().all():
                continue
            seq_ci = {"mean": seq_row.rho_mean.values[0], "lower": seq_row.rho_lower90.values[0], "upper": seq_row.rho_upper90.values[0]}
            film_ci = {"mean": film_row.rho_mean.values[0], "lower": film_row.rho_lower90.values[0], "upper": film_row.rho_upper90.values[0]} if len(film_row) else None
            for mech in ["concat", "perhost_heads_avg", "perhost_heads_nearest"]:
                mech_row = sub[sub.system == mech]
                if len(mech_row) == 0 or mech_row.rho_mean.isna().all():
                    continue
                mech_ci = {"mean": mech_row.rho_mean.values[0], "lower": mech_row.rho_lower90.values[0], "upper": mech_row.rho_upper90.values[0]}
                beats_seqonly = mech_ci["lower"] > seq_ci["upper"]
                beats_film = (film_ci is not None) and (mech_ci["lower"] > film_ci["upper"])
                verdict_rows.append({
                    "host": host, "readout": readout, "mechanism": mech,
                    "mech_mean": mech_ci["mean"], "seqonly_mean": seq_ci["mean"],
                    "film_mean": film_ci["mean"] if film_ci else None,
                    "beats_seqonly_distinguishably": bool(beats_seqonly),
                    "beats_film_distinguishably": bool(beats_film),
                    "beats_both_distinguishably": bool(beats_seqonly and beats_film),
                })
                film_mean_str = f"{film_ci['mean']:.3f}" if film_ci else "n/a"
                print(f"  {host} {readout} {mech}: mean={mech_ci['mean']:.3f} vs seqonly={seq_ci['mean']:.3f} "
                      f"(beats={beats_seqonly}) vs film={film_mean_str} (beats={beats_film})")

    vdf = pd.DataFrame(verdict_rows)
    vdf.to_csv(RESULTS / "gate8_5_conditioning_mechanisms_verdict.csv", index=False)
    n_wins = int(vdf["beats_both_distinguishably"].sum())
    n_total = len(vdf)
    print(f"\nOf {n_total} (host, readout, mechanism) comparisons: {n_wins} beat BOTH sequence-only AND FiLM distinguishably")
    winning = vdf[vdf["beats_both_distinguishably"]]
    if len(winning):
        print(winning.to_string(index=False))
    with open(RESULTS / "gate8_5_conditioning_mechanisms_verdict_summary.json", "w") as f:
        json.dump({"n_total": n_total, "n_wins_vs_both": n_wins,
                    "winning_cells": winning.to_dict("records")}, f, indent=2, default=str)
    print(f"Wrote {RESULTS / 'gate8_5_conditioning_mechanisms_verdict.csv'}")


if __name__ == "__main__":
    main()
