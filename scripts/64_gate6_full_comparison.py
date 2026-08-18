"""
GATE 6 - Task 3: the definitive comparison -- sequence-only / free-embedding /
genomic / physiology / DNABERT-2 / PromoGen2 [/ Evo 2, if available] -- across
every host, readout, and configuration, fold-resolved with 90% percentile
bootstrap intervals. Matches Gate 5/5.5 protocol exactly (scripts/51, /60).

CORRECTION TO A GATE 5.5 NUMBER, found and fixed while building this script:
out/results/gate5_5_per_host_ceiling.json's "actual_zeroshot_mean" field (used
to compute the widely-cited "49.2% / 50.7% / 44.8%" transcription
percent-of-ceiling figures in out/GATE5_5_MEMO.md and the Gate 6 prompt
itself) does NOT match the project's own officially bootstrapped sequence-
only zero-shot Spearman rho in 5 of 6 (host, readout) cells -- verified by
diffing against out/results/gate5_5_fourway_comparison.json, the file
GATE5_5_MEMO.md's own headline table is actually built from. No script for
gate5_5_per_host_ceiling.json exists anywhere in scripts/, meaning it was
produced by an unsaved ad-hoc computation -- exactly the failure mode the
project's own standing rules exist to catch (charter Part V rule 6: "verify
before trusting a stated number"). The cross-host measurement correlations
themselves (the ceiling values) were independently re-run via scripts/57 and
reproduce bit-for-bit -- only the "actual" side of the ratio was wrong.
This script recomputes percent-of-ceiling from ONLY verified, reproducible
sources (scripts/51's bootstrap CIs over the real per-fold spearman_rho
draws, and scripts/57's measurement correlations) and reports the corrected
figures. See out/GATE6_MEMO.md for the full writeup of this correction and
out/results/gate5_5_per_host_ceiling_CORRECTED.json for the fixed numbers.
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
FM_MODELS = ["dnabert2", "promogen2"]  # extend with "evo2" once/if available


# ---------- draw-extraction functions, one per system ----------

def seqonly_zeroshot_draws(host, readout, metric):
    d = json.load(open(OUT / "gate5_5_seqonly_loho_results.json"))
    out = []
    for f_ in range(5):
        val = d[host]["per_fold"][str(f_)]["eval"][readout][metric]
        out.append(np.array([val]) if val is not None else np.array([]))
    return out


def genomic_physiology_zeroshot_draws(host, variant, readout, metric):
    d = json.load(open(OUT / "gate4_loho_results.json"))
    out = []
    for f_ in range(5):
        val = d[host][variant]["per_fold"][str(f_)]["eval"][readout][metric]
        out.append(np.array([val]) if val is not None else np.array([]))
    return out


def b3_draws(host, readout, metric):
    d = json.load(open(OUT / "baselines" / "fold_rotation_all_baselines.json"))
    per_fold = d["b3_free_host_embedding"]["summary"][host][readout]["per_fold"]
    return [np.array([pf[metric]]) for pf in per_fold]


def fm_zeroshot_draws(model_tag, host, readout, metric):
    d = json.load(open(OUT / f"gate6_{model_tag}_loho_results.json"))
    out = []
    for f_ in range(5):
        val = d[host]["per_fold"][str(f_)]["eval"][readout][metric]
        out.append(np.array([val]) if val is not None else np.array([]))
    return out


def fm_n100_draws(model_tag, host, readout, metric):
    d = json.load(open(OUT / f"gate6_{model_tag}_calibration_curves.json"))
    out = []
    for f_ in range(5):
        draws = d[host][str(f_)]["N100_head_only"]
        out.append(np.array([x["eval"][readout][metric] for x in draws if x["eval"][readout][metric] is not None]))
    return out


def seqonly_n100_draws(host, readout, mechanism, metric):
    if mechanism == "head_only":
        d = json.load(open(OUT / "gate5_5_seqonly_calibration_curves.json"))
        out = []
        for f_ in range(5):
            draws = d[host][str(f_)]["N100_head_only"]
            out.append(np.array([x["eval"][readout][metric] for x in draws if x["eval"][readout][metric] is not None]))
        return out
    else:
        d = json.load(open(OUT / "gate5_5_seqonly_n100_topconv.json"))
        out = []
        for f_ in range(5):
            draws = d[host][str(f_)]
            out.append(np.array([x["eval"][readout][metric] for x in draws if x["eval"][readout][metric] is not None]))
        return out


def genomic_physiology_n100_draws(host, variant, readout, mechanism, metric):
    if mechanism == "head_only":
        d = json.load(open(OUT / "gate4_calibration_curves.json"))
        out = []
        for f_ in range(5):
            draws = d[host][variant][str(f_)]["N100_head_only"]
            out.append(np.array([x["eval"][readout][metric] for x in draws if x["eval"][readout][metric] is not None]))
        return out
    else:
        d = json.load(open(OUT / "gate4_5_n100_topconv_supplement.json"))
        out = []
        for f_ in range(5):
            draws = d[variant][host][str(f_)]
            out.append(np.array([x["eval"][readout][metric] for x in draws if x["eval"][readout][metric] is not None]))
        return out


SYSTEMS_ZEROSHOT = {
    "sequence_only": lambda h, r: seqonly_zeroshot_draws(h, r, "spearman_rho"),
    "free_embedding_B3": lambda h, r: b3_draws(h, r, "spearman_rho"),
    "genomic": lambda h, r: genomic_physiology_zeroshot_draws(h, "genomic", r, "spearman_rho"),
    "physiology": lambda h, r: genomic_physiology_zeroshot_draws(h, "physiology", r, "spearman_rho"),
    "dnabert2": lambda h, r: fm_zeroshot_draws("dnabert2", h, r, "spearman_rho"),
    "promogen2": lambda h, r: fm_zeroshot_draws("promogen2", h, r, "spearman_rho"),
}

SYSTEMS_N100 = {
    "sequence_only": lambda h, r, mech: seqonly_n100_draws(h, r, mech, "spearman_rho"),
    "genomic": lambda h, r, mech: genomic_physiology_n100_draws(h, "genomic", r, mech, "spearman_rho"),
    "physiology": lambda h, r, mech: genomic_physiology_n100_draws(h, "physiology", r, mech, "spearman_rho"),
    "dnabert2": lambda h, r, mech: (fm_n100_draws("dnabert2", h, r, "spearman_rho") if mech == "head_only" else None),
    "promogen2": lambda h, r, mech: (fm_n100_draws("promogen2", h, r, "spearman_rho") if mech == "head_only" else None),
}


def add_evo2_if_available():
    evo2_loho = OUT / "gate6_evo2_loho_results.json"
    evo2_calib = OUT / "gate6_evo2_calibration_curves.json"
    if evo2_loho.exists():
        FM_MODELS.append("evo2")
        SYSTEMS_ZEROSHOT["evo2"] = lambda h, r: fm_zeroshot_draws("evo2", h, r, "spearman_rho")
        if evo2_calib.exists():
            SYSTEMS_N100["evo2"] = lambda h, r, mech: (fm_n100_draws("evo2", h, r, "spearman_rho") if mech == "head_only" else None)
        return True
    return False


def build_comparison_table():
    rows = []
    for host in HOSTS:
        for readout in READOUTS:
            for system, fn in SYSTEMS_ZEROSHOT.items():
                ci = boot.bootstrap_ci_90(fn(host, readout))
                rows.append({"host": host, "readout": readout, "eval_point": "N0_zeroshot", "mechanism": "n/a",
                             "system": system, "rho_mean": ci["mean"], "rho_lower90": ci["lower"], "rho_upper90": ci["upper"]})
            for mechanism in ["head_only", "top_conv"]:
                for system, fn in SYSTEMS_N100.items():
                    draws = fn(host, readout, mechanism)
                    if draws is None:
                        continue
                    ci = boot.bootstrap_ci_90(draws)
                    rows.append({"host": host, "readout": readout, "eval_point": "N100", "mechanism": mechanism,
                                 "system": system, "rho_mean": ci["mean"], "rho_lower90": ci["lower"], "rho_upper90": ci["upper"]})
    return pd.DataFrame(rows)


def build_corrected_ceiling_table(df):
    """Percent-of-ceiling, recomputed from ONLY verified sources: scripts/57's
    measurement correlations (ceiling) and this script's own bootstrapped
    N0_zeroshot rho_mean (actual) -- see module docstring for why the
    previous out/results/gate5_5_per_host_ceiling.json is not trusted."""
    corr = json.load(open(OUT / "results" / "gate5_5_crosshost_measurement_correlation.json"))
    pair_map = {"EC": ["EC_BS", "EC_PA"], "BS": ["EC_BS", "BS_PA"], "PA": ["EC_PA", "BS_PA"]}
    systems_for_ceiling = ["sequence_only", "dnabert2", "promogen2"] + (["evo2"] if "evo2" in FM_MODELS else [])
    out = {}
    for host in HOSTS:
        out[host] = {}
        for readout in READOUTS:
            pairs = pair_map[host]
            ceiling_vals = [corr[readout][p]["rho"] for p in pairs]
            ceiling_max = max(ceiling_vals)
            ceiling_mean_pair = float(np.mean(ceiling_vals))
            out[host][readout] = {"ceiling_max": ceiling_max, "ceiling_mean_pair": ceiling_mean_pair, "systems": {}}
            for system in systems_for_ceiling:
                sub = df[(df.host == host) & (df.readout == readout) & (df.eval_point == "N0_zeroshot") & (df.system == system)]
                if len(sub) == 0 or sub.rho_mean.isna().all():
                    continue
                actual = float(sub.rho_mean.values[0])
                out[host][readout]["systems"][system] = {
                    "actual_zeroshot_rho": actual,
                    "pct_of_max_ceiling": 100 * actual / ceiling_max if ceiling_max else None,
                    "pct_of_mean_pair_ceiling": 100 * actual / ceiling_mean_pair if ceiling_mean_pair else None,
                }
    return out


RS241_CONFIGS = {"PRIMARY": ["EC", "BS", "PA"], "SECONDARY": ["EC", "PA"]}
RS241_HOSTS = ["SE", "VN", "CG"]


def seqonly_rs241_seed_draws(config, held_out, readout):
    d = json.load(open(OUT / "gate5_5_seqonly_rs241_results.json"))
    seed_draws = []
    for seed in range(5):
        entry = d[config]["per_seed"][str(seed)]["eval"][held_out][readout]
        if entry["spearman_rho"] is not None:
            seed_draws.append(np.array([entry["spearman_rho"]]))
    return seed_draws


def genomic_physiology_rs241_seed_draws(config, variant, held_out, readout):
    d = json.load(open(OUT / "gate5_rs241_results.json"))
    seed_draws = []
    for seed in range(5):
        entry = d[config][variant]["per_seed"][str(seed)]["eval"][held_out][readout]
        if entry["spearman_rho"] is not None:
            seed_draws.append(np.array([entry["spearman_rho"]]))
    return seed_draws


def fm_rs241_seed_draws(model_tag, config, held_out, readout):
    d = json.load(open(OUT / f"gate6_{model_tag}_rs241_results.json"))
    seed_draws = []
    for seed in range(5):
        entry = d[config]["per_seed"][str(seed)]["eval"][held_out][readout]
        if entry["spearman_rho"] is not None:
            seed_draws.append(np.array([entry["spearman_rho"]]))
    return seed_draws


def build_rs241_comparison_table():
    rows = []
    for config in RS241_CONFIGS:
        for held_out in RS241_HOSTS:
            for readout in READOUTS:
                systems = {
                    "sequence_only": seqonly_rs241_seed_draws(config, held_out, readout),
                    "genomic": genomic_physiology_rs241_seed_draws(config, "genomic", held_out, readout),
                    "physiology": genomic_physiology_rs241_seed_draws(config, "physiology", held_out, readout),
                    "dnabert2": fm_rs241_seed_draws("dnabert2", config, held_out, readout),
                    "promogen2": fm_rs241_seed_draws("promogen2", config, held_out, readout),
                }
                if "evo2" in FM_MODELS:
                    systems["evo2"] = fm_rs241_seed_draws("evo2", config, held_out, readout)
                for system, draws in systems.items():
                    ci = boot.bootstrap_ci_90(draws)
                    rows.append({"config": config, "held_out": held_out, "readout": readout, "system": system,
                                 "rho_mean": ci["mean"], "rho_lower90": ci["lower"], "rho_upper90": ci["upper"]})
    return pd.DataFrame(rows)


def main():
    had_evo2 = add_evo2_if_available()
    print(f"Evo2 results available: {had_evo2}. Systems: {list(SYSTEMS_ZEROSHOT.keys())}")

    df = build_comparison_table()
    df.to_csv(RESULTS / "gate6_full_comparison.csv", index=False)
    with open(RESULTS / "gate6_full_comparison.json", "w") as f:
        json.dump(df.to_dict(orient="records"), f, indent=2, default=str)
    print(f"Wrote gate6_full_comparison.{{csv,json}} ({len(df)} rows)")

    ceiling = build_corrected_ceiling_table(df)
    with open(RESULTS / "gate5_5_per_host_ceiling_CORRECTED.json", "w") as f:
        json.dump(ceiling, f, indent=2)
    print("\n=== CORRECTED percent-of-ceiling (supersedes gate5_5_per_host_ceiling.json) ===")
    for host in HOSTS:
        for readout in READOUTS:
            c = ceiling[host][readout]
            print(f"  {host} {readout}: ceiling_max={c['ceiling_max']:.3f}")
            for sysname, v in c["systems"].items():
                print(f"    {sysname}: rho={v['actual_zeroshot_rho']:.3f}, pct_of_max_ceiling={v['pct_of_max_ceiling']:.1f}%")
    print(f"\nWrote {RESULTS / 'gate5_5_per_host_ceiling_CORRECTED.json'}")

    # ---------- verdict: is each FM distinguishable from sequence-only? ----------
    print("\n" + "=" * 80)
    print("VERDICT: is each foundation model distinguishable from sequence-only?")
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
                    for fm in FM_MODELS:
                        fm_row = sub[sub.system == fm]
                        if len(fm_row) == 0 or fm_row.rho_mean.isna().all():
                            continue
                        fm_ci = {"mean": fm_row.rho_mean.values[0], "lower": fm_row.rho_lower90.values[0], "upper": fm_row.rho_upper90.values[0]}
                        no_overlap = (seq_ci["lower"] > fm_ci["upper"]) or (fm_ci["lower"] > seq_ci["upper"])
                        winner = "fm" if fm_ci["mean"] > seq_ci["mean"] else "sequence_only"
                        verdict_rows.append({"host": host, "readout": readout, "eval_point": eval_point, "mechanism": mechanism,
                                              "fm": fm, "fm_mean": fm_ci["mean"], "seqonly_mean": seq_ci["mean"],
                                              "distinguishable": bool(no_overlap), "point_winner": winner})
                        print(f"  {host} {readout} {eval_point} {mechanism} {fm}({fm_ci['mean']:.3f}) vs "
                              f"seqonly({seq_ci['mean']:.3f}): distinguishable={no_overlap}, winner={winner}")
    vdf = pd.DataFrame(verdict_rows)
    vdf.to_csv(RESULTS / "gate6_fm_verdict_table.csv", index=False)
    n_total = len(vdf)
    n_distinguishable = int(vdf["distinguishable"].sum())
    n_fm_wins = int(((vdf["distinguishable"]) & (vdf["point_winner"] == "fm")).sum())
    n_seqonly_wins = int(((vdf["distinguishable"]) & (vdf["point_winner"] == "sequence_only")).sum())
    print(f"\nOf {n_total} FM-vs-sequence-only comparisons: {n_distinguishable} distinguishable "
          f"({n_fm_wins} FM wins, {n_seqonly_wins} sequence-only wins)")
    with open(RESULTS / "gate6_fm_verdict_summary.json", "w") as f:
        json.dump({"n_total": n_total, "n_distinguishable": n_distinguishable,
                    "n_fm_wins": n_fm_wins, "n_seqonly_wins": n_seqonly_wins}, f, indent=2)
    print(f"Wrote {RESULTS / 'gate6_fm_verdict_table.csv'}")

    # ---------- RS241 ----------
    print("\n" + "=" * 80)
    print("RS241 zero-shot: sequence-only / genomic / physiology / FMs")
    print("=" * 80)
    rs_df = build_rs241_comparison_table()
    rs_df.to_csv(RESULTS / "gate6_rs241_comparison.csv", index=False)
    with open(RESULTS / "gate6_rs241_comparison.json", "w") as f:
        json.dump(rs_df.to_dict(orient="records"), f, indent=2, default=str)
    for _, r in rs_df.iterrows():
        print(f"  {r['config']} {r['held_out']} {r['readout']} {r['system']}: "
              f"rho={r['rho_mean']:.3f} [{r['rho_lower90']:.3f},{r['rho_upper90']:.3f}]" if r['rho_mean'] is not None
              else f"  {r['config']} {r['held_out']} {r['readout']} {r['system']}: n/a")
    print(f"Wrote {RESULTS / 'gate6_rs241_comparison.csv'}")

    print("\n=== RS241 VERDICT: FM vs sequence-only ===")
    rs_verdict_rows = []
    for config in RS241_CONFIGS:
        for held_out in RS241_HOSTS:
            for readout in READOUTS:
                sub = rs_df[(rs_df.config == config) & (rs_df.held_out == held_out) & (rs_df.readout == readout)]
                seq_row = sub[sub.system == "sequence_only"]
                if len(seq_row) == 0 or pd.isna(seq_row.rho_mean.values[0]):
                    continue
                seq_ci = {"mean": seq_row.rho_mean.values[0], "lower": seq_row.rho_lower90.values[0], "upper": seq_row.rho_upper90.values[0]}
                for fm in FM_MODELS:
                    fm_row = sub[sub.system == fm]
                    if len(fm_row) == 0 or pd.isna(fm_row.rho_mean.values[0]):
                        continue
                    fm_ci = {"mean": fm_row.rho_mean.values[0], "lower": fm_row.rho_lower90.values[0], "upper": fm_row.rho_upper90.values[0]}
                    no_overlap = (seq_ci["lower"] > fm_ci["upper"]) or (fm_ci["lower"] > seq_ci["upper"])
                    winner = "fm" if fm_ci["mean"] > seq_ci["mean"] else "sequence_only"
                    rs_verdict_rows.append({"config": config, "held_out": held_out, "readout": readout, "fm": fm,
                                             "fm_mean": fm_ci["mean"], "seqonly_mean": seq_ci["mean"],
                                             "distinguishable": bool(no_overlap), "point_winner": winner})
    rvdf = pd.DataFrame(rs_verdict_rows)
    rvdf.to_csv(RESULTS / "gate6_rs241_fm_verdict_table.csv", index=False)
    n_total_rs = len(rvdf)
    n_dist_rs = int(rvdf["distinguishable"].sum())
    n_fm_wins_rs = int(((rvdf["distinguishable"]) & (rvdf["point_winner"] == "fm")).sum())
    n_seq_wins_rs = int(((rvdf["distinguishable"]) & (rvdf["point_winner"] == "sequence_only")).sum())
    print(f"Of {n_total_rs} RS241 FM-vs-sequence-only comparisons: {n_dist_rs} distinguishable "
          f"({n_fm_wins_rs} FM wins, {n_seq_wins_rs} sequence-only wins)")
    with open(RESULTS / "gate6_rs241_fm_verdict_summary.json", "w") as f:
        json.dump({"n_total": n_total_rs, "n_distinguishable": n_dist_rs,
                    "n_fm_wins": n_fm_wins_rs, "n_seqonly_wins": n_seq_wins_rs}, f, indent=2)
    print(f"Wrote {RESULTS / 'gate6_rs241_fm_verdict_table.csv'}")


if __name__ == "__main__":
    main()
