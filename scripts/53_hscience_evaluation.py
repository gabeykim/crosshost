"""
GATE 5 - Task 2: H-SCIENCE, the Bernstein test. Genomic vs physiology host
features, full grid (3 hosts x 2 readouts x 2 mechanisms x 2 variants),
fold-resolved with 90% bootstrap CIs, plus the extrapolation-distance
covariate (Task 2.3).
"""
import sys
import json
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
m = import_module("40_film_cnn_model")
boot = import_module("51_bootstrap_utils")

OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

HOSTS = ["EC", "BS", "PA"]
READOUTS = ["transcription", "translation"]
VARIANTS = ["genomic", "physiology"]
MECHANISMS = {"head_only": "frozen-trunk (H-MAIN primary)", "top_conv": "unfrozen-conv4 (supplementary)"}
TRAIN_PAIRS = {"EC": ["BS", "PA"], "BS": ["EC", "PA"], "PA": ["EC", "BS"]}


def get_fold_draws(host, variant, readout, mechanism, metric):
    if mechanism == "head_only":
        with open(OUT / "gate4_calibration_curves.json") as f:
            cal = json.load(f)
        out = []
        for f in range(5):
            draws = cal[host][variant][str(f)]["N100_head_only"]
            out.append(np.array([x["eval"][readout][metric] for x in draws if x["eval"][readout][metric] is not None]))
        return out
    else:
        with open(OUT / "gate4_5_n100_topconv_supplement.json") as f:
            supp = json.load(f)
        out = []
        for f in range(5):
            draws = supp[variant][host][str(f)]
            out.append(np.array([x["eval"][readout][metric] for x in draws if x["eval"][readout][metric] is not None]))
        return out


def extrapolation_degree(z_df, held_out, train_hosts):
    held_vec = z_df.loc[held_out].values
    train_vecs = z_df.loc[train_hosts].values
    train_min, train_max = train_vecs.min(axis=0), train_vecs.max(axis=0)
    below = np.clip(train_min - held_vec, 0, None)
    above = np.clip(held_vec - train_max, 0, None)
    total = float((below + above).sum())
    per_dim = total / z_df.shape[1]
    return total, per_dim


def main():
    g_z, gcols, p_z, pcols, imputed = m.load_host_features()
    assert imputed == []

    print("=" * 70)
    print("H-SCIENCE full grid: 3 hosts x 2 readouts x 2 mechanisms x 2 variants")
    print("=" * 70)
    grid = {}
    for host in HOSTS:
        grid[host] = {}
        for readout in READOUTS:
            grid[host][readout] = {}
            for mech in MECHANISMS:
                grid[host][readout][mech] = {}
                for variant in VARIANTS:
                    cell = {}
                    for metric in ["spearman_rho", "mcc", "auc"]:
                        fd = get_fold_draws(host, variant, readout, mech, metric)
                        cell[metric] = boot.bootstrap_ci_90(fd)
                    grid[host][readout][mech][variant] = cell
                    print(f"  {host} {readout} {mech} {variant}: "
                          f"rho={cell['spearman_rho']['mean']:.3f} "
                          f"[{cell['spearman_rho']['lower']:.3f},{cell['spearman_rho']['upper']:.3f}]  "
                          f"auc={cell['auc']['mean']}")

    print()
    print("=" * 70)
    print("Overall verdict inputs: genomic vs physiology, per host, PRIMARY mechanism (rho)")
    print("=" * 70)
    per_host_verdict = {}
    for host in HOSTS:
        per_host_verdict[host] = {}
        for readout in READOUTS:
            g = grid[host][readout]["head_only"]["genomic"]["spearman_rho"]
            p = grid[host][readout]["head_only"]["physiology"]["spearman_rho"]
            g_wins = g["mean"] > p["mean"]
            overlap = boot.overlap_size(g, p)
            no_overlap = (g["lower"] > p["upper"]) or (p["lower"] > g["upper"])
            winner = "genomic" if g_wins else "physiology"
            distinguishable = no_overlap
            per_host_verdict[host][readout] = {
                "genomic_rho": g, "physiology_rho": p, "point_estimate_winner": winner,
                "intervals_distinguishable": distinguishable, "overlap": overlap,
            }
            print(f"  {host} {readout}: genomic={g['mean']:.3f} physiology={p['mean']:.3f} "
                  f"winner(point)={winner} distinguishable(90%CI)={distinguishable}")

    print()
    print("=" * 70)
    print("Extrapolation distance covariate (n=3 hosts -- NO statistical claim supported)")
    print("=" * 70)
    extrap = {"genomic": {}, "physiology": {}}
    for host in HOSTS:
        train = TRAIN_PAIRS[host]
        g_tot, g_per = extrapolation_degree(g_z, host, train)
        p_tot, p_per = extrapolation_degree(p_z, host, train)
        extrap["genomic"][host] = {"total": g_tot, "per_dim": g_per, "n_dims": len(gcols)}
        extrap["physiology"][host] = {"total": p_tot, "per_dim": p_per, "n_dims": len(pcols)}
        print(f"  {host}: genomic_extrap_total={g_tot:.3f} (per_dim={g_per:.4f}, n_dims={len(gcols)}), "
              f"physiology_extrap_total={p_tot:.3f} (per_dim={p_per:.4f}, n_dims={len(pcols)})")

    # performance vs extrapolation distance, primary mechanism, both readouts
    perf_vs_extrap = {}
    for variant in VARIANTS:
        perf_vs_extrap[variant] = {}
        for readout in READOUTS:
            rows = []
            for host in HOSTS:
                perf = grid[host][readout]["head_only"][variant]["spearman_rho"]["mean"]
                dist = extrap[variant][host]["total"]
                rows.append({"host": host, "extrapolation_distance": dist, "performance_rho": perf})
            rows_sorted = sorted(rows, key=lambda r: r["extrapolation_distance"])
            monotonic_decreasing = all(rows_sorted[i]["performance_rho"] >= rows_sorted[i+1]["performance_rho"]
                                        for i in range(len(rows_sorted)-1))
            perf_vs_extrap[variant][readout] = {"rows": rows, "sorted_by_distance": rows_sorted,
                                                  "monotonic_decreasing_performance": monotonic_decreasing}
            print(f"  {variant} {readout}: sorted_by_extrap={[(r['host'], round(r['extrapolation_distance'],2), round(r['performance_rho'],3)) for r in rows_sorted]} "
                  f"monotonic_decreasing={monotonic_decreasing}")

    output = {
        "full_grid": grid,
        "per_host_verdict_primary_mechanism": per_host_verdict,
        "extrapolation_distance": extrap,
        "performance_vs_extrapolation_distance": perf_vs_extrap,
        "n_hosts_caveat": "n=3 primary hosts -- no statistical claim about the extrapolation-distance relationship is supported; ordering and monotonicity are reported descriptively only",
    }
    with open(RESULTS / "gate5_hscience_results.json", "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nWrote {RESULTS / 'gate5_hscience_results.json'}")


if __name__ == "__main__":
    main()
