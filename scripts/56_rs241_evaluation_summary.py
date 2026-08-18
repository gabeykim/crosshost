"""
GATE 5 - Task 4: RS241 evaluation summary. Seed-resolved 90% bootstrap CIs
(5 seeds, treated as fold-equivalents in the same hierarchical bootstrap
utility as everywhere else in this gate -- disclosed as training-stochasticity
uncertainty, not partition-level fold uncertainty, since RS241 has no fold
structure). Extrapolation distance for SE/VN/CG using the same method as
scripts/53.
"""
import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
m = import_module("40_film_cnn_model")
boot = import_module("51_bootstrap_utils")

OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

CONFIGS = {"PRIMARY": ["EC", "BS", "PA"], "SECONDARY": ["EC", "PA"]}
VARIANTS = ["genomic", "physiology"]
RS241_HOSTS = ["SE", "VN", "CG"]
READOUTS = ["transcription", "translation"]


def extrapolation_degree(z_df, held_out, train_hosts):
    held_vec = z_df.loc[held_out].values
    train_vecs = z_df.loc[train_hosts].values
    train_min, train_max = train_vecs.min(axis=0), train_vecs.max(axis=0)
    below = np.clip(train_min - held_vec, 0, None)
    above = np.clip(held_vec - train_max, 0, None)
    total = float((below + above).sum())
    return total, total / z_df.shape[1]


def main():
    with open(OUT / "gate5_rs241_results.json") as f:
        raw = json.load(f)
    g_z, gcols, p_z, pcols, imputed = m.load_host_features()

    print("=" * 70)
    print("RS241 zero-shot (N=0) performance, seed-resolved 90% bootstrap CI")
    print("=" * 70)
    summary = {}
    rows = []
    for config, train_hosts in CONFIGS.items():
        summary[config] = {}
        for variant in VARIANTS:
            summary[config][variant] = {}
            for held_out in RS241_HOSTS:
                summary[config][variant][held_out] = {}
                for readout in READOUTS:
                    seed_draws = []
                    n_test_list = []
                    for seed in range(5):
                        entry = raw[config][variant]["per_seed"][str(seed)]["eval"][held_out][readout]
                        if entry["spearman_rho"] is not None:
                            seed_draws.append(np.array([entry["spearman_rho"]]))
                        n_test_list.append(entry["n_test"])
                    ci = boot.bootstrap_ci_90(seed_draws)
                    n_test = n_test_list[0] if n_test_list else None
                    summary[config][variant][held_out][readout] = {"rho_ci": ci, "n_test_corrected_207": n_test}
                    rows.append({"config": config, "variant": variant, "held_out": held_out, "readout": readout,
                                 "n_test": n_test, "rho_mean": ci["mean"], "rho_lower90": ci["lower"], "rho_upper90": ci["upper"]})
                    print(f"  {config} {variant} {held_out} {readout}: n={n_test}, rho={ci['mean']:.3f} "
                          f"[{ci['lower']:.3f},{ci['upper']:.3f}]")

    print()
    print("=" * 70)
    print("Extrapolation distance, RS241 hosts (SE/VN/CG) -- further out-of-distribution than any primary host?")
    print("=" * 70)
    extrap = {}
    for config, train_hosts in CONFIGS.items():
        extrap[config] = {"genomic": {}, "physiology": {}}
        for held_out in RS241_HOSTS:
            g_tot, g_per = extrapolation_degree(g_z, held_out, train_hosts)
            p_tot, p_per = extrapolation_degree(p_z, held_out, train_hosts)
            extrap[config]["genomic"][held_out] = {"total": g_tot, "per_dim": g_per}
            extrap[config]["physiology"][held_out] = {"total": p_tot, "per_dim": p_per}
            print(f"  {config} {held_out}: genomic_total={g_tot:.3f} (per_dim={g_per:.4f}), "
                  f"physiology_total={p_tot:.3f} (per_dim={p_per:.4f})")

    print("\n  For reference, primary-host genomic per_dim range was 0.309-0.806; physiology per_dim range was 0.170-0.615")

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "gate5_rs241_table.csv", index=False)
    print(f"\nWrote {RESULTS / 'gate5_rs241_table.csv'}")

    output = {"performance": summary, "extrapolation_distance": extrap,
              "n_seeds": 5, "seed_uncertainty_note": "5 random-init/shuffle seeds on IDENTICAL training data -- training-stochasticity uncertainty, NOT partition-level fold uncertainty (RS241 has no fold structure, it is a wholly reserved evaluation set)"}
    with open(RESULTS / "gate5_rs241_results_summary.json", "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"Wrote {RESULTS / 'gate5_rs241_results_summary.json'}")


if __name__ == "__main__":
    main()
