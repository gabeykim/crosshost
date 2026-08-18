"""
GATE 3.5 - Task 2 (compile step): update out/baselines/master_baseline_results.csv
with fold-resolved statistics (mean/std across the 5 folds, not just 10 draws
within fold 0) and regenerate the calibration-cost figures with fold-level
error bands layered on top of the original within-fold draw error bars.
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "out" / "baselines"
HOSTS = ["EC", "BS", "PA"]
HOST_NAMES = {"EC": "E. coli", "BS": "B. subtilis", "PA": "P. aeruginosa"}
N_VALUES = [0, 10, 30, 100, 300, 1000, 3000]


def load_all():
    fr = json.load(open(OUT / "fold_rotation_all_baselines.json"))
    rs241 = json.load(open(OUT / "baseline3_rs241_extension.json"))
    return fr, rs241


def build_master_table(fr):
    rows = []
    for host in HOSTS:
        for readout in ["transcription", "translation"]:
            # B1, fold-resolved
            b1 = fr["b1_mean_majority"]["summary"][host][readout]
            rows.append({"host": host, "readout": readout, "baseline": "B1_mean_majority", "config": "N/A",
                         "spearman_rho_mean": b1["spearman_rho"]["mean"], "spearman_rho_std_across_folds": b1["spearman_rho"]["std"],
                         "mcc_mean": b1["mcc"]["mean"], "mcc_std_across_folds": b1["mcc"]["std"],
                         "auc_mean": b1["auc"]["mean"], "auc_std_across_folds": b1["auc"]["std"],
                         "n_folds": b1["mcc"]["n_folds"]})
            # B2, fold-resolved per N
            for N in N_VALUES:
                b2 = fr["b2_calibration"]["summary"][host][readout][str(N)]
                rows.append({"host": host, "readout": readout, "baseline": "B2_per_host_calibration",
                             "config": f"N={N}",
                             "spearman_rho_mean": b2["rho_across_folds"]["mean"],
                             "spearman_rho_std_across_folds": b2["rho_across_folds"]["std"],
                             "mcc_mean": b2["mcc_across_folds"]["mean"],
                             "mcc_std_across_folds": b2["mcc_across_folds"]["std"],
                             "auc_mean": None, "auc_std_across_folds": None,
                             "n_folds": b2["rho_across_folds"]["n_folds"]})
            # B3, fold-resolved
            b3 = fr["b3_free_host_embedding"]["summary"][host][readout]
            rows.append({"host": host, "readout": readout, "baseline": "B3_free_host_embedding_diagnostic",
                         "config": "fallback=mean",
                         "spearman_rho_mean": b3["spearman_rho"]["mean"], "spearman_rho_std_across_folds": b3["spearman_rho"]["std"],
                         "mcc_mean": b3["mcc"]["mean"], "mcc_std_across_folds": b3["mcc"]["std"],
                         "auc_mean": b3["auc"]["mean"], "auc_std_across_folds": b3["auc"]["std"],
                         "n_folds": b3["spearman_rho"]["n_folds"]})
            # B4, fold-resolved
            b4 = fr["b4_biophysical"]["summary"][host][readout]
            rows.append({"host": host, "readout": readout, "baseline": "B4_biophysical", "config": "n/a",
                         "spearman_rho_mean": b4["spearman_rho"]["mean"], "spearman_rho_std_across_folds": b4["spearman_rho"]["std"],
                         "mcc_mean": None, "mcc_std_across_folds": None, "auc_mean": None, "auc_std_across_folds": None,
                         "n_folds": b4["spearman_rho"]["n_folds"]})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "master_baseline_results.csv", index=False)
    return df


def build_rs241_table(rs241):
    rows = []
    for config, hosts in rs241["results"].items():
        for held_out, readouts in hosts.items():
            for readout, r in readouts.items():
                rows.append({"config": config, "held_out_host": held_out, "readout": readout,
                             "n_test": r.get("n_test"), "mcc": r.get("mcc"), "auc": r.get("auc"),
                             "spearman_rho": r.get("spearman_rho")})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "baseline3_rs241_extension_table.csv", index=False)
    return df


def plot_calibration_curves_fold_resolved(fr):
    for readout in ["transcription", "translation"]:
        fig, axes = plt.subplots(1, 2, figsize=(13, 5))
        for ax, metric_key, metric_label in zip(axes, ["rho_across_folds", "mcc_across_folds"],
                                                   ["Spearman ρ (strength, on actives)", "MCC (active/inactive)"]):
            for host in HOSTS:
                means, fold_stds = [], []
                for N in N_VALUES:
                    s = fr["b2_calibration"]["summary"][host][readout][str(N)][metric_key]
                    means.append(s["mean"] if s["mean"] is not None else np.nan)
                    fold_stds.append(s["std"] if s["std"] is not None else 0)
                means, fold_stds = np.array(means), np.array(fold_stds)
                line = ax.errorbar(N_VALUES, means, yerr=fold_stds, marker="o", capsize=3,
                                    label=f"{HOST_NAMES[host]} (fold-to-fold band)")
                color = line.lines[0].get_color()
                ax.fill_between(N_VALUES, means - fold_stds, means + fold_stds, alpha=0.12, color=color)
            ax.set_xscale("symlog", linthresh=10)
            ax.set_xlabel("N calibration examples (log scale)")
            ax.set_ylabel(metric_label)
            ax.set_title(f"{readout.capitalize()}: {metric_label} vs N\n(error band = std across 5 held-out folds, mean of 10 draws per fold)")
            ax.legend(fontsize=8)
            ax.grid(alpha=0.3)
        plt.tight_layout()
        out_path = OUT / f"calibration_curve_{readout}_fold_resolved.png"
        plt.savefig(out_path, dpi=150)
        plt.close()
        print(f"Wrote {out_path}")


def print_variance_comparison_summary(fr):
    print("\n" + "=" * 78)
    print("FOLD-VARIANCE vs DRAW-VARIANCE (the key Task 2 methodological finding)")
    print("=" * 78)
    vc = fr["b2_calibration"]["variance_comparison"]
    dominates_count, total_count = 0, 0
    for host in HOSTS:
        for readout in ["transcription", "translation"]:
            for N in [10, 100, 1000, 3000]:
                v = vc[host][readout].get(str(N))
                if v and v["ratio_fold_std_over_draw_std"] is not None:
                    total_count += 1
                    if v["fold_variance_dominates"]:
                        dominates_count += 1
                    print(f"  {host} {readout} N={N}: fold_std={v['fold_to_fold_std_of_rho']:.4f}, "
                          f"draw_std={v['mean_within_fold_draw_to_draw_std_of_rho']:.4f}, "
                          f"ratio={v['ratio_fold_std_over_draw_std']:.2f}, "
                          f"fold_dominates={v['fold_variance_dominates']}")
    print(f"\n  Fold variance exceeds draw variance in {dominates_count}/{total_count} "
          f"(host, readout, N) combinations checked.")


def main():
    fr, rs241 = load_all()
    df = build_master_table(fr)
    print(df.to_string(index=False))
    print(f"\nWrote {OUT / 'master_baseline_results.csv'} (fold-resolved)")

    rs241_df = build_rs241_table(rs241)
    print(f"\n{rs241_df.to_string(index=False)}")
    print(f"Wrote {OUT / 'baseline3_rs241_extension_table.csv'}")

    plot_calibration_curves_fold_resolved(fr)
    print_variance_comparison_summary(fr)


if __name__ == "__main__":
    main()
