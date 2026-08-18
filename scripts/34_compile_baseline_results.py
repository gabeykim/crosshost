"""
GATE 3 - Baseline 5: compile all baseline results into one table plus the
calibration-cost curves (B2 performance vs N, with error bars over 10 draws,
overlaid with B1/B3/B4 as reference lines).
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
    b1 = json.load(open(OUT / "baseline1_mean_majority.json"))
    b2 = json.load(open(OUT / "baseline2_calibration_raw.json"))
    b3 = json.load(open(OUT / "baseline3_host_embedding.json"))
    b4 = json.load(open(OUT / "baseline4_biophysical.json"))
    return b1, b2, b3, b4


def build_master_table(b1, b2, b3, b4):
    rows = []
    for host in HOSTS:
        for readout_key, readout_label in [("transcription", "transcription"), ("translation", "translation")]:
            # B1
            r1 = b1[host][readout_label]
            rows.append({"host": host, "readout": readout_label, "baseline": "B1_mean_majority",
                         "config": "N/A", "n_test": r1["n_test_usable"],
                         "spearman_rho_mean": r1["spearman_rho"], "spearman_rho_std": None,
                         "mcc_mean": r1["mcc"], "mcc_std": None, "auc_mean": r1["auc"], "auc_std": None})
            # B2 per N
            for N in N_VALUES:
                draws = b2[host]["results"][readout_label][str(N)]
                rhos = [d["spearman_rho"] for d in draws if d["spearman_rho"] is not None]
                mccs = [d["mcc"] for d in draws if d["mcc"] is not None]
                aucs = [d["auc"] for d in draws if d["auc"] is not None]
                rows.append({"host": host, "readout": readout_label, "baseline": "B2_per_host_calibration",
                             "config": f"N={N}", "n_test": b2[host]["meta"][readout_label]["n_test"],
                             "spearman_rho_mean": float(np.mean(rhos)) if rhos else None,
                             "spearman_rho_std": float(np.std(rhos)) if rhos else None,
                             "mcc_mean": float(np.mean(mccs)) if mccs else None,
                             "mcc_std": float(np.std(mccs)) if mccs else None,
                             "auc_mean": float(np.mean(aucs)) if aucs else None,
                             "auc_std": float(np.std(aucs)) if aucs else None})
            # B3 (held-out host = this host)
            r3 = b3[host][readout_label]
            rows.append({"host": host, "readout": readout_label, "baseline": "B3_free_host_embedding_diagnostic",
                         "config": f"train={r3['train_hosts']}, fallback=mean", "n_test": r3["n_test"],
                         "spearman_rho_mean": r3["spearman_rho"], "spearman_rho_std": None,
                         "mcc_mean": r3["mcc"], "mcc_std": None, "auc_mean": r3["auc"], "auc_std": None})
            # B4
            r4 = b4[host][readout_label]
            rows.append({"host": host, "readout": readout_label, "baseline": "B4_biophysical",
                         "config": r4["feature"], "n_test": r4["n"],
                         "spearman_rho_mean": r4["spearman_rho"], "spearman_rho_std": None,
                         "mcc_mean": None, "mcc_std": None, "auc_mean": None, "auc_std": None})

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "master_baseline_results.csv", index=False)
    return df


def plot_calibration_curves(b1, b2, b3, b4):
    for readout in ["transcription", "translation"]:
        fig, axes = plt.subplots(1, 2, figsize=(13, 5))
        for ax, metric, metric_label in zip(axes, ["spearman_rho", "mcc"], ["Spearman ρ (strength, on actives)", "MCC (active/inactive)"]):
            for host in HOSTS:
                means, stds = [], []
                for N in N_VALUES:
                    draws = b2[host]["results"][readout][str(N)]
                    vals = [d[metric] for d in draws if d[metric] is not None]
                    means.append(np.mean(vals) if vals else np.nan)
                    stds.append(np.std(vals) if vals else 0)
                means, stds = np.array(means), np.array(stds)
                ax.errorbar(N_VALUES, means, yerr=stds, marker="o", capsize=3,
                            label=f"{HOST_NAMES[host]} (B2, sequence-only)")

                # B3 reference (as a horizontal dashed line at this host's held-out result)
                b3_val = b3[host][readout][metric]
                if b3_val is not None:
                    ax.axhline(b3_val, linestyle=":", alpha=0.4,
                               color=ax.lines[-1].get_color())

            ax.set_xscale("symlog", linthresh=10)
            ax.set_xlabel("N calibration examples (log scale)")
            ax.set_ylabel(metric_label)
            ax.set_title(f"{readout.capitalize()}: {metric_label} vs N")
            ax.legend(fontsize=8)
            ax.grid(alpha=0.3)
        plt.tight_layout()
        out_path = OUT / f"calibration_curve_{readout}.png"
        plt.savefig(out_path, dpi=150)
        plt.close()
        print(f"Wrote {out_path}")


def print_gate5_relevant_summary(b2):
    """H-MAIN context: per-host, does N=100 approach N>=3000 performance?"""
    print("\n" + "=" * 70)
    print("CALIBRATION-COST CONTEXT (for the eventual Gate 5 H-MAIN comparison)")
    print("=" * 70)
    for host in HOSTS:
        for readout in ["transcription", "translation"]:
            draws100 = b2[host]["results"][readout]["100"]
            draws3000 = b2[host]["results"][readout]["3000"]
            r100 = np.mean([d["spearman_rho"] for d in draws100 if d["spearman_rho"] is not None])
            r3000 = np.mean([d["spearman_rho"] for d in draws3000 if d["spearman_rho"] is not None])
            pct = 100 * r100 / r3000 if r3000 else float("nan")
            print(f"  {host} {readout}: rho(N=100)={r100:.3f}, rho(N=3000)={r3000:.3f}, "
                  f"N=100 reaches {pct:.1f}% of N=3000 performance")


def main():
    b1, b2, b3, b4 = load_all()
    df = build_master_table(b1, b2, b3, b4)
    print(df.to_string(index=False))
    print(f"\nWrote {OUT / 'master_baseline_results.csv'}")

    plot_calibration_curves(b1, b2, b3, b4)
    print_gate5_relevant_summary(b2)


if __name__ == "__main__":
    main()
