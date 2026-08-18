"""
GATE 4 - compile Task 2 (LOHO) and Task 3 (calibration curve) results into
machine-readable tables and figures. Does NOT compute any H-MAIN comparison
(Gate 5's job) -- reports the model's own performance only.
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "out"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

HOSTS = ["EC", "BS", "PA"]
HOST_NAMES = {"EC": "E. coli", "BS": "B. subtilis", "PA": "P. aeruginosa"}
VARIANTS = ["genomic", "physiology"]
READOUTS = ["transcription", "translation"]
CELLS = [(10, "head_only"), (30, "head_only"), (100, "head_only"),
         (300, "head_only"), (300, "top_conv"), (1000, "top_conv"), (3000, "top_conv")]


def build_loho_table(loho):
    rows = []
    for host in HOSTS:
        for variant in VARIANTS:
            for readout in READOUTS:
                rhos, mccs, aucs = [], [], []
                for fold in range(5):
                    e = loho[host][variant]["per_fold"][str(fold)]["eval"][readout]
                    if e["spearman_rho"] is not None:
                        rhos.append(e["spearman_rho"])
                    if e["mcc"] is not None:
                        mccs.append(e["mcc"])
                    if e["auc"] is not None:
                        aucs.append(e["auc"])
                rows.append({
                    "host": host, "variant": variant, "readout": readout, "N": "N=0 (LOHO base)",
                    "spearman_rho_mean": np.mean(rhos) if rhos else None,
                    "spearman_rho_std": np.std(rhos) if rhos else None,
                    "mcc_mean": np.mean(mccs) if mccs else None, "mcc_std": np.std(mccs) if mccs else None,
                    "auc_mean": np.mean(aucs) if aucs else None, "auc_std": np.std(aucs) if aucs else None,
                    "n_folds": len(rhos),
                })
    return pd.DataFrame(rows)


def build_calibration_table(cal, loho):
    rows = []
    for host in HOSTS:
        for variant in VARIANTS:
            for readout in READOUTS:
                rhos0 = [loho[host][variant]["per_fold"][str(f)]["eval"][readout]["spearman_rho"] for f in range(5)]
                rhos0 = [r for r in rhos0 if r is not None]
                rows.append({
                    "host": host, "variant": variant, "readout": readout, "N": 0, "mechanism": "none (LOHO zero-shot)",
                    "spearman_rho_mean": np.mean(rhos0) if rhos0 else None,
                    "spearman_rho_std_across_folds": np.std(rhos0) if rhos0 else None,
                    "mcc_mean": None, "auc_mean": None, "n_folds": len(rhos0),
                })
            for N, mech in CELLS:
                key = f"N{N}_{mech}"
                for readout in READOUTS:
                    fold_means_rho, fold_means_mcc, fold_means_auc = [], [], []
                    for fold in range(5):
                        draws = cal[host][variant][str(fold)][key]
                        rhos = [d["eval"][readout]["spearman_rho"] for d in draws if d["eval"][readout]["spearman_rho"] is not None]
                        mccs = [d["eval"][readout]["mcc"] for d in draws if d["eval"][readout]["mcc"] is not None]
                        aucs = [d["eval"][readout]["auc"] for d in draws if d["eval"][readout]["auc"] is not None]
                        if rhos:
                            fold_means_rho.append(np.mean(rhos))
                        if mccs:
                            fold_means_mcc.append(np.mean(mccs))
                        if aucs:
                            fold_means_auc.append(np.mean(aucs))
                    rows.append({
                        "host": host, "variant": variant, "readout": readout, "N": N, "mechanism": mech,
                        "spearman_rho_mean": np.mean(fold_means_rho) if fold_means_rho else None,
                        "spearman_rho_std_across_folds": np.std(fold_means_rho) if fold_means_rho else None,
                        "mcc_mean": np.mean(fold_means_mcc) if fold_means_mcc else None,
                        "auc_mean": np.mean(fold_means_auc) if fold_means_auc else None,
                        "n_folds": len(fold_means_rho),
                    })
    return pd.DataFrame(rows)


def plot_calibration_curves(cal, loho):
    for host in HOSTS:
        for readout in READOUTS:
            fig, ax = plt.subplots(figsize=(8, 5.5))
            for variant, color in [("genomic", "tab:blue"), ("physiology", "tab:orange")]:
                # N=0 point: the unmodified LOHO base model (Task 2), NOT recomputed by
                # Task 3 -- included here for a complete, honest picture (an early
                # draft of this figure omitted it and hid a real drop-then-recover
                # pattern for physiology at small N).
                rhos0 = [loho[host][variant]["per_fold"][str(f)]["eval"][readout]["spearman_rho"] for f in range(5)]
                rhos0 = [r for r in rhos0 if r is not None]
                ax.errorbar([0.5], [np.mean(rhos0)], yerr=[np.std(rhos0)], marker="*", markersize=14,
                            capsize=3, color=color, linestyle="none",
                            label=f"{variant} (N=0, LOHO zero-shot)")
                ns, means, stds, mechs = [], [], [], []
                for N, mech in CELLS:
                    key = f"N{N}_{mech}"
                    fold_means = []
                    for fold in range(5):
                        draws = cal[host][variant][str(fold)][key]
                        rhos = [d["eval"][readout]["spearman_rho"] for d in draws if d["eval"][readout]["spearman_rho"] is not None]
                        if rhos:
                            fold_means.append(np.mean(rhos))
                    ns.append(N)
                    means.append(np.mean(fold_means) if fold_means else np.nan)
                    stds.append(np.std(fold_means) if fold_means else 0)
                    mechs.append(mech)
                ns = np.array(ns)
                means = np.array(means)
                stds = np.array(stds)
                head_mask = np.array([m == "head_only" for m in mechs])
                top_mask = ~head_mask
                ax.errorbar(ns[head_mask], means[head_mask], yerr=stds[head_mask], marker="o",
                            capsize=3, color=color, linestyle="-", label=f"{variant} (frozen trunk)")
                ax.errorbar(ns[top_mask], means[top_mask], yerr=stds[top_mask], marker="s",
                            capsize=3, color=color, linestyle="--", label=f"{variant} (top conv unfrozen)")
            ax.axvline(300, color="gray", linestyle=":", alpha=0.5, label="N=300 crossover (both mechanisms shown)")
            ax.set_xscale("symlog", linthresh=10)
            ax.set_xlabel("N calibration examples from held-out host (log scale)")
            ax.set_ylabel("Spearman ρ (strength, on actives), fold-resolved mean ± std")
            ax.set_title(f"{HOST_NAMES[host]} held out, {readout}\nFiLM-CNN calibration curve, genomic vs physiology")
            ax.legend(fontsize=7.5, loc="best")
            ax.grid(alpha=0.3)
            plt.tight_layout()
            out_path = FIG / f"calibration_curve_{host}_{readout}.png"
            plt.savefig(out_path, dpi=150)
            plt.close()
            print(f"Wrote {out_path}")


def plot_hscience_comparison(loho):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    x = np.arange(len(HOSTS))
    width = 0.35
    for ax, readout in zip(axes, READOUTS):
        g_means, g_stds, p_means, p_stds = [], [], [], []
        for host in HOSTS:
            for variant, means, stds in [("genomic", g_means, g_stds), ("physiology", p_means, p_stds)]:
                rhos = [loho[host][variant]["per_fold"][str(f)]["eval"][readout]["spearman_rho"] for f in range(5)]
                rhos = [r for r in rhos if r is not None]
                means.append(np.mean(rhos))
                stds.append(np.std(rhos))
        ax.bar(x - width / 2, g_means, width, yerr=g_stds, capsize=4, label="genomic (37-D)", color="tab:blue")
        ax.bar(x + width / 2, p_means, width, yerr=p_stds, capsize=4, label="physiology (6-D)", color="tab:orange")
        ax.set_xticks(x)
        ax.set_xticklabels([HOST_NAMES[h] for h in HOSTS])
        ax.set_ylabel("Spearman ρ, LOHO base model (N=0), fold-resolved mean ± std")
        ax.set_title(f"{readout.capitalize()}: genomic vs physiology host features")
        ax.axhline(0, color="black", linewidth=0.8)
        ax.legend()
        ax.grid(alpha=0.3, axis="y")
    plt.tight_layout()
    out_path = FIG / "hscience_genomic_vs_physiology.png"
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Wrote {out_path}")


def main():
    with open(OUT / "gate4_loho_results.json") as f:
        loho = json.load(f)
    with open(OUT / "gate4_calibration_curves.json") as f:
        cal = json.load(f)

    loho_df = build_loho_table(loho)
    loho_df.to_csv(OUT / "gate4_loho_results_table.csv", index=False)
    print(loho_df.to_string(index=False))
    print(f"\nWrote {OUT / 'gate4_loho_results_table.csv'}")

    cal_df = build_calibration_table(cal, loho)
    cal_df.to_csv(OUT / "gate4_calibration_results_table.csv", index=False)
    print(f"\nWrote {OUT / 'gate4_calibration_results_table.csv'} ({len(cal_df)} rows)")

    plot_calibration_curves(cal, loho)
    plot_hscience_comparison(loho)


if __name__ == "__main__":
    main()
