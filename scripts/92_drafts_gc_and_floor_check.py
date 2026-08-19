"""
GATE 10 - Task 4.4/4.5: the GC/phylum confound and per-host floor-value
artifact check, on DRAFTS's own terms. Descriptive only, no models.
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

SPECIES = ["Ec", "Ef", "Se", "Ko", "Pa", "Pp", "Vn", "Bs", "Cg", "Ll"]
SPECIES_PHYLUM = {
    "Ec": "Proteobacteria", "Ef": "Proteobacteria", "Se": "Proteobacteria", "Ko": "Proteobacteria",
    "Pa": "Proteobacteria", "Pp": "Proteobacteria", "Vn": "Proteobacteria",
    "Bs": "Firmicutes", "Ll": "Firmicutes", "Cg": "Actinobacteria",
}


def gc_confound(df):
    print("=" * 80)
    print("TASK 4.4: source-genome GC content vs. activity, per recipient host")
    print("=" * 80)
    rows = []
    for sp in SPECIES:
        usable = df[f"{sp}_usable"]
        gc = df.loc[usable, "gc_pct"]
        tx = df.loc[usable, f"{sp}_tx"]
        log_tx = np.log10(tx.clip(lower=1e-6))
        rho = float(spearmanr(gc, log_tx).correlation)
        rows.append({"species": sp, "phylum": SPECIES_PHYLUM[sp], "n": int(usable.sum()),
                     "spearman_rho_gc_vs_log_tx": rho})
        print(f"  {sp} ({SPECIES_PHYLUM[sp]}): n={int(usable.sum())}, rho(GC%, log10 tx)={rho:.3f}")

    # source-phylum activity comparison: are Firmicutes-derived sequences more active,
    # as the paper reports, and does this hold in THIS project's re-derived data?
    print("\n--- Source-phylum (donor) activity comparison, per recipient host ---")
    phylum_rows = []
    for sp in SPECIES:
        usable = df[f"{sp}_usable"]
        sub = df.loc[usable, ["source_phylum", f"{sp}_tx"]].copy()
        sub["log_tx"] = np.log10(sub[f"{sp}_tx"].clip(lower=1e-6))
        by_phylum = sub.groupby("source_phylum")["log_tx"].agg(["mean", "count"])
        by_phylum = by_phylum[by_phylum["count"] >= 10].sort_values("mean", ascending=False)
        top_phylum = by_phylum.index[0] if len(by_phylum) else None
        firmicutes_mean = by_phylum.loc["Firmicutes", "mean"] if "Firmicutes" in by_phylum.index else None
        proteo_mean = by_phylum.loc["Proteobacteria", "mean"] if "Proteobacteria" in by_phylum.index else None
        phylum_rows.append({"recipient": sp, "highest_mean_activity_source_phylum": top_phylum,
                             "firmicutes_source_mean_log_tx": firmicutes_mean,
                             "proteobacteria_source_mean_log_tx": proteo_mean,
                             "firmicutes_gt_proteobacteria": (firmicutes_mean > proteo_mean)
                             if firmicutes_mean is not None and proteo_mean is not None else None})
        print(f"  recipient={sp}: highest-activity donor phylum={top_phylum}, "
              f"Firmicutes-donor mean log_tx={firmicutes_mean}, Proteobacteria-donor mean log_tx={proteo_mean}")

    return rows, phylum_rows


def floor_artifact(df):
    print("\n" + "=" * 80)
    print("TASK 4.5: per-host value distribution and floor/pile-up check")
    print("=" * 80)
    rows = []
    for sp in SPECIES:
        usable = df[f"{sp}_usable"]
        tx = df.loc[usable, f"{sp}_tx"]
        rounded = tx.round(4)
        mode_val = rounded.mode().iloc[0] if len(rounded) else None
        mode_frac = (rounded == mode_val).mean() if mode_val is not None else None
        n_unusable = int((~df[f"{sp}_usable"]).sum())
        unusable_frac = n_unusable / len(df)
        rows.append({
            "species": sp, "n_usable": int(usable.sum()), "unusable_fraction": round(unusable_frac, 4),
            "tx_min": float(tx.min()), "tx_median": float(tx.median()), "tx_max": float(tx.max()),
            "mode_value": float(mode_val) if mode_val is not None else None,
            "mode_fraction_of_usable": float(mode_frac) if mode_frac is not None else None,
            "floor_pileup_flag": bool(mode_frac is not None and mode_frac > 0.05),
        })
        flag = " <-- POSSIBLE PILE-UP" if (mode_frac is not None and mode_frac > 0.05) else ""
        print(f"  {sp}: unusable_frac={unusable_frac:.3f}, tx range=[{tx.min():.4f},{tx.max():.4f}], "
              f"median={tx.median():.4f}, mode={mode_val} (frac={mode_frac:.4f}){flag}")
    return rows


def main():
    df = pd.read_parquet(DATA / "drafts.parquet")
    gc_rows, phylum_rows = gc_confound(df)
    floor_rows = floor_artifact(df)

    with open(RESULTS / "gate10_gc_confound.json", "w") as f:
        json.dump({"gc_vs_activity": gc_rows, "source_phylum_activity": phylum_rows}, f, indent=2, default=str)
    with open(RESULTS / "gate10_floor_artifact_check.json", "w") as f:
        json.dump(floor_rows, f, indent=2, default=str)
    pd.DataFrame(gc_rows).to_csv(RESULTS / "gate10_gc_confound.csv", index=False)
    pd.DataFrame(floor_rows).to_csv(RESULTS / "gate10_floor_artifact_check.csv", index=False)
    print(f"\nWrote gate10_gc_confound.{{json,csv}}, gate10_floor_artifact_check.{{json,csv}}")


if __name__ == "__main__":
    main()
