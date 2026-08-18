"""
GATE 8.6 - Task 2a: FiLM's fold-to-fold instability, characterized properly
and compared against the two alternative mechanisms (concatenation,
per-host-heads) on the identical 5 folds. All from ALREADY-COMPUTED Gate
8.5 zero-shot LOHO results -- no retraining.

Reports fold-to-fold standard deviation of zero-shot Spearman rho, per
mechanism, per host, per readout -- the same statistic Gate 3.5
(out/baselines/fold_rotation_all_baselines.json) and Gate 4
(scripts/44_cnn_fold_variance_spotcheck.py) used for the fold-variance
comparison, extended here to cover all three Task 1A mechanisms side by
side rather than one model against one baseline.

NOTE ON "DRAW" VARIANCE: Gate 4's own fold-vs-draw comparison (scripts/44)
required a DRAW dimension -- repeated fine-tune runs at fixed data, from
Gate 4's N=100 calibration-curve protocol (10 draws/fold). Task 1A's
zero-shot LOHO protocol (this script's data source) trains exactly ONE
model per (host, fold) -- there is no repeated-draw dimension in the
already-computed data to reuse the way scripts/44 did. Computing a true
draw-level variance for the zero-shot cell would require re-training each
mechanism multiple times per fold with different weight-init seeds (data
held fixed) -- NOT part of Task 1A's original protocol and not run here,
disclosed as a scope limit rather than silently substituted. What IS
reported: fold-to-fold variance (real, from the actual completed runs) as
the primary comparison, which is itself the quantity Gate 3.5's own
headline finding (SR5) was about.
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

HOSTS = ["EC", "BS", "PA"]
READOUTS = ["transcription", "translation"]


def film_fold_values(host, readout):
    d = json.load(open(OUT / "gate4_loho_results.json"))
    return [d[host]["genomic"]["per_fold"][str(f)]["eval"][readout]["spearman_rho"] for f in range(5)]


def concat_fold_values(host, readout):
    d = json.load(open(OUT / "gate8_5_concat_loho_results.json"))
    return [d[host]["per_fold"][str(f)]["eval"][readout]["spearman_rho"] for f in range(5)]


def perhost_avg_fold_values(host, readout):
    d = json.load(open(OUT / "gate8_5_perhost_heads_loho_results.json"))
    return [d[host]["per_fold"][str(f)]["eval_avg"][readout]["spearman_rho"] for f in range(5)]


def seqonly_fold_values(host, readout):
    d = json.load(open(OUT / "gate5_5_seqonly_loho_results.json"))
    return [d[host]["per_fold"][str(f)]["eval"][readout]["spearman_rho"] for f in range(5)]


def main():
    rows = []
    for host in HOSTS:
        for readout in READOUTS:
            for mech_name, fn in [("sequence_only", seqonly_fold_values), ("film_genomic", film_fold_values),
                                    ("concat", concat_fold_values), ("perhost_heads_avg", perhost_avg_fold_values)]:
                vals = [v for v in fn(host, readout) if v is not None]
                if len(vals) < 2:
                    continue
                rows.append({
                    "host": host, "readout": readout, "mechanism": mech_name,
                    "fold_values": vals, "fold_mean": float(np.mean(vals)),
                    "fold_to_fold_std": float(np.std(vals)),
                    "fold_min": float(np.min(vals)), "fold_max": float(np.max(vals)),
                    "fold_range": float(np.max(vals) - np.min(vals)),
                })

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "gate8_6_mechanism_fold_variance.csv", index=False)
    pd.set_option("display.width", 200)
    print(df[["host", "readout", "mechanism", "fold_mean", "fold_to_fold_std", "fold_range"]].round(3).to_string(index=False))

    # headline comparison: EC transcription, all 4 mechanisms
    print("\n=== EC transcription: fold-to-fold std, all mechanisms ===")
    ec_tx = df[(df.host == "EC") & (df.readout == "transcription")]
    print(ec_tx[["mechanism", "fold_mean", "fold_to_fold_std", "fold_range"]].round(3).to_string(index=False))

    film_std = ec_tx[ec_tx.mechanism == "film_genomic"]["fold_to_fold_std"].values[0]
    concat_std = ec_tx[ec_tx.mechanism == "concat"]["fold_to_fold_std"].values[0]
    perhost_std = ec_tx[ec_tx.mechanism == "perhost_heads_avg"]["fold_to_fold_std"].values[0]
    seqonly_std = ec_tx[ec_tx.mechanism == "sequence_only"]["fold_to_fold_std"].values[0]

    summary = {
        "ec_transcription_fold_to_fold_std": {
            "film_genomic": float(film_std), "concat": float(concat_std),
            "perhost_heads_avg": float(perhost_std), "sequence_only": float(seqonly_std),
        },
        "film_std_over_concat_std": float(film_std / concat_std) if concat_std > 0 else None,
        "film_std_over_perhost_std": float(film_std / perhost_std) if perhost_std > 0 else None,
        "film_std_over_seqonly_std": float(film_std / seqonly_std) if seqonly_std > 0 else None,
        "note": ("Draw-level (weight-init) variance not computed -- Task 1A's zero-shot LOHO protocol trains one "
                 "model per (host,fold); a true fold-vs-draw ratio comparable to scripts/44 would require repeated "
                 "same-data reruns not part of the original protocol. Fold-to-fold std reported directly instead."),
    }
    with open(RESULTS / "gate8_6_mechanism_fold_variance_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\n{json.dumps(summary, indent=2)}")
    print(f"\nWrote {RESULTS / 'gate8_6_mechanism_fold_variance.csv'}, "
          f"{RESULTS / 'gate8_6_mechanism_fold_variance_summary.json'}")


if __name__ == "__main__":
    main()
