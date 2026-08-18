"""
GATE 5 - Task 3: H-DIAGNOSTIC. Host-biology representation (CNN, zero-shot
LOHO, N=0 -- no fine-tuning, matching B3's own "trained on pooled 2-host
data, evaluated zero-shot" structure) vs the free per-host lookup embedding
(Gate 3 B3) under its documented fallback (mean of training-host
embeddings). Diagnostic only -- evidence about biology vs identity learning,
not a gate.
"""
import sys
import json
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
boot = import_module("51_bootstrap_utils")

OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

HOSTS = ["EC", "BS", "PA"]
READOUTS = ["transcription", "translation"]
VARIANTS = ["genomic", "physiology"]


def get_cnn_zeroshot_fold_draws(host, variant, readout, metric):
    with open(OUT / "gate4_loho_results.json") as f:
        d = json.load(f)
    out = []
    for f in range(5):
        val = d[host][variant]["per_fold"][str(f)]["eval"][readout][metric]
        out.append(np.array([val]) if val is not None else np.array([]))
    return out


def get_b3_fold_draws(host, readout, metric):
    with open(OUT / "baselines" / "fold_rotation_all_baselines.json") as f:
        d = json.load(f)
    per_fold = d["b3_free_host_embedding"]["summary"][host][readout]["per_fold"]
    return [np.array([pf[metric]]) for pf in per_fold]


def main():
    print("H-DIAGNOSTIC: CNN zero-shot (N=0, LOHO base) vs B3 free-embedding")
    print("Fallback restated: mean of training-host one-hot embeddings (Gate 3), chosen because")
    print("with only 2 training hosts per LOHO fold, nearest-neighbor reduces to an arbitrary pick.")
    print("=" * 70)

    results = {}
    for host in HOSTS:
        results[host] = {}
        for readout in READOUTS:
            results[host][readout] = {}
            b3_rho = boot.bootstrap_ci_90(get_b3_fold_draws(host, readout, "spearman_rho"))
            b3_mcc = boot.bootstrap_ci_90(get_b3_fold_draws(host, readout, "mcc"))
            b3_auc = boot.bootstrap_ci_90(get_b3_fold_draws(host, readout, "auc"))
            results[host][readout]["b3_free_embedding"] = {"spearman_rho": b3_rho, "mcc": b3_mcc, "auc": b3_auc}
            for variant in VARIANTS:
                cnn_rho = boot.bootstrap_ci_90(get_cnn_zeroshot_fold_draws(host, variant, readout, "spearman_rho"))
                cnn_mcc = boot.bootstrap_ci_90(get_cnn_zeroshot_fold_draws(host, variant, readout, "mcc"))
                cnn_auc = boot.bootstrap_ci_90(get_cnn_zeroshot_fold_draws(host, variant, readout, "auc"))
                v = boot.verdict(cnn_rho, b3_rho)
                results[host][readout][variant] = {"spearman_rho": cnn_rho, "mcc": cnn_mcc, "auc": cnn_auc, "verdict_vs_b3": v}
                print(f"  {host} {readout} {variant}: rho={cnn_rho['mean']:.3f} [{cnn_rho['lower']:.3f},{cnn_rho['upper']:.3f}] "
                      f"vs B3={b3_rho['mean']:.3f} [{b3_rho['lower']:.3f},{b3_rho['upper']:.3f}] verdict={v}")

    with open(RESULTS / "gate5_hdiagnostic_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nWrote {RESULTS / 'gate5_hdiagnostic_results.json'}")


if __name__ == "__main__":
    main()
