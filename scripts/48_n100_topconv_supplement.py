"""
GATE 4.5 - Task 2: extend the N=300 both-mechanisms comparison down to N=100.
Gate 4's memo found frozen-trunk (head_only) fine-tuning at N<=300 inherits
and cannot fix a bad LOHO base model's calibration, and separately that
small-N physiology fine-tuning can transiently DEGRADE performance relative
to zero-shot -- both effects live inside N=10-300, which includes H-MAIN-TX's
own evaluation point (N=100). The pre-registration only required both
mechanisms to be compared at N=300; this task adds N=100.

REPORTING EXPANSION, NOT A HYPOTHESIS CHANGE: frozen-trunk (head_only)
remains the pre-registered PRIMARY mechanism for H-MAIN. top_conv at N=100 is
computed and reported as SUPPLEMENTARY context only -- see the dated
amendment appended to out/PREREGISTRATION.md, written before H-MAIN is
evaluated (Gate 5).

Runs for a given `variants` list (so genomic can run immediately, unaffected
by the physiology-vector fix in Task 1, while physiology waits for
scripts/46's retrain to complete first).
"""
import sys
import time
import json
import argparse
import numpy as np
import torch
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
m = import_module("40_film_cnn_model")
loho = import_module("42_loho_training")
cal43 = import_module("43_calibration_curves")

OUT = Path(__file__).resolve().parent.parent / "out"
MODELS_DIR = OUT / "models"
DEVICE = "mps" if torch.backends.mps.is_available() else "cpu"

HOSTS = ["EC", "BS", "PA"]
N_FOLDS = 5
N = 100
MECHANISM = "top_conv"
N_DRAWS = 10
FT_EPOCHS = cal43.FT_EPOCHS_LARGE
FT_PATIENCE = cal43.FT_PATIENCE


def run_variant(variant_name, feat_df, host_dim):
    results = {}
    for held_out in HOSTS:
        results[held_out] = {}
        for test_fold in range(N_FOLDS):
            seq_pool, hv_pool, targets_pool, masks_pool, _ = m.build_pooled_arrays(
                [held_out], feat_df, fold_filter_fn=lambda f, tf=test_fold: f != tf)
            seq_eval, hv_eval, targets_eval, masks_eval, _ = m.build_pooled_arrays(
                [held_out], feat_df, fold_filter_fn=lambda f, tf=test_fold: f == tf)
            pool_n = len(seq_pool)

            draws = []
            for draw_i in range(N_DRAWS):
                rng = np.random.default_rng(10_000 * N + 7 * draw_i + test_fold)
                n_sample = min(N, pool_n)
                sample_idx = rng.choice(pool_n, size=n_sample, replace=False)
                seq_ft, hv_ft = seq_pool[sample_idx], hv_pool[sample_idx]
                targets_ft = {k: v[sample_idx] for k, v in targets_pool.items()}
                masks_ft = {k: v[sample_idx] for k, v in masks_pool.items()}

                base_model = cal43.load_base_model(held_out, variant_name, test_fold, host_dim)
                ft_model, epochs_run = m.train_film_model(
                    seq_ft, hv_ft, targets_ft, masks_ft, host_dim=host_dim,
                    epochs=FT_EPOCHS, batch_size=min(1024, max(n_sample, 1)), seed=draw_i,
                    device=DEVICE, use_film=True, patience=FT_PATIENCE,
                    freeze_mode=MECHANISM, init_model=base_model)
                preds = m.predict_film(ft_model, seq_eval, hv_eval)
                eval_result = loho.evaluate(preds, targets_eval, masks_eval)
                draws.append({"eval": eval_result, "epochs_run": epochs_run, "n_actual": n_sample})
            results[held_out][test_fold] = draws
            tx_rhos = [d["eval"]["transcription"]["spearman_rho"] for d in draws if d["eval"]["transcription"]["spearman_rho"] is not None]
            print(f"  held_out={held_out} variant={variant_name} fold={test_fold} N=100 top_conv: "
                  f"tx_rho_mean={np.mean(tx_rhos) if tx_rhos else None}")
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--variants", nargs="+", default=["genomic", "physiology"])
    args = parser.parse_args()

    t0 = time.time()
    g_z, gcols, p_z, pcols, imputed = m.load_host_features()
    variant_map = {"genomic": (g_z, len(gcols)), "physiology": (p_z, len(pcols))}

    out_path = OUT / "gate4_5_n100_topconv_supplement.json"
    all_results = json.load(open(out_path)) if out_path.exists() else {}

    for variant_name in args.variants:
        feat_df, host_dim = variant_map[variant_name]
        print(f"=== Running N=100 top_conv supplement for variant={variant_name} ===")
        all_results[variant_name] = run_variant(variant_name, feat_df, host_dim)
        with open(out_path, "w") as f:
            json.dump(all_results, f, indent=2, default=str)

    print(f"\nTotal wall time: {time.time()-t0:.1f}s")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
