"""
GATE 4.5 - Task 5: regenerate the physiology-variant slice of Gate 4's
calibration-curve grid (scripts/43) using the retrained (scripts/46)
checkpoints. Genomic-variant entries in out/gate4_calibration_curves.json
are left untouched -- unaffected by the growth-rate vector fix.
"""
import sys
import time
import json
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
m = import_module("40_film_cnn_model")
loho = import_module("42_loho_training")
cal43 = import_module("43_calibration_curves")

OUT = Path(__file__).resolve().parent.parent / "out"
HOSTS = cal43.HOSTS
N_FOLDS = cal43.N_FOLDS
CELLS = cal43.CELLS
N_DRAWS = cal43.N_DRAWS
DEVICE = cal43.DEVICE


def main():
    t_start = time.time()
    g_z, gcols, p_z, pcols, imputed = m.load_host_features()
    assert imputed == [], "Expected no imputed physiology cells"

    with open(OUT / "gate4_calibration_curves.json") as f:
        all_results = json.load(f)

    n_cells_done = 0
    n_cells_total = len(HOSTS) * N_FOLDS * len(CELLS) * N_DRAWS
    for held_out in HOSTS:
        all_results[held_out]["physiology"] = {}
        for test_fold in range(N_FOLDS):
            seq_pool, hv_pool, targets_pool, masks_pool, _ = m.build_pooled_arrays(
                [held_out], p_z, fold_filter_fn=lambda f, tf=test_fold: f != tf)
            seq_eval, hv_eval, targets_eval, masks_eval, _ = m.build_pooled_arrays(
                [held_out], p_z, fold_filter_fn=lambda f, tf=test_fold: f == tf)
            pool_n = len(seq_pool)

            fold_cell_results = {}
            for N, mechanism in CELLS:
                draws = []
                epochs = cal43.FT_EPOCHS_SMALL if mechanism == "head_only" else cal43.FT_EPOCHS_LARGE
                for draw_i in range(N_DRAWS):
                    rng = np.random.default_rng(10_000 * N + 7 * draw_i + test_fold)
                    n_sample = min(N, pool_n)
                    sample_idx = rng.choice(pool_n, size=n_sample, replace=False)
                    seq_ft, hv_ft = seq_pool[sample_idx], hv_pool[sample_idx]
                    targets_ft = {k: v[sample_idx] for k, v in targets_pool.items()}
                    masks_ft = {k: v[sample_idx] for k, v in masks_pool.items()}

                    base_model = cal43.load_base_model(held_out, "physiology", test_fold, len(pcols))
                    ft_model, epochs_run = m.train_film_model(
                        seq_ft, hv_ft, targets_ft, masks_ft, host_dim=len(pcols),
                        epochs=epochs, batch_size=min(1024, max(n_sample, 1)), seed=draw_i,
                        device=DEVICE, use_film=True, patience=cal43.FT_PATIENCE,
                        freeze_mode=mechanism, init_model=base_model)
                    preds = m.predict_film(ft_model, seq_eval, hv_eval)
                    eval_result = loho.evaluate(preds, targets_eval, masks_eval)
                    draws.append({"eval": eval_result, "epochs_run": epochs_run, "n_actual": n_sample})
                    n_cells_done += 1
                fold_cell_results[cal43.cell_key(N, mechanism)] = draws
                tx_rhos = [d["eval"]["transcription"]["spearman_rho"] for d in draws if d["eval"]["transcription"]["spearman_rho"] is not None]
                print(f"[{n_cells_done}/{n_cells_total}] held_out={held_out} fold={test_fold} N={N} "
                      f"mech={mechanism}: tx_rho_mean={np.mean(tx_rhos) if tx_rhos else None}")

            all_results[held_out]["physiology"][test_fold] = fold_cell_results
            with open(OUT / "gate4_calibration_curves.json", "w") as f:
                json.dump(all_results, f, indent=2, default=str)

    print(f"\nTotal wall time: {time.time()-t_start:.1f}s")
    with open(OUT / "gate4_calibration_curves.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"Updated {OUT / 'gate4_calibration_curves.json'} (physiology entries replaced)")


if __name__ == "__main__":
    main()
