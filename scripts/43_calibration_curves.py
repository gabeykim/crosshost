"""
GATE 4 - Task 3: calibration curves for the cross-host model. For each
(held-out host, host-feature variant) LOHO base model trained in Task 2
(scripts/42), fine-tune on N examples drawn from the held-out host's own
data and evaluate on the remainder. N in {0, 10, 30, 100, 300, 1000, 3000},
>=10 draws per N (except N=0, a single deterministic evaluation of the
unmodified LOHO base model -- this is exactly Task 2's own held-out result,
reused rather than recomputed).

TRANSFER MECHANISM, per charter: freeze the CNN trunk (conv1-4/bn1-4),
refit only the FiLM generator + fc_shared + output heads for N<=100;
unfreeze the top conv layer (conv4/bn4) in addition for N>=300. Both
mechanisms are run at N=300 (the crossover point) so the difference is
visible, per instruction -- this yields 7 result "cells" per (base model,
readout... no, per base model, since both readouts train jointly): N=10,
N=30, N=100 (frozen/head_only), N=300-frozen, N=300-topconv, N=1000
(topconv), N=3000 (topconv).

Same N grid and >=10-draw protocol as Gate 3's B2 baseline (scripts/31) and
Gate 3.5's fold rotation (scripts/36), so the curves are directly
comparable to the k-mer/linear baseline's calibration-cost curves.
"""
import sys
import time
import json
import numpy as np
import torch
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
m = import_module("40_film_cnn_model")
loho = import_module("42_loho_training")

OUT = Path(__file__).resolve().parent.parent / "out"
MODELS_DIR = OUT / "models"
DEVICE = "mps" if torch.backends.mps.is_available() else "cpu"

HOSTS = ["EC", "BS", "PA"]
N_FOLDS = 5
N_VALUES_NONZERO = [10, 30, 100, 300, 1000, 3000]
N_DRAWS = 10
FT_EPOCHS_SMALL = 20   # N<=100, frozen trunk -- cheap, converges fast
FT_EPOCHS_LARGE = 15   # N>=300, top conv unfrozen -- more params, cap epochs for cost
FT_PATIENCE = 6


def cell_key(N, mechanism):
    return f"N{N}_{mechanism}"


CELLS = [(10, "head_only"), (30, "head_only"), (100, "head_only"),
         (300, "head_only"), (300, "top_conv"), (1000, "top_conv"), (3000, "top_conv")]


def load_base_model(held_out, variant_name, test_fold, host_dim):
    model = m.FiLMSequenceCNN(host_dim=host_dim, use_film=True)
    state = torch.load(MODELS_DIR / f"loho_{held_out}_{variant_name}_fold{test_fold}.pt",
                        map_location="cpu")
    model.load_state_dict(state)
    return model


def main():
    t_start = time.time()
    g_z, gcols, p_z, pcols, imputed = m.load_host_features()
    variants = {"genomic": (g_z, len(gcols)), "physiology": (p_z, len(pcols))}

    all_results = {}
    n_cells_done = 0
    n_cells_total = len(HOSTS) * len(variants) * N_FOLDS * len(CELLS) * N_DRAWS
    for held_out in HOSTS:
        all_results[held_out] = {}
        for variant_name, (feat_df, host_dim) in variants.items():
            all_results[held_out][variant_name] = {}
            host_vec_full = feat_df.loc[held_out].values.astype(np.float32)
            for test_fold in range(N_FOLDS):
                # held-out host's OWN data for this fold: fine-tune pool = fold!=test_fold
                # (mirrors Task 2's fold semantics -- the held-out host's fold==test_fold
                # rows are the ones Task 2 evaluated on with N=0; here we additionally draw
                # fine-tuning examples from the SAME held-out host, so those examples must
                # come from a DIFFERENT partition than the eval set to avoid testing on
                # what was fine-tuned on. Fine-tune pool = held-out host's fold!=test_fold;
                # eval = held-out host's fold==test_fold, held fixed -- identical eval set
                # to Task 2's N=0 result for this (host, variant, fold), enabling direct
                # comparison across the whole N=0..3000 curve.)
                seq_pool, hv_pool, targets_pool, masks_pool, _ = m.build_pooled_arrays(
                    [held_out], feat_df, fold_filter_fn=lambda f, tf=test_fold: f != tf)
                seq_eval, hv_eval, targets_eval, masks_eval, _ = m.build_pooled_arrays(
                    [held_out], feat_df, fold_filter_fn=lambda f, tf=test_fold: f == tf)
                pool_n = len(seq_pool)

                fold_cell_results = {}
                for N, mechanism in CELLS:
                    draws = []
                    n_draws_here = N_DRAWS
                    epochs = FT_EPOCHS_SMALL if mechanism == "head_only" else FT_EPOCHS_LARGE
                    for draw_i in range(n_draws_here):
                        rng = np.random.default_rng(10_000 * N + 7 * draw_i + test_fold)
                        n_sample = min(N, pool_n)
                        sample_idx = rng.choice(pool_n, size=n_sample, replace=False)
                        seq_ft = seq_pool[sample_idx]
                        hv_ft = hv_pool[sample_idx]
                        targets_ft = {k: v[sample_idx] for k, v in targets_pool.items()}
                        masks_ft = {k: v[sample_idx] for k, v in masks_pool.items()}

                        base_model = load_base_model(held_out, variant_name, test_fold, host_dim)
                        ft_model, epochs_run = m.train_film_model(
                            seq_ft, hv_ft, targets_ft, masks_ft, host_dim=host_dim,
                            epochs=epochs, batch_size=min(1024, max(n_sample, 1)), seed=draw_i,
                            device=DEVICE, use_film=True, patience=FT_PATIENCE,
                            freeze_mode=mechanism, init_model=base_model)
                        preds = m.predict_film(ft_model, seq_eval, hv_eval)
                        eval_result = loho.evaluate(preds, targets_eval, masks_eval)
                        draws.append({"eval": eval_result, "epochs_run": epochs_run, "n_actual": n_sample})
                        n_cells_done += 1
                    fold_cell_results[cell_key(N, mechanism)] = draws
                    print(f"[{n_cells_done}/{n_cells_total} draws] held_out={held_out} variant={variant_name} "
                          f"fold={test_fold} N={N} mech={mechanism}: "
                          f"tx_rho_mean={np.nanmean([d['eval']['transcription']['spearman_rho'] for d in draws if d['eval']['transcription']['spearman_rho'] is not None] or [np.nan]):.3f}")

                all_results[held_out][variant_name][test_fold] = fold_cell_results
                with open(OUT / "gate4_calibration_curves.json", "w") as f:
                    json.dump(all_results, f, indent=2, default=str)

    print(f"\nTotal wall time: {time.time()-t_start:.1f}s")
    with open(OUT / "gate4_calibration_curves.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"Wrote {OUT / 'gate4_calibration_curves.json'}")


if __name__ == "__main__":
    main()
