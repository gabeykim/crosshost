"""
GATE 5.5 - Task 1: full pipeline for the sequence-only ablation --
LOHO training (15 fits), calibration curves (1050 fine-tune draws), and
RS241 zero-shot (10 fits). One combined script so it can run unattended as
a single background job. Reuses scripts/40's data-loading (build_pooled_arrays,
build_target_arrays) and scripts/42/43/50's evaluation logic -- only the
model itself (scripts/58) differs, since the host_vec these functions return
is simply never passed to the sequence-only model.
"""
import sys
import time
import json
import numpy as np
import pandas as pd
import torch
from pathlib import Path
from scipy.stats import spearmanr
from sklearn.metrics import matthews_corrcoef, roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
m = import_module("40_film_cnn_model")
so = import_module("58_sequence_only_model")
loho = import_module("42_loho_training")
cal43 = import_module("43_calibration_curves")
rs241mod = import_module("50_rs241_train_and_eval")

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"
MODELS_DIR = OUT / "models"
MODELS_DIR.mkdir(parents=True, exist_ok=True)
DEVICE = "mps" if torch.backends.mps.is_available() else "cpu"

HOSTS = ["EC", "BS", "PA"]
N_FOLDS = 5
EPOCHS = 35
PATIENCE = 8
BATCH_SIZE = 1024
CELLS = cal43.CELLS
N_DRAWS = cal43.N_DRAWS
CONFIGS = {"PRIMARY": ["EC", "BS", "PA"], "SECONDARY": ["EC", "PA"]}
RS241_HOSTS = ["SE", "VN", "CG"]
N_SEEDS = 5


def stage1_loho(g_z):
    print("=" * 70)
    print("STAGE 1: LOHO training, sequence-only, 3 hosts x 5 folds = 15 fits")
    print("=" * 70)
    all_results = {}
    n_done = 0
    for held_out in HOSTS:
        train_hosts = [h for h in HOSTS if h != held_out]
        all_results[held_out] = {"per_fold": {}}
        for test_fold in range(N_FOLDS):
            t0 = time.time()
            seq_tr, _, targets_tr, masks_tr, _ = m.build_pooled_arrays(
                train_hosts, g_z, fold_filter_fn=lambda f, tf=test_fold: f != tf)
            model, epochs_run = so.train_seqonly_model(
                seq_tr, targets_tr, masks_tr, epochs=EPOCHS, batch_size=BATCH_SIZE,
                seed=test_fold, device=DEVICE, patience=PATIENCE)
            seq_te, _, targets_te, masks_te, _ = m.build_pooled_arrays(
                [held_out], g_z, fold_filter_fn=lambda f, tf=test_fold: f == tf)
            preds = so.predict_seqonly(model, seq_te)
            eval_result = loho.evaluate(preds, targets_te, masks_te)
            elapsed = time.time() - t0
            n_done += 1
            print(f"[{n_done}/15] held_out={held_out} fold={test_fold}: {elapsed:.1f}s, epochs_run={epochs_run}")
            for ro, r in eval_result.items():
                print(f"    {ro}: mcc={r['mcc']}, auc={r['auc']}, rho={r['spearman_rho']}")
            all_results[held_out]["per_fold"][test_fold] = {
                "eval": eval_result, "epochs_run": epochs_run, "wall_clock_sec": elapsed,
                "n_train": int(len(seq_tr)), "n_test": int(len(seq_te)),
            }
            torch.save(model.state_dict(), MODELS_DIR / f"seqonly_loho_{held_out}_fold{test_fold}.pt")
            with open(OUT / "gate5_5_seqonly_loho_results.json", "w") as f:
                json.dump(all_results, f, indent=2, default=str)
    print(f"Stage 1 wall time: {time.time()-t0:.1f}s (last fit)")
    return all_results


def load_seqonly_base(held_out, test_fold):
    model = so.SequenceOnlyCNN()
    state = torch.load(MODELS_DIR / f"seqonly_loho_{held_out}_fold{test_fold}.pt", map_location="cpu")
    model.load_state_dict(state)
    return model


def stage2_calibration(g_z):
    print("=" * 70)
    print("STAGE 2: calibration curves, sequence-only, 3 hosts x 5 folds x 7 cells x 10 draws")
    print("=" * 70)
    t_start = time.time()
    all_results = {}
    n_cells_done = 0
    n_cells_total = len(HOSTS) * N_FOLDS * len(CELLS) * N_DRAWS
    for held_out in HOSTS:
        all_results[held_out] = {}
        for test_fold in range(N_FOLDS):
            seq_pool, _, targets_pool, masks_pool, _ = m.build_pooled_arrays(
                [held_out], g_z, fold_filter_fn=lambda f, tf=test_fold: f != tf)
            seq_eval, _, targets_eval, masks_eval, _ = m.build_pooled_arrays(
                [held_out], g_z, fold_filter_fn=lambda f, tf=test_fold: f == tf)
            pool_n = len(seq_pool)
            fold_cell_results = {}
            for N, mechanism in CELLS:
                draws = []
                epochs = cal43.FT_EPOCHS_SMALL if mechanism == "head_only" else cal43.FT_EPOCHS_LARGE
                for draw_i in range(N_DRAWS):
                    rng = np.random.default_rng(10_000 * N + 7 * draw_i + test_fold)
                    n_sample = min(N, pool_n)
                    sample_idx = rng.choice(pool_n, size=n_sample, replace=False)
                    seq_ft = seq_pool[sample_idx]
                    targets_ft = {k: v[sample_idx] for k, v in targets_pool.items()}
                    masks_ft = {k: v[sample_idx] for k, v in masks_pool.items()}
                    base_model = load_seqonly_base(held_out, test_fold)
                    ft_model, epochs_run = so.train_seqonly_model(
                        seq_ft, targets_ft, masks_ft, epochs=epochs,
                        batch_size=min(1024, max(n_sample, 1)), seed=draw_i, device=DEVICE,
                        patience=cal43.FT_PATIENCE, freeze_mode=mechanism, init_model=base_model)
                    preds = so.predict_seqonly(ft_model, seq_eval)
                    eval_result = loho.evaluate(preds, targets_eval, masks_eval)
                    draws.append({"eval": eval_result, "epochs_run": epochs_run, "n_actual": n_sample})
                    n_cells_done += 1
                fold_cell_results[cal43.cell_key(N, mechanism)] = draws
                tx_rhos = [d["eval"]["transcription"]["spearman_rho"] for d in draws if d["eval"]["transcription"]["spearman_rho"] is not None]
                print(f"[{n_cells_done}/{n_cells_total}] held_out={held_out} fold={test_fold} N={N} mech={mechanism}: "
                      f"tx_rho_mean={np.mean(tx_rhos) if tx_rhos else None}")
            all_results[held_out][test_fold] = fold_cell_results
            with open(OUT / "gate5_5_seqonly_calibration_curves.json", "w") as f:
                json.dump(all_results, f, indent=2, default=str)
    print(f"Stage 2 wall time: {time.time()-t_start:.1f}s")
    return all_results


def stage2b_n100_topconv(g_z):
    print("=" * 70)
    print("STAGE 2b: N=100 top_conv supplement, sequence-only")
    print("=" * 70)
    t0 = time.time()
    results = {}
    for held_out in HOSTS:
        results[held_out] = {}
        for test_fold in range(N_FOLDS):
            seq_pool, _, targets_pool, masks_pool, _ = m.build_pooled_arrays(
                [held_out], g_z, fold_filter_fn=lambda f, tf=test_fold: f != tf)
            seq_eval, _, targets_eval, masks_eval, _ = m.build_pooled_arrays(
                [held_out], g_z, fold_filter_fn=lambda f, tf=test_fold: f == tf)
            pool_n = len(seq_pool)
            draws = []
            N = 100
            for draw_i in range(N_DRAWS):
                rng = np.random.default_rng(10_000 * N + 7 * draw_i + test_fold)
                n_sample = min(N, pool_n)
                sample_idx = rng.choice(pool_n, size=n_sample, replace=False)
                seq_ft = seq_pool[sample_idx]
                targets_ft = {k: v[sample_idx] for k, v in targets_pool.items()}
                masks_ft = {k: v[sample_idx] for k, v in masks_pool.items()}
                base_model = load_seqonly_base(held_out, test_fold)
                ft_model, epochs_run = so.train_seqonly_model(
                    seq_ft, targets_ft, masks_ft, epochs=cal43.FT_EPOCHS_LARGE,
                    batch_size=min(1024, max(n_sample, 1)), seed=draw_i, device=DEVICE,
                    patience=cal43.FT_PATIENCE, freeze_mode="top_conv", init_model=base_model)
                preds = so.predict_seqonly(ft_model, seq_eval)
                eval_result = loho.evaluate(preds, targets_eval, masks_eval)
                draws.append({"eval": eval_result, "epochs_run": epochs_run, "n_actual": n_sample})
            results[held_out][test_fold] = draws
            print(f"  held_out={held_out} fold={test_fold} N=100 top_conv done")
    with open(OUT / "gate5_5_seqonly_n100_topconv.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"Stage 2b wall time: {time.time()-t0:.1f}s")
    return results


def stage3_rs241(g_z):
    print("=" * 70)
    print("STAGE 3: RS241 zero-shot, sequence-only, 2 configs x 5 seeds = 10 fits")
    print("=" * 70)
    t_start = time.time()
    rs241_seq, onehot = rs241mod.load_rs241_eval_arrays()
    print(f"RS241 sequence-available rows: {len(rs241_seq)} / 241")

    all_results = {}
    n_done = 0
    for config_name, train_hosts in CONFIGS.items():
        config_tag = "_".join(train_hosts)
        all_results[config_name] = {"per_seed": {}}
        seq_tr, _, targets_tr, masks_tr, _ = m.build_pooled_arrays(train_hosts, g_z, fold_filter_fn=None)
        print(f"\n=== config={config_name}: pool N={len(seq_tr)} ===")
        for seed in range(N_SEEDS):
            t0 = time.time()
            model, epochs_run = so.train_seqonly_model(
                seq_tr, targets_tr, masks_tr, epochs=EPOCHS, batch_size=BATCH_SIZE,
                seed=seed, device=DEVICE, patience=PATIENCE)
            elapsed = time.time() - t0
            n_done += 1

            preds = so.predict_seqonly(model, onehot)
            seed_result = {}
            for held_out in RS241_HOSTS:
                res_by_readout = {}
                for readout, active_prefix, strength_col in [
                    ("transcription", "tx", f"tx_log2_{held_out}"),
                    ("translation", "tl", f"tl_log10_{held_out}"),
                ]:
                    usable_suffix = "tx" if readout == "transcription" else "tl"
                    usable_col = f"usable_heldout_{held_out}_train_{config_tag}_{usable_suffix}"
                    usable_mask = rs241_seq[usable_col].fillna(False).values.astype(bool)
                    n_test = int(usable_mask.sum())
                    if n_test < 3:
                        res_by_readout[readout] = {"n_test": n_test, "mcc": None, "auc": None, "spearman_rho": None}
                        continue
                    y_strength = rs241_seq.loc[usable_mask, strength_col].values.astype(np.float64)
                    pred_strength = preds[f"{active_prefix}_strength"][usable_mask]
                    rho = spearmanr(y_strength, pred_strength).correlation if len(set(np.round(y_strength, 6))) > 1 else None
                    if readout == "transcription":
                        y_active = (y_strength > 0).astype(int)
                        prob = 1 / (1 + np.exp(-preds["tx_active_logit"][usable_mask]))
                        pred_class = (prob > 0.5).astype(int)
                        mcc = matthews_corrcoef(y_active, pred_class) if len(set(y_active)) > 1 else None
                        try:
                            auc = roc_auc_score(y_active, prob) if len(set(y_active)) > 1 else None
                        except Exception:
                            auc = None
                    else:
                        mcc, auc = None, None
                    res_by_readout[readout] = {"n_test": n_test,
                                                "mcc": float(mcc) if mcc is not None and mcc == mcc else None,
                                                "auc": float(auc) if auc is not None and auc == auc else None,
                                                "spearman_rho": float(rho) if rho is not None and rho == rho else None}
                seed_result[held_out] = res_by_readout

            all_results[config_name]["per_seed"][seed] = {
                "eval": seed_result, "epochs_run": epochs_run, "wall_clock_sec": elapsed, "n_train": int(len(seq_tr))
            }
            print(f"[{n_done}/10] config={config_name} seed={seed}: {elapsed:.1f}s, epochs_run={epochs_run}")
            for h in RS241_HOSTS:
                print(f"    {h}: {seed_result[h]}")
            torch.save(model.state_dict(), MODELS_DIR / f"seqonly_rs241_{config_name}_seed{seed}.pt")
            with open(OUT / "gate5_5_seqonly_rs241_results.json", "w") as f:
                json.dump(all_results, f, indent=2, default=str)
    print(f"Stage 3 wall time: {time.time()-t_start:.1f}s")
    return all_results


def main():
    t_all = time.time()
    g_z, gcols, p_z, pcols, imputed = m.load_host_features()
    assert imputed == []

    stage1_loho(g_z)
    stage2_calibration(g_z)
    stage2b_n100_topconv(g_z)
    stage3_rs241(g_z)

    print(f"\n=== TOTAL PIPELINE WALL TIME: {time.time()-t_all:.1f}s ===")


if __name__ == "__main__":
    main()
