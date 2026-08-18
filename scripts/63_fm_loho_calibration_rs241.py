"""
GATE 6 - Task 1/2/3: full evaluation of a frozen-FM-embedding + lightweight-
head system, matching Gate 5's protocol exactly: LOHO zero-shot (scripts/42
pattern), N-calibration curve (scripts/43 pattern), RS241 zero-shot
(scripts/50 pattern). Same folds, same held-out hosts, same metrics
(scripts/42's evaluate()), same fold semantics.

ONE MECHANISM, not two, and this is a real protocol difference from the CNN
disclosed here rather than silently matched: the CNN's calibration curve had
two transfer mechanisms (frozen-trunk head_only vs unfrozen-conv4 top_conv)
because the CNN's OWN trunk was being adapted. Here there is no trunk to
adapt -- the FM embedding is fixed by construction (that is the entire
"frozen embeddings" premise Task 1 specifies) -- only the lightweight head
(scripts/62's FMHeadMLP) is ever trained, at every N. Calling this
"head_only" throughout for naming consistency with the rest of the project's
result tables, but there is no "top_conv"-equivalent second mechanism for a
frozen-embedding system and none is reported.

No host_vec / variant split: FM embeddings carry no host information (same
structural position as the sequence-only ablation, scripts/58-60), so there
is a single result per (host, readout, fold) -- not a genomic/physiology
pair.
"""
import sys
import time
import json
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
m40 = import_module("40_film_cnn_model")
loho = import_module("42_loho_training")
fm = import_module("62_fm_head_model")

OUT = Path(__file__).resolve().parent.parent / "out"
MODELS_DIR = OUT / "models"
MODELS_DIR.mkdir(parents=True, exist_ok=True)
DATA = Path(__file__).resolve().parent.parent / "data"

HOSTS = ["EC", "BS", "PA"]
N_FOLDS = 5
DEVICE = None  # resolved per-call inside fm.train_fm_head via torch default; set explicitly below

N_VALUES_NONZERO = [10, 30, 100, 300, 1000, 3000]
N_DRAWS = 10
FT_EPOCHS = 60
FT_PATIENCE = 10

RS241_CONFIGS = {"PRIMARY": ["EC", "BS", "PA"], "SECONDARY": ["EC", "PA"]}
RS241_HOSTS = ["SE", "VN", "CG"]
N_SEEDS = 5


def build_targets_masks_for_rows(host, row_idx, all_targets_masks):
    t, msk = all_targets_masks[host]
    targets = {k: v[row_idx] for k, v in t.items()}
    masks = {k: v[row_idx] for k, v in msk.items()}
    return targets, masks


def load_all_host_target_arrays():
    """Per host: (targets, masks) dicts row-aligned to the shared 29042-row
    library order (identical across hosts, verified in scripts/61)."""
    out = {}
    for h in HOSTS:
        d = dict(np.load(DATA / "baseline_cache" / f"{h}_baseline_data.npz", allow_pickle=True))
        t, msk = m40.build_target_arrays(d)
        out[h] = (t, msk)
    return out


def load_rs241_targets():
    import pandas as pd
    rs241 = pd.read_parquet(DATA / "rs241.parquet")
    rs241["id"] = rs241["id"].astype(str)
    lib = pd.read_parquet(DATA / "three_host_library.parquet")
    lib["OLIGO ID"] = lib["OLIGO ID"].astype(str)
    seq_map = dict(zip(lib["OLIGO ID"], lib["Regulatory Sequence"]))
    rs241["has_sequence"] = rs241["id"].isin(seq_map)
    rs241_seq = rs241[rs241["has_sequence"]].copy().reset_index(drop=True)
    return rs241_seq


def evaluate_rs241_fm(preds, rs241_seq, held_out, config_tag):
    from scipy.stats import spearmanr
    from sklearn.metrics import matthews_corrcoef, roc_auc_score
    results = {}
    for readout, active_col_prefix, strength_col, usable_suffix in [
        ("transcription", "tx", f"tx_log2_{held_out}", "tx"),
        ("translation", "tl", f"tl_log10_{held_out}", "tl"),
    ]:
        usable_col = f"usable_heldout_{held_out}_train_{config_tag}_{usable_suffix}"
        usable_mask = rs241_seq[usable_col].fillna(False).values.astype(bool)
        n_test = int(usable_mask.sum())
        if n_test < 3:
            results[readout] = {"n_test": n_test, "mcc": None, "auc": None, "spearman_rho": None}
            continue
        y_strength = rs241_seq.loc[usable_mask, strength_col].values.astype(np.float64)
        pred_strength = preds[f"{active_col_prefix}_strength"][usable_mask]
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
        results[readout] = {
            "n_test": n_test,
            "mcc": float(mcc) if mcc is not None and mcc == mcc else None,
            "auc": float(auc) if auc is not None and auc == auc else None,
            "spearman_rho": float(rho) if rho is not None and rho == rho else None,
        }
    return results


def main(model_tag, device=None):
    t_start = time.time()
    device = device or ("mps" if __import__("torch").backends.mps.is_available() else "cpu")
    print(f"=== FM head evaluation: model_tag={model_tag}, device={device} ===")

    oligo_ids, emb, fold = fm.load_library_embeddings(model_tag)
    emb_dim = emb.shape[1]
    all_tm = load_all_host_target_arrays()
    print(f"Library: {len(oligo_ids)} rows, emb_dim={emb_dim}")

    # ---------- Task: LOHO zero-shot (N=0), matches scripts/42 fold semantics ----------
    loho_results = {}
    n_fits = 0
    n_fits_total = len(HOSTS) * N_FOLDS
    for held_out in HOSTS:
        train_hosts = [h for h in HOSTS if h != held_out]
        loho_results[held_out] = {"per_fold": {}}
        for test_fold in range(N_FOLDS):
            t0 = time.time()
            train_idx = np.where(fold != test_fold)[0]
            emb_tr_parts, targets_tr_parts, masks_tr_parts = [], {k: [] for k in m40.TARGET_KEYS}, {k: [] for k in m40.TARGET_KEYS}
            for h in train_hosts:
                t, msk = all_tm[h]
                emb_tr_parts.append(emb[train_idx])
                for k in m40.TARGET_KEYS:
                    targets_tr_parts[k].append(t[k][train_idx])
                    masks_tr_parts[k].append(msk[k][train_idx])
            emb_tr = np.concatenate(emb_tr_parts)
            targets_tr = {k: np.concatenate(v) for k, v in targets_tr_parts.items()}
            masks_tr = {k: np.concatenate(v) for k, v in masks_tr_parts.items()}

            model, epochs_run = fm.train_fm_head(emb_tr, targets_tr, masks_tr, emb_dim=emb_dim,
                                                  epochs=FT_EPOCHS, seed=test_fold, device=device,
                                                  patience=FT_PATIENCE)

            test_idx = np.where(fold == test_fold)[0]
            t_te, msk_te = all_tm[held_out]
            targets_te = {k: v[test_idx] for k, v in t_te.items()}
            masks_te = {k: v[test_idx] for k, v in msk_te.items()}
            preds = fm.predict_fm_head(model, emb[test_idx])
            eval_result = loho.evaluate(preds, targets_te, masks_te)

            elapsed = time.time() - t0
            n_fits += 1
            print(f"[{n_fits}/{n_fits_total}] LOHO held_out={held_out} fold={test_fold}: "
                  f"epochs_run={epochs_run}, {elapsed:.2f}s, n_train={len(emb_tr)}, n_test={len(test_idx)}")
            for ro, r in eval_result.items():
                print(f"    {ro}: mcc={r['mcc']}, auc={r['auc']}, rho={r['spearman_rho']}")

            loho_results[held_out]["per_fold"][test_fold] = {
                "eval": eval_result, "epochs_run": epochs_run,
                "wall_clock_sec": elapsed, "n_train": int(len(emb_tr)), "n_test": int(len(test_idx)),
            }
            import torch
            torch.save(model.state_dict(), MODELS_DIR / f"fm_{model_tag}_loho_{held_out}_fold{test_fold}.pt")
    with open(OUT / f"gate6_{model_tag}_loho_results.json", "w") as f:
        json.dump(loho_results, f, indent=2, default=str)
    print(f"Wrote gate6_{model_tag}_loho_results.json")

    # ---------- Task: N-calibration curve, matches scripts/43 fold semantics ----------
    calib_results = {}
    n_cells = 0
    n_cells_total = len(HOSTS) * N_FOLDS * len(N_VALUES_NONZERO) * N_DRAWS
    for held_out in HOSTS:
        calib_results[held_out] = {}
        t_te_full, msk_te_full = all_tm[held_out]
        for test_fold in range(N_FOLDS):
            pool_idx = np.where(fold != test_fold)[0]
            eval_idx = np.where(fold == test_fold)[0]
            emb_pool = emb[pool_idx]
            targets_pool = {k: v[pool_idx] for k, v in t_te_full.items()}
            masks_pool = {k: v[pool_idx] for k, v in msk_te_full.items()}
            emb_eval = emb[eval_idx]
            targets_eval = {k: v[eval_idx] for k, v in t_te_full.items()}
            masks_eval = {k: v[eval_idx] for k, v in msk_te_full.items()}
            pool_n = len(emb_pool)

            base_model_path = MODELS_DIR / f"fm_{model_tag}_loho_{held_out}_fold{test_fold}.pt"
            fold_cell_results = {}
            for N in N_VALUES_NONZERO:
                draws = []
                for draw_i in range(N_DRAWS):
                    rng = np.random.default_rng(10_000 * N + 7 * draw_i + test_fold)
                    n_sample = min(N, pool_n)
                    sample_idx = rng.choice(pool_n, size=n_sample, replace=False)
                    emb_ft = emb_pool[sample_idx]
                    targets_ft = {k: v[sample_idx] for k, v in targets_pool.items()}
                    masks_ft = {k: v[sample_idx] for k, v in masks_pool.items()}

                    import torch
                    base_model = fm.FMHeadMLP(emb_dim)
                    base_model.load_state_dict(torch.load(base_model_path, map_location="cpu"))
                    ft_model, epochs_run = fm.train_fm_head(
                        emb_ft, targets_ft, masks_ft, emb_dim=emb_dim, epochs=FT_EPOCHS,
                        batch_size=min(1024, max(n_sample, 1)), seed=draw_i, device=device,
                        patience=FT_PATIENCE, init_model=base_model)
                    preds = fm.predict_fm_head(ft_model, emb_eval)
                    eval_result = loho.evaluate(preds, targets_eval, masks_eval)
                    draws.append({"eval": eval_result, "epochs_run": epochs_run, "n_actual": n_sample})
                    n_cells += 1
                fold_cell_results[f"N{N}_head_only"] = draws
                tx_rho = np.nanmean([d["eval"]["transcription"]["spearman_rho"] for d in draws
                                      if d["eval"]["transcription"]["spearman_rho"] is not None] or [np.nan])
                print(f"[{n_cells}/{n_cells_total} draws] held_out={held_out} fold={test_fold} N={N}: "
                      f"tx_rho_mean={tx_rho:.3f}")
            calib_results[held_out][test_fold] = fold_cell_results
            with open(OUT / f"gate6_{model_tag}_calibration_curves.json", "w") as f:
                json.dump(calib_results, f, indent=2, default=str)
    with open(OUT / f"gate6_{model_tag}_calibration_curves.json", "w") as f:
        json.dump(calib_results, f, indent=2, default=str)
    print(f"Wrote gate6_{model_tag}_calibration_curves.json")

    # ---------- Task: RS241 zero-shot, matches scripts/50 ----------
    rs241_ids, rs241_emb = fm.load_rs241_embeddings(model_tag)
    rs241_seq = load_rs241_targets()
    assert list(rs241_seq["id"]) == list(rs241_ids), "RS241 id order mismatch between targets and embeddings"

    rs241_results = {}
    n_rs = 0
    n_rs_total = len(RS241_CONFIGS) * N_SEEDS
    for config_name, train_hosts in RS241_CONFIGS.items():
        config_tag = "_".join(train_hosts)
        rs241_results[config_name] = {"per_seed": {}}
        emb_tr_parts, targets_tr_parts, masks_tr_parts = [], {k: [] for k in m40.TARGET_KEYS}, {k: [] for k in m40.TARGET_KEYS}
        for h in train_hosts:
            t, msk = all_tm[h]
            emb_tr_parts.append(emb)  # all folds -- no exclusion, matches scripts/50
            for k in m40.TARGET_KEYS:
                targets_tr_parts[k].append(t[k])
                masks_tr_parts[k].append(msk[k])
        emb_tr = np.concatenate(emb_tr_parts)
        targets_tr = {k: np.concatenate(v) for k, v in targets_tr_parts.items()}
        masks_tr = {k: np.concatenate(v) for k, v in masks_tr_parts.items()}
        print(f"\n=== RS241 config={config_name}: pool N={len(emb_tr)} ===")

        for seed in range(N_SEEDS):
            t0 = time.time()
            model, epochs_run = fm.train_fm_head(emb_tr, targets_tr, masks_tr, emb_dim=emb_dim,
                                                  epochs=FT_EPOCHS, seed=seed, device=device,
                                                  patience=FT_PATIENCE)
            preds = fm.predict_fm_head(model, rs241_emb)
            seed_result = {}
            for held_out in RS241_HOSTS:
                seed_result[held_out] = evaluate_rs241_fm(preds, rs241_seq, held_out, config_tag)
            elapsed = time.time() - t0
            n_rs += 1
            print(f"[{n_rs}/{n_rs_total}] config={config_name} seed={seed}: {elapsed:.2f}s, epochs_run={epochs_run}")
            for h in RS241_HOSTS:
                print(f"    {h}: {seed_result[h]}")
            rs241_results[config_name]["per_seed"][seed] = {
                "eval": seed_result, "epochs_run": epochs_run, "wall_clock_sec": elapsed, "n_train": int(len(emb_tr))
            }
        with open(OUT / f"gate6_{model_tag}_rs241_results.json", "w") as f:
            json.dump(rs241_results, f, indent=2, default=str)
    with open(OUT / f"gate6_{model_tag}_rs241_results.json", "w") as f:
        json.dump(rs241_results, f, indent=2, default=str)
    print(f"Wrote gate6_{model_tag}_rs241_results.json")

    print(f"\nTotal wall time for {model_tag}: {time.time()-t_start:.1f}s")


if __name__ == "__main__":
    tag = sys.argv[1] if len(sys.argv) > 1 else "dnabert2"
    main(tag)
