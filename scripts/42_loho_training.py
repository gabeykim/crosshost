"""
GATE 4 - Task 2: leave-one-host-out training, both host-feature variants.

COMPUTE-EFFICIENT DESIGN, disclosed: the task grid is specified as "3 held-out
hosts x 2 readouts x 2 host-feature variants = 12 configurations, each across
5 folds." Because Task 1's architecture trains transcription and translation
JOINTLY from one shared trunk (two output heads on one model, per the
charter's own spec -- "Two output heads: transcription and translation as
separate tasks" -- NOT two separate models), a single trained model produces
BOTH readouts' results at once. The grid therefore requires 3 held-out hosts
x 2 host-feature variants x 5 folds = 30 actual model fits, not 60 -- each
fit is read off into 2 of the 12 reported (host, readout, variant) rows.
This is a direct, faithful consequence of the joint-head architecture Task 1
specifies, not a scope reduction -- 60 independent single-readout fits would
in fact CONTRADICT the joint-architecture spec.

FOLD SEMANTICS, matching Gate 3.5's B3 LOHO convention exactly (scripts/32,
36): for held-out host H and test_fold f, evaluate on H's fold==f rows;
train on the two OTHER hosts' pooled fold!=f rows (not all their data) --
this keeps fold f "quarantined" consistently across the whole pipeline for
that rotation round, matching every other LOHO baseline already built.

Same seeds (=fold index, for reproducibility) and same folds are used across
the genomic/physiology variant pairing for a given (held-out host, fold) --
only the host feature vector fed to the FiLM generator differs -- so
H-SCIENCE's genomic-vs-physiology comparison is not confounded by anything
else.
"""
import sys
import time
import json
import numpy as np
import torch
from pathlib import Path
from scipy.stats import spearmanr
from sklearn.metrics import matthews_corrcoef, roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
m = import_module("40_film_cnn_model")

OUT = Path(__file__).resolve().parent.parent / "out"
MODELS_DIR = OUT / "models"
MODELS_DIR.mkdir(parents=True, exist_ok=True)
DEVICE = "mps" if torch.backends.mps.is_available() else "cpu"

HOSTS = ["EC", "BS", "PA"]
N_FOLDS = 5
EPOCHS = 35
PATIENCE = 8
BATCH_SIZE = 1024


def evaluate(preds, targets, masks):
    out = {}
    for readout, active_key, strength_key in [("transcription", "tx_active", "tx_strength"),
                                                 ("translation", "tl_active", "tl_strength")]:
        active_mask = masks[active_key]
        y_active = targets[active_key][active_mask].astype(int)
        prob = 1 / (1 + np.exp(-preds[f"{active_key}_logit"][active_mask]))
        pred_class = (prob > 0.5).astype(int)
        mcc = matthews_corrcoef(y_active, pred_class) if len(set(y_active)) > 1 else None
        try:
            auc = roc_auc_score(y_active, prob) if len(set(y_active)) > 1 else None
        except Exception:
            auc = None

        reg_mask = masks[strength_key]
        rho = None
        if reg_mask.sum() > 1:
            y_true = targets[strength_key][reg_mask]
            y_pred = preds[strength_key][reg_mask]
            if len(set(np.round(y_true, 6))) > 1:
                rho = spearmanr(y_true, y_pred).correlation
        out[readout] = {
            "n_active_eval": int(active_mask.sum()), "n_strength_eval": int(reg_mask.sum()),
            "mcc": float(mcc) if mcc is not None and mcc == mcc else None,
            "auc": float(auc) if auc is not None and auc == auc else None,
            "spearman_rho": float(rho) if rho is not None and rho == rho else None,
        }
    return out


def main():
    t_start = time.time()
    g_z, gcols, p_z, pcols, imputed = m.load_host_features()
    variants = {"genomic": (g_z, len(gcols)), "physiology": (p_z, len(pcols))}
    print(f"Device={DEVICE}, imputed cells={imputed}")

    all_results = {}
    n_fits_done = 0
    n_fits_total = len(HOSTS) * len(variants) * N_FOLDS
    for held_out in HOSTS:
        train_hosts = [h for h in HOSTS if h != held_out]
        all_results[held_out] = {}
        for variant_name, (feat_df, host_dim) in variants.items():
            all_results[held_out][variant_name] = {"per_fold": {}}
            for test_fold in range(N_FOLDS):
                t0 = time.time()
                seq_tr, hv_tr, targets_tr, masks_tr, _ = m.build_pooled_arrays(
                    train_hosts, feat_df, fold_filter_fn=lambda f, tf=test_fold: f != tf)
                model, epochs_run = m.train_film_model(
                    seq_tr, hv_tr, targets_tr, masks_tr, host_dim=host_dim,
                    epochs=EPOCHS, batch_size=BATCH_SIZE, seed=test_fold, device=DEVICE,
                    use_film=True, patience=PATIENCE, verbose=False)

                seq_te, hv_te, targets_te, masks_te, _ = m.build_pooled_arrays(
                    [held_out], feat_df, fold_filter_fn=lambda f, tf=test_fold: f == tf)
                preds = m.predict_film(model, seq_te, hv_te)
                eval_result = evaluate(preds, targets_te, masks_te)

                elapsed = time.time() - t0
                n_fits_done += 1
                print(f"[{n_fits_done}/{n_fits_total}] held_out={held_out} variant={variant_name} "
                      f"fold={test_fold}: epochs_run={epochs_run}, {elapsed:.1f}s, "
                      f"n_train={len(seq_tr)}, n_test={len(seq_te)}")
                for ro, r in eval_result.items():
                    print(f"    {ro}: mcc={r['mcc']}, auc={r['auc']}, rho={r['spearman_rho']}")

                all_results[held_out][variant_name]["per_fold"][test_fold] = {
                    "eval": eval_result, "epochs_run": epochs_run, "wall_clock_sec": elapsed,
                    "n_train": int(len(seq_tr)), "n_test": int(len(seq_te)),
                }
                model_path = MODELS_DIR / f"loho_{held_out}_{variant_name}_fold{test_fold}.pt"
                torch.save(model.state_dict(), model_path)

                with open(OUT / "gate4_loho_results.json", "w") as f:
                    json.dump(all_results, f, indent=2, default=str)

    print(f"\nTotal wall time: {time.time()-t_start:.1f}s")
    with open(OUT / "gate4_loho_results.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"Wrote {OUT / 'gate4_loho_results.json'}")


if __name__ == "__main__":
    main()
