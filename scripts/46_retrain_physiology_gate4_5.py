"""
GATE 4.5 - Task 5: retrain the 15 physiology-variant LOHO fits (3 held-out
hosts x 5 folds; both readouts share one model per fit, per Gate 4's joint
architecture) after Task 1's growth-rate vector correction (B. subtilis
growth_rate_mu_h: imputed mean-of-other-5-hosts 1.852 -> independently
sourced, condition-matched 1.7, see scripts/12).

Genomic-variant fits are NOT retrained -- unaffected by the physiology vector
change, and Task 4's investigation found the PA sign-inversion issue to be
genuine variance in the physiology extrapolation, not a bug requiring wider
retraining (see out/gate4_5_pa_investigation.json).

Identical seeds, folds, architecture, and hyperparameters as the original
Gate 4 Task 2 run (scripts/42) -- only the physiology feature vector's BS row
differs. This is deliberate: keeps genomic-vs-physiology a clean
single-variable comparison per the task's own instruction.
"""
import sys
import time
import json
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
EPOCHS = loho.EPOCHS
PATIENCE = loho.PATIENCE
BATCH_SIZE = loho.BATCH_SIZE


def main():
    t_start = time.time()
    g_z, gcols, p_z, pcols, imputed = m.load_host_features()
    print(f"Device={DEVICE}, imputed cells={imputed} (should be [] now -- confirms Task 1 fix applied)")
    assert imputed == [], "Expected no imputed physiology cells after Gate 4.5 Task 1 fix"

    all_results = {}
    n_fits_done = 0
    for held_out in HOSTS:
        train_hosts = [h for h in HOSTS if h != held_out]
        all_results[held_out] = {"per_fold": {}}
        for test_fold in range(N_FOLDS):
            t0 = time.time()
            seq_tr, hv_tr, targets_tr, masks_tr, _ = m.build_pooled_arrays(
                train_hosts, p_z, fold_filter_fn=lambda f, tf=test_fold: f != tf)
            model, epochs_run = m.train_film_model(
                seq_tr, hv_tr, targets_tr, masks_tr, host_dim=len(pcols),
                epochs=EPOCHS, batch_size=BATCH_SIZE, seed=test_fold, device=DEVICE,
                use_film=True, patience=PATIENCE, verbose=False)

            seq_te, hv_te, targets_te, masks_te, _ = m.build_pooled_arrays(
                [held_out], p_z, fold_filter_fn=lambda f, tf=test_fold: f == tf)
            preds = m.predict_film(model, seq_te, hv_te)
            eval_result = loho.evaluate(preds, targets_te, masks_te)

            elapsed = time.time() - t0
            n_fits_done += 1
            print(f"[{n_fits_done}/15] held_out={held_out} variant=physiology fold={test_fold}: "
                  f"epochs_run={epochs_run}, {elapsed:.1f}s, n_train={len(seq_tr)}, n_test={len(seq_te)}")
            for ro, r in eval_result.items():
                print(f"    {ro}: mcc={r['mcc']}, auc={r['auc']}, rho={r['spearman_rho']}")

            all_results[held_out]["per_fold"][test_fold] = {
                "eval": eval_result, "epochs_run": epochs_run, "wall_clock_sec": elapsed,
                "n_train": int(len(seq_tr)), "n_test": int(len(seq_te)),
            }
            # overwrite the ORIGINAL LOHO checkpoint -- this IS the corrected model now
            model_path = MODELS_DIR / f"loho_{held_out}_physiology_fold{test_fold}.pt"
            torch.save(model.state_dict(), model_path)

    print(f"\nTotal wall time: {time.time()-t_start:.1f}s")

    # merge into the existing gate4_loho_results.json, replacing only the physiology entries
    with open(OUT / "gate4_loho_results.json") as f:
        full_results = json.load(f)
    for held_out in HOSTS:
        full_results[held_out]["physiology"] = all_results[held_out]
    with open(OUT / "gate4_loho_results.json", "w") as f:
        json.dump(full_results, f, indent=2, default=str)
    print(f"Updated {OUT / 'gate4_loho_results.json'} (physiology entries replaced)")

    with open(OUT / "gate4_5_physiology_retrain.json", "w") as f:
        json.dump({"results": all_results, "wall_clock_sec": time.time() - t_start}, f, indent=2, default=str)
    print(f"Wrote {OUT / 'gate4_5_physiology_retrain.json'}")


if __name__ == "__main__":
    main()
