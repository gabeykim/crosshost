"""
GATE 4 - Task 1 validation: does FiLM conditioning actually condition on the
host, or does it silently learn gamma=1/beta=0 and become a no-op? This is
the single most important check before any real training grid is run --
per the task's own framing, a FiLM layer that's structurally present but
functionally inert would invalidate the whole experiment (H-SCIENCE and
H-DIAGNOSTIC both assume the model actually uses the host vector).

Trains a FiLM-enabled and a FiLM-disabled (use_film=False, the exact same
architecture with the host branch removed) model on the SAME pooled EC+PA
data (both readouts, genomic host features), then:
  1. Reports parameter counts for both.
  2. Direct sensitivity test: for a FIXED batch of held-out sequences, swaps
     the host vector between EC's and PA's genomic feature vector and
     measures how much each model's predictions change. This is the cleanest
     possible check -- for the disabled model this MUST be exactly 0 by
     construction (host_vec is never read); for the enabled model a
     meaningfully large difference is direct proof FiLM is conditioning.
  3. Reports gamma/beta statistics (mean, std, min, max per layer) evaluated
     at the 3 primary hosts' genomic feature vectors, after training.
  4. Reports MPS seconds/epoch at this realistic scale.
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

OUT = Path(__file__).resolve().parent.parent / "out"
DEVICE = "mps" if torch.backends.mps.is_available() else "cpu"
EPOCHS = 15
SEED = 0


def main():
    print(f"Device: {DEVICE}")
    g_z, gcols, p_z, pcols, imputed = m.load_host_features()
    print(f"genomic host_dim={len(gcols)}, physiology host_dim={len(pcols)}, imputed_cells={imputed}")

    print("\n=== Building pooled EC+PA training data (all folds, both readouts) ===")
    seq, host_vec, targets, masks, host_of_row = m.build_pooled_arrays(["EC", "PA"], g_z)
    print(f"  N={len(seq)}, host counts: EC={sum(host_of_row=='EC')}, PA={sum(host_of_row=='PA')}")
    for k in m.TARGET_KEYS:
        print(f"  mask[{k}].sum()={masks[k].sum()}")

    print(f"\n=== Training FiLM-ENABLED model, {EPOCHS} epochs ===")
    t0 = time.time()
    model_on, epochs_run_on = m.train_film_model(seq, host_vec, targets, masks, host_dim=len(gcols),
                                                   epochs=EPOCHS, batch_size=1024, seed=SEED,
                                                   device=DEVICE, use_film=True, verbose=True)
    t_on = time.time() - t0
    print(f"  done in {t_on:.1f}s ({epochs_run_on} epochs run, {t_on/epochs_run_on:.2f}s/epoch)")

    print(f"\n=== Training FiLM-DISABLED model (control), {EPOCHS} epochs ===")
    t0 = time.time()
    model_off, epochs_run_off = m.train_film_model(seq, host_vec, targets, masks, host_dim=len(gcols),
                                                     epochs=EPOCHS, batch_size=1024, seed=SEED,
                                                     device=DEVICE, use_film=False, verbose=True)
    t_off = time.time() - t0
    print(f"  done in {t_off:.1f}s ({epochs_run_off} epochs run)")

    print("\n=== Direct sensitivity test: swap host vector, same sequences, measure prediction change ===")
    test_seq = seq[:500]
    ec_vec = np.tile(g_z.loc["EC"].values.astype(np.float32), (500, 1))
    pa_vec = np.tile(g_z.loc["PA"].values.astype(np.float32), (500, 1))

    pred_on_ec = m.predict_film(model_on, test_seq, ec_vec)
    pred_on_pa = m.predict_film(model_on, test_seq, pa_vec)
    pred_off_ec = m.predict_film(model_off, test_seq, ec_vec)
    pred_off_pa = m.predict_film(model_off, test_seq, pa_vec)

    sensitivity = {}
    for k in ["tx_active_logit", "tx_strength", "tl_active_logit", "tl_strength"]:
        diff_on = float(np.mean(np.abs(pred_on_ec[k] - pred_on_pa[k])))
        diff_off = float(np.mean(np.abs(pred_off_ec[k] - pred_off_pa[k])))
        sensitivity[k] = {"film_enabled_mean_abs_diff": diff_on, "film_disabled_mean_abs_diff": diff_off}
        print(f"  {k}: FiLM-enabled |diff|={diff_on:.4f}   FiLM-disabled |diff|={diff_off:.6f} (should be ~0)")

    print("\n=== FiLM gamma/beta statistics after training (evaluated at 3 primary hosts' genomic vectors) ===")
    host_vecs_3 = torch.tensor(g_z.loc[["EC", "BS", "PA"]].values.astype(np.float32))
    gb_stats = m.get_film_gamma_beta_stats(model_on, host_vecs_3)
    for layer, stats in gb_stats.items():
        flag = ""
        if layer.startswith("gamma") and abs(stats["mean"] - 1.0) < 0.02 and stats["std"] < 0.02:
            flag = "  <<< WARNING: gamma looks uniformly ~1, possible inert FiLM"
        if layer.startswith("beta") and abs(stats["mean"]) < 0.02 and stats["std"] < 0.02:
            flag = "  <<< WARNING: beta looks uniformly ~0, possible inert FiLM"
        print(f"  {layer}: mean={stats['mean']:.4f}, std={stats['std']:.4f}, "
              f"min={stats['min']:.4f}, max={stats['max']:.4f}{flag}")

    n_params_on = m.count_params(model_on)
    n_params_off = m.count_params(model_off)
    print(f"\nParam counts: FiLM-enabled={n_params_on}, FiLM-disabled={n_params_off}")

    verdict = "FiLM IS CONDITIONING" if all(sensitivity[k]["film_enabled_mean_abs_diff"] > 10 * max(sensitivity[k]["film_disabled_mean_abs_diff"], 1e-6) for k in sensitivity) else "INCONCLUSIVE OR INERT -- INVESTIGATE"
    print(f"\nVERDICT: {verdict}")

    output = {
        "device": DEVICE, "epochs_requested": EPOCHS,
        "epochs_run_film_on": epochs_run_on, "epochs_run_film_off": epochs_run_off,
        "wall_clock_sec_film_on": t_on, "wall_clock_sec_film_off": t_off,
        "sec_per_epoch_film_on": t_on / epochs_run_on,
        "n_train_examples": int(len(seq)),
        "param_counts": {"film_enabled": n_params_on, "film_disabled": n_params_off},
        "host_vector_swap_sensitivity": sensitivity,
        "film_gamma_beta_stats_after_training": gb_stats,
        "verdict": verdict,
        "imputed_host_feature_cells": imputed,
    }
    with open(OUT / "gate4_film_sanity_check.json", "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nWrote {OUT / 'gate4_film_sanity_check.json'}")

    import os
    os.makedirs(Path(__file__).resolve().parent.parent / "out" / "models", exist_ok=True)
    torch.save(model_on.state_dict(), Path(__file__).resolve().parent.parent / "out" / "models" / "film_sanity_check_model_on.pt")


if __name__ == "__main__":
    main()
