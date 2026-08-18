"""
GATE 8.5 - Task 1B: predict the shift, not the level.

Reframe: instead of predicting activity(target host) from sequence alone,
give the model the REFERENCE host's measured value and ask it to predict
the DIFFERENCE target - reference. Rationale (task's own): the absolute
level may be dominated by sequence properties the model already captures,
while the cross-host SHIFT is the quantity that might encode host-specific
biology, and handing the model the reference value for free may isolate a
learnable signal currently buried in the level-prediction framing.

DATA: for each ordered host pair (ref, target) among the 3 primary hosts (6
ordered pairs) and each readout, restricted to sequences with a usable
STRENGTH value in BOTH hosts -- i.e. exactly the "co-active" mask already
established in scripts/57 and 80 (strength usability already implies
usable-for-regression & active; there is no broader "pooled" set for a
continuous shift target the way there was for the binary active/inactive
question in Task 1C). Target = target_strength - ref_strength, on the SAME
transformed scale scripts/40's build_target_arrays already uses for each
readout (log1p(tx_norm) for transcription, protein_log10 for translation)
-- not a separately-invented scale.

SPLIT: the standard 5-fold sequence-cluster assignment (same folds as every
other model in this project), NOT leave-one-host-out -- Task 1B is not
asking the model to generalize to an unseen host (the reference host's true
value is handed to it at both train and test time); it is asking whether
the shift is predictable from sequence + reference value, evaluated on
held-out SEQUENCES.

TWO VARIANTS: "no_condition" (input = pooled seq features + reference
value only) and "with_condition" (+ target host's genomic z-features,
37-D, same vector as everywhere else in this project) -- per the task,
conditioning is most likely to help HERE specifically, since the target is
now an explicit host-difference rather than an absolute level.

BASELINES (no training; computed once per fold from TRAIN-fold statistics
only, to avoid leakage into the test-fold comparison):
  - zero-shift: predict 0 for every test row ("assume full portability")
  - mean-shift: predict the TRAIN-fold's mean true shift for every test row

METRICS: Spearman rho of the model's predictions (meaningful only for the
model, since the two baselines are constants with zero prediction variance
and therefore an undefined Spearman correlation with anything) AND fraction
of shift-variance explained relative to each baseline,
frac_var_explained = 1 - MSE(model)/MSE(baseline) -- this is the head-to-
head number the task's "can the model beat 'assume it transfers unchanged'"
question actually requires, since Spearman cannot compare a model against a
constant predictor.

DISCLOSED COMPUTE-BOUNDED CHOICE: epochs=25/patience=6 (vs 35/8 used by
Gate 4's from-scratch LOHO fits) and a smaller architecture width are NOT
used here -- the trunk is identical to scripts/78's Trunk. The only
reduction is epoch budget, justified by (a) each pair's training pool is
already known to be far smaller than Gate 4's ~46-58K-row LOHO pools (the
co-active pools this reuses range ~300-9,700 rows per PAPER_FRAMING.md
finding #1 and scripts/80's Part 1 output), so convergence is expected
faster per-epoch already, and (b) the task grid here (6 pairs x 2 readouts
x 2 conditioning variants x 5 folds = 120 fits) is 4x Task 1A's fit count;
running Gate 4's full epoch budget on all 120 would cost proportionally
more wall-clock than the rest of this gate combined. This is a "same
folds, same intervals" task per its own text, not an "identical protocol"
one (contrast Task 1A, which IS run at Gate 4's exact epoch/patience
budget) -- disclosed here and in out/GATE8_5_MEMO.md, not silently applied.
"""
import sys
import time
import json
import numpy as np
import torch
import torch.nn as nn
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
m40 = import_module("40_film_cnn_model")
c78 = import_module("78_conditioning_mechanisms")

DATA = Path(__file__).resolve().parent.parent / "data"
CACHE = DATA / "baseline_cache"
OUT = Path(__file__).resolve().parent.parent / "out"
MODELS_DIR = OUT / "models"
MODELS_DIR.mkdir(parents=True, exist_ok=True)
DEVICE = "mps" if torch.backends.mps.is_available() else "cpu"

HOSTS = ["EC", "BS", "PA"]
ORDERED_PAIRS = [(r, t) for r in HOSTS for t in HOSTS if r != t]
READOUTS = [("transcription", "tx_norm", "tx_usable", "tx_active", "log1p"),
            ("translation", "protein_log10", "tl_usable_for_regression", "tl_active", "identity")]
N_FOLDS = 5
EPOCHS = 25
PATIENCE = 6
BATCH_SIZE = 1024


def strength_transform(raw, kind):
    if kind == "log1p":
        return np.log1p(np.clip(raw, 0, None)).astype(np.float32)
    return raw.astype(np.float32)


def load_pair_data(ref, target, usable_col, active_col, raw_col, transform):
    d_ref = np.load(CACHE / f"{ref}_baseline_data.npz", allow_pickle=True)
    d_tgt = np.load(CACHE / f"{target}_baseline_data.npz", allow_pickle=True)
    both_mask = (d_ref[usable_col] & d_ref[active_col] & d_tgt[usable_col] & d_tgt[active_col])
    idx = np.where(both_mask)[0]
    seq = d_ref["onehot"][idx]  # same onehot regardless of host (sequence identity, not host-specific)
    fold = d_ref["fold"][idx]
    ref_strength = strength_transform(d_ref[raw_col][idx], transform)
    tgt_strength = strength_transform(d_tgt[raw_col][idx], transform)
    shift = tgt_strength - ref_strength
    return seq, fold, ref_strength, shift, int(len(idx))


# ---------------------------------------------------------------------------
# Model: trunk -> concat[pooled, ref_value, (host_vec)] -> hidden -> scalar shift
# ---------------------------------------------------------------------------
class ShiftCNN(nn.Module):
    def __init__(self, host_dim=0):
        super().__init__()
        self.trunk = c78.Trunk()
        in_dim = 64 + 1 + host_dim
        self.fc1 = nn.Linear(in_dim, 32)
        self.out = nn.Linear(32, 1)

    def forward(self, seq, ref_value, host_vec=None):
        pooled = self.trunk(seq)
        parts = [pooled, ref_value.unsqueeze(-1)]
        if host_vec is not None:
            parts.append(host_vec)
        x = torch.cat(parts, dim=1)
        h = torch.relu(self.fc1(x))
        return self.out(h).squeeze(-1)


def train_shift_model(seq, ref_value, host_vec, shift, host_dim, epochs=EPOCHS, lr=1e-3, weight_decay=1e-4,
                       batch_size=BATCH_SIZE, seed=0, device="cpu", patience=PATIENCE):
    torch.manual_seed(seed)
    dev = torch.device(device)
    model = ShiftCNN(host_dim=host_dim).to(dev)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=max(epochs, 1))
    mse = nn.MSELoss()

    X = torch.tensor(seq, dtype=torch.float32).to(dev)
    R = torch.tensor(ref_value, dtype=torch.float32).to(dev)
    Y = torch.tensor(shift, dtype=torch.float32).to(dev)
    H = torch.tensor(host_vec, dtype=torch.float32).to(dev) if host_vec is not None else None

    n = len(X)
    best_loss = float("inf")
    patience_ctr = 0
    epochs_run = 0
    for epoch in range(epochs):
        model.train()
        perm = torch.randperm(n)
        epoch_loss, n_batches = 0.0, 0
        for i in range(0, n, batch_size):
            idx = perm[i:i + batch_size]
            opt.zero_grad()
            hv = H[idx] if H is not None else None
            pred = model(X[idx], R[idx], hv)
            loss = mse(pred, Y[idx])
            loss.backward()
            opt.step()
            epoch_loss += loss.item()
            n_batches += 1
        sched.step()
        epochs_run = epoch + 1
        if n_batches == 0:
            break
        epoch_loss /= n_batches
        if epoch_loss < best_loss - 1e-6:
            best_loss = epoch_loss
            patience_ctr = 0
        else:
            patience_ctr += 1
            if patience_ctr >= patience and epoch > 8:
                break
    model.eval()
    return model, epochs_run


def predict_shift(model, seq, ref_value, host_vec=None, batch_size=1024):
    model.eval()
    dev = next(model.parameters()).device
    n = len(seq)
    out = []
    with torch.no_grad():
        for i in range(0, n, batch_size):
            X = torch.tensor(seq[i:i + batch_size], dtype=torch.float32).to(dev)
            R = torch.tensor(ref_value[i:i + batch_size], dtype=torch.float32).to(dev)
            HV = torch.tensor(host_vec[i:i + batch_size], dtype=torch.float32).to(dev) if host_vec is not None else None
            pred = model(X, R, HV)
            out.append(pred.cpu().numpy())
    return np.concatenate(out)


def main():
    t_all = time.time()
    g_z, gcols, p_z, pcols, imputed = m40.load_host_features()
    assert imputed == []
    host_dim = len(gcols)
    print(f"Device={DEVICE}, host_dim={host_dim}")

    all_results = {}
    n_fits_done, n_fits_total = 0, len(ORDERED_PAIRS) * len(READOUTS) * 2 * N_FOLDS
    for ref, target in ORDERED_PAIRS:
        pair_key = f"{ref}_to_{target}"
        all_results[pair_key] = {}
        for readout_name, raw_col, usable_col, active_col, transform in READOUTS:
            seq, fold, ref_strength, shift, n_total = load_pair_data(ref, target, usable_col, active_col, raw_col, transform)
            print(f"\n=== {pair_key} {readout_name}: n_usable_both={n_total} ===")
            all_results[pair_key][readout_name] = {"n_total": n_total, "per_fold": {}}
            if n_total < 20:
                print("  SKIPPED (n<20, insufficient for a 5-fold split)")
                continue
            target_hv = g_z.loc[target].values.astype(np.float32)

            for test_fold in range(N_FOLDS):
                train_mask = fold != test_fold
                test_mask = fold == test_fold
                n_test = int(test_mask.sum())
                if n_test < 5 or train_mask.sum() < 20:
                    continue

                seq_tr, ref_tr, shift_tr = seq[train_mask], ref_strength[train_mask], shift[train_mask]
                seq_te, ref_te, shift_te = seq[test_mask], ref_strength[test_mask], shift[test_mask]
                hv_tr = np.tile(target_hv, (len(seq_tr), 1))
                hv_te = np.tile(target_hv, (len(seq_te), 1))

                mean_shift_train = float(shift_tr.mean())
                mse_zero = float(np.mean(shift_te ** 2))
                mse_mean = float(np.mean((shift_te - mean_shift_train) ** 2))

                fold_result = {"n_train": int(train_mask.sum()), "n_test": n_test,
                               "mean_shift_train": mean_shift_train, "mse_zero_baseline": mse_zero,
                               "mse_mean_shift_baseline": mse_mean}

                t0 = time.time()
                for variant, use_cond in [("no_condition", False), ("with_condition", True)]:
                    model, epochs_run = train_shift_model(
                        seq_tr, ref_tr, hv_tr if use_cond else None, shift_tr,
                        host_dim=host_dim if use_cond else 0, seed=test_fold, device=DEVICE)
                    pred = predict_shift(model, seq_te, ref_te, hv_te if use_cond else None)
                    mse_model = float(np.mean((pred - shift_te) ** 2))
                    rho = spearmanr(shift_te, pred).correlation if len(set(np.round(shift_te, 6))) > 1 else None
                    fold_result[variant] = {
                        "epochs_run": epochs_run, "mse_model": mse_model,
                        "spearman_rho": float(rho) if rho is not None and rho == rho else None,
                        "frac_var_explained_vs_zero": 1.0 - mse_model / mse_zero if mse_zero > 0 else None,
                        "frac_var_explained_vs_mean": 1.0 - mse_model / mse_mean if mse_mean > 0 else None,
                    }
                    n_fits_done += 1
                elapsed = time.time() - t0
                fold_result["wall_clock_sec"] = elapsed
                all_results[pair_key][readout_name]["per_fold"][test_fold] = fold_result
                print(f"[{n_fits_done}/{n_fits_total}] {pair_key} {readout_name} fold={test_fold} "
                      f"({elapsed:.1f}s): no_cond rho={fold_result['no_condition']['spearman_rho']} "
                      f"fve_zero={fold_result['no_condition']['frac_var_explained_vs_zero']:.3f} | "
                      f"with_cond rho={fold_result['with_condition']['spearman_rho']} "
                      f"fve_zero={fold_result['with_condition']['frac_var_explained_vs_zero']:.3f}")

                with open(OUT / "results" / "gate8_5_shift_prediction.json", "w") as f:
                    json.dump(all_results, f, indent=2, default=str)

    print(f"\n=== TOTAL WALL TIME: {time.time()-t_all:.1f}s ===")


if __name__ == "__main__":
    (OUT / "results").mkdir(parents=True, exist_ok=True)
    main()
