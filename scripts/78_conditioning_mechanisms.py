"""
GATE 8.5 - Task 1A: alternative conditioning mechanisms.

The objection this answers: Gate 4/5's negative result rests on ONE
conditioning mechanism (FiLM at conv layers 2-3), fit from only 2 training
hosts per LOHO fold -- a regime FiLM's own authors note it can underperform
in. This script implements two architecturally distinct alternatives under
the IDENTICAL protocol as Gate 4's scripts/42_loho_training.py: same 5 frozen
folds, same seeds (=fold index), same conv trunk (4 layers, 128/128/64/64,
k=15/9/5/3, BatchNorm, residual wrapping conv4 -- byte-identical to
SequenceOnlyCNN / FiLMSequenceCNN's trunk), same training hyperparameters
(epochs=35, patience=8, batch_size=1024, AdamW lr=1e-3/wd=1e-4, cosine
schedule), same masked joint 4-head loss, same evaluate() function
(scripts/42), same 90% percentile bootstrap CI (scripts/51). Only the
conditioning MECHANISM differs. Both readouts, all 3 primary hosts LOHO,
zero-shot (N=0) evaluation point -- matching Gate 4's own primary evaluation
point, the one H-MAIN and this project's central claim are built on.

HOST-FEATURE VARIANT: genomic only (37-D, z-scored across all 6 hosts, same
as scripts/40's load_host_features()). DISCLOSED SCOPE CHOICE, following the
precedent set by Gate 7's feature-group ablation ("genomic variant, one
representative config"): running both genomic AND physiology variants across
3 mechanisms would double an already compute-heavy task (Gate 4's own 30
zero-shot LOHO fits took a measured 27,135s/7.5hrs total). Genomic is chosen
because it is the stronger of the two host-feature vectors in every prior
gate's comparison (H-SCIENCE's one distinguishable cell favors genomic; nothing
favors physiology) -- if a conditioning mechanism is going to recover signal
anywhere, the richer/better-performing feature vector is the fairer test.

MECHANISM 1 -- Concatenation. Host vector appended to the pooled (64-d)
representation before fc_shared; fc_shared's input dim becomes 64+host_dim.
This is architecturally "FiLM with gamma fixed at 1" -- a pure additive
conditional bias, cheap, single shared head set (same as FiLM/sequence-only).

MECHANISM 2 -- Separate per-host output heads. Shared conv trunk (identical
to Mechanism 1's, MINUS the concatenation -- trunk never sees host
information at all), but each of the 2 TRAINING hosts in a given LOHO fold
gets its own private set of 4 output heads (tx_active/tx_strength/
tl_active/tl_strength), selected per-row by that row's true host during
training (vectorized: all heads score every row, then gathered by host
index -- see PerHostHeadsCNN.forward). At held-out-host eval time there is,
by construction, no head trained for that host. RULE (stated in advance,
per the task's requirement): report predict-space averaging (apply both
trained hosts' heads to the held-out host's pooled features, then average
the two resulting predictions/logits) as the PRIMARY number, because
prediction-space averaging does not assume the two heads' weight spaces are
linearly interpolable (which is not guaranteed, especially for classifier
logits) -- weight-space averaging was considered and rejected for this
reason. Nearest-host selection (apply only the single training host closest
to the held-out host in genomic-feature Euclidean distance) is ALSO computed
and reported as a secondary column, at zero marginal training cost (same
trained model, different combination rule at eval time).

MECHANISM 3 -- hypernetwork: NOT RUN. The task marks it "optional if cheap."
Given the two mandatory mechanisms alone already require 2 x 3 x 5 = 30
from-scratch LOHO fits at Gate 4's own measured ~904s/fit average (~7.5hrs
total, this script), a third mechanism trained to the same standard would add
a comparable amount of additional wall-clock for a question the two mandatory
mechanisms are already positioned to answer conclusively (does ANY
alternative conditioning mechanism beat sequence-only). Skipped, disclosed,
not silently dropped -- see out/GATE8_5_MEMO.md.
"""
import sys
import time
import json
import numpy as np
import torch
import torch.nn as nn
from pathlib import Path
from scipy.stats import spearmanr
from sklearn.metrics import matthews_corrcoef, roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
m40 = import_module("40_film_cnn_model")
loho = import_module("42_loho_training")

OUT = Path(__file__).resolve().parent.parent / "out"
MODELS_DIR = OUT / "models"
MODELS_DIR.mkdir(parents=True, exist_ok=True)
DEVICE = "mps" if torch.backends.mps.is_available() else "cpu"

HOSTS = ["EC", "BS", "PA"]
N_FOLDS = 5
EPOCHS = 35
PATIENCE = 8
BATCH_SIZE = 1024

evaluate = loho.evaluate


# ---------------------------------------------------------------------------
# Shared trunk (identical to SequenceOnlyCNN / FiLMSequenceCNN's conv stack)
# ---------------------------------------------------------------------------
class Trunk(nn.Module):
    def __init__(self):
        super().__init__()
        self.conv1 = nn.Conv1d(4, 128, kernel_size=15, padding=7)
        self.bn1 = nn.BatchNorm1d(128)
        self.conv2 = nn.Conv1d(128, 128, kernel_size=9, padding=4)
        self.bn2 = nn.BatchNorm1d(128)
        self.conv3 = nn.Conv1d(128, 64, kernel_size=5, padding=2)
        self.bn3 = nn.BatchNorm1d(64)
        self.conv4 = nn.Conv1d(64, 64, kernel_size=3, padding=1)
        self.bn4 = nn.BatchNorm1d(64)
        self.dropout = nn.Dropout(0.2)
        self.pool = nn.AdaptiveAvgPool1d(1)

    def forward(self, seq):
        h = torch.relu(self.bn1(self.conv1(seq)))
        h = torch.relu(self.bn2(self.conv2(h)))
        h = self.dropout(h)
        h3 = torch.relu(self.bn3(self.conv3(h)))
        h4 = self.bn4(self.conv4(h3))
        h4 = torch.relu(h4 + h3)
        return self.pool(h4).squeeze(-1)  # (batch, 64)


def count_params(model):
    return sum(p.numel() for p in model.parameters() if p.requires_grad)


# ---------------------------------------------------------------------------
# Mechanism 1: concatenation
# ---------------------------------------------------------------------------
class ConcatCNN(nn.Module):
    def __init__(self, host_dim):
        super().__init__()
        self.trunk = Trunk()
        self.fc_shared = nn.Linear(64 + host_dim, 64)
        self.head_tx_active = nn.Linear(64, 1)
        self.head_tx_strength = nn.Linear(64, 1)
        self.head_tl_active = nn.Linear(64, 1)
        self.head_tl_strength = nn.Linear(64, 1)

    def forward(self, seq, host_vec):
        pooled = self.trunk(seq)
        cat = torch.cat([pooled, host_vec], dim=1)
        shared = torch.relu(self.fc_shared(cat))
        return {
            "tx_active_logit": self.head_tx_active(shared).squeeze(-1),
            "tx_strength": self.head_tx_strength(shared).squeeze(-1),
            "tl_active_logit": self.head_tl_active(shared).squeeze(-1),
            "tl_strength": self.head_tl_strength(shared).squeeze(-1),
        }


def train_concat_model(seq, host_vec, targets, masks, host_dim, epochs=EPOCHS, lr=1e-3, weight_decay=1e-4,
                        batch_size=BATCH_SIZE, seed=0, device="cpu", patience=PATIENCE):
    torch.manual_seed(seed)
    dev = torch.device(device)
    model = ConcatCNN(host_dim).to(dev)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=max(epochs, 1))
    bce = nn.BCEWithLogitsLoss(reduction="none")
    mse = nn.MSELoss(reduction="none")

    X = torch.tensor(seq, dtype=torch.float32).to(dev)
    H = torch.tensor(host_vec, dtype=torch.float32).to(dev)
    T = {k: torch.tensor(v, dtype=torch.float32).to(dev) for k, v in targets.items()}
    M = {k: torch.tensor(v, dtype=torch.bool).to(dev) for k, v in masks.items()}

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
            preds = model(X[idx], H[idx])
            targets_b = {k: v[idx] for k, v in T.items()}
            masks_b = {k: v[idx] for k, v in M.items()}
            loss, _ = m40.compute_loss(preds, targets_b, masks_b, bce, mse)
            loss.backward()
            opt.step()
            epoch_loss += loss.item()
            n_batches += 1
        sched.step()
        epochs_run = epoch + 1
        if n_batches == 0:
            break
        epoch_loss /= n_batches
        if epoch_loss < best_loss - 1e-5:
            best_loss = epoch_loss
            patience_ctr = 0
        else:
            patience_ctr += 1
            if patience_ctr >= patience and epoch > 10:
                break
    model.eval()
    return model, epochs_run


def predict_concat(model, seq, host_vec, batch_size=1024):
    """Chunked per SR7 (out/CHARTER_AMENDMENTS.md) -- never call a raw
    unbatched forward pass over a pooled multi-fold/host array."""
    model.eval()
    dev = next(model.parameters()).device
    n = len(seq)
    out = None
    with torch.no_grad():
        for i in range(0, n, batch_size):
            X = torch.tensor(seq[i:i + batch_size], dtype=torch.float32).to(dev)
            H = torch.tensor(host_vec[i:i + batch_size], dtype=torch.float32).to(dev)
            preds = model(X, H)
            if out is None:
                out = {k: [] for k in preds}
            for k, v in preds.items():
                out[k].append(v.cpu().numpy())
    return {k: np.concatenate(v) for k, v in out.items()}


# ---------------------------------------------------------------------------
# Mechanism 2: separate per-host output heads
# ---------------------------------------------------------------------------
class PerHostHeadsCNN(nn.Module):
    def __init__(self, train_hosts):
        super().__init__()
        self.trunk = Trunk()
        self.train_hosts = list(train_hosts)
        self.fc_shared = nn.ModuleDict({h: nn.Linear(64, 64) for h in self.train_hosts})
        self.heads = nn.ModuleDict({
            h: nn.ModuleDict({
                "tx_active": nn.Linear(64, 1), "tx_strength": nn.Linear(64, 1),
                "tl_active": nn.Linear(64, 1), "tl_strength": nn.Linear(64, 1),
            }) for h in self.train_hosts
        })

    def forward_all_hosts(self, seq):
        """Returns dict {host: {output_key: (batch,) tensor}} -- every
        training host's private head-set applied to every row. Cheap (heads
        are small linear layers); used both for training (select by true
        host) and eval (combine across hosts)."""
        pooled = self.trunk(seq)
        out = {}
        for h in self.train_hosts:
            shared = torch.relu(self.fc_shared[h](pooled))
            out[h] = {
                "tx_active_logit": self.heads[h]["tx_active"](shared).squeeze(-1),
                "tx_strength": self.heads[h]["tx_strength"](shared).squeeze(-1),
                "tl_active_logit": self.heads[h]["tl_active"](shared).squeeze(-1),
                "tl_strength": self.heads[h]["tl_strength"](shared).squeeze(-1),
            }
        return out

    def forward_selected(self, seq, host_of_row_idx):
        """Training path: host_of_row_idx (batch,) long tensor, values index
        into self.train_hosts. Selects each row's OWN host's head output."""
        all_out = self.forward_all_hosts(seq)
        keys = ["tx_active_logit", "tx_strength", "tl_active_logit", "tl_strength"]
        stacked = {k: torch.stack([all_out[h][k] for h in self.train_hosts], dim=0) for k in keys}  # (n_hosts, batch)
        selected = {k: stacked[k].gather(0, host_of_row_idx.unsqueeze(0)).squeeze(0) for k in keys}
        return selected


def train_perhost_heads_model(seq, host_of_row, targets, masks, train_hosts, epochs=EPOCHS, lr=1e-3,
                               weight_decay=1e-4, batch_size=BATCH_SIZE, seed=0, device="cpu", patience=PATIENCE):
    torch.manual_seed(seed)
    dev = torch.device(device)
    model = PerHostHeadsCNN(train_hosts).to(dev)
    host_idx_map = {h: i for i, h in enumerate(model.train_hosts)}
    host_idx = np.array([host_idx_map[h] for h in host_of_row], dtype=np.int64)

    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=max(epochs, 1))
    bce = nn.BCEWithLogitsLoss(reduction="none")
    mse = nn.MSELoss(reduction="none")

    X = torch.tensor(seq, dtype=torch.float32).to(dev)
    HI = torch.tensor(host_idx, dtype=torch.long).to(dev)
    T = {k: torch.tensor(v, dtype=torch.float32).to(dev) for k, v in targets.items()}
    M = {k: torch.tensor(v, dtype=torch.bool).to(dev) for k, v in masks.items()}

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
            selected = model.forward_selected(X[idx], HI[idx])
            targets_b = {k: v[idx] for k, v in T.items()}
            masks_b = {k: v[idx] for k, v in M.items()}
            loss, _ = m40.compute_loss(selected, targets_b, masks_b, bce, mse)
            loss.backward()
            opt.step()
            epoch_loss += loss.item()
            n_batches += 1
        sched.step()
        epochs_run = epoch + 1
        if n_batches == 0:
            break
        epoch_loss /= n_batches
        if epoch_loss < best_loss - 1e-5:
            best_loss = epoch_loss
            patience_ctr = 0
        else:
            patience_ctr += 1
            if patience_ctr >= patience and epoch > 10:
                break
    model.eval()
    return model, epochs_run


def predict_perhost_heads(model, seq, batch_size=1024, nearest_host=None):
    """Returns (preds_avg, preds_nearest). preds_avg: mean prediction across
    all trained hosts' heads (logits averaged pre-sigmoid, matching standard
    ensembling convention). preds_nearest: single nearest training host's
    head only (nearest_host: str, must be in model.train_hosts)."""
    model.eval()
    dev = next(model.parameters()).device
    n = len(seq)
    keys = ["tx_active_logit", "tx_strength", "tl_active_logit", "tl_strength"]
    avg_parts = {k: [] for k in keys}
    nearest_parts = {k: [] for k in keys}
    with torch.no_grad():
        for i in range(0, n, batch_size):
            X = torch.tensor(seq[i:i + batch_size], dtype=torch.float32).to(dev)
            all_out = model.forward_all_hosts(X)
            for k in keys:
                stacked = torch.stack([all_out[h][k] for h in model.train_hosts], dim=0)
                avg_parts[k].append(stacked.mean(dim=0).cpu().numpy())
                if nearest_host is not None:
                    nearest_parts[k].append(all_out[nearest_host][k].cpu().numpy())
    preds_avg = {k: np.concatenate(v) for k, v in avg_parts.items()}
    preds_nearest = {k: np.concatenate(v) for k, v in nearest_parts.items()} if nearest_host is not None else None
    return preds_avg, preds_nearest


def nearest_training_host(held_out, train_hosts, g_z):
    dists = {h: float(np.linalg.norm(g_z.loc[held_out].values - g_z.loc[h].values)) for h in train_hosts}
    return min(dists, key=dists.get), dists


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------
def main():
    t_all = time.time()
    g_z, gcols, p_z, pcols, imputed = m40.load_host_features()
    assert imputed == []
    host_dim = len(gcols)
    print(f"Device={DEVICE}, host_dim={host_dim}")

    concat_results = {}
    perhost_results = {}
    n_fits_done = 0
    n_fits_total = len(HOSTS) * N_FOLDS * 2

    for held_out in HOSTS:
        train_hosts = [h for h in HOSTS if h != held_out]
        nearest_h, dists = nearest_training_host(held_out, train_hosts, g_z)
        print(f"held_out={held_out}: train_hosts={train_hosts}, nearest={nearest_h}, dists={dists}")
        concat_results[held_out] = {"per_fold": {}}
        perhost_results[held_out] = {"train_hosts": train_hosts, "nearest_host": nearest_h,
                                      "host_distances": dists, "per_fold": {}}

        for test_fold in range(N_FOLDS):
            seq_tr, hv_tr, targets_tr, masks_tr, host_row_tr = m40.build_pooled_arrays(
                train_hosts, g_z, fold_filter_fn=lambda f, tf=test_fold: f != tf)
            seq_te, hv_te, targets_te, masks_te, _ = m40.build_pooled_arrays(
                [held_out], g_z, fold_filter_fn=lambda f, tf=test_fold: f == tf)

            # --- Mechanism 1: concatenation ---
            t0 = time.time()
            model_c, epochs_c = train_concat_model(seq_tr, hv_tr, targets_tr, masks_tr, host_dim=host_dim,
                                                     seed=test_fold, device=DEVICE)
            preds_c = predict_concat(model_c, seq_te, hv_te)
            eval_c = evaluate(preds_c, targets_te, masks_te)
            elapsed_c = time.time() - t0
            concat_results[held_out]["per_fold"][test_fold] = {
                "eval": eval_c, "epochs_run": epochs_c, "wall_clock_sec": elapsed_c,
                "n_train": int(len(seq_tr)), "n_test": int(len(seq_te)),
            }
            torch.save(model_c.state_dict(), MODELS_DIR / f"concat_{held_out}_fold{test_fold}.pt")
            n_fits_done += 1
            print(f"[{n_fits_done}/{n_fits_total}] CONCAT held_out={held_out} fold={test_fold}: "
                  f"{elapsed_c:.1f}s, epochs={epochs_c}")
            for ro, r in eval_c.items():
                print(f"    {ro}: mcc={r['mcc']}, auc={r['auc']}, rho={r['spearman_rho']}")

            # --- Mechanism 2: per-host heads ---
            t0 = time.time()
            model_p, epochs_p = train_perhost_heads_model(seq_tr, host_row_tr, targets_tr, masks_tr,
                                                            train_hosts=train_hosts, seed=test_fold, device=DEVICE)
            preds_avg, preds_nearest = predict_perhost_heads(model_p, seq_te, nearest_host=nearest_h)
            eval_avg = evaluate(preds_avg, targets_te, masks_te)
            eval_nearest = evaluate(preds_nearest, targets_te, masks_te)
            elapsed_p = time.time() - t0
            perhost_results[held_out]["per_fold"][test_fold] = {
                "eval_avg": eval_avg, "eval_nearest": eval_nearest, "epochs_run": epochs_p,
                "wall_clock_sec": elapsed_p, "n_train": int(len(seq_tr)), "n_test": int(len(seq_te)),
            }
            torch.save(model_p.state_dict(), MODELS_DIR / f"perhost_{held_out}_fold{test_fold}.pt")
            n_fits_done += 1
            print(f"[{n_fits_done}/{n_fits_total}] PERHOST held_out={held_out} fold={test_fold}: "
                  f"{elapsed_p:.1f}s, epochs={epochs_p}")
            for ro, r in eval_avg.items():
                print(f"    {ro} (avg): mcc={r['mcc']}, auc={r['auc']}, rho={r['spearman_rho']}")

            with open(OUT / "gate8_5_concat_loho_results.json", "w") as f:
                json.dump(concat_results, f, indent=2, default=str)
            with open(OUT / "gate8_5_perhost_heads_loho_results.json", "w") as f:
                json.dump(perhost_results, f, indent=2, default=str)

    print(f"\n=== TOTAL WALL TIME: {time.time()-t_all:.1f}s ===")


if __name__ == "__main__":
    main()
