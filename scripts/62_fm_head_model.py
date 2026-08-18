"""
GATE 6 - Task 1/2/3: the "lightweight supervised head" trained on top of a
FROZEN foundation-model embedding (charter's preferred protocol: "Extract
[FM] embeddings for the 165bp sequences, train a small supervised head
(ridge, or a shallow MLP) on the same folds"). A shallow MLP is used (not
ridge) so the SAME two-stage active/strength loss and evaluate() function
from scripts/40 and scripts/42 apply unchanged -- identical metrics to every
other system in this project (mcc, auc, spearman_rho via
`42_loho_training.evaluate`), not a bespoke metric requiring re-justification.

Structurally this is the same role scripts/58's SequenceOnlyCNN plays for the
"sequence-only" ablation, with one difference: there the trunk (conv1-4) is
TRAINED per LOHO fit; here the trunk (the FM) is fixed and pre-computed
(scripts/61) -- only this small head is ever trained. No host_vec input
anywhere, matching sequence-only/no-host-conditioning by construction (an FM
embedding of the sequence alone carries no host information to condition on).
"""
import torch
import torch.nn as nn
import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
m40 = import_module("40_film_cnn_model")

DATA = Path(__file__).resolve().parent.parent / "data"
EMB_DIR = DATA / "fm_embeddings"


class FMHeadMLP(nn.Module):
    def __init__(self, emb_dim, hidden=64):
        super().__init__()
        self.fc1 = nn.Linear(emb_dim, hidden)
        self.head_tx_active = nn.Linear(hidden, 1)
        self.head_tx_strength = nn.Linear(hidden, 1)
        self.head_tl_active = nn.Linear(hidden, 1)
        self.head_tl_strength = nn.Linear(hidden, 1)

    def forward(self, emb):
        h = torch.relu(self.fc1(emb))
        return {
            "tx_active_logit": self.head_tx_active(h).squeeze(-1),
            "tx_strength": self.head_tx_strength(h).squeeze(-1),
            "tl_active_logit": self.head_tl_active(h).squeeze(-1),
            "tl_strength": self.head_tl_strength(h).squeeze(-1),
        }


def count_params(model):
    return sum(p.numel() for p in model.parameters() if p.requires_grad)


def train_fm_head(emb, targets, masks, emb_dim, epochs=60, lr=2e-3, weight_decay=1e-4,
                   batch_size=1024, seed=0, device="cpu", patience=10, init_model=None):
    torch.manual_seed(seed)
    dev = torch.device(device)
    model = init_model.to(dev) if init_model is not None else FMHeadMLP(emb_dim).to(dev)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=max(epochs, 1))
    bce = nn.BCEWithLogitsLoss(reduction="none")
    mse = nn.MSELoss(reduction="none")

    X = torch.tensor(emb, dtype=torch.float32).to(dev)
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
            preds = model(X[idx])
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
        if epoch_loss < best_loss - 1e-6:
            best_loss = epoch_loss
            patience_ctr = 0
        else:
            patience_ctr += 1
            if patience_ctr >= patience and epoch > 15:
                break
    model.eval()
    return model, epochs_run


def predict_fm_head(model, emb):
    model.eval()
    dev = next(model.parameters()).device
    with torch.no_grad():
        X = torch.tensor(emb, dtype=torch.float32).to(dev)
        preds = model(X)
        return {k: v.cpu().numpy() for k, v in preds.items()}


def load_library_embeddings(model_tag):
    """Returns (oligo_ids, embeddings (N,D), fold) row-aligned to
    data/baseline_cache/{host}_baseline_data.npz (verified byte-identical
    oligo_id order across hosts in scripts/61)."""
    d = dict(np.load(EMB_DIR / f"{model_tag}_library.npz", allow_pickle=True))
    return d["oligo_ids"], d["embeddings"], d["fold"]


def load_rs241_embeddings(model_tag):
    d = dict(np.load(EMB_DIR / f"{model_tag}_rs241.npz", allow_pickle=True))
    return d["ids"], d["embeddings"]


if __name__ == "__main__":
    model = FMHeadMLP(emb_dim=768)
    print(f"FMHeadMLP(emb_dim=768) params: {count_params(model)}")
