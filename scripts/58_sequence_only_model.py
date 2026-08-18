"""
GATE 5.5 - Task 1: the sequence-only ablation architecture. NOT Gate 4's
FiLM-disabled control -- that model still carried a `use_film` branch and
accepted (and ignored) a host_vec argument; its identity-transform path
(gamma=1, beta=0) was numerically a no-op but the code branch existed. This
class has NO host_vec parameter anywhere, no FiLMGenerator import, no
conditioning pathway of any kind -- the architecture literally cannot be
handed host information. Same trunk otherwise: 4 conv layers (128/128/64/64,
k=15/9/5/3), BatchNorm, residual wrapping conv4, two jointly-trained output
heads (transcription, translation) each two-stage (active + strength).

213,956 params -- IDENTICAL to Gate 4's FiLM-disabled control's parameter
count (confirms that control had zero conditioning parameters even though
its code path retained the (unused) argument).
"""
import torch
import torch.nn as nn
import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
m40 = import_module("40_film_cnn_model")


class SequenceOnlyCNN(nn.Module):
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
        self.fc_shared = nn.Linear(64, 64)
        self.head_tx_active = nn.Linear(64, 1)
        self.head_tx_strength = nn.Linear(64, 1)
        self.head_tl_active = nn.Linear(64, 1)
        self.head_tl_strength = nn.Linear(64, 1)

    def forward(self, seq):
        h = torch.relu(self.bn1(self.conv1(seq)))
        h = torch.relu(self.bn2(self.conv2(h)))
        h = self.dropout(h)
        h3 = torch.relu(self.bn3(self.conv3(h)))
        h4 = self.bn4(self.conv4(h3))
        h4 = torch.relu(h4 + h3)
        pooled = self.pool(h4).squeeze(-1)
        shared = torch.relu(self.fc_shared(pooled))
        return {
            "tx_active_logit": self.head_tx_active(shared).squeeze(-1),
            "tx_strength": self.head_tx_strength(shared).squeeze(-1),
            "tl_active_logit": self.head_tl_active(shared).squeeze(-1),
            "tl_strength": self.head_tl_strength(shared).squeeze(-1),
        }


def count_params(model):
    return sum(p.numel() for p in model.parameters() if p.requires_grad)


FREEZE_MODES = (None, "head_only", "top_conv")


def set_freeze_mode(model, mode):
    """Same semantics as scripts/40's version, minus the FiLM generator
    (which does not exist here). head_only: freeze conv1-4/bn1-4, only
    fc_shared+heads train. top_conv: also unfreeze conv4/bn4."""
    assert mode in FREEZE_MODES
    always_trainable = [model.fc_shared, model.head_tx_active, model.head_tx_strength,
                         model.head_tl_active, model.head_tl_strength]
    frozen_modules = [model.conv1, model.bn1, model.conv2, model.bn2, model.conv3, model.bn3,
                       model.conv4, model.bn4]
    if mode is None:
        for mod in frozen_modules:
            for p in mod.parameters():
                p.requires_grad = True
            mod.train()
        return
    for mod in frozen_modules:
        for p in mod.parameters():
            p.requires_grad = False
        mod.eval()
    if mode == "top_conv":
        for p in model.conv4.parameters():
            p.requires_grad = True
        for p in model.bn4.parameters():
            p.requires_grad = True
        model.conv4.train()
        model.bn4.train()
    for mod in always_trainable:
        for p in mod.parameters():
            p.requires_grad = True


def compute_loss(preds, targets, masks, bce_none, mse_none):
    """Identical vectorized-masked-loss logic to scripts/40 (same bugs already
    fixed there: no boolean indexing, targets pre-sanitized for NaN/Inf by
    build_target_arrays)."""
    total = 0.0
    for name, kind in m40.LOSS_HEADS:
        mask = masks[name].float()
        denom = mask.sum().clamp(min=1.0)
        p = preds[f"{name}_logit" if kind == "bce" else name]
        y = targets[name]
        elementwise = bce_none(p, y) if kind == "bce" else mse_none(p, y)
        loss = (elementwise * mask).sum() / denom
        total = total + loss
    return total


def train_seqonly_model(seq, targets, masks, epochs=40, lr=1e-3, weight_decay=1e-4,
                         batch_size=1024, seed=0, device="cpu", patience=8,
                         freeze_mode=None, init_model=None):
    torch.manual_seed(seed)
    dev = torch.device(device)
    model = init_model.to(dev) if init_model is not None else SequenceOnlyCNN().to(dev)
    if freeze_mode is not None:
        set_freeze_mode(model, freeze_mode)
    trainable = [p for p in model.parameters() if p.requires_grad]
    opt = torch.optim.AdamW(trainable, lr=lr, weight_decay=weight_decay)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=max(epochs, 1))
    bce = nn.BCEWithLogitsLoss(reduction="none")
    mse = nn.MSELoss(reduction="none")

    X = torch.tensor(seq, dtype=torch.float32).to(dev)
    T = {k: torch.tensor(v, dtype=torch.float32).to(dev) for k, v in targets.items()}
    M = {k: torch.tensor(v, dtype=torch.bool).to(dev) for k, v in masks.items()}

    n = len(X)
    best_loss = float("inf")
    patience_ctr = 0
    epochs_run = 0
    for epoch in range(epochs):
        model.train()
        if freeze_mode is not None:
            set_freeze_mode(model, freeze_mode)
        perm = torch.randperm(n)
        epoch_loss, n_batches = 0.0, 0
        for i in range(0, n, batch_size):
            idx = perm[i:i + batch_size]
            opt.zero_grad()
            preds = model(X[idx])
            targets_b = {k: v[idx] for k, v in T.items()}
            masks_b = {k: v[idx] for k, v in M.items()}
            loss = compute_loss(preds, targets_b, masks_b, bce, mse)
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


def predict_seqonly(model, seq, batch_size=1024):
    """GATE 7 FIX: a single unbatched forward pass over the whole input array
    hits a severe non-linear performance cliff (measured: 5,810 rows in
    28.5s, 23,232 rows -- 4x the data -- in 8,454s, ~300x the cost, on CPU;
    Gate 4 documented the same shape of cliff on MPS above batch~2048). This
    is a missing loop, not a slow machine. Chunking at batch_size=1024
    (this project's already-established safe size) fixes it; see
    scripts/audit_batching_fix.py for the numerical-identity verification
    and before/after benchmark this fix was validated against."""
    model.eval()
    dev = next(model.parameters()).device
    n = len(seq)
    out = None
    with torch.no_grad():
        for i in range(0, n, batch_size):
            X = torch.tensor(seq[i:i + batch_size], dtype=torch.float32).to(dev)
            preds = model(X)
            if out is None:
                out = {k: [] for k in preds}
            for k, v in preds.items():
                out[k].append(v.cpu().numpy())
    return {k: np.concatenate(v) for k, v in out.items()}


if __name__ == "__main__":
    model = SequenceOnlyCNN()
    print(f"SequenceOnlyCNN params: {count_params(model)}")
    from importlib import import_module
    m40check = import_module("40_film_cnn_model")
    control = m40check.FiLMSequenceCNN(host_dim=37, use_film=False)
    print(f"Gate 4 FiLM-disabled control params: {m40check.count_params(control)}")
    assert count_params(model) == m40check.count_params(control), "param count mismatch -- investigate"
    print("Confirmed: identical parameter count to Gate 4's FiLM-disabled control.")
