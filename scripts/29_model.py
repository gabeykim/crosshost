"""
GATE 3 - Shared sequence-only CNN architecture for baselines B2/B3.

Architecture family matches the charter's Gate 4 spec (4 conv layers, kernel
widths 15/9/5/3, channels 128/128/64/64) MINUS host conditioning -- this is
what isolates the host-conditioning contribution when Gate 4's model is
built and compared against this baseline later.

SIMPLIFICATION, disclosed: the charter's Gate 4 design is one model with two
output heads (transcription + translation) trained jointly. For these
baselines, transcription and translation are trained as SEPARATE model
instances of the same architecture (one single-readout output: active-logit
+ strength-regression-value), because B2 requires many independent per-host,
per-N, per-draw fits and a shared multi-task model would entangle the two
readouts' calibration curves in a way that complicates the N-sweep. Same
architecture family, different training protocol.

No BatchNorm: N can be as low as 10 in the calibration sweep, where batch
statistics would be unusable. Dropout + weight decay used for regularization
instead.
"""
import torch
import torch.nn as nn
import numpy as np


class SequenceCNN(nn.Module):
    """~220K parameters. Input: (batch, 4, 165) one-hot. Outputs: (active_logit, strength_pred)."""

    def __init__(self, host_embed_dim=0):
        super().__init__()
        self.host_embed_dim = host_embed_dim
        self.conv1 = nn.Conv1d(4, 128, kernel_size=15, padding=7)
        self.conv2 = nn.Conv1d(128, 128, kernel_size=9, padding=4)
        self.conv3 = nn.Conv1d(128, 64, kernel_size=5, padding=2)
        self.conv4 = nn.Conv1d(64, 64, kernel_size=3, padding=1)
        self.dropout = nn.Dropout(0.2)
        self.pool = nn.AdaptiveAvgPool1d(1)
        fc_in = 64 + host_embed_dim
        self.fc_shared = nn.Linear(fc_in, 64)
        self.fc_active = nn.Linear(64, 1)
        self.fc_strength = nn.Linear(64, 1)

    def forward(self, x, host_embed=None):
        h = torch.relu(self.conv1(x))
        h = torch.relu(self.conv2(h))
        h = self.dropout(h)
        h = torch.relu(self.conv3(h))
        h = torch.relu(self.conv4(h))
        h = self.pool(h).squeeze(-1)  # (batch, 64)
        if self.host_embed_dim > 0 and host_embed is not None:
            h = torch.cat([h, host_embed], dim=1)
        h = torch.relu(self.fc_shared(h))
        active_logit = self.fc_active(h).squeeze(-1)
        strength = self.fc_strength(h).squeeze(-1)
        return active_logit, strength


def count_params(model):
    return sum(p.numel() for p in model.parameters() if p.requires_grad)


def train_model(X_train, y_active_train, y_strength_train, active_mask_train,
                 X_val=None, epochs=60, lr=1e-3, weight_decay=1e-4, batch_size=64,
                 seed=0, verbose=False, host_embed_train=None, host_embed_dim=0,
                 patience=10, device="cpu"):
    """Trains a SequenceCNN. y_strength_train / active_mask_train select which
    rows have a usable strength-regression target (only actives with a
    non-floor measurement). Returns the trained model.
    Early stopping on training loss plateau (no held-out val split inside
    this function -- N can be as small as 10, too small to further split).
    device: 'cpu' (default, unchanged behavior) or 'mps' (Apple Silicon GPU
    -- added Gate 3.5, see scripts/38_cnn_compute_benchmark.py; measured
    ~4-20x speedup over CPU in this environment depending on batch size)."""
    torch.manual_seed(seed)
    dev = torch.device(device)
    model = SequenceCNN(host_embed_dim=host_embed_dim).to(dev)
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
    bce = nn.BCEWithLogitsLoss()
    mse = nn.MSELoss()

    X = torch.tensor(X_train, dtype=torch.float32).to(dev)
    y_active = torch.tensor(y_active_train, dtype=torch.float32).to(dev)
    y_strength = torch.tensor(y_strength_train, dtype=torch.float32).to(dev)
    active_mask = torch.tensor(active_mask_train, dtype=torch.bool).to(dev)
    he = torch.tensor(host_embed_train, dtype=torch.float32).to(dev) if host_embed_train is not None else None

    n = len(X)
    best_loss = float("inf")
    patience_ctr = 0
    for epoch in range(epochs):
        model.train()
        perm = torch.randperm(n)
        epoch_loss = 0.0
        for i in range(0, n, batch_size):
            idx = perm[i:i + batch_size]
            opt.zero_grad()
            xb = X[idx]
            heb = he[idx] if he is not None else None
            active_logit, strength = model(xb, heb)
            loss = bce(active_logit, y_active[idx])
            am = active_mask[idx]
            if am.sum() > 0:
                loss = loss + mse(strength[am], y_strength[idx][am])
            loss.backward()
            opt.step()
            epoch_loss += loss.item() * len(idx)
        epoch_loss /= n
        if verbose and epoch % 10 == 0:
            print(f"    epoch {epoch}: loss={epoch_loss:.4f}")
        if epoch_loss < best_loss - 1e-4:
            best_loss = epoch_loss
            patience_ctr = 0
        else:
            patience_ctr += 1
            if patience_ctr >= patience and epoch > 15:
                break
    return model


def predict(model, X, host_embed=None):
    model.eval()
    dev = next(model.parameters()).device
    with torch.no_grad():
        X = torch.tensor(X, dtype=torch.float32).to(dev)
        he = torch.tensor(host_embed, dtype=torch.float32).to(dev) if host_embed is not None else None
        active_logit, strength = model(X, he)
        active_prob = torch.sigmoid(active_logit).cpu().numpy()
        strength = strength.cpu().numpy()
    return active_prob, strength
