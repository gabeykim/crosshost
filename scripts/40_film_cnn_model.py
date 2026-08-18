"""
GATE 4 - FiLM-conditioned CNN. Builds on scripts/29_model.py's validated
plain-CNN architecture (213K params, kernel widths 15/9/5/3 matching the
165bp part length verified in Gate 1) and adds:
  - BatchNorm after every conv layer
  - FiLM conditioning at layers 2 and 3: a small MLP maps the host feature
    vector to per-channel (gamma, beta), applied as gamma*x + beta
    immediately after BatchNorm, before the ReLU (standard FiLM placement)
  - A residual connection wrapping layer 4: conv3's post-FiLM, post-ReLU
    output (64 channels) is added to conv4's BN output (64 channels) before
    the final ReLU. DESIGN CHOICE, disclosed: the charter text ("residual
    connection on layers 3-4") is ambiguous about exactly which tensors are
    added; layers 2 and 3 have mismatched channel counts (128 vs 64) so a
    skip from layer-2-output to layer-4-output is not shape-compatible
    without a projection. The chosen interpretation -- a residual block
    wrapping conv4, with conv3's output as the skip -- is shape-compatible
    without an extra projection layer and is the standard "residual block"
    reading of the phrase.
  - Two output heads (transcription, translation), each two-stage
    (active-logit classifier + strength regressor) = 4 scalar outputs from
    ONE shared trunk, trained JOINTLY -- unlike Gate 3's B2 CNN validation,
    which trained separate model instances per readout as a disclosed
    simplification. Gate 4 implements the true joint architecture the
    charter specifies.

MPS PERFORMANCE NOTES (profiled during Gate 4 prep, both load-bearing for
every training script that follows):
  1. Boolean-mask loss indexing (preds[...][mask]) cost ~875ms/batch of pure
     MPS gather overhead for this model -- fixed by computing an elementwise
     (reduction='none') loss and multiplying by a float mask instead (see
     compute_loss below). Cut per-batch cost roughly in half.
  2. Batch size has a sharp, non-monotonic cliff on this backend for this
     architecture: batch=1024 measured at 0.52ms/example; batch=2048 at
     1.00ms/example (2x worse); batch=4096 and 8192 at 16-18ms/example
     (30-35x worse than batch=1024). This is NOT the smooth "larger batch
     amortizes overhead" curve seen on CUDA/CPU -- something in Apple's
     Metal backend (likely memory-pressure-triggered kernel fallback, given
     this machine's 8GB unified memory) falls off a cliff above ~2048 for
     this model. batch_size=1024 is used everywhere as a result; DO NOT
     increase it without re-profiling first.

Host feature normalization: z-scored using statistics computed across ALL
6 hosts (not fold-local training-hosts-only). DISCLOSED DEVIATION from SR3's
literal fold-local rule, with reasoning: SR3 governs statistics fit from the
sequence-level TARGET/label data (e.g. the translation floor value), where
using test-fold information would leak label information. Host genomic/
physiology features are PUBLIC COVARIATES describing host identity, not
labels -- a practitioner deploying this model on a genuinely new host would
know that host's genome sequence and proteome composition regardless of
whether any regulatory-activity data exists for it yet. Standardizing by the
full 6-host distribution (computed once, fixed, not re-fit per fold) does
not leak any regulatory-activity information about the held-out host; it is
the domain-generalization-literature convention of treating domain
descriptors as known. With only 2 training hosts per LOHO fold, computing
z-score statistics from n=2 would in any case be numerically unstable.
"""
import torch
import torch.nn as nn
import numpy as np
import pandas as pd
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data"
CACHE = DATA / "baseline_cache"
TARGET_KEYS = ["tx_active", "tx_strength", "tl_active", "tl_strength"]

GENOMIC_FEATURE_COLS = None  # populated by load_host_features()
PHYSIOLOGY_FEATURE_COLS = ["ribosomal_proteins_fraction", "rnap_core_fraction", "sigma_factors_fraction",
                           "chaperones_fraction", "elongation_factors_fraction", "growth_rate_mu_h"]


def load_host_features():
    """Returns (genomic_df, genomic_cols, physiology_df, physiology_cols),
    each df indexed by host code, z-scored across all 6 hosts."""
    global GENOMIC_FEATURE_COLS
    g = pd.read_parquet(DATA / "hosts_genomic.parquet").set_index("host")
    genomic_cols = [c for c in g.columns if g[c].dtype != object and g[c].dtype != "string"]
    GENOMIC_FEATURE_COLS = genomic_cols
    g_num = g[genomic_cols].astype(np.float64)
    g_z = (g_num - g_num.mean()) / g_num.std().replace(0, 1.0)

    p = pd.read_parquet(DATA / "hosts_physiology.parquet").set_index("host")
    p_num = p[PHYSIOLOGY_FEATURE_COLS].astype(np.float64)

    # DATA GAP, found while wiring up Gate 4 (not previously blocking anything
    # downstream of Gate 1, which only flagged it as "unverified/blocked"):
    # B. subtilis's growth_rate_mu_h is NaN -- the Zhu et al. 2025 PNAS
    # supplementary table value was never extractable (Gate 1 finding,
    # unresolved since). BS is both a training host (EC/PA held out) and
    # THE primary held-out host (H-MAIN) for the physiology variant, so a
    # missing feature here cannot be silently dropped or left as NaN.
    # FIX, disclosed: impute as the mean of the other 5 hosts' raw
    # growth_rate_mu_h values -- a neutral, conservative fallback, not a
    # fabricated literature citation. Flagged in the returned imputed-cells
    # list; also reported in out/GATE4_MEMO.md WHAT I COULD NOT DO.
    imputed_cells = []
    for col in p_num.columns:
        nan_hosts = p_num.index[p_num[col].isna()]
        for h in nan_hosts:
            fill_val = p_num[col].drop(index=h).mean()
            p_num.loc[h, col] = fill_val
            imputed_cells.append({"host": h, "feature": col, "imputed_value": float(fill_val),
                                   "method": "mean_of_other_5_hosts"})
            print(f"  [IMPUTED] {h}.{col} = {fill_val:.4f} (mean of other 5 hosts; "
                  f"original source value unavailable, see scripts/40 docstring)")

    p_z = (p_num - p_num.mean()) / p_num.std().replace(0, 1.0)

    return g_z, genomic_cols, p_z, PHYSIOLOGY_FEATURE_COLS, imputed_cells


class FiLMGenerator(nn.Module):
    def __init__(self, host_dim, hidden=32, film2_channels=128, film3_channels=64):
        super().__init__()
        self.film2_channels = film2_channels
        self.film3_channels = film3_channels
        out_dim = 2 * film2_channels + 2 * film3_channels
        self.net = nn.Sequential(nn.Linear(host_dim, hidden), nn.ReLU(), nn.Linear(hidden, out_dim))

    def forward(self, host_vec):
        out = self.net(host_vec)
        c2, c3 = self.film2_channels, self.film3_channels
        gamma2, beta2, gamma3, beta3 = torch.split(out, [c2, c2, c3, c3], dim=1)
        return gamma2, beta2, gamma3, beta3


class FiLMSequenceCNN(nn.Module):
    """~227K params (genomic variant) / ~226K (physiology variant).
    Input: seq (batch, 4, 165) one-hot, host_vec (batch, host_dim).
    Outputs: dict with tx_active_logit, tx_strength, tl_active_logit, tl_strength.
    use_film=False disables conditioning entirely (host_vec ignored, gamma
    fixed at 1, beta fixed at 0) -- the sanity-check control condition."""

    def __init__(self, host_dim, use_film=True):
        super().__init__()
        self.use_film = use_film
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
        self.film_gen = FiLMGenerator(host_dim, film2_channels=128, film3_channels=64) if use_film else None

    def forward(self, seq, host_vec=None):
        b = seq.shape[0]
        if self.use_film and host_vec is not None:
            gamma2, beta2, gamma3, beta3 = self.film_gen(host_vec)
        else:
            gamma2 = torch.ones(b, 128, device=seq.device)
            beta2 = torch.zeros(b, 128, device=seq.device)
            gamma3 = torch.ones(b, 64, device=seq.device)
            beta3 = torch.zeros(b, 64, device=seq.device)

        h = torch.relu(self.bn1(self.conv1(seq)))
        h = self.bn2(self.conv2(h))
        h = gamma2.unsqueeze(-1) * h + beta2.unsqueeze(-1)
        h = torch.relu(h)
        h = self.dropout(h)

        h3 = self.bn3(self.conv3(h))
        h3 = gamma3.unsqueeze(-1) * h3 + beta3.unsqueeze(-1)
        h3 = torch.relu(h3)

        h4 = self.bn4(self.conv4(h3))
        h4 = torch.relu(h4 + h3)  # residual block wrapping conv4, skip = conv3's output

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


def get_film_gamma_beta_stats(model, host_vecs):
    """host_vecs: (n_hosts, host_dim) tensor (any device). Returns per-layer
    gamma/beta stats (mean, std, min, max) across the given hosts, for the
    sanity-check report."""
    model.eval()
    dev = next(model.parameters()).device
    with torch.no_grad():
        gamma2, beta2, gamma3, beta3 = model.film_gen(host_vecs.to(dev))
    stats = {}
    for name, t in [("gamma2", gamma2), ("beta2", beta2), ("gamma3", gamma3), ("beta3", beta3)]:
        arr = t.cpu().numpy()
        stats[name] = {"mean": float(arr.mean()), "std": float(arr.std()),
                       "min": float(arr.min()), "max": float(arr.max())}
    return stats


FREEZE_MODES = (None, "head_only", "top_conv")


def set_freeze_mode(model, mode):
    """None: everything trainable (Task 2, from-scratch LOHO training).
    'head_only': freeze conv1-4/bn1-4 entirely; only FiLM generator,
    fc_shared, and the 4 output heads train -- Task 3's N<=100 transfer
    mechanism (frozen BN running stats too, which is what keeps small-batch
    fine-tuning numerically sane at N as low as 10).
    'top_conv': as head_only, but also unfreeze conv4/bn4 -- Task 3's
    N>=300 transfer mechanism."""
    assert mode in FREEZE_MODES
    always_trainable = [model.film_gen, model.fc_shared, model.head_tx_active,
                         model.head_tx_strength, model.head_tl_active, model.head_tl_strength]
    frozen_modules = [model.conv1, model.bn1, model.conv2, model.bn2, model.conv3, model.bn3,
                      model.conv4, model.bn4]
    if mode is None:
        for m in frozen_modules:
            for p in m.parameters():
                p.requires_grad = True
        for m in frozen_modules:
            m.train()
        return
    for m in frozen_modules:
        for p in m.parameters():
            p.requires_grad = False
        m.eval()  # freezes BN running stats too
    if mode == "top_conv":
        for p in model.conv4.parameters():
            p.requires_grad = True
        for p in model.bn4.parameters():
            p.requires_grad = True
        model.conv4.train()
        model.bn4.train()
    for m in always_trainable:
        if m is not None:
            for p in m.parameters():
                p.requires_grad = True


LOSS_HEADS = [("tx_active", "bce"), ("tx_strength", "mse"), ("tl_active", "bce"), ("tl_strength", "mse")]


def compute_loss(preds, targets, masks, bce_none, mse_none):
    """VECTORIZED masked loss -- no boolean indexing. Profiled during Gate 4
    Task 2 prep: boolean-mask gather (preds[...][mask]) cost ~875ms/batch of
    pure overhead on MPS for this model (measured: 1607ms/batch with
    boolean indexing vs 732ms/batch on the identical unmasked forward+
    backward), evidently because MPS's gather/compaction kernels are not
    well optimized for this shape. Fix: multiply the elementwise
    (reduction='none') loss by a float mask and normalize by the mask sum,
    entirely avoiding dynamic-shape gather ops. Cut the masking overhead to
    ~200ms/batch (a remaining, accepted cost). bce_none/mse_none must be
    constructed with reduction='none'."""
    total = 0.0
    for name, kind in LOSS_HEADS:
        mask = masks[name].float()
        denom = mask.sum().clamp(min=1.0)
        p = preds[f"{name}_logit" if kind == "bce" else name]
        y = targets[name]
        elementwise = bce_none(p, y) if kind == "bce" else mse_none(p, y)
        loss = (elementwise * mask).sum() / denom
        total = total + loss
    return total, {}


def train_film_model(seq, host_vec, targets, masks, host_dim, epochs=40, lr=1e-3, weight_decay=1e-4,
                      batch_size=1024, seed=0, device="cpu", use_film=True, patience=8,
                      verbose=False, freeze_mode=None, init_model=None):
    """seq: (N,4,165) float32. host_vec: (N,host_dim) float32 (per-example,
    repeated per that example's host). targets/masks: dicts of (N,) float32 /
    bool arrays for tx_active, tx_strength, tl_active, tl_strength.
    init_model: an existing FiLMSequenceCNN to continue training (fine-tuning,
    Task 3) instead of a fresh random init (Task 2)."""
    torch.manual_seed(seed)
    dev = torch.device(device)
    model = init_model.to(dev) if init_model is not None else FiLMSequenceCNN(host_dim, use_film=use_film).to(dev)
    if freeze_mode is not None:
        set_freeze_mode(model, freeze_mode)
    trainable = [p for p in model.parameters() if p.requires_grad]
    opt = torch.optim.AdamW(trainable, lr=lr, weight_decay=weight_decay)
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
        if freeze_mode is not None:
            set_freeze_mode(model, freeze_mode)  # re-assert eval() on frozen BN after model.train()
        perm = torch.randperm(n)
        epoch_loss, n_batches = 0.0, 0
        for i in range(0, n, batch_size):
            idx = perm[i:i + batch_size]
            opt.zero_grad()
            preds = model(X[idx], H[idx])
            targets_b = {k: v[idx] for k, v in T.items()}
            masks_b = {k: v[idx] for k, v in M.items()}
            loss, _ = compute_loss(preds, targets_b, masks_b, bce, mse)
            loss.backward()
            opt.step()
            epoch_loss += loss.item()
            n_batches += 1
        sched.step()
        epochs_run = epoch + 1
        if n_batches == 0:
            break
        epoch_loss /= n_batches
        if verbose and epoch % 5 == 0:
            print(f"    epoch {epoch}: loss={epoch_loss:.4f}")
        if epoch_loss < best_loss - 1e-5:
            best_loss = epoch_loss
            patience_ctr = 0
        else:
            patience_ctr += 1
            if patience_ctr >= patience and epoch > 10:
                break
    model.eval()
    return model, epochs_run


def build_target_arrays(d):
    """From a baseline_cache dict (one host), returns (targets, masks) dicts
    row-aligned to d['onehot'], ready for train_film_model. Uses the
    floor-corrected translation definitions (Gate 3): tl_strength target is
    raw protein_log10, masked by tl_usable_for_regression & tl_active (i.e.
    excludes floor-pinned rows from the regression loss); tx_strength is
    log1p(tx_norm), matching scripts/31's convention exactly.

    NaN/INF SANITIZATION, required by the vectorized masked loss
    (compute_loss): protein_log10 is NaN, and tx_norm is literal +Inf
    (398 EC rows / 613 PA rows, verified during Gate 4 debugging), at
    exactly the rows where their masks are False -- these are rows with no
    measurement at all (or, for tx_norm, a division/log artifact of the
    paper's own normalization on an unusable row), not floor-pinned values,
    which are numeric. compute_loss multiplies an elementwise loss by a 0/1
    mask BEFORE summing; 0 * NaN = NaN in IEEE float, so an unmasked NaN
    silently poisons the entire batch loss -- caught immediately (loss was
    NaN from epoch 0). FIRST FIX ATTEMPT (nan_to_num with only nan=0.0
    specified) was INCOMPLETE and re-broke silently: nan_to_num's default
    posinf replacement is ~3.4e38 (float32 max), not 0 -- so the +Inf rows
    became a finite-but-enormous number instead of NaN, which still blew up
    MSE to non-finite (verified: elementwise MSE loss was non-finite even
    though every input tensor read as "finite" by torch.isfinite, because
    3.4e38 squared overflows). Fix: explicitly zero out nan, posinf, AND
    neginf. The mask still correctly excludes these rows from contributing
    any gradient, since 0 * finite_loss_at_that_row = 0."""
    targets = {
        "tx_active": d["tx_active"].astype(np.float32),
        "tx_strength": np.nan_to_num(np.log1p(np.clip(d["tx_norm"], 0, None)).astype(np.float32),
                                      nan=0.0, posinf=0.0, neginf=0.0),
        "tl_active": d["tl_active"].astype(np.float32),
        "tl_strength": np.nan_to_num(d["protein_log10"].astype(np.float32), nan=0.0, posinf=0.0, neginf=0.0),
    }
    masks = {
        "tx_active": d["tx_usable"].astype(bool),
        "tx_strength": (d["tx_usable"].astype(bool) & d["tx_active"].astype(bool)),
        "tl_active": d["tl_usable"].astype(bool),
        "tl_strength": (d["tl_usable_for_regression"].astype(bool) & d["tl_active"].astype(bool)),
    }
    return targets, masks


def build_pooled_arrays(hosts, host_feature_df, fold_filter_fn=None):
    """hosts: list of host codes to pool. host_feature_df: z-scored df
    indexed by host (genomic or physiology, from load_host_features()).
    fold_filter_fn(fold_array)->bool row mask, or None for all rows of that
    host. Returns seq, host_vec (per-row, that row's host's fixed feature
    vector), targets, masks, host_of_row (diagnostic)."""
    seq_parts, hv_parts, host_row_parts = [], [], []
    targets_parts = {k: [] for k in TARGET_KEYS}
    masks_parts = {k: [] for k in TARGET_KEYS}
    for h in hosts:
        d = dict(np.load(CACHE / f"{h}_baseline_data.npz", allow_pickle=True))
        fold = d["fold"]
        row_mask = fold_filter_fn(fold) if fold_filter_fn is not None else np.ones(len(fold), dtype=bool)
        idx = np.where(row_mask)[0]
        t, m = build_target_arrays(d)
        seq_parts.append(d["onehot"][idx])
        hv = np.tile(host_feature_df.loc[h].values.astype(np.float32), (len(idx), 1))
        hv_parts.append(hv)
        host_row_parts.append(np.array([h] * len(idx)))
        for k in TARGET_KEYS:
            targets_parts[k].append(t[k][idx])
            masks_parts[k].append(m[k][idx])
    seq = np.concatenate(seq_parts)
    host_vec = np.concatenate(hv_parts)
    host_of_row = np.concatenate(host_row_parts)
    targets = {k: np.concatenate(v) for k, v in targets_parts.items()}
    masks = {k: np.concatenate(v) for k, v in masks_parts.items()}
    return seq, host_vec, targets, masks, host_of_row


def predict_film(model, seq, host_vec, batch_size=1024):
    """GATE 7 FIX: same unbatched-forward-pass performance cliff as
    predict_seqonly (scripts/58) -- see that function's docstring for the
    measured numbers. Chunked at batch_size=1024, numerically verified
    identical to the unbatched path (scripts/audit_batching_fix.py)."""
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
