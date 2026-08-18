"""Evaluation API: give it predictions, get back the full metric suite with
fold-resolved 90% percentile bootstrap intervals -- the identical protocol
used throughout the CROSSHOST project (Gates 5-8), ported verbatim from
scripts/51_bootstrap_utils.py and scripts/42_loho_training.py's evaluate().
"""
import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import matthews_corrcoef, roc_auc_score

N_BOOT = 10_000
BOOT_SEED = 0


def bootstrap_ci_90(fold_draws, n_boot=N_BOOT, seed=BOOT_SEED):
    """fold_draws: list of arrays (one per fold) of raw draw-level values
    (e.g. Spearman rho per draw, or a length-1 array for a deterministic
    single evaluation). Hierarchical bootstrap: resample fold identities
    with replacement, then resample within each selected fold's own draws
    with replacement; the statistic per iteration is the mean of the
    resulting fold-means. 10,000 iterations, 90% CI = 5th/95th percentile.

    Returns dict: mean, lower, upper, n_folds, n_boot."""
    fold_draws = [np.asarray(fd, dtype=np.float64) for fd in fold_draws if len(fd) > 0]
    n_folds = len(fold_draws)
    if n_folds == 0:
        return {"mean": None, "lower": None, "upper": None, "n_folds": 0}
    rng = np.random.default_rng(seed)
    boot_means = np.empty(n_boot, dtype=np.float64)
    for b in range(n_boot):
        fold_idx = rng.integers(0, n_folds, size=n_folds)
        iter_fold_means = np.empty(n_folds, dtype=np.float64)
        for j, fi in enumerate(fold_idx):
            draws = fold_draws[fi]
            resampled = rng.choice(draws, size=len(draws), replace=True)
            iter_fold_means[j] = resampled.mean()
        boot_means[b] = iter_fold_means.mean()
    lower, upper = np.percentile(boot_means, [5, 95])
    point_mean = float(np.mean([fd.mean() for fd in fold_draws]))
    return {"mean": point_mean, "lower": float(lower), "upper": float(upper),
            "n_folds": n_folds, "n_boot": n_boot}


def verdict(model_ci, baseline_ci):
    """MET / NOT MET -- INTERVALS OVERLAP / NOT MET -- POINT ESTIMATE FAVORS BASELINE."""
    if model_ci["mean"] is None or baseline_ci["mean"] is None:
        return "UNDETERMINED (missing data)"
    if model_ci["lower"] > baseline_ci["upper"]:
        return "MET"
    if model_ci["mean"] > baseline_ci["mean"]:
        return "NOT MET -- INTERVALS OVERLAP"
    return "NOT MET -- POINT ESTIMATE FAVORS BASELINE"


def _evaluate_one_readout(active_prob, active_true, strength_pred, strength_true):
    """active_prob: predicted P(active) in [0,1]. active_true: bool/0-1.
    strength_pred/strength_true: continuous strength predictions/targets,
    evaluated only where strength_true is not NaN (i.e. the active subset)."""
    active_prob = np.asarray(active_prob)
    active_true = np.asarray(active_true).astype(int)
    mcc = matthews_corrcoef(active_true, (active_prob > 0.5).astype(int)) if len(set(active_true)) > 1 else None
    try:
        auc = roc_auc_score(active_true, active_prob) if len(set(active_true)) > 1 else None
    except Exception:
        auc = None

    strength_pred = np.asarray(strength_pred)
    strength_true = np.asarray(strength_true, dtype=float)
    mask = ~np.isnan(strength_true)
    rho = None
    if mask.sum() > 1 and len(set(np.round(strength_true[mask], 6))) > 1:
        rho = spearmanr(strength_true[mask], strength_pred[mask]).correlation

    return {"n_active_eval": int(len(active_true)), "n_strength_eval": int(mask.sum()),
            "mcc": float(mcc) if mcc is not None and mcc == mcc else None,
            "auc": float(auc) if auc is not None and auc == auc else None,
            "spearman_rho": float(rho) if rho is not None and rho == rho else None}


def evaluate(predictions_per_fold, targets_per_fold, readouts=("transcription", "translation")):
    """The public evaluation entry point.

    predictions_per_fold / targets_per_fold: each a list (one entry per
    fold) of dicts, one dict per readout, each dict with keys
    'active_prob' (predicted P(active), or None if not scoring classification),
    'active_true', 'strength_pred' (or None if not scoring regression),
    'strength_true'. Missing readouts for a fold may be omitted.

    Returns: {readout: {'mcc': {...ci...}, 'auc': {...ci...}, 'spearman_rho': {...ci...},
                          'per_fold': [...]}}

    This is the SAME metric suite and SAME 90% bootstrap protocol used
    throughout this project's own gates (H-MAIN, H-SCIENCE, H-DIAGNOSTIC,
    the sequence-only ablation, the foundation-model comparison, and Gate 7's
    conformal-calibration work) -- a submission scored with this function is
    directly comparable to every number in out/PAPER_FRAMING.md."""
    n_folds = len(predictions_per_fold)
    out = {}
    for readout in readouts:
        per_fold_results = []
        mcc_draws, auc_draws, rho_draws = [], [], []
        for f in range(n_folds):
            pred = predictions_per_fold[f].get(readout)
            targ = targets_per_fold[f].get(readout)
            if pred is None or targ is None:
                continue
            r = _evaluate_one_readout(pred.get("active_prob"), targ["active_true"],
                                       pred.get("strength_pred"), targ["strength_true"])
            per_fold_results.append(r)
            if r["mcc"] is not None:
                mcc_draws.append(np.array([r["mcc"]]))
            if r["auc"] is not None:
                auc_draws.append(np.array([r["auc"]]))
            if r["spearman_rho"] is not None:
                rho_draws.append(np.array([r["spearman_rho"]]))
        out[readout] = {
            "mcc": bootstrap_ci_90(mcc_draws),
            "auc": bootstrap_ci_90(auc_draws),
            "spearman_rho": bootstrap_ci_90(rho_draws),
            "per_fold": per_fold_results,
        }
    return out
