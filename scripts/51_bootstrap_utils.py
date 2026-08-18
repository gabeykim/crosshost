"""
GATE 5 - shared 90% percentile bootstrap CI utility, per
out/PREREGISTRATION.md Amendment 2.

Mechanically: resample the fold identities with replacement (5 draws from 5
folds), and for each selected fold, resample from that fold's own raw draws
with replacement too (hierarchical / finer-resolution bootstrap, the option
Amendment 2 names as preferred when per-draw data is available -- it is,
throughout this gate). The statistic per bootstrap iteration is the mean of
the resulting 5 fold-means. 10,000 iterations. 90% CI = 5th/95th percentile
of the resulting distribution.

random seed is fixed for reproducibility -- NOT for tuning the result. The
same seed (0) is used for every bootstrap call in this gate.
"""
import numpy as np

N_BOOT = 10_000
BOOT_SEED = 0


def bootstrap_ci_90(fold_draws, n_boot=N_BOOT, seed=BOOT_SEED):
    """fold_draws: list of arrays (one per fold) of raw draw-level values
    (e.g. Spearman rho per draw). Folds with a single deterministic value
    (e.g. N=0 zero-shot, or N=max-pool with only 1 draw) should pass a
    length-1 array -- the inner resample is then a no-op (always returns
    that same value), which is correct.
    Returns dict: mean, lower (5th pct), upper (95th pct), n_folds, boot_means (for plotting)."""
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
    return {
        "mean": point_mean,
        "lower": float(lower),
        "upper": float(upper),
        "n_folds": n_folds,
        "n_boot": n_boot,
    }


def verdict(model_ci, baseline_ci):
    """Returns one of: MET, NOT MET -- INTERVALS OVERLAP, NOT MET -- POINT ESTIMATE FAVORS BASELINE."""
    if model_ci["mean"] is None or baseline_ci["mean"] is None:
        return "UNDETERMINED (missing data)"
    if model_ci["lower"] > baseline_ci["upper"]:
        return "MET"
    if model_ci["mean"] > baseline_ci["mean"]:
        return "NOT MET -- INTERVALS OVERLAP"
    return "NOT MET -- POINT ESTIMATE FAVORS BASELINE"


def overlap_size(model_ci, baseline_ci):
    """Size of interval overlap (0 if none, positive = extent of overlap)."""
    if model_ci["lower"] is None or baseline_ci["lower"] is None:
        return None
    lo = max(model_ci["lower"], baseline_ci["lower"])
    hi = min(model_ci["upper"], baseline_ci["upper"])
    return max(0.0, hi - lo)
