"""
GATE 3 - Baseline 1: mean / majority predictor (the floor).

Held-out evaluation set = fold 0 of the frozen split (fixed across all
baselines in this gate for consistency and speed -- full 5-fold rotation
was not run for the per-host baselines given this environment's compute
constraints, see out/GATE3_MEMO.md "WHAT I COULD NOT DO"). Training pool
= folds 1-4 for that host.

For each host and readout:
  - classifier: predict the majority class (active/inactive) rate from the
    training pool; report MCC (degenerate/0 for a constant predictor -- MCC
    is undefined for a classifier that never varies, reported as such) and
    AUC (0.5 by construction for a constant score).
  - regressor: predict the training-pool mean of the (log1p) strength value
    for all actives in the eval set; report Spearman rho (== 0 by
    construction, since a constant prediction has no rank information --
    reported explicitly rather than omitted, since "the floor is exactly 0"
    is itself the informative fact for the calibration curve's y-intercept).
"""
import numpy as np
import json
from pathlib import Path
from scipy.stats import spearmanr
from sklearn.metrics import matthews_corrcoef, roc_auc_score, mean_squared_error

DATA = Path(__file__).resolve().parent.parent / "data"
CACHE = DATA / "baseline_cache"
OUT = Path(__file__).resolve().parent.parent / "out" / "baselines"
OUT.mkdir(parents=True, exist_ok=True)

HOSTS = ["EC", "BS", "PA"]
TEST_FOLD = 0


def evaluate_readout(d, usable_key, active_key, strength_key, usable_reg_key=None):
    usable = d[usable_key]
    active = d[active_key]
    strength = d[strength_key]
    fold = d["fold"]
    usable_reg = d[usable_reg_key] if usable_reg_key else usable

    train_mask = usable & (fold != TEST_FOLD)
    test_mask = usable & (fold == TEST_FOLD)

    majority_rate = active[train_mask].mean()
    majority_class = int(majority_rate > 0.5)

    y_test_active = active[test_mask].astype(int)
    # constant-probability predictor: predicts majority_rate for everyone
    pred_prob = np.full(y_test_active.shape, majority_rate)
    pred_class = np.full(y_test_active.shape, majority_class)

    mcc = matthews_corrcoef(y_test_active, pred_class) if len(set(y_test_active)) > 1 else float("nan")
    try:
        auc = roc_auc_score(y_test_active, pred_prob) if len(set(y_test_active)) > 1 else float("nan")
    except Exception:
        auc = float("nan")

    # regression on actives (train mean of log1p strength, applied to test actives)
    train_reg_mask = usable_reg & (fold != TEST_FOLD) & active
    test_reg_mask = usable_reg & (fold == TEST_FOLD) & active
    train_strength = np.log1p(np.clip(strength[train_reg_mask], 0, None)) if strength_key == "tx_norm" \
        else strength[train_reg_mask]
    test_strength = np.log1p(np.clip(strength[test_reg_mask], 0, None)) if strength_key == "tx_norm" \
        else strength[test_reg_mask]

    mean_pred = train_strength.mean() if len(train_strength) else float("nan")
    preds = np.full(test_strength.shape, mean_pred)
    rho = spearmanr(test_strength, preds).correlation if len(test_strength) > 1 else float("nan")
    rmse = float(np.sqrt(mean_squared_error(test_strength, preds))) if len(test_strength) else float("nan")

    return {
        "n_train_usable": int(train_mask.sum()), "n_test_usable": int(test_mask.sum()),
        "majority_active_rate_train": float(majority_rate), "mcc": float(mcc) if mcc == mcc else None,
        "auc": float(auc) if auc == auc else None,
        "n_train_actives_for_regression": int(train_reg_mask.sum()),
        "n_test_actives_for_regression": int(test_reg_mask.sum()),
        "spearman_rho": float(rho) if rho == rho else None,
        "rmse_log_strength": rmse if rmse == rmse else None,
    }


def main():
    results = {}
    for host in HOSTS:
        d = dict(np.load(CACHE / f"{host}_baseline_data.npz", allow_pickle=True))
        print(f"\n=== {host} ===")

        tx_result = evaluate_readout(d, "tx_usable", "tx_active", "tx_norm")
        print(f"  transcription: {tx_result}")

        tl_result = evaluate_readout(d, "tl_usable", "tl_active", "protein_log10",
                                      usable_reg_key="tl_usable_for_regression")
        print(f"  translation: {tl_result}")

        results[host] = {"transcription": tx_result, "translation": tl_result}

    with open(OUT / "baseline1_mean_majority.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {OUT / 'baseline1_mean_majority.json'}")


if __name__ == "__main__":
    main()
