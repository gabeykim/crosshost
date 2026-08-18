"""
GATE 3.5 - Task 1: is B. subtilis translation-regression data sufficient to
support H-MAIN as currently specified?

H-MAIN (charter Gate 5): cross-host model at N=100 matches/beats the
per-host-only baseline at N>=3000, on B. subtilis held out.

The Gate 3 floor-value fix (scripts/28) left B. subtilis with only 1,101
usable translation-regression examples FOLD-WIDE (all 5 folds combined), not
the ~11,564 raw "usable" count. This script establishes, with numbers, exact
per-fold usable N for both the regression task and the classification task
(which is less affected -- floor-pinned rows are still informative as
"inactive"), whether the N=3000 training-pool arm is actually reachable for
ANY fold-rotation, extends the calibration curve to the real ceiling using
the same k-mer/linear substitute as Gate 3's B2, and lays out (without
choosing) the H-MAIN respecification options this implies.
"""
import numpy as np
import json
import sys
from pathlib import Path
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import matthews_corrcoef, roc_auc_score
import warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
b2mod = import_module("31_baseline2_calibration")

DATA = Path(__file__).resolve().parent.parent / "data"
CACHE = DATA / "baseline_cache"
OUT = Path(__file__).resolve().parent.parent / "out"
OUT.mkdir(parents=True, exist_ok=True)

N_FOLDS = 5
LOGREG_C = 1.0
RIDGE_ALPHA = 10.0
CURVE_N_VALUES_BASE = [10, 30, 100, 300]  # + max-available, computed below
N_DRAWS = 10


def per_fold_counts(d):
    fold = d["fold"]
    rows = []
    for f in range(N_FOLDS):
        mask = fold == f
        rows.append({
            "fold": f,
            "n_total_rows": int(mask.sum()),
            "tl_usable": int((mask & d["tl_usable"]).sum()),
            "tl_usable_for_regression": int((mask & d["tl_usable_for_regression"]).sum()),
            "tl_active_of_usable": int((mask & d["tl_active"]).sum()),
            "tl_inactive_of_usable": int((mask & d["tl_usable"] & ~d["tl_active"]).sum()),
        })
    return rows


def n3000_arm_check(d):
    """For each fold held out as TEST, is train-pool >= 3000 for regression?
    Also report the classification-task train-pool (larger, less affected)."""
    fold = d["fold"]
    rows = []
    for test_fold in range(N_FOLDS):
        train_mask = fold != test_fold
        reg_pool = int((train_mask & d["tl_usable_for_regression"]).sum())
        clf_pool = int((train_mask & d["tl_usable"]).sum())
        test_reg = int((~train_mask & d["tl_usable_for_regression"]).sum())
        test_clf = int((~train_mask & d["tl_usable"]).sum())
        rows.append({
            "held_out_fold": test_fold,
            "regression_train_pool": reg_pool,
            "regression_test_n": test_reg,
            "n3000_arm_computable_regression": reg_pool >= 3000,
            "classification_train_pool": clf_pool,
            "classification_test_n": test_clf,
            "n3000_arm_computable_classification": clf_pool >= 3000,
        })
    return rows


def learning_curve(d, kmer_features, test_fold, max_n, task):
    """task: 'regression' (Ridge on real, non-floor protein_log10 values) or
    'classification' (LogisticRegression, active-vs-inactive, floor rows
    count as negative/inactive class)."""
    fold = d["fold"]
    if task == "regression":
        usable_reg = d["tl_usable_for_regression"]
        train_pool_idx = np.where(usable_reg & (fold != test_fold))[0]
        test_idx = np.where(usable_reg & (fold == test_fold))[0]
        y_all = d["protein_log10"]
    else:
        usable = d["tl_usable"]
        train_pool_idx = np.where(usable & (fold != test_fold))[0]
        test_idx = np.where(usable & (fold == test_fold))[0]
        y_active_all = d["tl_active"].astype(int)

    n_values = [n for n in CURVE_N_VALUES_BASE if n < len(train_pool_idx)] + [len(train_pool_idx)]
    n_values = sorted(set(n_values))

    curve = {}
    for N in n_values:
        draws = []
        n_draws_here = 1 if N == len(train_pool_idx) else N_DRAWS
        for draw_i in range(n_draws_here):
            rng = np.random.default_rng(2000 * N + draw_i)
            sample_idx = train_pool_idx if N >= len(train_pool_idx) else \
                rng.choice(train_pool_idx, size=N, replace=False)

            if task == "regression":
                X_train = kmer_features[sample_idx]
                y_train = y_all[sample_idx]
                X_test = kmer_features[test_idx]
                y_test = y_all[test_idx]
                if len(set(np.round(y_train, 6))) > 1 and len(sample_idx) >= 3:
                    reg = Ridge(alpha=RIDGE_ALPHA).fit(X_train, y_train)
                    y_pred = reg.predict(X_test)
                    rho = spearmanr(y_test, y_pred).correlation if len(set(np.round(y_test, 6))) > 1 else None
                else:
                    rho = None
                draws.append({"spearman_rho": float(rho) if rho is not None and rho == rho else None,
                               "n_actual_train": int(len(sample_idx))})
            else:
                X_train = kmer_features[sample_idx]
                y_train = y_active_all[sample_idx]
                X_test = kmer_features[test_idx]
                y_test = y_active_all[test_idx]
                if len(set(y_train)) > 1:
                    clf = LogisticRegression(C=LOGREG_C, max_iter=300).fit(X_train, y_train)
                    prob = clf.predict_proba(X_test)[:, 1]
                    pred = (prob > 0.5).astype(int)
                    mcc = matthews_corrcoef(y_test, pred) if len(set(y_test)) > 1 else None
                    try:
                        auc = roc_auc_score(y_test, prob) if len(set(y_test)) > 1 else None
                    except Exception:
                        auc = None
                else:
                    mcc, auc = None, None
                draws.append({"mcc": float(mcc) if mcc is not None and mcc == mcc else None,
                               "auc": float(auc) if auc is not None and auc == auc else None,
                               "n_actual_train": int(len(sample_idx))})
        curve[N] = draws
    return curve, {"train_pool_total": int(len(train_pool_idx)), "test_n": int(len(test_idx))}


def summarize_curve(curve, metric):
    out = {}
    for N, draws in curve.items():
        vals = [dr[metric] for dr in draws if dr.get(metric) is not None]
        out[N] = {"mean": float(np.mean(vals)) if vals else None,
                   "std": float(np.std(vals)) if vals else None,
                   "n_draws": len(draws)}
    return out


def reconcile_with_gate3_nominal_scheme(d, kmer_features, test_fold, n_values=(10, 30, 100, 300, 3000)):
    """Gate 3's original B2 script draws N from the FULL tl_usable pool
    (floor-pinned + active mixed) and fits the regressor only on whichever
    sampled rows happen to be non-floor -- i.e. nominal N != actual
    regression-training N. This reproduces that scheme exactly (same
    sampling logic as scripts/31) to report the EFFECTIVE regression N per
    nominal N, explaining why this script's N=100 (rho=0.28, drawn directly
    from regression-usable rows) looks nothing like Gate 3's reported
    N=100 (rho=0.071, drawn from the full usable pool)."""
    fold = d["fold"]
    usable = d["tl_usable"]
    active = d["tl_active"]
    usable_reg = d["tl_usable_for_regression"]
    strength = d["protein_log10"]
    train_pool_idx = np.where(usable & (fold != test_fold))[0]
    test_idx = np.where(usable & (fold == test_fold))[0]
    X_test = kmer_features[test_idx]
    active_mask_test = active[test_idx]

    rows = {}
    for N in n_values:
        effective_ns, rhos = [], []
        n_draws_here = 10
        for draw_i in range(n_draws_here):
            rng = np.random.default_rng(3000 * N + draw_i)
            sample_idx = train_pool_idx if N >= len(train_pool_idx) else \
                rng.choice(train_pool_idx, size=N, replace=False)
            active_mask_train = active[sample_idx] & usable_reg[sample_idx]
            effective_ns.append(int(active_mask_train.sum()))
            if active_mask_train.sum() >= 3:
                Xa = kmer_features[sample_idx][active_mask_train]
                ya = strength[sample_idx][active_mask_train]
                if len(set(np.round(ya, 6))) > 1:
                    reg = Ridge(alpha=RIDGE_ALPHA).fit(Xa, ya)
                    if active_mask_test.sum() > 1:
                        y_pred = reg.predict(X_test[active_mask_test])
                        y_true = strength[test_idx][active_mask_test]
                        if len(set(np.round(y_true, 6))) > 1:
                            rhos.append(spearmanr(y_true, y_pred).correlation)
        rows[N] = {
            "nominal_n": N,
            "mean_effective_regression_n": float(np.mean(effective_ns)),
            "mean_spearman_rho": float(np.mean(rhos)) if rhos else None,
            "n_valid_draws": len(rhos),
        }
    return rows


def interpret_headroom(reg_curve_summary, max_n_reached, floor_frac):
    """Decide: genuine headroom / data floor / undetermined, from the shape
    of the curve (still rising at max N?) and the noise level at max N."""
    ns = sorted(reg_curve_summary.keys())
    if len(ns) < 3:
        return "UNDETERMINED", "Fewer than 3 usable N points on the curve -- not enough resolution to judge curve shape."
    last3 = ns[-3:]
    vals = [reg_curve_summary[n]["mean"] for n in last3 if reg_curve_summary[n]["mean"] is not None]
    if len(vals) < 3:
        return "UNDETERMINED", "Regression rho is undefined/None at one or more of the top-3 N points (too few actives or degenerate fit)."
    slope_last_segment = vals[-1] - vals[-2]
    slope_first_segment = vals[-2] - vals[-3]
    still_rising = slope_last_segment > 0.01
    std_at_max = reg_curve_summary[ns[-1]]["std"]
    mean_at_max = reg_curve_summary[ns[-1]]["mean"]
    noisy = (std_at_max is not None and mean_at_max is not None and mean_at_max > 0
             and std_at_max / max(mean_at_max, 1e-6) > 0.5)
    if still_rising and not noisy:
        return "GENUINE HEADROOM", (
            f"rho is still increasing from N={ns[-2]} (rho={vals[-2]:.3f}) to N={ns[-1]} "
            f"(rho={vals[-1]:.3f}), a rise of {slope_last_segment:+.3f} over the last step, "
            f"and is not yet dominated by noise at the largest N reached. The per-host baseline "
            f"has not exhausted its own achievable ceiling -- it is plausible a cross-host model "
            f"with more effective N (via pooling) could do better still.")
    if noisy:
        return "DATA FLOOR / UNDETERMINED", (
            f"At the largest N reached ({ns[-1]}, n_draws={reg_curve_summary[ns[-1]]['n_draws']}), "
            f"draw-to-draw std ({std_at_max:.3f}) is large relative to the mean ({mean_at_max:.3f}) -- "
            f"the estimate itself is not stable enough to confidently call this either headroom or a "
            f"floor. With only {max_n_reached} total usable examples ({floor_frac:.1%} of raw 'usable' "
            f"rows were floor-pinned and excluded), a single-draw estimate at max N has no error bar at all.")
    return "DATA FLOOR", (
        f"rho has plateaued or is declining by N={ns[-1]} (rho={vals[-1]:.3f} vs {vals[-2]:.3f} at "
        f"N={ns[-2]}) -- more data at the SAME sequence-only, no-host-conditioning setup would not "
        f"obviously help further; the constraint looks like task difficulty, not data starvation.")


def premise_correction():
    """Gate 3's B2 (scripts/31) draws nominal N from the FULL tl_usable pool
    (9,300/fold for BS -- floor-pinned + active mixed), not from the
    943-per-fold tl_usable_for_regression pool. This means the N=3,000 arm
    WAS already computable and WAS already run in Gate 3 -- pulled directly
    from the archived out/baselines/baseline2_calibration_raw.json (source
    of truth, not re-derived) to correct the premise that motivated this
    task with numbers, not assertion."""
    with open(Path(__file__).resolve().parent.parent / "out" / "baselines" / "baseline2_calibration_raw.json") as f:
        b2 = json.load(f)
    bs_tl = b2["BS"]["results"]["translation"]
    meta = b2["BS"]["meta"]["translation"]
    rows = {}
    for N in ["10", "30", "100", "300", "1000", "3000"]:
        draws = bs_tl[N]
        rhos = [x["spearman_rho"] for x in draws if x["spearman_rho"] is not None]
        rows[N] = {
            "n_draws_with_valid_rho": len(rhos), "n_draws_total": len(draws),
            "mean_rho": float(np.mean(rhos)) if rhos else None,
            "std_rho": float(np.std(rhos)) if rhos else None,
            "n_actual_train_drawn": draws[0]["n_actual_train"] if draws else None,
        }
    return {"train_pool_used_by_gate3_b2": meta["n_train_pool"],
            "note": "This IS the tl_usable pool (floor+active mixed), not tl_usable_for_regression -- "
                    "confirms N=3000 was a genuine 3000-example draw from a 9300-example pool, not capped.",
            "archived_gate3_numbers_by_nominal_n": rows}


def main():
    d = dict(np.load(CACHE / "BS_baseline_data.npz", allow_pickle=True))
    kmer_features = b2mod.build_kmer_features(d["onehot"])
    floor_frac = 0.899  # from Gate 3 (out/state.json gate_3.translation_floor_value_finding)

    print("=== Per-fold usable-N breakdown (B. subtilis, translation) ===")
    fold_counts = per_fold_counts(d)
    for r in fold_counts:
        print(f"  fold{r['fold']}: n_rows={r['n_total_rows']}, tl_usable={r['tl_usable']}, "
              f"tl_usable_for_regression={r['tl_usable_for_regression']}, "
              f"active={r['tl_active_of_usable']}, inactive={r['tl_inactive_of_usable']}")
    total_usable = sum(r["tl_usable"] for r in fold_counts)
    total_usable_reg = sum(r["tl_usable_for_regression"] for r in fold_counts)
    total_active = sum(r["tl_active_of_usable"] for r in fold_counts)
    total_inactive = sum(r["tl_inactive_of_usable"] for r in fold_counts)
    print(f"  TOTAL: tl_usable={total_usable}, tl_usable_for_regression={total_usable_reg}, "
          f"active={total_active}, inactive={total_inactive}")

    print("\n=== N=3000 baseline-arm computability, per held-out-fold scenario ===")
    arm_check = n3000_arm_check(d)
    for r in arm_check:
        print(f"  held_out_fold={r['held_out_fold']}: regression_train_pool={r['regression_train_pool']} "
              f"(N=3000 computable: {r['n3000_arm_computable_regression']}), "
              f"classification_train_pool={r['classification_train_pool']} "
              f"(N=3000 computable: {r['n3000_arm_computable_classification']})")
    max_reg_pool = max(r["regression_train_pool"] for r in arm_check)
    max_clf_pool = max(r["classification_train_pool"] for r in arm_check)
    any_n3000_reg = any(r["n3000_arm_computable_regression"] for r in arm_check)
    any_n3000_clf = any(r["n3000_arm_computable_classification"] for r in arm_check)
    print(f"  Max regression train pool across all 5 fold-rotations: {max_reg_pool}")
    print(f"  Max classification train pool across all 5 fold-rotations: {max_clf_pool}")
    print(f"  N=3000 regression arm computable on ANY fold: {any_n3000_reg}")
    print(f"  N=3000 classification arm computable on ANY fold: {any_n3000_clf}")

    TEST_FOLD = 0  # match every other Gate 3 baseline for direct comparability
    max_n_reg = [r for r in arm_check if r["held_out_fold"] == TEST_FOLD][0]["regression_train_pool"]
    max_n_clf = [r for r in arm_check if r["held_out_fold"] == TEST_FOLD][0]["classification_train_pool"]

    print(f"\n=== Learning curve (TEST_FOLD={TEST_FOLD}, matching rest of Gate 3) ===")
    print(f"Regression: max available train pool = {max_n_reg}")
    reg_curve, reg_meta = learning_curve(d, kmer_features, TEST_FOLD, max_n_reg, "regression")
    reg_summary = summarize_curve(reg_curve, "spearman_rho")
    for N, s in reg_summary.items():
        print(f"    N={N:>5d}: rho={s['mean']}, std={s['std']}, n_draws={s['n_draws']}")

    print(f"Classification: max available train pool = {max_n_clf}")
    clf_curve, clf_meta = learning_curve(d, kmer_features, TEST_FOLD, max_n_clf, "classification")
    clf_summary_mcc = summarize_curve(clf_curve, "mcc")
    clf_summary_auc = summarize_curve(clf_curve, "auc")
    for N in clf_summary_mcc:
        print(f"    N={N:>5d}: mcc={clf_summary_mcc[N]['mean']}, auc={clf_summary_auc[N]['mean']}")

    verdict, verdict_reason = interpret_headroom(reg_summary, max_n_reg, floor_frac)
    print(f"\n=== Interpretation of the 24.5% figure ===")
    print(f"  VERDICT: {verdict}")
    print(f"  {verdict_reason}")

    print(f"\n=== PREMISE CORRECTION: was the N=3000 arm actually computable? ===")
    pc = premise_correction()
    print(f"  Gate 3's B2 train pool for BS translation: {pc['train_pool_used_by_gate3_b2']} "
          f"(= tl_usable, NOT tl_usable_for_regression)")
    print(f"  {pc['note']}")
    for N, r in pc["archived_gate3_numbers_by_nominal_n"].items():
        print(f"    nominal_N={N:>5s}: n_actual_train_drawn={r['n_actual_train_drawn']}, "
              f"mean_rho={r['mean_rho']}, std_rho={r['std_rho']}")
    print(f"  CONCLUSION: the N=3000 arm WAS computable and WAS already run (9,300-example pool >> 3000). "
          f"The '1,101 usable regression examples' figure describes a DIFFERENT quantity -- rows that turn "
          f"out non-floor/regression-informative -- not a cap on how many parts can be drawn/measured.")

    print(f"\n=== Reconciliation with Gate 3's original nominal-N sampling scheme ===")
    print(f"  (Gate 3's B2 drew N from the FULL tl_usable pool, mixing floor+active rows;")
    print(f"   this script's curve above draws N directly from tl_usable_for_regression only.)")
    reconciliation = reconcile_with_gate3_nominal_scheme(d, kmer_features, TEST_FOLD)
    for N, r in reconciliation.items():
        print(f"    nominal_N={N:>5d}: mean_effective_regression_n={r['mean_effective_regression_n']:.1f}, "
              f"mean_rho={r['mean_spearman_rho']}")

    # H-MAIN respecification options -- numbers only, no choice made
    bs_tx_total_usable = 15697  # from Gate 3 cache (tx_usable count, not floor-limited)
    print(f"\n=== H-MAIN respecification options (numbers only -- NOT chosen here) ===")
    print(f"  Option A (transcription primary): BS transcription tx_usable={bs_tx_total_usable} total, "
          f"train pool per fold ~{int(bs_tx_total_usable*4/5)} -- N=3000 arm fully computable, "
          f"no respecification needed for this readout.")
    print(f"  Option B (lower translation arm to data ceiling): max regression train pool={max_reg_pool} "
          f"-> use N={max_reg_pool} in place of N=3000; naive multiplier {max_reg_pool}/100="
          f"{max_reg_pool/100:.1f}x instead of 30x.")
    print(f"  Option C (dual-arm): N=3000 kept for transcription, N={max_reg_pool} (or nearest round "
          f"number below it) for translation, multipliers reported per-readout rather than one number.")

    results = {
        "per_fold_counts": fold_counts,
        "totals": {"tl_usable": total_usable, "tl_usable_for_regression": total_usable_reg,
                   "tl_active": total_active, "tl_inactive": total_inactive},
        "n3000_arm_check": arm_check,
        "max_regression_train_pool_any_fold": max_reg_pool,
        "max_classification_train_pool_any_fold": max_clf_pool,
        "n3000_regression_arm_computable_any_fold": any_n3000_reg,
        "n3000_classification_arm_computable_any_fold": any_n3000_clf,
        "learning_curve_test_fold": TEST_FOLD,
        "regression_curve": {"summary": reg_summary, "meta": reg_meta, "raw": reg_curve},
        "classification_curve": {"summary_mcc": clf_summary_mcc, "summary_auc": clf_summary_auc,
                                  "meta": clf_meta, "raw": clf_curve},
        "headroom_interpretation": {"verdict": verdict, "reason": verdict_reason},
        "reconciliation_with_gate3_nominal_scheme": reconciliation,
        "h_main_respecification_options": {
            "option_a_transcription_primary": {
                "bs_transcription_tx_usable_total": bs_tx_total_usable,
                "n3000_arm_computable": True,
                "description": "Make transcription H-MAIN's primary readout (BS has ample N there); "
                                "report translation as secondary with its N limitation stated plainly."
            },
            "option_b_lower_translation_arm": {
                "max_regression_train_pool_any_fold": max_reg_pool,
                "implied_multiplier_vs_n100": round(max_reg_pool / 100, 1),
                "description": "Replace the N>=3000 translation comparison arm with the data's actual "
                                "ceiling; report the resulting (smaller) data-efficiency multiplier honestly."
            },
            "option_c_dual_arm_per_readout": {
                "transcription_arm": 3000,
                "translation_arm": max_reg_pool,
                "description": "Keep N=3000 for transcription; use the translation ceiling for "
                                "translation; report two multipliers instead of one project-wide number."
            },
        },
    }
    with open(OUT / "gate3_5_bs_translation_feasibility.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nWrote {OUT / 'gate3_5_bs_translation_feasibility.json'}")


if __name__ == "__main__":
    main()
