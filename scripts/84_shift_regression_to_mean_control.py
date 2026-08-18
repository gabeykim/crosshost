"""
GATE 8.6 - Task 1: regression-to-the-mean control on Task 1B (LOAD-BEARING).

THE CONCERN: Task 1B (scripts/81) handed the model the reference host's
measured value as an input and asked it to predict the shift to the target
host. If shift correlates with the reference value itself (regression to
the mean -- a high-activity sequence in the reference host has more room to
fall than to rise), a model could predict shift from the reference value
ALONE, with zero sequence information, and look exactly like recovered
cross-host signal.

FOUR CHECKS, run for every (pair, readout, variant) cell -- reference-only
and raw-correlation for all 24 cells (cheap); shuffled-sequence for the 11
cells scripts/81 found DISTINGUISHABLY beat the mean-shift baseline only
(the ones actually at stake):

1. REFERENCE-ONLY BASELINE: predict shift from the reference host's
   measured value alone (ordinary least squares, 1 feature, fit on
   TRAIN-fold data, evaluated on TEST-fold data -- same 5 folds as
   scripts/81). No sequence input at all.

2. DECISIVE COMPARISON: does scripts/81's original sequence+reference model
   beat this reference-only model, distinguishably (90% CI, no overlap)?

3. RAW RELATIONSHIP: Spearman(shift, reference_value) directly on the full
   co-usable population per (pair, readout) -- descriptive, no fold split
   needed. If strongly negative, that IS the regression-to-the-mean
   mechanism and must be reported regardless of how (2) resolves.

4. SHUFFLED-SEQUENCE CONTROL: re-run scripts/81's exact model (same
   architecture, same conditioning variant, same folds) with the sequence
   array randomly permuted relative to (reference_value, shift, host_vec,
   fold) -- each row keeps its own true reference value, shift target, and
   fold assignment, but is handed a RANDOM OTHER row's sequence. This
   isolates how much of the original model's performance the ARCHITECTURE
   itself can extract from reference_value alone (any residual signal here
   is an upper bound on what sequence could possibly be contributing, since
   the sequence pathway now carries pure noise) and gives a second,
   empirical (not just analytic-OLS) reference-only-equivalent number.
"""
import sys
import time
import json
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
s81 = import_module("81_shift_prediction")
boot = import_module("51_bootstrap_utils")
m40 = import_module("40_film_cnn_model")

OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)
DEVICE = s81.DEVICE
N_FOLDS = 5

# The 11 (pair, readout, variant) cells scripts/81 found distinguishably beat
# the mean-shift baseline (fve_mean_lower90 > 0) -- read directly from the
# committed summary, not hand-copied, so this list can't drift from the data.
WINNING_CELLS_SOURCE = RESULTS / "gate8_5_shift_prediction_summary.csv"


def get_winning_cells():
    df = pd.read_csv(WINNING_CELLS_SOURCE)
    winners = df[df["fve_mean_lower90"] > 0]
    return [(r["pair"], r["readout"], r["variant"]) for _, r in winners.iterrows()]


def reference_only_ols(ref_tr, shift_tr, ref_te, shift_te):
    """Closed-form OLS, 1 feature (reference value) -> shift. Returns
    (pred, rho, mse)."""
    A = np.vstack([ref_tr, np.ones_like(ref_tr)]).T
    coef, *_ = np.linalg.lstsq(A, shift_tr, rcond=None)
    slope, intercept = coef
    pred = slope * ref_te + intercept
    mse = float(np.mean((pred - shift_te) ** 2))
    rho = spearmanr(shift_te, pred).correlation if len(set(np.round(shift_te, 6))) > 1 else None
    return pred, (float(rho) if rho is not None and rho == rho else None), mse, float(slope), float(intercept)


def run_reference_only_and_raw_corr(g_z):
    rows_ref = []
    rows_raw = []
    for ref, target in s81.ORDERED_PAIRS:
        pair_key = f"{ref}_to_{target}"
        for readout_name, raw_col, usable_col, active_col, transform in s81.READOUTS:
            seq, fold, ref_strength, shift, n_total = s81.load_pair_data(
                ref, target, usable_col, active_col, raw_col, transform)
            if n_total < 20:
                continue
            raw_rho = spearmanr(shift, ref_strength).correlation if len(set(np.round(shift, 6))) > 1 else None
            rows_raw.append({"pair": pair_key, "readout": readout_name, "n_total": n_total,
                              "spearman_shift_vs_reference": float(raw_rho) if raw_rho == raw_rho else None})

            per_fold = {}
            for test_fold in range(N_FOLDS):
                train_mask = fold != test_fold
                test_mask = fold == test_fold
                if test_mask.sum() < 5 or train_mask.sum() < 20:
                    continue
                ref_tr, shift_tr = ref_strength[train_mask], shift[train_mask]
                ref_te, shift_te = ref_strength[test_mask], shift[test_mask]
                pred, rho, mse, slope, intercept = reference_only_ols(ref_tr, shift_tr, ref_te, shift_te)
                mean_shift_train = float(shift_tr.mean())
                mse_zero = float(np.mean(shift_te ** 2))
                mse_mean = float(np.mean((shift_te - mean_shift_train) ** 2))
                per_fold[test_fold] = {
                    "spearman_rho": rho, "mse": mse, "slope": slope, "intercept": intercept,
                    "frac_var_explained_vs_zero": 1.0 - mse / mse_zero if mse_zero > 0 else None,
                    "frac_var_explained_vs_mean": 1.0 - mse / mse_mean if mse_mean > 0 else None,
                }
            for f_, r in per_fold.items():
                rows_ref.append({"pair": pair_key, "readout": readout_name, "fold": f_, **r})
            print(f"  reference-only {pair_key} {readout_name}: raw_corr(shift,ref)={raw_rho}, "
                  f"mean fold rho={np.mean([r['spearman_rho'] for r in per_fold.values() if r['spearman_rho'] is not None]):.3f}")

    df_ref = pd.DataFrame(rows_ref)
    df_raw = pd.DataFrame(rows_raw)
    df_ref.to_csv(RESULTS / "gate8_6_reference_only_baseline.csv", index=False)
    df_raw.to_csv(RESULTS / "gate8_6_shift_vs_reference_correlation.csv", index=False)
    return df_ref, df_raw


def run_shuffled_sequence_control(g_z, winning_cells):
    host_dim = len(g_z.columns)
    rows = []
    n_done, n_total = 0, sum(5 for _ in winning_cells)
    for pair, readout_name, variant in winning_cells:
        ref, target = pair.split("_to_")
        raw_col, usable_col, active_col, transform = next(
            (rc, uc, ac, tr) for rn, rc, uc, ac, tr in s81.READOUTS if rn == readout_name)
        seq, fold, ref_strength, shift, n_total_rows = s81.load_pair_data(
            ref, target, usable_col, active_col, raw_col, transform)
        target_hv = g_z.loc[target].values.astype(np.float32)
        use_cond = (variant == "with_condition")

        rng = np.random.default_rng(hash((pair, readout_name)) % (2**31))
        perm = rng.permutation(len(seq))
        seq_shuffled = seq[perm]  # each row now paired with a RANDOM OTHER row's sequence

        for test_fold in range(N_FOLDS):
            train_mask = fold != test_fold
            test_mask = fold == test_fold
            n_test = int(test_mask.sum())
            if n_test < 5 or train_mask.sum() < 20:
                continue
            seq_tr, ref_tr, shift_tr = seq_shuffled[train_mask], ref_strength[train_mask], shift[train_mask]
            seq_te, ref_te, shift_te = seq_shuffled[test_mask], ref_strength[test_mask], shift[test_mask]
            hv_tr = np.tile(target_hv, (len(seq_tr), 1)) if use_cond else None
            hv_te = np.tile(target_hv, (len(seq_te), 1)) if use_cond else None

            t0 = time.time()
            model, epochs_run = s81.train_shift_model(
                seq_tr, ref_tr, hv_tr, shift_tr, host_dim=host_dim if use_cond else 0,
                seed=test_fold, device=DEVICE)
            pred = s81.predict_shift(model, seq_te, ref_te, hv_te)
            mse_model = float(np.mean((pred - shift_te) ** 2))
            mean_shift_train = float(shift_tr.mean())
            mse_zero = float(np.mean(shift_te ** 2))
            mse_mean = float(np.mean((shift_te - mean_shift_train) ** 2))
            rho = spearmanr(shift_te, pred).correlation if len(set(np.round(shift_te, 6))) > 1 else None
            elapsed = time.time() - t0
            n_done += 1
            rows.append({
                "pair": pair, "readout": readout_name, "variant": variant, "fold": test_fold,
                "epochs_run": epochs_run, "wall_clock_sec": elapsed,
                "spearman_rho": float(rho) if rho is not None and rho == rho else None,
                "mse_model": mse_model,
                "frac_var_explained_vs_zero": 1.0 - mse_model / mse_zero if mse_zero > 0 else None,
                "frac_var_explained_vs_mean": 1.0 - mse_model / mse_mean if mse_mean > 0 else None,
            })
            print(f"[{n_done}/{n_total}] SHUFFLED {pair} {readout_name} {variant} fold={test_fold} "
                  f"({elapsed:.1f}s): rho={rho}, fve_zero={rows[-1]['frac_var_explained_vs_zero']}")
            with open(RESULTS / "gate8_6_shuffled_sequence_control.json", "w") as f:
                json.dump(rows, f, indent=2, default=str)

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "gate8_6_shuffled_sequence_control.csv", index=False)
    return df


def build_decisive_comparison(winning_cells):
    """Bootstrap-CI'd: original (seq+reference) vs reference-only OLS vs
    shuffled-sequence, per winning cell."""
    with open(OUT / "results" / "gate8_5_shift_prediction.json") as f:
        original = json.load(f)
    df_ref = pd.read_csv(RESULTS / "gate8_6_reference_only_baseline.csv")
    df_shuf = pd.read_csv(RESULTS / "gate8_6_shuffled_sequence_control.csv")

    rows = []
    for pair, readout_name, variant in winning_cells:
        orig_node = original[pair][readout_name]["per_fold"]
        orig_draws = [np.array([orig_node[str(f_)][variant]["spearman_rho"]])
                      for f_ in range(N_FOLDS) if str(f_) in orig_node]
        orig_fve_draws = [np.array([orig_node[str(f_)][variant]["frac_var_explained_vs_mean"]])
                           for f_ in range(N_FOLDS) if str(f_) in orig_node]

        ref_sub = df_ref[(df_ref.pair == pair) & (df_ref.readout == readout_name)]
        ref_draws = [np.array([v]) for v in ref_sub["spearman_rho"].dropna()]
        ref_fve_draws = [np.array([v]) for v in ref_sub["frac_var_explained_vs_mean"].dropna()]

        shuf_sub = df_shuf[(df_shuf.pair == pair) & (df_shuf.readout == readout_name) & (df_shuf.variant == variant)]
        shuf_draws = [np.array([v]) for v in shuf_sub["spearman_rho"].dropna()]
        shuf_fve_draws = [np.array([v]) for v in shuf_sub["frac_var_explained_vs_mean"].dropna()]

        ci_orig = boot.bootstrap_ci_90(orig_draws)
        ci_ref = boot.bootstrap_ci_90(ref_draws)
        ci_shuf = boot.bootstrap_ci_90(shuf_draws)
        ci_orig_fve = boot.bootstrap_ci_90(orig_fve_draws)
        ci_ref_fve = boot.bootstrap_ci_90(ref_fve_draws)
        ci_shuf_fve = boot.bootstrap_ci_90(shuf_fve_draws)

        beats_ref_distinguishably = (ci_orig["lower"] is not None and ci_ref["upper"] is not None
                                      and ci_orig["lower"] > ci_ref["upper"])
        effect_size_rho = (ci_orig["mean"] - ci_ref["mean"]) if ci_orig["mean"] is not None and ci_ref["mean"] is not None else None
        effect_size_fve = (ci_orig_fve["mean"] - ci_ref_fve["mean"]) if ci_orig_fve["mean"] is not None and ci_ref_fve["mean"] is not None else None

        rows.append({
            "pair": pair, "readout": readout_name, "variant": variant,
            "orig_rho_mean": ci_orig["mean"], "orig_rho_lower90": ci_orig["lower"], "orig_rho_upper90": ci_orig["upper"],
            "refonly_rho_mean": ci_ref["mean"], "refonly_rho_lower90": ci_ref["lower"], "refonly_rho_upper90": ci_ref["upper"],
            "shuffled_rho_mean": ci_shuf["mean"], "shuffled_rho_lower90": ci_shuf["lower"], "shuffled_rho_upper90": ci_shuf["upper"],
            "orig_fve_mean_mean": ci_orig_fve["mean"], "refonly_fve_mean_mean": ci_ref_fve["mean"], "shuffled_fve_mean_mean": ci_shuf_fve["mean"],
            "beats_refonly_distinguishably": bool(beats_ref_distinguishably),
            "effect_size_rho_orig_minus_refonly": effect_size_rho,
            "effect_size_fve_mean_orig_minus_refonly": effect_size_fve,
        })
        print(f"  {pair} {readout_name} {variant}: orig={ci_orig['mean']:.3f} [{ci_orig['lower']:.3f},{ci_orig['upper']:.3f}] "
              f"vs ref-only={ci_ref['mean']:.3f} [{ci_ref['lower']:.3f},{ci_ref['upper']:.3f}] "
              f"vs shuffled-seq={ci_shuf['mean']:.3f} [{ci_shuf['lower']:.3f},{ci_shuf['upper']:.3f}] "
              f"-- beats ref-only distinguishably: {beats_ref_distinguishably}")

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "gate8_6_shift_controls_decisive_comparison.csv", index=False)
    with open(RESULTS / "gate8_6_shift_controls.json", "w") as f:
        json.dump(rows, f, indent=2, default=str)
    return df


def main():
    g_z, gcols, p_z, pcols, imputed = m40.load_host_features()
    assert imputed == []

    print("=" * 80)
    print("STEP 1+3: reference-only OLS baseline + raw shift-vs-reference correlation (all 24 cells)")
    print("=" * 80)
    run_reference_only_and_raw_corr(g_z)

    winning_cells = get_winning_cells()
    print(f"\n{len(winning_cells)} winning cells (beat mean-shift distinguishably in scripts/81): {winning_cells}")

    print("\n" + "=" * 80)
    print(f"STEP 4: shuffled-sequence control ({len(winning_cells)} winning cells x 5 folds)")
    print("=" * 80)
    run_shuffled_sequence_control(g_z, winning_cells)

    print("\n" + "=" * 80)
    print("STEP 2+5: decisive comparison + effect size, per winning cell")
    print("=" * 80)
    build_decisive_comparison(winning_cells)

    print("\n=== TASK 1 CONTROLS COMPLETE ===")


if __name__ == "__main__":
    main()
