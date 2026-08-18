"""
GATE 8.5 - Task 1C: restrict to the co-active subset.

PART 1 -- pooled vs. co-active-restricted cross-host measurement correlation.

IMPORTANT CORRECTION TO THE TASK'S OWN PREMISE, found while implementing
this: the task asks to compute "cross-host measurement correlation
restricted to sequences active in both hosts... alongside the pooled
figures" as if this were new. It is NOT new for the "restricted" half --
scripts/57_crosshost_measurement_correlation.py (Gate 5.5), which produced
the ρ≈0.75 (EC-PA) / ρ≈0.16-0.26 (BS-pairs) numbers that PAPER_FRAMING.md's
finding #1 leads with, ALREADY restricts to `{h1}_tx_usable & {h1}_tx_active
& {h2}_tx_usable & {h2}_tx_active` -- i.e. active in BOTH hosts. That
number IS the co-active-restricted figure, not a pooled one. This script
reproduces it here (verified bit-identical below) specifically so it can
be reported side-by-side with the genuinely-missing comparison: the POOLED
figure (usable in both hosts, active status ignored, so zero/inactive rows
are included). That pooled-vs-restricted contrast is the real new
information this task adds.

PART 2 -- does predictability improve on the co-active subset? Reuses the
ALREADY-TRAINED, already-saved LOHO checkpoints for sequence-only
(out/models/seqonly_loho_{host}_fold{f}.pt) and genomic FiLM
(out/models/loho_{host}_genomic_fold{f}.pt, the stronger of the two
host-feature variants per H-SCIENCE and Task 1A's own choice) -- NO
retraining. For each (held-out host, training-partner host) pair, the
held-out host's test-fold rows are restricted to those ALSO active in the
training-partner host (using the same row alignment guaranteed by
scripts/40's baseline_cache convention: all 3 hosts' .npz caches share
identical oligo_ids/fold ordering, verified inline below), and the SAME
saved model is re-scored on that restricted subset. This is compared
against the model's already-known full-test-set number.

HONESTY, per the task's explicit instruction: this is NOT a fair predictor
comparison. Knowing a sequence is active in a TRAINING host is information
a real held-out-host deployment would not have in advance (you would not
yet know which of your new host's sequences are "co-active" with anything
without already having measured them, which defeats the point of
predicting). Any rho improvement here is a statement about WHERE cross-host
conservation lives in the data, not an improved predictor -- exactly the
framing the task requires and PAPER_FRAMING.md must preserve.
"""
import sys
import json
import numpy as np
import pandas as pd
import torch
from pathlib import Path
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
m40 = import_module("40_film_cnn_model")
so = import_module("58_sequence_only_model")
loho = import_module("42_loho_training")
boot = import_module("51_bootstrap_utils")

DATA = Path(__file__).resolve().parent.parent / "data"
CACHE = DATA / "baseline_cache"
OUT = Path(__file__).resolve().parent.parent / "out"
MODELS_DIR = OUT / "models"
RESULTS = OUT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

HOSTS = ["EC", "BS", "PA"]
PAIRS = [("EC", "BS"), ("EC", "PA"), ("BS", "PA")]
N_FOLDS = 5


# ---------------------------------------------------------------------------
# PART 1
# ---------------------------------------------------------------------------
def part1_pooled_vs_coactive():
    df = pd.read_parquet(DATA / "three_host.parquet")

    tl_active_col, tl_usable_reg_col = {}, {}
    for host in HOSTS:
        usable = df[f"{host}_tl_usable"].values.astype(bool)
        protein = df[f"{host}_protein_log10"].values.astype(np.float32)
        vals = protein[usable]
        rounded = np.round(vals, 2)
        u, c = np.unique(rounded, return_counts=True)
        mode_val, mode_frac = u[np.argmax(c)], c[np.argmax(c)] / len(vals)
        is_floor = mode_frac > 0.05
        if is_floor:
            at_floor = usable & (np.round(protein, 2) == mode_val)
            tl_usable_reg_col[host] = usable & ~at_floor
            tl_active_col[host] = usable & (np.round(protein, 2) > mode_val)
        else:
            tl_usable_reg_col[host] = usable.copy()
            tl_active_col[host] = usable & (protein > np.median(vals))

    rows = []
    print("=" * 80)
    print("PART 1: pooled (usable-in-both) vs. co-active-restricted (active-in-both)")
    print("=" * 80)
    for h1, h2 in PAIRS:
        # TRANSCRIPTION
        usable_both = (df[f"{h1}_tx_usable"] & df[f"{h2}_tx_usable"]).values
        active_both = (df[f"{h1}_tx_usable"] & df[f"{h1}_tx_active"] &
                        df[f"{h2}_tx_usable"] & df[f"{h2}_tx_active"]).values
        n_pooled, n_active = int(usable_both.sum()), int(active_both.sum())
        rho_pooled = spearmanr(df.loc[usable_both, f"{h1}_tx_norm"], df.loc[usable_both, f"{h2}_tx_norm"]).correlation if n_pooled > 1 else None
        rho_active = spearmanr(df.loc[active_both, f"{h1}_tx_norm"], df.loc[active_both, f"{h2}_tx_norm"]).correlation if n_active > 1 else None
        rows.append({"pair": f"{h1}_{h2}", "readout": "transcription", "n_pooled": n_pooled, "rho_pooled": rho_pooled,
                      "n_coactive": n_active, "rho_coactive": rho_active})
        print(f"  tx {h1}-{h2}: pooled n={n_pooled} rho={rho_pooled:.4f} | co-active n={n_active} rho={rho_active:.4f}")

        # TRANSLATION
        usable_both_tl = (tl_usable_reg_col[h1] & tl_usable_reg_col[h2])
        active_both_tl = (tl_active_col[h1] & tl_active_col[h2] & tl_usable_reg_col[h1] & tl_usable_reg_col[h2])
        n_pooled_tl, n_active_tl = int(usable_both_tl.sum()), int(active_both_tl.sum())
        rho_pooled_tl = spearmanr(df.loc[usable_both_tl, f"{h1}_protein_log10"], df.loc[usable_both_tl, f"{h2}_protein_log10"]).correlation if n_pooled_tl > 1 else None
        rho_active_tl = spearmanr(df.loc[active_both_tl, f"{h1}_protein_log10"], df.loc[active_both_tl, f"{h2}_protein_log10"]).correlation if n_active_tl > 1 else None
        rows.append({"pair": f"{h1}_{h2}", "readout": "translation", "n_pooled": n_pooled_tl, "rho_pooled": rho_pooled_tl,
                      "n_coactive": n_active_tl, "rho_coactive": rho_active_tl})
        print(f"  tl {h1}-{h2}: pooled n={n_pooled_tl} rho={rho_pooled_tl:.4f} | co-active n={n_active_tl} rho={rho_active_tl:.4f}")

    df_out = pd.DataFrame(rows)
    df_out.to_csv(RESULTS / "gate8_5_coactive_correlation.csv", index=False)

    # sanity check: verify the "active_both" numbers reproduce script 57's committed output bit-for-bit
    with open(RESULTS / "gate5_5_crosshost_measurement_correlation.json") as f:
        prior = json.load(f)
    mismatches = []
    for _, r in df_out.iterrows():
        prior_val = prior[r["readout"]][r["pair"]]
        if prior_val["n"] != r["n_coactive"] or (prior_val["rho"] is not None and r["rho_coactive"] is not None
                                                   and abs(prior_val["rho"] - r["rho_coactive"]) > 1e-9):
            mismatches.append({"pair": r["pair"], "readout": r["readout"], "prior": prior_val,
                                "recomputed": {"n": r["n_coactive"], "rho": r["rho_coactive"]}})
    print(f"\nBit-identical check against scripts/57's committed output: "
          f"{'PASS, 0 mismatches' if not mismatches else f'{len(mismatches)} MISMATCHES: ' + str(mismatches)}")

    with open(RESULTS / "gate8_5_coactive_correlation_sanity_check.json", "w") as f:
        json.dump({"mismatches": mismatches, "n_checked": len(df_out)}, f, indent=2)

    print(f"Wrote {RESULTS / 'gate8_5_coactive_correlation.csv'}")
    return df_out


# ---------------------------------------------------------------------------
# PART 2
# ---------------------------------------------------------------------------
def load_host_active_masks():
    masks = {}
    for h in HOSTS:
        d = np.load(CACHE / f"{h}_baseline_data.npz", allow_pickle=True)
        masks[h] = {"tx_active": d["tx_active"].astype(bool), "tl_active": d["tl_active"].astype(bool),
                     "tx_usable": d["tx_usable"].astype(bool),
                     "tl_usable_for_regression": d["tl_usable_for_regression"].astype(bool),
                     "fold": d["fold"]}
    # verify row alignment across hosts, as documented
    fold_ec = masks["EC"]["fold"]
    for h in ["BS", "PA"]:
        assert (masks[h]["fold"] == fold_ec).all(), f"row misalignment detected for {h}"
    return masks


def load_seqonly_model(held_out, fold):
    model = so.SequenceOnlyCNN()
    state = torch.load(MODELS_DIR / f"seqonly_loho_{held_out}_fold{fold}.pt", map_location="cpu")
    model.load_state_dict(state)
    return model


def load_film_genomic_model(held_out, fold, host_dim):
    model = m40.FiLMSequenceCNN(host_dim=host_dim, use_film=True)
    state = torch.load(MODELS_DIR / f"loho_{held_out}_genomic_fold{fold}.pt", map_location="cpu")
    model.load_state_dict(state)
    return model


def part2_rerun_on_coactive(active_masks, g_z, gcols):
    host_dim = len(gcols)
    rows = []
    print("\n" + "=" * 80)
    print("PART 2: re-scoring EXISTING saved LOHO checkpoints on the co-active-restricted test subset")
    print("=" * 80)

    for held_out in HOSTS:
        train_hosts = [h for h in HOSTS if h != held_out]
        for partner in train_hosts:
            for readout, active_key, strength_key in [("transcription", "tx_active", "tx_strength"),
                                                         ("translation", "tl_active", "tl_strength")]:
                usable_key = "tx_usable" if readout == "transcription" else "tl_usable_for_regression"
                for test_fold in range(N_FOLDS):
                    fold_mask = active_masks[held_out]["fold"] == test_fold
                    idx_full = np.where(fold_mask)[0]

                    # full test set (already-known numbers, recomputed here for a clean paired comparison)
                    seq_full, hv_full, targets_full, masks_full, _ = m40.build_pooled_arrays(
                        [held_out], g_z, fold_filter_fn=lambda f, tf=test_fold: f == tf)

                    # co-active-restricted: held-out host active/usable AND training-partner active/usable, same row idx
                    coactive_mask_within_fold = (
                        active_masks[held_out][active_key][idx_full] & active_masks[held_out][usable_key][idx_full] &
                        active_masks[partner][active_key][idx_full] & active_masks[partner][usable_key][idx_full]
                    )
                    n_coactive = int(coactive_mask_within_fold.sum())

                    seqonly_model = load_seqonly_model(held_out, test_fold)
                    film_model = load_film_genomic_model(held_out, test_fold, host_dim)

                    preds_seq_full = so.predict_seqonly(seqonly_model, seq_full)
                    preds_film_full = m40.predict_film(film_model, seq_full, hv_full)
                    eval_full = {"sequence_only": loho.evaluate(preds_seq_full, targets_full, masks_full)[readout],
                                  "film_genomic": loho.evaluate(preds_film_full, targets_full, masks_full)[readout]}

                    if n_coactive > 1:
                        seq_ca = seq_full[coactive_mask_within_fold]
                        hv_ca = hv_full[coactive_mask_within_fold]
                        targets_ca = {k: v[coactive_mask_within_fold] for k, v in targets_full.items()}
                        masks_ca = {k: v[coactive_mask_within_fold] for k, v in masks_full.items()}
                        preds_seq_ca = so.predict_seqonly(seqonly_model, seq_ca)
                        preds_film_ca = m40.predict_film(film_model, seq_ca, hv_ca)
                        eval_ca = {"sequence_only": loho.evaluate(preds_seq_ca, targets_ca, masks_ca)[readout],
                                   "film_genomic": loho.evaluate(preds_film_ca, targets_ca, masks_ca)[readout]}
                    else:
                        eval_ca = {"sequence_only": {"spearman_rho": None, "n_strength_eval": 0},
                                   "film_genomic": {"spearman_rho": None, "n_strength_eval": 0}}

                    for system in ["sequence_only", "film_genomic"]:
                        rows.append({
                            "held_out": held_out, "training_partner": partner, "readout": readout, "fold": test_fold,
                            "system": system, "n_coactive": n_coactive,
                            "rho_full_test": eval_full[system]["spearman_rho"],
                            "rho_coactive_subset": eval_ca[system]["spearman_rho"],
                        })
                    print(f"  held_out={held_out} partner={partner} {readout} fold={test_fold}: "
                          f"n_coactive={n_coactive}, seqonly rho full={eval_full['sequence_only']['spearman_rho']} "
                          f"-> coactive={eval_ca['sequence_only']['spearman_rho']}; "
                          f"film rho full={eval_full['film_genomic']['spearman_rho']} "
                          f"-> coactive={eval_ca['film_genomic']['spearman_rho']}")

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "gate8_5_coactive_model_rerun.csv", index=False)

    # bootstrap CI summary per (held_out, partner, readout, system): full vs coactive
    summary_rows = []
    for (held_out, partner, readout, system), grp in df.groupby(["held_out", "training_partner", "readout", "system"]):
        full_draws = [np.array([v]) if v is not None and v == v else np.array([]) for v in grp["rho_full_test"]]
        ca_draws = [np.array([v]) if v is not None and v == v else np.array([]) for v in grp["rho_coactive_subset"]]
        ci_full = boot.bootstrap_ci_90(full_draws)
        ci_ca = boot.bootstrap_ci_90(ca_draws)
        summary_rows.append({
            "held_out": held_out, "training_partner": partner, "readout": readout, "system": system,
            "n_coactive_mean": float(grp["n_coactive"].mean()),
            "rho_full_mean": ci_full["mean"], "rho_full_lower90": ci_full["lower"], "rho_full_upper90": ci_full["upper"],
            "rho_coactive_mean": ci_ca["mean"], "rho_coactive_lower90": ci_ca["lower"], "rho_coactive_upper90": ci_ca["upper"],
        })
    sdf = pd.DataFrame(summary_rows)
    sdf.to_csv(RESULTS / "gate8_5_coactive_model_rerun_summary.csv", index=False)
    print(f"\nWrote {RESULTS / 'gate8_5_coactive_model_rerun.csv'}, {RESULTS / 'gate8_5_coactive_model_rerun_summary.csv'}")
    return df, sdf


def main():
    df1 = part1_pooled_vs_coactive()
    material_diff = (df1["rho_coactive"] - df1["rho_pooled"]).abs().max()
    print(f"\nMax |rho_coactive - rho_pooled| across all pairs/readouts: {material_diff:.4f}")

    g_z, gcols, p_z, pcols, imputed = m40.load_host_features()
    assert imputed == []
    active_masks = load_host_active_masks()
    df2, sdf2 = part2_rerun_on_coactive(active_masks, g_z, gcols)

    print("\n=== TASK 1C COMPLETE ===")


if __name__ == "__main__":
    main()
