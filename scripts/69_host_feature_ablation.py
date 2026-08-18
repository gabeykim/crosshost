"""
GATE 7 - Task 3: bounded host-feature ablation, genomic variant, ONE
representative configuration (B. subtilis held out -- the charter's official
H-MAIN stringent test host), both readouts, all 5 folds (cheap: no
retraining, see method below). One working session, not a full ablation
program, per the task's explicit scope instruction.

METHOD -- inference-time feature-group knockout (mean-ablation), not
retraining: the already-trained, already-saved LOHO checkpoints
(out/models/loho_BS_genomic_fold{0-4}.pt) are reused unchanged. Each
feature group's z-scored columns are set to 0 (the population mean, since
these are z-scored across all 6 hosts -- scripts/40's load_host_features)
in the host vector fed to the model at prediction time, and the resulting
drop in Spearman rho vs the unablated baseline is reported. This is a
standard, well-established lightweight ablation method (feature knockout)
that avoids the cost of retraining 25 additional models (5 groups x 5
folds) from scratch, appropriate given the task's explicit "keep this
tightly bounded" instruction.

FEATURE GROUPS, mapped from data/hosts_genomic.parquet's 37 columns to the
5 named groups in the task text (5 columns are NOT part of any named group
and are always retained: genome_size_bp, gc_content, n_16s_rrna_copies,
n_trna_genes_total, n_ribosomal_proteins -- general genome-architecture
features the task's 5-group list does not name):
  sigma-factor:    n_sigma_factors_total, n_sigma70_primary, n_sigma_ecf, n_sigma54
  anti-SD:         anti_sd_rnaduplex_mfe_kcalmol
  tAI/codon:       tai_genome_wide_mean, gc3_content, enc_nc, trna_n_* (21 columns)
  RNAP:            n_rnap_core_subunits
  chaperone/heme:  n_chaperones, n_heme_biosynthesis_genes
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
loho = import_module("42_loho_training")

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"
MODELS_DIR = OUT / "models"

HELD_OUT = "BS"
N_FOLDS = 5

FEATURE_GROUPS = {
    "sigma_factor": ["n_sigma_factors_total", "n_sigma70_primary", "n_sigma_ecf", "n_sigma54"],
    "anti_SD": ["anti_sd_rnaduplex_mfe_kcalmol"],
    "tAI_codon": ["tai_genome_wide_mean", "gc3_content", "enc_nc"] + [
        f"trna_n_{aa}" for aa in ["Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
                                    "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "fMet"]],
    "RNAP": ["n_rnap_core_subunits"],
    "chaperone_heme": ["n_chaperones", "n_heme_biosynthesis_genes"],
}


def main():
    g_z, gcols, p_z, pcols, imputed = m40.load_host_features()
    for grp, cols in FEATURE_GROUPS.items():
        missing = [c for c in cols if c not in gcols]
        assert not missing, f"columns not found in genomic feature set: {missing}"
    col_idx = {c: i for i, c in enumerate(gcols)}
    print(f"37 genomic columns; {sum(len(v) for v in FEATURE_GROUPS.values())} covered by named groups, "
          f"{37 - sum(len(v) for v in FEATURE_GROUPS.values())} always-retained (ungrouped)")

    results = {"baseline": {}, "ablated": {g: {} for g in FEATURE_GROUPS}}
    for readout in ["transcription", "translation"]:
        results["baseline"][readout] = []
        for g in FEATURE_GROUPS:
            results["ablated"][g][readout] = []

    for f in range(N_FOLDS):
        seq_te, hv_te, targets_te, masks_te, _ = m40.build_pooled_arrays(
            [HELD_OUT], g_z, fold_filter_fn=lambda fold_arr, tf=f: fold_arr == tf)
        model = m40.FiLMSequenceCNN(host_dim=len(gcols), use_film=True)
        model.load_state_dict(torch.load(MODELS_DIR / f"loho_{HELD_OUT}_genomic_fold{f}.pt", map_location="cpu"))

        preds_base = m40.predict_film(model, seq_te, hv_te)
        eval_base = loho.evaluate(preds_base, targets_te, masks_te)
        for readout in ["transcription", "translation"]:
            results["baseline"][readout].append(eval_base[readout]["spearman_rho"])

        for grp, cols in FEATURE_GROUPS.items():
            hv_ablated = hv_te.copy()
            for c in cols:
                hv_ablated[:, col_idx[c]] = 0.0  # z-scored -> 0 is the population mean
            preds_ab = m40.predict_film(model, seq_te, hv_ablated)
            eval_ab = loho.evaluate(preds_ab, targets_te, masks_te)
            for readout in ["transcription", "translation"]:
                results["ablated"][grp][readout].append(eval_ab[readout]["spearman_rho"])

        print(f"fold {f}: baseline tx={eval_base['transcription']['spearman_rho']:.3f} "
              f"tl={eval_base['translation']['spearman_rho']:.3f}")

    print("\n" + "=" * 78)
    print(f"FEATURE-GROUP ABLATION, genomic variant, {HELD_OUT} held out, 5 folds")
    print("=" * 78)
    summary = {"baseline": {}, "ablation_delta": {}}
    for readout in ["transcription", "translation"]:
        base_vals = [v for v in results["baseline"][readout] if v is not None]
        base_mean = float(np.mean(base_vals))
        summary["baseline"][readout] = {"mean": base_mean, "std": float(np.std(base_vals)), "per_fold": base_vals}
        print(f"\n{readout}: baseline rho = {base_mean:.4f} (std={np.std(base_vals):.4f})")
        summary["ablation_delta"][readout] = {}
        for grp in FEATURE_GROUPS:
            ab_vals = [v for v in results["ablated"][grp][readout] if v is not None]
            ab_mean = float(np.mean(ab_vals))
            delta = ab_mean - base_mean
            summary["ablation_delta"][readout][grp] = {"ablated_mean": ab_mean, "delta": delta,
                                                         "delta_within_1_baseline_std": abs(delta) < np.std(base_vals)}
            flag = "" if abs(delta) >= np.std(base_vals) else "  (within 1 baseline-fold-std -- likely noise)"
            print(f"  drop {grp:16s}: rho={ab_mean:.4f}  delta={delta:+.4f}{flag}")

    all_deltas = [abs(summary["ablation_delta"][r][g]["delta"]) for r in ["transcription", "translation"] for g in FEATURE_GROUPS]
    all_within_noise = all(summary["ablation_delta"][r][g]["delta_within_1_baseline_std"]
                            for r in ["transcription", "translation"] for g in FEATURE_GROUPS)
    verdict = ("NO FEATURE GROUP CARRIES DETECTABLE SIGNAL -- every ablation's rho change is within 1 "
               "baseline fold-to-fold standard deviation, indistinguishable from noise" if all_within_noise else
               "AT LEAST ONE FEATURE GROUP SHOWS A CHANGE EXCEEDING FOLD-TO-FOLD NOISE -- see per-group deltas")
    summary["verdict"] = verdict
    summary["max_abs_delta"] = float(max(all_deltas))
    print(f"\nVERDICT: {verdict}")
    print(f"Max |delta| across all groups/readouts: {max(all_deltas):.4f}")

    with open(RESULTS / "gate7_host_feature_ablation.json", "w") as f:
        json.dump({"raw": results, "summary": summary}, f, indent=2, default=str)
    print(f"\nWrote {RESULTS / 'gate7_host_feature_ablation.json'}")


if __name__ == "__main__":
    main()
