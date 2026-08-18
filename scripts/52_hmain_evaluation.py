"""
GATE 5 - Task 1: evaluate H-MAIN, both arms, per out/PREREGISTRATION.md
(including both dated amendments). Applies the pre-registered rules
mechanically -- no threshold/interval/mechanism is adjusted here.
"""
import sys
import json
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
boot = import_module("51_bootstrap_utils")

OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

VARIANTS = ["genomic", "physiology"]
MECHANISMS = {"head_only": "frozen-trunk (PRE-REGISTERED PRIMARY)", "top_conv": "unfrozen-conv4 (SUPPLEMENTARY)"}


def get_baseline_fold_draws(readout, N, host="BS"):
    """k-mer/linear baseline, from Gate 3.5's fold-rotated raw data."""
    with open(OUT / "baselines" / "fold_rotation_b2_raw.json") as f:
        d = json.load(f)
    entry = d[host][readout][str(N)]
    return [np.array([x["spearman_rho"] for x in entry[str(f)]["draws"] if x["spearman_rho"] is not None]) for f in range(5)]


def get_tl_ceiling_fold_draws():
    with open(OUT / "gate5_bs_translation_ceiling_raw.json") as f:
        d = json.load(f)
    out = []
    for f in range(5):
        curve = d[str(f)]["curve_raw"]
        draws = curve.get("300", [])
        out.append(np.array([x["spearman_rho"] for x in draws if x["spearman_rho"] is not None]))
    return out


def get_model_fold_draws(variant, readout, mechanism, host="BS"):
    if mechanism == "head_only":
        with open(OUT / "gate4_calibration_curves.json") as f:
            cal = json.load(f)
        out = []
        for f in range(5):
            draws = cal[host][variant][str(f)]["N100_head_only"]
            out.append(np.array([x["eval"][readout]["spearman_rho"] for x in draws if x["eval"][readout]["spearman_rho"] is not None]))
        return out
    else:  # top_conv
        with open(OUT / "gate4_5_n100_topconv_supplement.json") as f:
            supp = json.load(f)
        out = []
        for f in range(5):
            draws = supp[variant][host][str(f)]
            out.append(np.array([x["eval"][readout]["spearman_rho"] for x in draws if x["eval"][readout]["spearman_rho"] is not None]))
        return out


def multiplier_tx(model_point, baseline_curve_folds_by_n):
    """What N does the baseline need to match the model's point estimate?
    baseline_curve_folds_by_n: dict N -> fold_draws (list of arrays)."""
    ns = sorted(baseline_curve_folds_by_n.keys())
    for n in ns:
        fd = baseline_curve_folds_by_n[n]
        pt = np.mean([d.mean() for d in fd if len(d) > 0])
        if pt >= model_point:
            return n, pt
    return None, None  # baseline never reaches it at any tested N


def main():
    results = {"H_MAIN_TX": {}, "H_MAIN_TL": {}}

    print("=" * 70)
    print("H-MAIN-TX (transcription): model@N=100 vs baseline@N=3000")
    print("=" * 70)
    baseline_100_folds = get_baseline_fold_draws("transcription", 100)
    baseline_3000_folds = get_baseline_fold_draws("transcription", 3000)
    baseline_100_ci = boot.bootstrap_ci_90(baseline_100_folds)
    baseline_3000_ci = boot.bootstrap_ci_90(baseline_3000_folds)
    print(f"baseline N=100:  {baseline_100_ci}")
    print(f"baseline N=3000: {baseline_3000_ci}")

    # full baseline curve for multiplier computation
    baseline_curve_by_n = {n: get_baseline_fold_draws("transcription", n) for n in [0, 10, 30, 100, 300, 1000, 3000]}

    results["H_MAIN_TX"]["baseline"] = {"N100": baseline_100_ci, "N3000": baseline_3000_ci}
    results["H_MAIN_TX"]["model"] = {}
    for variant in VARIANTS:
        results["H_MAIN_TX"]["model"][variant] = {}
        for mech in MECHANISMS:
            model_folds = get_model_fold_draws(variant, "transcription", mech)
            model_ci = boot.bootstrap_ci_90(model_folds)
            v = boot.verdict(model_ci, baseline_3000_ci)
            ov = boot.overlap_size(model_ci, baseline_3000_ci)
            mult_n, mult_pt = multiplier_tx(model_ci["mean"], baseline_curve_by_n)
            print(f"  variant={variant} mechanism={MECHANISMS[mech]}: model={model_ci} "
                  f"verdict={v} overlap={ov} multiplier_N={mult_n}")
            results["H_MAIN_TX"]["model"][variant][mech] = {
                "mechanism_label": MECHANISMS[mech], "model_ci": model_ci,
                "verdict": v, "overlap_with_baseline_upper": ov,
                "data_efficiency_multiplier_baseline_N_to_match": mult_n,
                "baseline_value_at_multiplier_N": mult_pt,
            }

    print()
    print("=" * 70)
    print("H-MAIN-TL (translation): model@N=100 vs baseline SATURATION CEILING (N~300)")
    print("=" * 70)
    tl_ceiling_folds = get_tl_ceiling_fold_draws()
    tl_ceiling_ci = boot.bootstrap_ci_90(tl_ceiling_folds)
    print(f"baseline ceiling (N=300, direct-regression-only draws): {tl_ceiling_ci}")
    results["H_MAIN_TL"]["baseline_ceiling"] = tl_ceiling_ci
    results["H_MAIN_TL"]["baseline_ceiling_note"] = (
        "Ceiling the per-host baseline CANNOT cross at any available N (flat from N=300 to max pool, "
        "Gate 3.5/Gate 4 finding). No 'multiplier N' exists in the usual sense -- if the model beats this, "
        "the baseline cannot match it at ANY N, reported as such rather than a finite multiplier.")

    results["H_MAIN_TL"]["model"] = {}
    for variant in VARIANTS:
        results["H_MAIN_TL"]["model"][variant] = {}
        for mech in MECHANISMS:
            model_folds = get_model_fold_draws(variant, "translation", mech)
            model_ci = boot.bootstrap_ci_90(model_folds)
            v = boot.verdict(model_ci, tl_ceiling_ci)
            ov = boot.overlap_size(model_ci, tl_ceiling_ci)
            print(f"  variant={variant} mechanism={MECHANISMS[mech]}: model={model_ci} verdict={v} overlap={ov}")
            results["H_MAIN_TL"]["model"][variant][mech] = {
                "mechanism_label": MECHANISMS[mech], "model_ci": model_ci,
                "verdict": v, "overlap_with_baseline_upper": ov,
                "data_efficiency_multiplier": ("N/A -- baseline ceiling, model exceeded it at N=100" if v == "MET"
                                                 else "N/A -- baseline ceiling not exceeded"),
            }

    print()
    print("=" * 70)
    print("CONTEXT (charter Part V / Gate 5 rule 5): P. aeruginosa and E. coli, host-level")
    print("(same N=100-vs-N=3000 structure as H-MAIN-TX for BOTH readouts -- neither PA nor EC")
    print(" has a documented translation saturation-ceiling issue like BS's, so no special TL")
    print(" ceiling arm is defined for them; this is host-level CONTEXT for the charter's PARTIAL")
    print(" PASS clause, not a second pre-registered hypothesis.)")
    print("=" * 70)
    results["host_level_context_PA_EC"] = {}
    for host in ["PA", "EC"]:
        results["host_level_context_PA_EC"][host] = {}
        for readout in ["transcription", "translation"]:
            baseline_3000_h = boot.bootstrap_ci_90(get_baseline_fold_draws(readout, 3000, host=host))
            results["host_level_context_PA_EC"][host][readout] = {"baseline_N3000": baseline_3000_h, "model": {}}
            for variant in VARIANTS:
                results["host_level_context_PA_EC"][host][readout]["model"][variant] = {}
                for mech in MECHANISMS:
                    model_folds_h = get_model_fold_draws(variant, readout, mech, host=host)
                    model_ci_h = boot.bootstrap_ci_90(model_folds_h)
                    v_h = boot.verdict(model_ci_h, baseline_3000_h)
                    print(f"  {host} {readout} variant={variant} mech={MECHANISMS[mech]}: "
                          f"model_mean={model_ci_h['mean']:.3f} baseline_N3000_mean={baseline_3000_h['mean']:.3f} verdict={v_h}")
                    results["host_level_context_PA_EC"][host][readout]["model"][variant][mech] = {
                        "model_ci": model_ci_h, "verdict": v_h,
                    }

    with open(RESULTS / "gate5_hmain_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nWrote {RESULTS / 'gate5_hmain_results.json'}")


if __name__ == "__main__":
    main()
