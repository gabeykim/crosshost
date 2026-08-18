"""
GATE 8 - Task 1: measurement-reliability attenuation check. Does the
EC-PA (rho~0.75) vs. BS-pairs (rho~0.16-0.26) cross-host measurement
correlation gap survive correction for measurement unreliability, or does
it substantially close?

===========================================================================
REPLICATE STRUCTURE AVAILABLE -- and what is NOT available, stated first
===========================================================================
Gate 1 already established [VERIFIED] that the manuscript-described
barcode-replicate QC set (4,778 RSs, technical replicates) is NOT
reconstructable from the released supplementary tables (CHARTER_AMENDMENTS.md
C1). That path is closed for all three hosts.

What IS available: raw/NIHMS945382-supplement-5.xlsx, "Supplementary Data
Table 3", contains E. coli transcription measured across FIVE growth
conditions (LB_exp, NaCl_exp, Fe_exp, LB-stat, M9-exp) for the SAME 29,249
sequences, same RNA_count/DNA_count/tx_norm structure as the main library.
This is EC-TRANSCRIPTION-ONLY -- no equivalent condition-series exists in
the released tables for B. subtilis, P. aeruginosa, or for translation in
any host (verified: the "Robust RSs" sheet's "Protein (a.u.)" column is a
single value per sequence, not a repeated-measures series).

METHODOLOGICAL CAVEAT, stated before any number is computed: growth
CONDITIONS are not technical replicates. Some of the sequence-to-sequence
variation in tx_norm across LB/NaCl/Fe/LB-stat/M9 is genuine condition-
dependent regulation (real biology), not measurement noise. Treating
cross-condition correlation as a reliability estimate therefore BUNDLES
true condition-effects into what classical test theory calls "noise" --
this UNDERESTIMATES true single-condition measurement reliability, which
means any disattenuation correction computed from it is an UPPER BOUND on
how much of the EC-PA/BS gap could plausibly be explained by measurement
error alone. Reported as such throughout, not as a precise point estimate.
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr
from itertools import combinations

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

RAW = Path(__file__).resolve().parent.parent / "raw"
CONDITIONS = ["LB_exp", "NaCl_exp", "Fe_exp", "LB-stat", "M9-exp"]
SENSITIVITY_RELIABILITIES = [0.5, 0.7, 0.9]


def load_condition_sheet(sheet):
    df = pd.read_excel(RAW / "NIHMS945382-supplement-5.xlsx", sheet_name=sheet)
    dna = df["DNA_count"]
    rna = df["RNA_count"]
    tx_usable = (dna != 0) & ((rna + dna) >= 15)   # identical to scripts/13's tx_usable definition
    tx_active = tx_usable & (rna > 0)               # identical to scripts/13's tx_active definition
    df["tx_usable"] = tx_usable
    df["tx_active"] = tx_active
    return df[["OLIGO ID", "tx_norm", "tx_usable", "tx_active"]]


def ec_cross_condition_reliability():
    print("=" * 78)
    print("EC cross-condition transcription correlation (5 growth conditions, Table 3)")
    print("=" * 78)
    sheets = {c: load_condition_sheet(c) for c in CONDITIONS}
    pairwise = {}
    for c1, c2 in combinations(CONDITIONS, 2):
        d1, d2 = sheets[c1], sheets[c2]
        merged = d1.merge(d2, on="OLIGO ID", suffixes=("_1", "_2"))
        both_active = merged["tx_active_1"] & merged["tx_active_2"]
        n = int(both_active.sum())
        v1 = merged.loc[both_active, "tx_norm_1"].values
        v2 = merged.loc[both_active, "tx_norm_2"].values
        rho = spearmanr(v1, v2).correlation if n > 1 else None
        pairwise[f"{c1}_vs_{c2}"] = {"n": n, "rho": float(rho) if rho == rho else None}
        print(f"  {c1} vs {c2}: n={n}, rho={rho:.4f}" if rho == rho else f"  {c1} vs {c2}: n={n}, undefined")

    rhos = [v["rho"] for v in pairwise.values() if v["rho"] is not None]
    mean_pairwise_rho = float(np.mean(rhos))
    print(f"\n  Mean pairwise cross-condition rho (10 pairs): {mean_pairwise_rho:.4f}")
    print(f"  Range: [{min(rhos):.4f}, {max(rhos):.4f}]")
    print("  Under classical test theory (each condition = independent noisy draw of the same\n"
          "  underlying quantity, equal reliability per condition), corr(condition_i, condition_j) = reliability.\n"
          "  This mean IS the implied reliability estimate directly -- no further transformation needed.")

    return {"pairwise": pairwise, "mean_pairwise_rho": mean_pairwise_rho,
            "implied_reliability_ec_tx": mean_pairwise_rho,
            "caveat": "growth conditions are not technical replicates; this bundles real condition-effects into "
                      "the 'noise' term and therefore UNDERESTIMATES true single-condition reliability -- treat "
                      "the implied reliability as a LOWER bound and any resulting disattenuation as an UPPER bound "
                      "on how much of the cross-host gap measurement error alone could explain"}


def robust_rs_icc_crosscheck():
    """Cross-check via the 'Robust RSs' sheet's stdev(log2) -- a DIFFERENT
    estimation method (intraclass-correlation-style) on a DIFFERENT (and
    selection-biased -- chosen for being condition-robust) 100-sequence
    subset. Reported as a secondary sanity check, not the primary estimate,
    because selecting for low within-sequence variance mechanically inflates
    the reliability this method will find."""
    df = pd.read_excel(RAW / "NIHMS945382-supplement-5.xlsx", sheet_name="Robust RSs")
    within_var = float((df["stdev (log2)"] ** 2).mean())
    between_var = float(df["log2 mean activity"].var())
    icc = between_var / (between_var + within_var)
    print(f"\n  Robust-RSs-subset ICC cross-check (n={len(df)}, SELECTION-BIASED toward low variance): "
          f"between_var={between_var:.3f}, mean_within_var={within_var:.3f}, ICC={icc:.4f}")
    return {"n": len(df), "between_var": between_var, "mean_within_var": within_var, "icc": icc,
            "caveat": "this 100-sequence subset was selected FOR cross-condition robustness (low variance) -- "
                      "reliability estimate from it is biased upward relative to the full population, reported "
                      "as a secondary cross-check only"}


def disattenuate(rho_obs, rel_a, rel_b):
    denom = np.sqrt(rel_a * rel_b)
    if denom <= 0:
        return None
    corrected = rho_obs / denom
    return float(np.clip(corrected, -1.0, 1.0)), corrected > 1.0  # (clipped value, was_clipped)


def main():
    ec_result = ec_cross_condition_reliability()
    robust_result = robust_rs_icc_crosscheck()
    implied_ec_reliability = ec_result["implied_reliability_ec_tx"]

    corr = json.load(open(RESULTS / "gate5_5_crosshost_measurement_correlation.json"))

    print("\n" + "=" * 78)
    print("DISATTENUATION -- Scenario 1: EC uses its own empirically-derived reliability;")
    print("BS/PA use the sensitivity grid (no direct estimate available for either)")
    print("=" * 78)
    scenario1 = {}
    for readout in ["transcription", "translation"]:
        scenario1[readout] = {}
        for pair, hosts in [("EC_BS", ("EC", "BS")), ("EC_PA", ("EC", "PA")), ("BS_PA", ("BS", "PA"))]:
            rho_obs = corr[readout][pair]["rho"]
            scenario1[readout][pair] = {}
            for rel_other in SENSITIVITY_RELIABILITIES:
                if "EC" in hosts:
                    rel_a = implied_ec_reliability if readout == "transcription" else None
                else:
                    rel_a = None
                # EC gets its own estimate (tx only); the non-EC host(s) get the sensitivity value;
                # for translation, no empirical estimate exists for EC either -- both get the sensitivity value
                if readout == "transcription" and "EC" in hosts:
                    rel_ec = implied_ec_reliability
                    rel_partner = rel_other
                    corrected, clipped = disattenuate(rho_obs, rel_ec, rel_partner)
                else:
                    corrected, clipped = disattenuate(rho_obs, rel_other, rel_other)
                scenario1[readout][pair][f"other_host_reliability_{rel_other}"] = {
                    "rho_observed": rho_obs, "rho_corrected": corrected, "clipped_at_1": clipped}
            print(f"  {readout} {pair}: rho_obs={rho_obs:.3f}, corrected @ assumed-reliability "
                  f"{SENSITIVITY_RELIABILITIES} = "
                  f"{[round(scenario1[readout][pair][f'other_host_reliability_{r}']['rho_corrected'], 3) for r in SENSITIVITY_RELIABILITIES]}"
                  + (f"  [EC uses its own empirical reliability={implied_ec_reliability:.3f}]" if ("EC" in hosts and readout == "transcription") else ""))

    print("\n" + "=" * 78)
    print("DISATTENUATION -- Scenario 2: fully symmetric sensitivity grid, BOTH hosts in every")
    print("pair assumed at the same reliability (0.5 / 0.7 / 0.9), all pairs, both readouts")
    print("=" * 78)
    scenario2 = {}
    for readout in ["transcription", "translation"]:
        scenario2[readout] = {}
        for pair in ["EC_BS", "EC_PA", "BS_PA"]:
            rho_obs = corr[readout][pair]["rho"]
            scenario2[readout][pair] = {"rho_observed": rho_obs}
            for rel in SENSITIVITY_RELIABILITIES:
                corrected, clipped = disattenuate(rho_obs, rel, rel)
                scenario2[readout][pair][f"reliability_{rel}"] = {"rho_corrected": corrected, "clipped_at_1": clipped}
            print(f"  {readout} {pair}: rho_obs={rho_obs:.3f}, corrected @ "
                  f"{SENSITIVITY_RELIABILITIES} = "
                  f"{[round(scenario2[readout][pair][f'reliability_{r}']['rho_corrected'], 3) for r in SENSITIVITY_RELIABILITIES]}")

    print("\n" + "=" * 78)
    print("GAP-SURVIVAL CHECK: does corrected BS-pair correlation approach corrected EC-PA?")
    print("=" * 78)
    gap_check = {}
    for readout in ["transcription", "translation"]:
        ec_pa_obs = corr[readout]["EC_PA"]["rho"]
        bs_pairs_obs = [corr[readout]["EC_BS"]["rho"], corr[readout]["BS_PA"]["rho"]]
        gap_check[readout] = {"ec_pa_observed": ec_pa_obs, "bs_pairs_observed_range": [min(bs_pairs_obs), max(bs_pairs_obs)]}
        for rel in SENSITIVITY_RELIABILITIES:
            ec_pa_corrected = disattenuate(ec_pa_obs, rel, rel)[0]
            bs_pairs_corrected = [disattenuate(v, rel, rel)[0] for v in bs_pairs_obs]
            gap_ratio = np.mean(bs_pairs_corrected) / ec_pa_corrected if ec_pa_corrected else None
            gap_check[readout][f"reliability_{rel}"] = {
                "ec_pa_corrected": ec_pa_corrected, "bs_pairs_corrected": bs_pairs_corrected,
                "bs_over_ecpa_ratio": gap_ratio}
            print(f"  {readout} @ reliability={rel}: EC-PA corrected={ec_pa_corrected:.3f}, "
                  f"BS-pairs corrected={[round(x,3) for x in bs_pairs_corrected]}, "
                  f"BS/EC-PA ratio={gap_ratio:.3f} (observed ratio was {np.mean(bs_pairs_obs)/ec_pa_obs:.3f})")

    output = {
        "ec_cross_condition_reliability": ec_result,
        "robust_rs_icc_crosscheck": robust_result,
        "scenario1_ec_empirical_plus_sensitivity": scenario1,
        "scenario2_symmetric_sensitivity_grid": scenario2,
        "gap_survival_check": gap_check,
    }
    with open(RESULTS / "gate8_attenuation_analysis.json", "w") as f:
        json.dump(output, f, indent=2, default=str)

    rows = []
    for scenario_name, scenario in [("symmetric_grid", scenario2)]:
        for readout, pairs in scenario.items():
            for pair, vals in pairs.items():
                for rel in SENSITIVITY_RELIABILITIES:
                    rows.append({"scenario": scenario_name, "readout": readout, "pair": pair,
                                 "assumed_reliability": rel, "rho_observed": vals["rho_observed"],
                                 "rho_corrected": vals[f"reliability_{rel}"]["rho_corrected"],
                                 "clipped_at_1": vals[f"reliability_{rel}"]["clipped_at_1"]})
    pd.DataFrame(rows).to_csv(RESULTS / "gate8_attenuation_analysis.csv", index=False)
    print(f"\nWrote {RESULTS / 'gate8_attenuation_analysis.json'} and .csv")


if __name__ == "__main__":
    main()
