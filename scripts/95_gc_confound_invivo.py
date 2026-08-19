"""
GATE 10.5 - Task 2 [LOAD-BEARING]: does the in-vivo EC-PA vs. BS-pairs
cross-host correlation contrast survive controlling for source-genome GC
content?

MOTIVATION: Gate 10 found GC%-vs-activity Spearman correlations of -0.49
to -0.74 in ALL TEN DRAFTS cell-free hosts -- a strong, universal
covariate. If the SAME relationship holds in vivo, part of the apparent
cross-host correlation structure could be a shared property of the
sequences (their source-genome GC content) rather than of the hosts
themselves -- two hosts could correlate simply because both respond to GC
the same way, independent of anything host-specific.

METHOD for the GC-controlled ("partial") correlation: rank-transform each
host's activity value and GC% (Spearman's own definition operates on
ranks), regress ranked-activity ~ ranked-GC per host (ordinary least
squares), take the residuals, then Pearson-correlate the two hosts'
residuals. This is the standard rank-residual approach to a partial
Spearman correlation and is reported ALONGSIDE the closed-form partial-
correlation formula (pcorr = (r_xy - r_xz*r_yz) / sqrt((1-r_xz^2)(1-r_yz^2)))
as a cross-check -- if the two methods disagree substantially, that itself
is worth flagging.

POPULATION: matches scripts/57's "active in both hosts of the pair"
definition exactly (the same restriction underlying this project's
headline EC-PA rho~0.75 / BS-pairs rho~0.16-0.26 numbers), so the
GC-controlled figures are directly comparable to the numbers already in
PAPER_FRAMING.md, not a different population.
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr, rankdata

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

HOSTS = ["EC", "BS", "PA"]
PAIRS = [("EC", "BS"), ("EC", "PA"), ("BS", "PA")]


def gc_content(seq):
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    return 100.0 * gc / len(seq)


def rank_residual(activity, gc):
    """Rank-transform both, OLS-regress ranked activity on ranked GC,
    return residuals -- the input to the residual-based partial correlation."""
    r_act = rankdata(activity)
    r_gc = rankdata(gc)
    slope, intercept = np.polyfit(r_gc, r_act, 1)
    resid = r_act - (slope * r_gc + intercept)
    return resid


def partial_spearman_formula(r_xy, r_xz, r_yz):
    denom = np.sqrt((1 - r_xz ** 2) * (1 - r_yz ** 2))
    if denom < 1e-9:
        return None
    return (r_xy - r_xz * r_yz) / denom


def load_data():
    lib = pd.read_parquet(DATA / "three_host_library.parquet")
    three_host = pd.read_parquet(DATA / "three_host.parquet")
    assert (lib["OLIGO ID"].values == three_host["OLIGO ID"].values).all()
    df = three_host.copy()
    df["source_phylum"] = lib["phylum"].values
    df["gc_pct"] = lib["Regulatory Sequence"].apply(gc_content).values
    return df


def get_floor_and_active(df, host):
    """Reproduces scripts/57's translation floor-detection logic exactly."""
    usable = df[f"{host}_tl_usable"].values.astype(bool)
    protein = df[f"{host}_protein_log10"].values.astype(np.float32)
    vals = protein[usable]
    rounded = np.round(vals, 2)
    u, c = np.unique(rounded, return_counts=True)
    mode_val, mode_frac = u[np.argmax(c)], c[np.argmax(c)] / len(vals)
    is_floor = mode_frac > 0.05
    if is_floor:
        at_floor = usable & (np.round(protein, 2) == mode_val)
        tl_usable_reg = usable & ~at_floor
        tl_active = usable & (np.round(protein, 2) > mode_val)
    else:
        tl_usable_reg = usable.copy()
        tl_active = usable & (protein > np.median(vals))
    return tl_usable_reg, tl_active


def task1_gc_vs_activity_invivo(df):
    print("=" * 80)
    print("TASK 2.1: GC% vs. activity, in-vivo, per host, both readouts")
    print("=" * 80)
    rows = []
    for host in HOSTS:
        tx_mask = df[f"{host}_tx_usable"] & df[f"{host}_tx_active"]
        rho_tx = float(spearmanr(df.loc[tx_mask, "gc_pct"], df.loc[tx_mask, f"{host}_tx_norm"]).correlation)
        rows.append({"host": host, "readout": "transcription", "n": int(tx_mask.sum()), "spearman_rho_gc_vs_activity": rho_tx})
        print(f"  {host} transcription: n={int(tx_mask.sum())}, rho(GC%, tx_norm)={rho_tx:.3f}")

        tl_usable_reg, tl_active = get_floor_and_active(df, host)
        tl_mask = tl_usable_reg & tl_active
        rho_tl = float(spearmanr(df.loc[tl_mask, "gc_pct"], df.loc[tl_mask, f"{host}_protein_log10"]).correlation)
        rows.append({"host": host, "readout": "translation", "n": int(tl_mask.sum()), "spearman_rho_gc_vs_activity": rho_tl})
        print(f"  {host} translation:   n={int(tl_mask.sum())}, rho(GC%, protein_log10)={rho_tl:.3f}")

    print(f"\nCOMPARE TO DRAFTS (cell-free): range was -0.74 to -0.49 across all 10 hosts.")
    return rows


def task2_partial_correlations(df):
    print("\n" + "=" * 80)
    print("TASK 2.2: raw vs. GC-controlled cross-host correlations, all 3 pairs, both readouts")
    print("=" * 80)
    results = {}
    for readout in ["transcription", "translation"]:
        results[readout] = {}
        for h1, h2 in PAIRS:
            if readout == "transcription":
                both = (df[f"{h1}_tx_usable"] & df[f"{h1}_tx_active"] &
                        df[f"{h2}_tx_usable"] & df[f"{h2}_tx_active"])
                v1_col, v2_col = f"{h1}_tx_norm", f"{h2}_tx_norm"
            else:
                u1, a1 = get_floor_and_active(df, h1)
                u2, a2 = get_floor_and_active(df, h2)
                both = u1 & a1 & u2 & a2
                v1_col, v2_col = f"{h1}_protein_log10", f"{h2}_protein_log10"

            n = int(both.sum())
            v1 = df.loc[both, v1_col].values
            v2 = df.loc[both, v2_col].values
            gc = df.loc[both, "gc_pct"].values

            r_xy = spearmanr(v1, v2).correlation
            r_xz = spearmanr(v1, gc).correlation
            r_yz = spearmanr(v2, gc).correlation
            r_partial_formula = partial_spearman_formula(r_xy, r_xz, r_yz)

            resid1 = rank_residual(v1, gc)
            resid2 = rank_residual(v2, gc)
            r_partial_residual = float(np.corrcoef(resid1, resid2)[0, 1])

            attenuation_pct = (100 * (r_xy - r_partial_residual) / r_xy) if r_xy else None

            results[readout][f"{h1}_{h2}"] = {
                "n": n, "rho_raw": float(r_xy),
                "rho_gc_controlled_formula": float(r_partial_formula) if r_partial_formula is not None else None,
                "rho_gc_controlled_residual": r_partial_residual,
                "pct_attenuation_residual_method": float(attenuation_pct) if attenuation_pct is not None else None,
            }
            print(f"  {readout} {h1}-{h2}: n={n}, raw rho={r_xy:.3f}, "
                  f"GC-controlled (formula)={r_partial_formula:.3f}, (residual)={r_partial_residual:.3f} "
                  f"[{attenuation_pct:+.1f}% change]")
    return results


def task3_verdict(partial_results):
    print("\n" + "=" * 80)
    print("TASK 2.3: DECISIVE QUESTION -- does the EC-PA vs. BS-pairs contrast survive?")
    print("=" * 80)
    verdicts = {}
    for readout in ["transcription", "translation"]:
        r = partial_results[readout]
        ec_pa_raw = r["EC_PA"]["rho_raw"]
        ec_pa_gc = r["EC_PA"]["rho_gc_controlled_residual"]
        bs_pairs_raw = [r["EC_BS"]["rho_raw"], r["BS_PA"]["rho_raw"]]
        bs_pairs_gc = [r["EC_BS"]["rho_gc_controlled_residual"], r["BS_PA"]["rho_gc_controlled_residual"]]

        gap_raw = np.mean(bs_pairs_raw) / ec_pa_raw if ec_pa_raw else None
        gap_gc = np.mean(bs_pairs_gc) / ec_pa_gc if ec_pa_gc else None
        gap_ratio_change_pct = 100 * (gap_gc - gap_raw) / gap_raw if gap_raw else None

        # verdict thresholds: "survives intact" if gap ratio changes by <10% relative;
        # "attenuated" if it changes by 10-50%; "does not survive" if BS-pairs approach
        # EC-PA (gap ratio moves substantially toward 1.0) or contrast direction reverses
        if ec_pa_gc <= 0 or (min(bs_pairs_gc) >= ec_pa_gc):
            verdict = "DOES NOT SURVIVE"
        elif abs(gap_ratio_change_pct) < 10:
            verdict = "SURVIVES INTACT"
        elif gap_gc > 0.5:
            verdict = "DOES NOT SURVIVE"
        else:
            verdict = "SURVIVES, ATTENUATED"

        verdicts[readout] = {
            "ec_pa_raw": ec_pa_raw, "ec_pa_gc_controlled": ec_pa_gc,
            "bs_pairs_raw_mean": float(np.mean(bs_pairs_raw)), "bs_pairs_gc_controlled_mean": float(np.mean(bs_pairs_gc)),
            "gap_ratio_raw": float(gap_raw), "gap_ratio_gc_controlled": float(gap_gc),
            "gap_ratio_pct_change": float(gap_ratio_change_pct),
            "verdict": verdict,
        }
        print(f"  {readout}: EC-PA raw={ec_pa_raw:.3f} -> GC-controlled={ec_pa_gc:.3f}; "
              f"BS-pairs mean raw={np.mean(bs_pairs_raw):.3f} -> GC-controlled={np.mean(bs_pairs_gc):.3f}")
        print(f"    gap ratio (BS-pairs/EC-PA): raw={gap_raw:.3f} -> GC-controlled={gap_gc:.3f} "
              f"({gap_ratio_change_pct:+.1f}% change)")
        print(f"    VERDICT: {verdict}")
    return verdicts


def task4_phylum_stratified(df):
    print("\n" + "=" * 80)
    print("TASK 2.4: cross-host correlation within source-phylum strata")
    print("=" * 80)
    phylum_counts = df["source_phylum"].value_counts()
    major_phyla = phylum_counts[phylum_counts >= 500].index.tolist()
    print(f"  Phyla with >=500 sequences: {dict(phylum_counts[phylum_counts >= 500])}")

    results = {}
    for readout in ["transcription", "translation"]:
        results[readout] = {}
        for h1, h2 in PAIRS:
            results[readout][f"{h1}_{h2}"] = {}
            for phylum in major_phyla:
                phylum_mask = df["source_phylum"] == phylum
                if readout == "transcription":
                    both = (df[f"{h1}_tx_usable"] & df[f"{h1}_tx_active"] &
                            df[f"{h2}_tx_usable"] & df[f"{h2}_tx_active"] & phylum_mask)
                    v1_col, v2_col = f"{h1}_tx_norm", f"{h2}_tx_norm"
                else:
                    u1, a1 = get_floor_and_active(df, h1)
                    u2, a2 = get_floor_and_active(df, h2)
                    both = u1 & a1 & u2 & a2 & phylum_mask.values
                    v1_col, v2_col = f"{h1}_protein_log10", f"{h2}_protein_log10"
                n = int(both.sum())
                if n > 20:
                    rho = float(spearmanr(df.loc[both, v1_col], df.loc[both, v2_col]).correlation)
                else:
                    rho = None
                results[readout][f"{h1}_{h2}"][phylum] = {"n": n, "spearman_rho": rho}
                print(f"  {readout} {h1}-{h2} | source={phylum}: n={n}, rho={rho}")
    return results


def disattenuate(rho_obs, rel_a, rel_b):
    denom = np.sqrt(rel_a * rel_b)
    corrected = rho_obs / denom if denom > 0 else None
    clipped = corrected is not None and corrected > 1.0
    return (min(corrected, 1.0) if clipped else corrected), clipped


def task5_ceiling_interaction(partial_results):
    print("\n" + "=" * 80)
    print("TASK 2.5: does GC control change the gap ratio Gate 8 reported (0.341 @ reliability=0.9)?")
    print("=" * 80)
    # Gate 8's own reliability anchor: EC transcription = 0.912
    results = {}
    for readout in ["transcription", "translation"]:
        r = partial_results[readout]
        ec_pa_raw = r["EC_PA"]["rho_raw"]
        bs_pairs_raw = [r["EC_BS"]["rho_raw"], r["BS_PA"]["rho_raw"]]
        ec_pa_gc = r["EC_PA"]["rho_gc_controlled_residual"]
        bs_pairs_gc = [r["EC_BS"]["rho_gc_controlled_residual"], r["BS_PA"]["rho_gc_controlled_residual"]]

        row = {}
        for rel in [0.9]:
            ec_pa_disatt_raw = disattenuate(ec_pa_raw, rel, rel)[0]
            bs_pairs_disatt_raw = [disattenuate(v, rel, rel)[0] for v in bs_pairs_raw]
            gap_disatt_raw = np.mean(bs_pairs_disatt_raw) / ec_pa_disatt_raw

            ec_pa_disatt_gc = disattenuate(ec_pa_gc, rel, rel)[0]
            bs_pairs_disatt_gc = [disattenuate(v, rel, rel)[0] for v in bs_pairs_gc]
            gap_disatt_gc = np.mean(bs_pairs_disatt_gc) / ec_pa_disatt_gc

            row[f"reliability_{rel}"] = {
                "gap_ratio_disattenuated_only_raw_correlations": float(gap_disatt_raw),
                "gap_ratio_disattenuated_AND_gc_controlled": float(gap_disatt_gc),
            }
            print(f"  {readout} @ reliability={rel}: gap ratio (disattenuation only) = {gap_disatt_raw:.3f}; "
                  f"gap ratio (disattenuation + GC control) = {gap_disatt_gc:.3f}")
        results[readout] = row
    return results


def main():
    df = load_data()
    print(f"Loaded {len(df)} sequences, GC content computed directly from 165bp sequence text.\n")

    t1 = task1_gc_vs_activity_invivo(df)
    t2 = task2_partial_correlations(df)
    t3 = task3_verdict(t2)
    t4 = task4_phylum_stratified(df)
    t5 = task5_ceiling_interaction(t2)

    output = {
        "task1_gc_vs_activity_invivo": t1,
        "task2_raw_vs_gc_controlled_correlations": t2,
        "task3_verdict": t3,
        "task4_phylum_stratified": t4,
        "task5_gap_ratio_with_gc_control": t5,
    }
    with open(RESULTS / "gate10_5_gc_confound.json", "w") as f:
        json.dump(output, f, indent=2, default=str)

    pd.DataFrame(t1).to_csv(RESULTS / "gate10_5_gc_vs_activity_invivo.csv", index=False)
    rows2 = [{"readout": ro, "pair": p, **v} for ro in t2 for p, v in t2[ro].items()]
    pd.DataFrame(rows2).to_csv(RESULTS / "gate10_5_gc_confound.csv", index=False)

    print(f"\nWrote {RESULTS / 'gate10_5_gc_confound.json'}, {RESULTS / 'gate10_5_gc_confound.csv'}")

    overall_verdict = set(v["verdict"] for v in t3.values())
    print(f"\n=== OVERALL VERDICT (both readouts): {overall_verdict} ===")


if __name__ == "__main__":
    main()
