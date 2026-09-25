"""
GATE 15 - Is the cell-free/in-vivo modality contrast a restriction artifact?

THE CHALLENGE (raised externally, and correct on its face):
  Section 3.2 compares cell-free EC-BS rho=0.597 against in-vivo EC-BS
  rho=0.258. The in-vivo 0.258 is the CO-ACTIVE value (pooled is 0.655).
  Every DRAFTS correlation in gate10_crosshost_correlations.csv is
  computed over `n_shared_usable`. If DRAFTS-"usable" means the same
  thing as Johns-"usable", then the published comparison puts a pooled
  cell-free number against a co-active in-vivo number, and the like-for-
  like pooled comparison (0.677 vs 0.655) shows no modality difference
  at all.

WHAT THIS SCRIPT ESTABLISHES FIRST (Task 1), because it decides everything:
  The two datasets do NOT use "usable" to mean the same thing.

    Johns (in vivo):  usable = DNA template present (rna may be 0)
                      active = usable AND rna > 0
                      -> inactive rows are KEPT, with tx_norm exactly 0
                      -> 38.9% (EC) / 72.0% (BS) of usable rows are inactive

    DRAFTS (cell-free): usable = numeric tx present, which REQUIRES rna>0
                      rows with template but no RNA are marked
                      'no_RNA_counts' and EXCLUDED from usable
                      -> DRAFTS-usable already applies the rna>0 filter
                      -> i.e. DRAFTS-usable is the analogue of Johns-ACTIVE,
                         not of Johns-usable

  So the regime-matched reconstruction is:
    in-vivo POOLED      == DRAFTS usable UNION no_RNA_counts (tx:=0)
    in-vivo CO-ACTIVE   == DRAFTS usable (both species)

  This is a naming collision between two datasets from the same lab, of
  exactly the same kind as the "Pa" collision caught in Gate 10.

WHAT IT THEN COMPUTES (Task 2), reporting ALL thresholds, not a chosen one:
  For every DRAFTS species pair, on both the full shared set and the
  82-sequence Johns-overlap set:
    pooled_equiv    both species DNA-adequate; no_RNA rows carried as 0
    co_detected     both species RNA-detected  (what the manuscript used)
    both_above_p25  both above own 25th percentile
    both_above_med  both above own median
    matched_26pct   both in own top 26.04%, matching the in-vivo co-active
                    retention rate (3668/14088) -- a QUANTILE match, which
                    is a deliberately harsher and NOT definitionally
                    equivalent test; see the memo.
  Plus IQR range-restriction diagnostics, mirroring Gate 8.6's in-vivo check.

No training. No manuscript edits. Reads only existing frozen tables.
"""
import json
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "out" / "results"

SPECIES = ["Ec", "Ef", "Se", "Ko", "Pa", "Pp", "Vn", "Bs", "Cg", "Ll"]
SPECIES_FULL = {
    "Ec": "E. coli", "Ef": "E. fergusonii", "Se": "S. enterica", "Ko": "K. oxytoca",
    "Pa": "P. agglomerans", "Pp": "P. putida", "Vn": "V. natriegens",
    "Bs": "B. subtilis", "Cg": "C. glutamicum", "Ll": "L. lactis",
}
# in-vivo EC-BS co-active retention: 3668 / 14088
INVIVO_COACTIVE_RETENTION = 3668 / 14088


def load_drafts(path):
    df = pd.read_parquet(path)
    out = {}
    for s in SPECIES:
        tx = pd.to_numeric(df[f"{s}_tx"], errors="coerce")
        usable = df[f"{s}_usable"].astype(bool)
        reason = df[f"{s}_unusable_reason"].astype("string")
        # 'no_RNA_counts' = template present, zero RNA detected = INACTIVE
        # (the in-vivo rna==0 class). 'no_DNA_counts'/'low_DNA_counts' =
        # template absent/insufficient = NOT MEASURABLE (in-vivo not-usable).
        no_rna = reason.fillna("") == "no_RNA_counts"
        activity = tx.copy()
        activity[no_rna] = 0.0                 # carry inactives as true zeros
        dna_adequate = usable | no_rna          # in-vivo "usable" analogue
        out[s] = {"tx": tx, "usable": usable, "no_rna": no_rna,
                  "activity": activity, "dna_adequate": dna_adequate}
    return df, out


def describe_activity(cols, label):
    """Task 1: is 'active' even definable for DRAFTS transcription?"""
    rows = []
    for s in SPECIES:
        c = cols[s]
        tx, usable, no_rna, dna_ok = c["tx"], c["usable"], c["no_rna"], c["dna_adequate"]
        v = tx[usable].dropna()
        logv = np.log10(v[v > 0])
        # detection-floor / pile-up check: most common rounded value
        top_frac = (v.round(4).value_counts().iloc[0] / len(v)) if len(v) else np.nan
        rows.append({
            "dataset": label, "species": s, "species_full": SPECIES_FULL[s],
            "n_dna_adequate": int(dna_ok.sum()),
            "n_rna_detected_usable": int(usable.sum()),
            "n_no_rna_inactive": int(no_rna.sum()),
            "inactive_fraction_of_dna_adequate": float(no_rna.sum() / dna_ok.sum()) if dna_ok.sum() else np.nan,
            "tx_min": float(v.min()) if len(v) else np.nan,
            "tx_p25": float(v.quantile(.25)) if len(v) else np.nan,
            "tx_median": float(v.median()) if len(v) else np.nan,
            "tx_p75": float(v.quantile(.75)) if len(v) else np.nan,
            "tx_max": float(v.max()) if len(v) else np.nan,
            "n_exact_zero_among_usable": int((v == 0).sum()),
            "max_single_rounded_value_fraction": float(top_frac),
            "log10_tx_skew": float(pd.Series(logv).skew()) if len(logv) > 2 else np.nan,
            "log10_tx_kurtosis": float(pd.Series(logv).kurtosis()) if len(logv) > 2 else np.nan,
        })
    return pd.DataFrame(rows)


def pair_regimes(cols, a, b):
    """Task 2: rho and n for each restriction regime, for one species pair."""
    A, B = cols[a], cols[b]
    res = {}

    # 1. pooled-equivalent: DNA adequate in both, inactives carried as 0
    m = A["dna_adequate"] & B["dna_adequate"]
    xa, xb = A["activity"][m], B["activity"][m]
    ok = xa.notna() & xb.notna()
    res["pooled_equiv"] = _rho(xa[ok], xb[ok])

    # 2. co-detected: RNA detected in both (== in-vivo co-active analogue)
    m = A["usable"] & B["usable"]
    xa, xb = A["tx"][m], B["tx"][m]
    ok = xa.notna() & xb.notna()
    codet_a, codet_b = xa[ok], xb[ok]
    res["co_detected"] = _rho(codet_a, codet_b)

    # thresholds are computed on each species' own usable distribution
    for name, q in [("both_above_p25", .25), ("both_above_med", .50)]:
        ta = A["tx"][A["usable"]].quantile(q)
        tb = B["tx"][B["usable"]].quantile(q)
        m = A["usable"] & B["usable"] & (A["tx"] > ta) & (B["tx"] > tb)
        xa, xb = A["tx"][m], B["tx"][m]
        ok = xa.notna() & xb.notna()
        res[name] = _rho(xa[ok], xb[ok])

    # 5. quantile-matched to the in-vivo INTERSECTION retention rate (26.0%),
    #    applied to each species' own DNA-adequate distribution. This is the
    #    harshest variant: in vivo the 26% is the *intersection* of two
    #    per-host filters that individually retain 61% (EC) and 28% (BS),
    #    so demanding top-26% of BOTH cuts far more deeply than in vivo did.
    qa = A["activity"][A["dna_adequate"]].quantile(1 - INVIVO_COACTIVE_RETENTION)
    qb = B["activity"][B["dna_adequate"]].quantile(1 - INVIVO_COACTIVE_RETENTION)
    m26 = A["dna_adequate"] & B["dna_adequate"] & (A["activity"] >= qa) & (B["activity"] >= qb)
    xa, xb = A["activity"][m26], B["activity"][m26]
    ok = xa.notna() & xb.notna()
    res["matched_26pct"] = _rho(xa[ok], xb[ok])

    # 6. matched to each host's OWN in-vivo active fraction, then intersected
    #    -- the faithful analogue of what the in-vivo co-active filter does.
    #    EC retains 61.1% of its usable rows, BS 28.0%; applying those same
    #    per-species retention rates here and intersecting reproduces the
    #    in-vivo filter's *structure*, not just its final retained fraction.
    per_host = {"Ec": 15040 / 24613, "Bs": 4435 / 15848, "Pa": 17886 / 21473}
    fa, fb = per_host.get(a), per_host.get(b)
    if fa is not None and fb is not None:
        qa2 = A["activity"][A["dna_adequate"]].quantile(1 - fa)
        qb2 = B["activity"][B["dna_adequate"]].quantile(1 - fb)
        m2 = A["dna_adequate"] & B["dna_adequate"] & (A["activity"] >= qa2) & (B["activity"] >= qb2)
        xa, xb = A["activity"][m2], B["activity"][m2]
        ok = xa.notna() & xb.notna()
        res["matched_per_host_active_frac"] = _rho(xa[ok], xb[ok])
    else:
        res["matched_per_host_active_frac"] = {"n": 0, "rho": None}
        m2 = None

    # range-restriction diagnostic (mirrors Gate 8.6's in-vivo IQR check):
    # if a restricted subset's IQR is WIDER, mechanical range restriction
    # cannot explain a correlation drop in that subset.
    def iqr(x):
        x = pd.Series(x).dropna()
        return float(x.quantile(.75) - x.quantile(.25)) if len(x) > 3 else np.nan
    pooled_m = A["dna_adequate"] & B["dna_adequate"]
    res["_diag"] = {
        "iqr_pooled_a": iqr(A["activity"][pooled_m]),
        "iqr_codetected_a": iqr(codet_a),
        "iqr_matched26_a": iqr(A["activity"][m26]),
        "iqr_pooled_b": iqr(B["activity"][pooled_m]),
        "iqr_codetected_b": iqr(codet_b),
        "iqr_matched26_b": iqr(B["activity"][m26]),
    }
    return res


def _rho(x, y):
    n = int(len(x))
    if n < 3:
        return {"n": n, "rho": None}
    r = spearmanr(x, y).correlation
    return {"n": n, "rho": float(r) if r == r else None}


def main():
    RESULTS.mkdir(parents=True, exist_ok=True)
    full_df, full = load_drafts(DATA / "drafts.parquet")
    ovl_df, ovl = load_drafts(DATA / "drafts_johns_overlap.parquet")

    # ---- Task 1 -------------------------------------------------------
    d1 = describe_activity(full, "drafts_full_1047")
    print("=" * 78)
    print("TASK 1: DRAFTS cell-free activity distribution, per species")
    print("=" * 78)
    print(d1[["species", "n_dna_adequate", "n_rna_detected_usable", "n_no_rna_inactive",
              "inactive_fraction_of_dna_adequate", "n_exact_zero_among_usable",
              "max_single_rounded_value_fraction", "tx_min", "tx_median"]].to_string(index=False))

    invivo_inactive = {"EC": 1 - 15040 / 24613, "BS": 1 - 4435 / 15848, "PA": 1 - 17886 / 21473}
    print("\nIN-VIVO comparison (inactive fraction of usable, from three_host.parquet):")
    for k, v in invivo_inactive.items():
        print(f"  {k}: {v*100:.1f}% inactive (tx_norm exactly 0)")
    print(f"\nDRAFTS inactive fraction range: "
          f"{d1['inactive_fraction_of_dna_adequate'].min()*100:.1f}% - "
          f"{d1['inactive_fraction_of_dna_adequate'].max()*100:.1f}%")

    # ---- Task 2 -------------------------------------------------------
    regimes = ["pooled_equiv", "co_detected", "both_above_p25", "both_above_med",
               "matched_per_host_active_frac", "matched_26pct"]
    rows, diags = [], []
    for setname, cols in [("full_shared_1047", full), ("johns_overlap_112", ovl)]:
        for a, b in itertools.combinations(SPECIES, 2):
            r = pair_regimes(cols, a, b)
            rec = {"set": setname, "species_a": a, "species_b": b,
                   "pair": f"{a}-{b}", "is_EC_BS": (a, b) == ("Ec", "Bs")}
            for g in regimes:
                rec[f"{g}_n"] = r[g]["n"]
                rec[f"{g}_rho"] = r[g]["rho"]
            rows.append(rec)
            if (a, b) == ("Ec", "Bs"):
                diags.append({"set": setname, **r["_diag"]})
    pairs = pd.DataFrame(rows)

    print("\n" + "=" * 78)
    print("TASK 2: EC-BS under every restriction regime")
    print("=" * 78)
    ecbs = pairs[pairs.is_EC_BS]
    for _, r in ecbs.iterrows():
        print(f"\n  [{r['set']}]")
        for g in regimes:
            rho = r[f"{g}_rho"]
            print(f"    {g:16} n={int(r[f'{g}_n']):5d}  rho={rho:.4f}" if rho is not None
                  else f"    {g:16} n={int(r[f'{g}_n']):5d}  rho=None")

    print("\n  IN-VIVO EC-BS reference (gate8_5_coactive_correlation.csv):")
    print(f"    pooled            n=14088  rho=0.6551")
    print(f"    co-active         n= 3668  rho=0.2575   <-- the manuscript's 0.258")

    print("\n  Range-restriction diagnostic (IQR: pooled-equiv -> co-detected -> matched26):")
    for d in diags:
        print(f"    [{d['set']}] Ec: {d['iqr_pooled_a']:.4g} -> {d['iqr_codetected_a']:.4g} -> {d['iqr_matched26_a']:.4g}")
        print(f"    [{d['set']}] Bs: {d['iqr_pooled_b']:.4g} -> {d['iqr_codetected_b']:.4g} -> {d['iqr_matched26_b']:.4g}")

    print("\n" + "=" * 78)
    print("TASK 2b: all 45 pairs, rho range under each regime (full shared set)")
    print("=" * 78)
    fs = pairs[pairs.set == "full_shared_1047"]
    summary45 = {}
    for g in regimes:
        v = fs[f"{g}_rho"].dropna()
        summary45[g] = {"n_pairs": int(len(v)), "min": float(v.min()), "max": float(v.max()),
                        "mean": float(v.mean()), "median": float(v.median())}
        print(f"  {g:16} {len(v):2d} pairs  range [{v.min():.3f}, {v.max():.3f}]  mean {v.mean():.3f}")

    # ---- outputs ------------------------------------------------------
    pairs.to_csv(RESULTS / "gate15_coactive_modality.csv", index=False)
    d1.to_csv(RESULTS / "gate15_drafts_activity_distribution.csv", index=False)

    ecbs_full = ecbs[ecbs.set == "full_shared_1047"].iloc[0]
    ecbs_ovl = ecbs[ecbs.set == "johns_overlap_112"].iloc[0]
    payload = {
        "task1_activity_distribution": {
            "drafts_inactive_fraction_of_dna_adequate": {
                r["species"]: r["inactive_fraction_of_dna_adequate"] for _, r in d1.iterrows()},
            "invivo_inactive_fraction_of_usable": invivo_inactive,
            "drafts_exact_zeros_among_usable": int(d1["n_exact_zero_among_usable"].sum()),
            "drafts_max_single_rounded_value_fraction": float(d1["max_single_rounded_value_fraction"].max()),
            "naming_collision": ("DRAFTS-'usable' requires rna>0 and is therefore the analogue of "
                                  "Johns-'active', NOT of Johns-'usable'. Johns keeps rna==0 rows as "
                                  "usable-but-inactive with tx_norm exactly 0; DRAFTS excludes them as "
                                  "'no_RNA_counts'."),
        },
        "task2_ec_bs": {
            "full_shared_1047": {g: {"n": int(ecbs_full[f"{g}_n"]), "rho": ecbs_full[f"{g}_rho"]} for g in regimes},
            "johns_overlap_112": {g: {"n": int(ecbs_ovl[f"{g}_n"]), "rho": ecbs_ovl[f"{g}_rho"]} for g in regimes},
            "invivo_reference": {"pooled": {"n": 14088, "rho": 0.6550978207525179},
                                  "co_active": {"n": 3668, "rho": 0.2575175805660613}},
            "range_restriction_diagnostic": diags,
        },
        "task2b_all_45_pairs_full_set": summary45,
        "invivo_coactive_retention_used_for_quantile_match": INVIVO_COACTIVE_RETENTION,
    }
    with open(RESULTS / "gate15_coactive_modality.json", "w") as f:
        json.dump(payload, f, indent=2)
    print(f"\nWrote {RESULTS/'gate15_coactive_modality.json'}, .csv, "
          f"{RESULTS/'gate15_drafts_activity_distribution.csv'}")


if __name__ == "__main__":
    main()
