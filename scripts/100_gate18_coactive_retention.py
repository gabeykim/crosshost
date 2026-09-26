"""
GATE 18 - Provenance for the co-active retention rates, and verification of
the two matched-restriction regime definitions.

WHY THIS EXISTS: the retention figures (93.6% cell-free vs 26.0% in vivo)
were derived by hand from gate15_coactive_modality and
gate8_5_coactive_correlation.csv. This project has twice had unsourced
numbers propagate into written documents (the percent-of-ceiling metric,
the shift-prediction result), so a figure carrying the paper's mechanism
gets a producing script.

METHOD: every number is RECOMPUTED FROM PRIMARY DATA
(data/three_host.parquet, data/drafts.parquet), replicating the exact
filter definitions used by scripts/80 (in vivo) and scripts/98 (cell-free),
then CROSS-CHECKED against those scripts' own committed outputs. A
disagreement between the recomputation and the stored file is reported as
a failure, not reconciled silently.

TASK 1 output: out/results/gate18_coactive_retention.csv
TASK 2 output: printed regime verification (thresholds, populations,
               comparison operators, achieved retention)
"""
import json
import itertools
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "out" / "results"

HOSTS = ["EC", "BS", "PA"]
INVIVO_PAIRS = [("EC", "BS"), ("EC", "PA"), ("BS", "PA")]
SPECIES = ["Ec", "Ef", "Se", "Ko", "Pa", "Pp", "Vn", "Bs", "Cg", "Ll"]
INVIVO_COACTIVE_RETENTION = 3668 / 14088   # the constant scripts/98 uses

problems = []


def invivo_masks(df):
    """Replicate scripts/80_coactive_subset.py exactly, including the
    translation floor logic (which makes translation's 'pooled' already
    floor-corrected -- NOT the analogue of transcription's pooled)."""
    tl_active, tl_usable_reg, floor_info = {}, {}, {}
    for h in HOSTS:
        usable = df[f"{h}_tl_usable"].values.astype(bool)
        protein = df[f"{h}_protein_log10"].values.astype(np.float32)
        vals = protein[usable]
        rounded = np.round(vals, 2)
        u, c = np.unique(rounded, return_counts=True)
        mode_val, mode_frac = u[np.argmax(c)], c[np.argmax(c)] / len(vals)
        is_floor = mode_frac > 0.05
        floor_info[h] = {"floor_value": float(mode_val), "floor_fraction": float(mode_frac),
                          "floor_detected": bool(is_floor)}
        if is_floor:
            at_floor = usable & (np.round(protein, 2) == mode_val)
            tl_usable_reg[h] = usable & ~at_floor
            tl_active[h] = usable & (np.round(protein, 2) > mode_val)
        else:
            tl_usable_reg[h] = usable.copy()
            tl_active[h] = usable & (protein > np.median(vals))
    return tl_active, tl_usable_reg, floor_info


def invivo_rows(df):
    tl_active, tl_usable_reg, floor_info = invivo_masks(df)
    rows = []
    for h1, h2 in INVIVO_PAIRS:
        # TRANSCRIPTION: pooled = usable both (rna==0 rows KEPT, tx_norm 0)
        #                co-active = usable AND rna>0, both hosts
        pooled = (df[f"{h1}_tx_usable"] & df[f"{h2}_tx_usable"]).values
        coact = (df[f"{h1}_tx_usable"] & df[f"{h1}_tx_active"] &
                 df[f"{h2}_tx_usable"] & df[f"{h2}_tx_active"]).values
        rows.append({"dataset": "in_vivo_johns", "set": "full_library_29249",
                     "pair": f"{h1}-{h2}", "readout": "transcription",
                     "n_pooled": int(pooled.sum()), "n_coactive": int(coact.sum()),
                     "pooled_definition": "tx_usable in both (DNA present; rna==0 rows kept as tx_norm=0)",
                     "coactive_definition": "tx_usable AND tx_active (rna>0) in both"})
        # TRANSLATION: pooled is ALREADY floor-corrected -- see note column
        pooled_tl = (tl_usable_reg[h1] & tl_usable_reg[h2])
        coact_tl = (tl_active[h1] & tl_active[h2] & tl_usable_reg[h1] & tl_usable_reg[h2])
        rows.append({"dataset": "in_vivo_johns", "set": "full_library_29249",
                     "pair": f"{h1}-{h2}", "readout": "translation",
                     "n_pooled": int(pooled_tl.sum()), "n_coactive": int(coact_tl.sum()),
                     "pooled_definition": "tl_usable AND NOT floor-pinned, both hosts (ALREADY floor-corrected)",
                     "coactive_definition": "above floor value in both hosts"})
    return rows, floor_info


def load_cellfree(path):
    df = pd.read_parquet(path)
    out = {}
    for s in SPECIES:
        tx = pd.to_numeric(df[f"{s}_tx"], errors="coerce")
        usable = df[f"{s}_usable"].astype(bool)
        reason = df[f"{s}_unusable_reason"].astype("string").fillna("")
        no_rna = reason == "no_RNA_counts"
        activity = tx.copy()
        activity[no_rna] = 0.0
        out[s] = {"tx": tx, "usable": usable, "no_rna": no_rna,
                  "activity": activity, "dna_adequate": usable | no_rna}
    return out


def cellfree_rows(cols, setname):
    rows = []
    for a, b in itertools.combinations(SPECIES, 2):
        A, B = cols[a], cols[b]
        pooled = (A["dna_adequate"] & B["dna_adequate"])
        codet = (A["usable"] & B["usable"])
        rows.append({"dataset": "cell_free_drafts", "set": setname,
                     "pair": f"{a}-{b}", "readout": "transcription",
                     "n_pooled": int(pooled.sum()), "n_coactive": int(codet.sum()),
                     "pooled_definition": "DNA-adequate in both (usable OR no_RNA_counts); no_RNA carried as activity 0",
                     "coactive_definition": "DRAFTS-usable in both (= RNA detected in both)"})
    return rows


def verify_regimes(cols):
    """TASK 2: what do matched_per_host_active_frac and matched_26pct
    ACTUALLY do? Recomputed independently of scripts/98."""
    A, B = cols["Ec"], cols["Bs"]
    per_host = {"Ec": 15040 / 24613, "Bs": 4435 / 15848}
    out = {}

    # matched_per_host_active_frac
    fa, fb = per_host["Ec"], per_host["Bs"]
    qa = A["activity"][A["dna_adequate"]].quantile(1 - fa)
    qb = B["activity"][B["dna_adequate"]].quantile(1 - fb)
    m = A["dna_adequate"] & B["dna_adequate"] & (A["activity"] >= qa) & (B["activity"] >= qb)
    n_pool = int((A["dna_adequate"] & B["dna_adequate"]).sum())
    out["matched_per_host_active_frac"] = {
        "population_thresholds_computed_over": "each species' DNA-adequate (pooled-equivalent) activity distribution",
        "Ec_fraction_retained_target": fa, "Ec_threshold": float(qa),
        "Bs_fraction_retained_target": fb, "Bs_threshold": float(qb),
        "comparison_operator": ">= (inclusive)",
        "n_after_intersection": int(m.sum()),
        "achieved_retention_of_pooled": int(m.sum()) / n_pool,
        "source_of_fractions": "EC 15040/24613 and BS 4435/15848 = each host's own in-vivo tx_active/tx_usable",
    }

    # matched_26pct
    qa2 = A["activity"][A["dna_adequate"]].quantile(1 - INVIVO_COACTIVE_RETENTION)
    qb2 = B["activity"][B["dna_adequate"]].quantile(1 - INVIVO_COACTIVE_RETENTION)
    m2 = A["dna_adequate"] & B["dna_adequate"] & (A["activity"] >= qa2) & (B["activity"] >= qb2)
    out["matched_26pct"] = {
        "population_thresholds_computed_over": "each species' DNA-adequate (pooled-equivalent) activity distribution",
        "uniform_fraction_applied_to_both": INVIVO_COACTIVE_RETENTION,
        "Ec_threshold": float(qa2), "Bs_threshold": float(qb2),
        "comparison_operator": ">= (inclusive)",
        "n_after_intersection": int(m2.sum()),
        "achieved_retention_of_pooled": int(m2.sum()) / n_pool,
        "source_of_fraction": "in-vivo EC-BS transcription co-active retention, 3668/14088",
    }

    # the other two regimes, for the inconsistency note
    ta = A["tx"][A["usable"]].quantile(.25)
    out["_note_p25_med_regimes"] = {
        "population_thresholds_computed_over": "each species' USABLE (co-detected-eligible) distribution -- NOT DNA-adequate",
        "comparison_operator": "> (strict)",
        "why_it_matters": ("the two matched regimes use DNA-adequate + '>=', the p25/median regimes use "
                            "usable + '>' -- defensible per-regime choices, but the six rows are not all "
                            "on identical footing and a reader comparing them should know"),
        "example_Ec_p25_threshold_on_usable": float(ta),
    }
    return out


def main():
    RESULTS.mkdir(parents=True, exist_ok=True)
    three = pd.read_parquet(DATA / "three_host.parquet")
    rows, floor_info = invivo_rows(three)
    full = load_cellfree(DATA / "drafts.parquet")
    ovl = load_cellfree(DATA / "drafts_johns_overlap.parquet")
    rows += cellfree_rows(full, "drafts_full_1047")
    rows += cellfree_rows(ovl, "drafts_johns_overlap_112")

    df = pd.DataFrame(rows)
    df["coactive_retention_rate"] = df.n_coactive / df.n_pooled
    df["coactive_retention_pct"] = (df.coactive_retention_rate * 100).round(2)

    # ---- cross-check recomputation against the committed outputs --------
    print("=" * 78)
    print("CROSS-CHECK: recomputed-from-primary vs. committed script outputs")
    print("=" * 78)
    stored_iv = pd.read_csv(RESULTS / "gate8_5_coactive_correlation.csv")
    for _, s in stored_iv.iterrows():
        pair = s["pair"].replace("_", "-")
        mine = df[(df.dataset == "in_vivo_johns") & (df.pair == pair) & (df.readout == s["readout"])]
        m = mine.iloc[0]
        ok = (int(m.n_pooled) == int(s.n_pooled)) and (int(m.n_coactive) == int(s.n_coactive))
        print(f"  in-vivo {pair:6} {s['readout']:14} stored({s.n_pooled},{s.n_coactive}) "
              f"recomputed({m.n_pooled},{m.n_coactive})  {'MATCH' if ok else 'MISMATCH'}")
        if not ok:
            problems.append(f"in-vivo {pair} {s['readout']}: stored vs recomputed n disagree")

    stored_cf = pd.read_csv(RESULTS / "gate15_coactive_modality.csv")
    sc = stored_cf[stored_cf.set == "full_shared_1047"]
    nmis = 0
    for _, s in sc.iterrows():
        mine = df[(df.dataset == "cell_free_drafts") & (df.set == "drafts_full_1047") & (df.pair == s["pair"])]
        m = mine.iloc[0]
        if int(m.n_pooled) != int(s.pooled_equiv_n) or int(m.n_coactive) != int(s.co_detected_n):
            nmis += 1
            problems.append(f"cell-free {s['pair']}: stored vs recomputed n disagree")
    print(f"  cell-free 45 pairs (full set): {45 - nmis}/45 MATCH, {nmis} mismatches")

    # ---- TASK 1: the two headline figures -------------------------------
    print("\n" + "=" * 78)
    print("TASK 1: the two headline retention figures")
    print("=" * 78)
    cf = df[(df.dataset == "cell_free_drafts") & (df.set == "drafts_full_1047") & (df.pair == "Ec-Bs")].iloc[0]
    iv = df[(df.dataset == "in_vivo_johns") & (df.pair == "EC-BS") & (df.readout == "transcription")].iloc[0]
    print(f"  cell-free EC-BS: {cf.n_coactive} / {cf.n_pooled} = {cf.coactive_retention_pct:.2f}%   "
          f"(claimed 93.6%, counts 807 of 862)")
    print(f"  in-vivo   EC-BS: {iv.n_coactive} / {iv.n_pooled} = {iv.coactive_retention_pct:.2f}%   "
          f"(claimed 26.0%, counts 3,668 of 14,088)")
    claim_ok = {
        "cell_free_93_6": bool(int(cf.n_coactive) == 807 and int(cf.n_pooled) == 862
                                and abs(cf.coactive_retention_pct - 93.6) < 0.05),
        "in_vivo_26_0": bool(int(iv.n_coactive) == 3668 and int(iv.n_pooled) == 14088
                              and abs(iv.coactive_retention_pct - 26.0) < 0.05),
    }
    for k, v in claim_ok.items():
        print(f"    {k}: {'CONFIRMED' if v else 'WRONG'}")
        if not v:
            problems.append(f"headline retention figure {k} does not hold")

    print("\n  All pairs/readouts, co-active retention rate:")
    for ds, setname in [("in_vivo_johns", "full_library_29249"), ("cell_free_drafts", "drafts_full_1047")]:
        sub = df[(df.dataset == ds) & (df.set == setname)]
        if ds == "cell_free_drafts":
            print(f"    cell-free (45 pairs): min {sub.coactive_retention_pct.min():.1f}%  "
                  f"median {sub.coactive_retention_pct.median():.1f}%  max {sub.coactive_retention_pct.max():.1f}%")
            lo = sub.nsmallest(3, "coactive_retention_pct")
            for _, r in lo.iterrows():
                print(f"       lowest: {r['pair']:8} {r.n_coactive}/{r.n_pooled} = {r.coactive_retention_pct:.1f}%")
        else:
            for _, r in sub.iterrows():
                print(f"    in-vivo {r['pair']:6} {r['readout']:14} {r.n_coactive:6}/{r.n_pooled:6} = "
                      f"{r.coactive_retention_pct:6.2f}%")

    print("\n  Translation floor detection (why translation 'pooled' is already restricted):")
    for h, fi in floor_info.items():
        print(f"    {h}: floor={fi['floor_value']:.2f} pinning {fi['floor_fraction']*100:.1f}% "
              f"-> floor_detected={fi['floor_detected']}")

    # ---- TASK 2 ---------------------------------------------------------
    print("\n" + "=" * 78)
    print("TASK 2: what the two matched regimes actually do (recomputed)")
    print("=" * 78)
    reg = verify_regimes(full)
    for name in ["matched_per_host_active_frac", "matched_26pct"]:
        print(f"\n  {name}:")
        for k, v in reg[name].items():
            print(f"    {k}: {v}")
    print(f"\n  {'_note_p25_med_regimes'}:")
    for k, v in reg["_note_p25_med_regimes"].items():
        print(f"    {k}: {v}")

    df.to_csv(RESULTS / "gate18_coactive_retention.csv", index=False)
    with open(RESULTS / "gate18_regime_verification.json", "w") as f:
        json.dump({"headline_claims": claim_ok, "regimes": reg, "translation_floor": floor_info,
                    "problems": problems}, f, indent=2)
    print(f"\nWrote {RESULTS/'gate18_coactive_retention.csv'} ({len(df)} rows), "
          f"{RESULTS/'gate18_regime_verification.json'}")
    print(f"\nPROBLEMS: {len(problems)}" + ("" if not problems else " -> " + "; ".join(problems)))


if __name__ == "__main__":
    main()
