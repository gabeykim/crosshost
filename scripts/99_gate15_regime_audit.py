"""
GATE 15 - Task 4: restriction-regime audit.

For every correlation figure that appears in out/PREPRINT/MANUSCRIPT.md,
record which restriction regime produced it and its n, so the manuscript
can be rewritten knowing which numbers are which. The manuscript currently
reports pooled and co-active figures side by side without labelling them.

Values are READ FROM THE RESULT FILES, not retyped, wherever a result file
holds them -- the point of this table is to be checkable.

REGIME VOCABULARY (the two datasets collide on the word "usable"):
  pooled        in-vivo: DNA template present in both hosts; rna==0 rows
                KEPT with tx_norm exactly 0
  co-active     in-vivo: pooled AND rna>0 in both hosts
  co-detected   cell-free: DRAFTS-"usable" in both species, which already
                REQUIRES rna>0 -- definitionally the co-active operation,
                despite the column being named n_shared_usable
  n/a           model-performance correlations (prediction vs. truth), not
                cross-host measurement agreement -- listed so the table is
                exhaustive, not because a restriction regime applies
"""
import csv
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "out" / "results"

co = pd.read_csv(RESULTS / "gate8_5_coactive_correlation.csv").set_index(["pair", "readout"])
gc = json.load(open(RESULTS / "gate10_5_gc_confound.json"))
dr = pd.read_csv(RESULTS / "gate10_crosshost_correlations.csv")
mod = json.load(open(RESULTS / "gate10_modality_comparison.json"))
g15 = json.load(open(RESULTS / "gate15_coactive_modality.json"))

rows = []


def add(section, figure, value, regime, n, source, note=""):
    rows.append({"manuscript_section": section, "figure": figure, "value": value,
                 "restriction_regime": regime, "n": n, "source_file": source, "note": note})


# ---- Section 3.1: in-vivo cross-host measurement agreement ----------------
for pair, rd, label in [("EC_PA", "transcription", "EC-PA transcription"),
                        ("EC_BS", "transcription", "EC-BS transcription"),
                        ("BS_PA", "transcription", "BS-PA transcription"),
                        ("EC_PA", "translation", "EC-PA translation"),
                        ("EC_BS", "translation", "EC-BS translation"),
                        ("BS_PA", "translation", "BS-PA translation")]:
    r = co.loc[(pair, rd)]
    add("3.1", f"{label} (headline cross-host correlation)", round(float(r.rho_coactive), 4),
        "co-active", int(r.n_coactive), "gate8_5_coactive_correlation.csv")
    add("3.1 / 3.6", f"{label} POOLED counterpart", round(float(r.rho_pooled), 4),
        "pooled", int(r.n_pooled), "gate8_5_coactive_correlation.csv",
        "reported in 3.6 as the pre-restriction value; NOT reported in 3.1")

# GC-controlled versions (Gate 10.5) -- computed on the same co-active masks
for rd in ["transcription", "translation"]:
    for pair in ["EC_BS", "EC_PA", "BS_PA"]:
        d = gc["task2_raw_vs_gc_controlled_correlations"][rd][pair]
        add("3.1", f"{pair.replace('_','-')} {rd} GC-CONTROLLED",
            round(d["rho_gc_controlled_formula"], 4), "co-active", d["n"],
            "gate10_5_gc_confound.json",
            "GC control applied to the co-active set; n identical to the co-active n above")

for rd in ["transcription", "translation"]:
    v = gc["task3_verdict"][rd]
    add("3.1", f"gap ratio {rd} (raw)", round(v["gap_ratio_raw"], 4), "co-active (derived)",
        "-", "gate10_5_gc_confound.json", "ratio of co-active correlations")
    add("3.1", f"gap ratio {rd} (GC-controlled)", round(v["gap_ratio_gc_controlled"], 4),
        "co-active (derived)", "-", "gate10_5_gc_confound.json", "")

for pair, phy in [("EC_PA", "Proteobacteria"), ("EC_PA", "Firmicutes")]:
    d = gc["task4_phylum_stratified"]["transcription"][pair][phy]
    add("3.1", f"{pair.replace('_','-')} transcription, {phy} stratum",
        round(d["spearman_rho"], 4), "co-active", d["n"], "gate10_5_gc_confound.json", "")

add("3.1", "disattenuated gap ratio @ reliability 0.9", 0.341, "co-active (derived)", "-",
    "gate8_attenuation_analysis.json", "disattenuation applied to co-active correlations")

# ---- Section 3.2: DRAFTS cell-free ---------------------------------------
add("3.2", "DRAFTS 45-pair range (min)", round(float(dr.spearman_rho.min()), 4),
    "co-detected", int(dr.loc[dr.spearman_rho.idxmin(), "n_shared_usable"]),
    "gate10_crosshost_correlations.csv",
    "column is named n_shared_usable, but DRAFTS-usable REQUIRES rna>0 -- this is the co-active operation")
add("3.2", "DRAFTS 45-pair range (max)", round(float(dr.spearman_rho.max()), 4),
    "co-detected", int(dr.loc[dr.spearman_rho.idxmax(), "n_shared_usable"]),
    "gate10_crosshost_correlations.csv", "")
add("3.2", "DRAFTS same-phylum mean", 0.852, "co-detected", "17 pairs",
    "gate10_crosshost_correlations_summary.json", "")
add("3.2", "DRAFTS cross-phylum mean", 0.769, "co-detected", "28 pairs",
    "gate10_crosshost_correlations_summary.json", "")
add("3.2", "B. subtilis mean correlation with other 9 species", 0.783, "co-detected", "9 pairs",
    "gate10_crosshost_correlations_summary.json", "")

a = mod["comparison_a_own_computation"]
add("3.2", "HEADLINE cell-free EC-BS", round(a["cell_free_DRAFTS_EC_BS"]["spearman_rho"], 4),
    "co-detected", a["cell_free_DRAFTS_EC_BS"]["n"], "gate10_modality_comparison.json",
    "Johns-overlap subset; script 91 filter is Ec_usable & Bs_usable = RNA detected in both")
add("3.2", "HEADLINE in-vivo EC-BS (full-library reference)",
    round(a["in_vivo_Johns_EC_BS_full_library_reference"]["spearman_rho"], 4), "co-active",
    a["in_vivo_Johns_EC_BS_full_library_reference"]["n"], "gate10_modality_comparison.json",
    "script 91 filter is tx_usable & tx_active in both hosts")
add("3.2", "in-vivo EC-BS on the 15-seq overlap (reported, not used as evidence)",
    round(a["in_vivo_Johns_EC_BS_on_same_112seq_overlap"]["spearman_rho"], 4), "co-active",
    a["in_vivo_Johns_EC_BS_on_same_112seq_overlap"]["n"], "gate10_modality_comparison.json", "")

b2 = mod["comparison_b_published_and_recomputed_rs234"]["b2_recomputed_from_drafts_own_released_data_NOT_a_published_quote"]["per_species"]
lo = min(v["spearman_rho_invitro_vs_invivo"] for v in b2.values())
hi = max(v["spearman_rho_invitro_vs_invivo"] for v in b2.values())
add("3.2", "RS234 within-species in-vitro vs in-vivo (range across 7 species)",
    f"{lo:.3f}-{hi:.3f}", "co-detected in BOTH modalities",
    f"{min(v['n'] for v in b2.values())}-{max(v['n'] for v in b2.values())}",
    "gate10_modality_comparison.json",
    "RS234 uses the same DRAFTS sentinel convention on both arms, so both sides exclude no-RNA rows -- internally regime-matched")

# ---- Gate 15's newly computed regime-matched figures ----------------------
for regime, d in g15["task2_ec_bs"]["full_shared_1047"].items():
    if d["rho"] is not None:
        add("NEW (Gate 15, not yet in manuscript)", f"cell-free EC-BS, {regime}, full shared set",
            round(d["rho"], 4), regime, d["n"], "gate15_coactive_modality.json", "")

# ---- Section 3.6: restriction effect, explicitly reported as a pair -------
add("3.6", "EC-BS transcription pooled -> co-active", "0.655 -> 0.258", "both, explicitly labelled",
    "14088 -> 3668", "gate8_5_coactive_correlation.csv", "the one place the manuscript does distinguish")
add("3.6", "BS-PA transcription pooled -> co-active", "0.508 -> 0.257", "both, explicitly labelled",
    "11969 -> 2099", "gate8_5_coactive_correlation.csv", "")
add("3.6", "EC-PA transcription pooled -> co-active", "0.621 -> 0.754", "both, explicitly labelled",
    "19643 -> 9741", "gate8_5_coactive_correlation.csv", "")

# ---- model-performance correlations: listed for completeness -------------
for fig, val in [("H-MAIN closest call, model", 0.213), ("H-MAIN closest call, baseline", 0.218),
                 ("sequence-only EC transcription", 0.367), ("concat EC transcription", 0.555),
                 ("per-host-heads-avg EC transcription", 0.615),
                 ("shift-prediction original EC->BS", 0.432),
                 ("shift-prediction shuffled-sequence control", 0.557)]:
    add("3.4 / 3.7", fig, val, "n/a - model prediction vs truth", "-", "various",
        "not a cross-host measurement correlation; no restriction regime applies")

out = RESULTS / "gate15_restriction_regime_audit.csv"
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["manuscript_section", "figure", "value",
                                       "restriction_regime", "n", "source_file", "note"])
    w.writeheader()
    w.writerows(rows)

df = pd.DataFrame(rows)
print(df[["manuscript_section", "figure", "value", "restriction_regime", "n"]].to_string(index=False))
print(f"\n{len(rows)} figures audited -> {out}")
print("\nRegime counts:")
print(df.restriction_regime.value_counts().to_string())
