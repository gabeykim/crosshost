"""
GATE 10 - Task 4.3 [LOAD-BEARING]: is the cell-free host effect the same
object as the in-vivo host effect?

TWO INDEPENDENT COMPARISONS:

(A) This project's own computation: on the 112 DRAFTS-Johns overlapping
sequences (Ec/Bs only -- P. aeruginosa is NOT in DRAFTS, see scripts/89),
compute the cell-free (DRAFTS) Ec-Bs cross-host correlation and compare it
directly against Johns's own in-vivo Ec-Bs correlation, ON THE SAME
SEQUENCE SET where possible.

(B) The paper's OWN published RS234 in-vitro-vs-in-vivo comparison
(Appendix Fig. S5 / Appendix Table S6) -- reported as published values,
clearly labeled as not this project's computation. This is DRAFTS's own
in-vivo re-measurement for 7 species (Ec, Se, Pp, Vn, Ko, Bs, Cg), run by
the DRAFTS authors specifically to validate their cell-free system --
NOT a reuse of Johns et al. 2018's original in-vivo numbers. This
distinction matters: comparison (A) checks DRAFTS-cell-free against
JOHNS'S ORIGINAL in-vivo measurement; comparison (B) checks DRAFTS's OWN
cell-free against DRAFTS'S OWN in-vivo re-measurement, a same-lab,
same-timepoint pair that is a cleaner apples-to-apples check but does not
by itself validate against Johns 2018's numbers.
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)


def comparison_a_own_computation():
    """DRAFTS cell-free Ec-Bs correlation vs Johns in-vivo Ec-Bs correlation,
    on the DRAFTS-Johns overlap (the only pair with both hosts present in
    both datasets -- P. aeruginosa is absent from DRAFTS)."""
    joined = pd.read_parquet(DATA / "drafts_johns_overlap.parquet")
    three_host = pd.read_parquet(DATA / "three_host.parquet")

    # cell-free (DRAFTS), on the 112-sequence overlap
    both_usable_cf = joined["Ec_usable"] & joined["Bs_usable"]
    n_cf = int(both_usable_cf.sum())
    rho_cf = float(spearmanr(joined.loc[both_usable_cf, "Ec_tx"], joined.loc[both_usable_cf, "Bs_tx"]).correlation) if n_cf > 2 else None

    # in-vivo (Johns), on the SAME 112-sequence overlap, matched via johns_oligo_id
    johns_indexed = three_host.set_index(three_host["OLIGO ID"].astype(int))
    matched = johns_indexed.reindex(joined["johns_oligo_id"].dropna().astype(int))
    both_usable_vivo = matched["EC_tx_usable"] & matched["EC_tx_active"] & matched["BS_tx_usable"] & matched["BS_tx_active"]
    n_vivo = int(both_usable_vivo.sum())
    rho_vivo = float(spearmanr(matched.loc[both_usable_vivo, "EC_tx_norm"],
                                matched.loc[both_usable_vivo, "BS_tx_norm"]).correlation) if n_vivo > 2 else None

    # in-vivo (Johns), on the FULL three-host library (already-established project number, for reference)
    full_both = (three_host["EC_tx_usable"] & three_host["EC_tx_active"] &
                 three_host["BS_tx_usable"] & three_host["BS_tx_active"])
    n_vivo_full = int(full_both.sum())
    rho_vivo_full = float(spearmanr(three_host.loc[full_both, "EC_tx_norm"],
                                     three_host.loc[full_both, "BS_tx_norm"]).correlation)

    result = {
        "comparison": "A: own computation, DRAFTS cell-free vs Johns in-vivo, EC-BS pair (only pair with both hosts in both datasets)",
        "cell_free_DRAFTS_EC_BS": {"n": n_cf, "spearman_rho": rho_cf},
        "in_vivo_Johns_EC_BS_on_same_112seq_overlap": {"n": n_vivo, "spearman_rho": rho_vivo},
        "in_vivo_Johns_EC_BS_full_library_reference": {"n": n_vivo_full, "spearman_rho": rho_vivo_full},
        "caveat": "n is small (the DRAFTS-Johns sequence overlap is only 112 total, ~80-95 usable per pair) "
                  "-- read the comparison directionally, not as a precise estimate.",
    }
    print("=== COMPARISON A: own computation ===")
    print(json.dumps(result, indent=2, default=str))
    return result


RS234_SPECIES = ["Ec", "Se", "Pp", "Vn", "Ko", "Bs", "Cg"]


def comparison_b_published_and_recomputed_rs234():
    """Two related but distinct numbers, both reported, clearly labeled:

    (b1) DIRECTLY TRANSCRIBED published values: Appendix Table S6
    "Reproducibility between biological replicates" gives Pearson r
    between biological replicate 1 and 2 for E. coli in vivo (RS29249,
    Figure 1D row) = 0.8932, and for DRAFTS E. coli cell-free (same row)
    = 0.9333. This is REPLICATE reproducibility (rep1 vs rep2 within the
    same modality), NOT an in-vitro-vs-in-vivo cross-modality correlation
    -- worth having regardless, since it is a second, independent E. coli
    reliability estimate (Gate 8's own estimate, from 5 growth conditions,
    was 0.912/0.929) that arrives via true biological replicates instead.

    (b2) NOT found as a clean extractable table: the actual RS234
    in-vitro-vs-in-vivo correlation coefficients (as opposed to replicate
    reproducibility) appear only as annotations on the Appendix Figure S5
    scatter-plot images themselves -- pdftotext cannot extract text baked
    into a plot image, and no separate results-text table states these
    numbers directly for all 7 species. RECOMPUTED INSTEAD, directly from
    DRAFTS's own released Appendix Fig. S5 source-data sheet
    (data/drafts_rs234_invivo.parquet, parsed in scripts/88) -- this is
    computed by this project, not a literal quote from the paper, but it
    uses exactly the data DRAFTS released for this specific comparison, so
    it should closely track whatever value the figure's on-image
    annotation shows. Labeled honestly as recomputed-from-their-data, not
    published-value, throughout.
    """
    b1_published = {
        "source": "raw/drafts/msb198875_Appendix.pdf, Appendix Table S6, 'Figure 1D' rows -- TRANSCRIBED DIRECTLY, not recomputed",
        "metric": "Pearson r between biological replicate 1 and replicate 2 (NOT in-vitro-vs-in-vivo)",
        "E_coli_in_vivo_RS29249": 0.8932,
        "E_coli_DRAFTS_cellfree_RS29249": 0.9333,
        "relevance": "A second, independent E. coli reliability estimate via true biological replicates, "
                     "comparable in spirit to Gate 8's growth-condition-based estimate (0.912/0.929) but not the "
                     "same measurement -- both point to E. coli transcription reliability in the low-to-mid 0.9s.",
    }

    df = pd.read_parquet(DATA / "drafts_rs234_invivo.parquet")
    b2_recomputed = {}
    for sp in RS234_SPECIES:
        both = df[f"{sp}_invitro_usable"] & df[f"{sp}_invivo_usable"]
        n = int(both.sum())
        rho = float(spearmanr(df.loc[both, f"{sp}_invitro_tx"], df.loc[both, f"{sp}_invivo_tx"]).correlation) if n > 2 else None
        b2_recomputed[sp] = {"n": n, "spearman_rho_invitro_vs_invivo": rho}
        print(f"  RS234 {sp}: in-vitro vs in-vivo, n={n}, Spearman rho={rho:.3f}" if rho is not None else f"  RS234 {sp}: insufficient data")

    result = {
        "b1_published_replicate_reproducibility": b1_published,
        "b2_recomputed_from_drafts_own_released_data_NOT_a_published_quote": {
            "source": "data/drafts_rs234_invivo.parquet (parsed from msb198875_SourceData_Fig2.xlsx, sheet 'Appendix Fig. S5')",
            "metric": "Spearman rho, DRAFTS in-vitro (cell-free) tx vs DRAFTS in-vivo tx, same 234-sequence RS234 library, per species",
            "per_species": b2_recomputed,
            "note_pa_not_tested": "P. aeruginosa is not in DRAFTS at all (see scripts/89); it is also not one of "
                                    "RS234's 7 in-vivo comparison species regardless.",
        },
    }
    print("\n=== COMPARISON B: published replicate reliability + recomputed RS234 in-vitro-vs-in-vivo ===")
    print(json.dumps(result, indent=2, default=str))
    return result


def main():
    a = comparison_a_own_computation()
    b = comparison_b_published_and_recomputed_rs234()
    with open(RESULTS / "gate10_modality_comparison.json", "w") as f:
        json.dump({"comparison_a_own_computation": a, "comparison_b_published_and_recomputed_rs234": b}, f, indent=2, default=str)
    print(f"\nWrote {RESULTS / 'gate10_modality_comparison.json'}")


if __name__ == "__main__":
    main()
