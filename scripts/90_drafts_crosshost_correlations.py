"""
GATE 10 - Task 4.1/4.2: cross-host measurement correlation, all 45 DRAFTS
species pairs, on the shared-usable subset per pair (Spearman) -- the
direct cell-free analog of Johns et al.'s in-vivo EC-PA/BS-pairs contrast
this project's whole central claim is built on.
"""
import json
import itertools
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

SPECIES = ["Ec", "Ef", "Se", "Ko", "Pa", "Pp", "Vn", "Bs", "Cg", "Ll"]
SPECIES_PHYLUM = {
    "Ec": "Proteobacteria", "Ef": "Proteobacteria", "Se": "Proteobacteria", "Ko": "Proteobacteria",
    "Pa": "Proteobacteria", "Pp": "Proteobacteria", "Vn": "Proteobacteria",
    "Bs": "Firmicutes", "Ll": "Firmicutes", "Cg": "Actinobacteria",
}
# Rough phylogenetic ordering within Proteobacteria (Enterobacterales cluster first,
# then more distant Gammaproteobacteria/others), then Firmicutes, then Actinobacteria --
# for a phylogeny-ordered correlation-matrix figure. Not a rigorous tree, just an
# ordering convenient for visually checking clade structure.
PHYLO_ORDER = ["Ec", "Ef", "Se", "Ko", "Pa", "Pp", "Vn", "Bs", "Ll", "Cg"]


def main():
    df = pd.read_parquet(DATA / "drafts.parquet")

    rows = []
    matrix = pd.DataFrame(index=SPECIES, columns=SPECIES, dtype=float)
    n_matrix = pd.DataFrame(index=SPECIES, columns=SPECIES, dtype=int)

    for s1, s2 in itertools.combinations(SPECIES, 2):
        both_usable = df[f"{s1}_usable"] & df[f"{s2}_usable"]
        n = int(both_usable.sum())
        v1 = df.loc[both_usable, f"{s1}_tx"]
        v2 = df.loc[both_usable, f"{s2}_tx"]
        rho = float(spearmanr(v1, v2).correlation) if n > 2 else None
        same_phylum = SPECIES_PHYLUM[s1] == SPECIES_PHYLUM[s2]
        rows.append({"species_a": s1, "species_b": s2, "phylum_a": SPECIES_PHYLUM[s1],
                     "phylum_b": SPECIES_PHYLUM[s2], "same_phylum": same_phylum,
                     "n_shared_usable": n, "spearman_rho": rho})
        matrix.loc[s1, s2] = rho
        matrix.loc[s2, s1] = rho
        n_matrix.loc[s1, s2] = n
        n_matrix.loc[s2, s1] = n
    for s in SPECIES:
        matrix.loc[s, s] = 1.0
        n_matrix.loc[s, s] = int(df[f"{s}_usable"].sum())

    assert len(rows) == 45, f"expected 45 pairs, got {len(rows)}"
    result_df = pd.DataFrame(rows).sort_values("spearman_rho", ascending=False)
    result_df.to_csv(RESULTS / "gate10_crosshost_correlations.csv", index=False)
    matrix.to_csv(RESULTS / "gate10_crosshost_correlation_matrix.csv")
    matrix_ordered = matrix.loc[PHYLO_ORDER, PHYLO_ORDER]
    matrix_ordered.to_csv(RESULTS / "gate10_crosshost_correlation_matrix_phylo_ordered.csv")

    print("=" * 80)
    print("ALL 45 PAIRWISE CROSS-HOST CORRELATIONS (Spearman, shared-usable subset)")
    print("=" * 80)
    print(result_df.round(3).to_string(index=False))

    same_phylum_rhos = result_df[result_df.same_phylum]["spearman_rho"].dropna()
    diff_phylum_rhos = result_df[~result_df.same_phylum]["spearman_rho"].dropna()
    print(f"\nSame-phylum pairs: n={len(same_phylum_rhos)}, mean rho={same_phylum_rhos.mean():.3f}, "
          f"range [{same_phylum_rhos.min():.3f}, {same_phylum_rhos.max():.3f}]")
    print(f"Cross-phylum pairs: n={len(diff_phylum_rhos)}, mean rho={diff_phylum_rhos.mean():.3f}, "
          f"range [{diff_phylum_rhos.min():.3f}, {diff_phylum_rhos.max():.3f}]")

    # Outlier check: are Bs, Ll, Cg (non-Proteobacteria) outliers the way BS is in vivo?
    print("\n=== Per-species mean correlation with all OTHER species (outlier check) ===")
    outlier_rows = []
    for s in SPECIES:
        others = result_df[(result_df.species_a == s) | (result_df.species_b == s)]
        mean_rho = others["spearman_rho"].mean()
        outlier_rows.append({"species": s, "phylum": SPECIES_PHYLUM[s], "mean_rho_with_others": mean_rho})
        print(f"  {s} ({SPECIES_PHYLUM[s]}): mean rho with other 9 species = {mean_rho:.3f}")
    outlier_df = pd.DataFrame(outlier_rows).sort_values("mean_rho_with_others")

    summary = {
        "n_pairs": 45,
        "same_phylum_mean_rho": float(same_phylum_rhos.mean()),
        "cross_phylum_mean_rho": float(diff_phylum_rhos.mean()),
        "highest_pair": {"species": [result_df.iloc[0]["species_a"], result_df.iloc[0]["species_b"]],
                          "rho": float(result_df.iloc[0]["spearman_rho"])},
        "lowest_pair": {"species": [result_df.iloc[-1]["species_a"], result_df.iloc[-1]["species_b"]],
                         "rho": float(result_df.iloc[-1]["spearman_rho"])},
        "per_species_mean_rho_ranked_lowest_first": outlier_df.to_dict("records"),
        "phylogenetic_pattern_holds": bool(same_phylum_rhos.mean() > diff_phylum_rhos.mean()),
    }
    with open(RESULTS / "gate10_crosshost_correlations_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\nWrote {RESULTS / 'gate10_crosshost_correlations.csv'}, "
          f"{RESULTS / 'gate10_crosshost_correlation_matrix.csv'}, "
          f"{RESULTS / 'gate10_crosshost_correlations_summary.json'}")


if __name__ == "__main__":
    main()
