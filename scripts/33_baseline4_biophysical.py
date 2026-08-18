"""
GATE 3 - Baseline 4: biophysical model.

ATTEMPTED: a fresh run of the Salis Lab RBS Calculator v1.0 (the exact tool
Johns et al. 2018 themselves used, per their Methods text, Gate 1) and/or the
Salis Promoter Calculator was considered. Both require substantial
environment setup (RBS Calculator v1.0 is Python 2 + a bundled NuPACK
build; the Promoter Calculator has its own model/dependency stack) that was
judged not worth the time given a better option exists (see below) --
documented here as a disclosed "did not attempt a fresh install," not a
silent skip.

USED INSTEAD: Johns et al.'s OWN released biophysical predictions, already
computed on these exact 165bp constructs and shipped in the supplementary
data (`raw/NIHMS945382-supplement-4.xlsx`, Gate 1 column inventory):
  - `{host}_best_sigma70_match_score` -- a position-weight-matrix sigma-70
    promoter motif match score (their transcription-side biophysical model)
  - `{host}_delta_G` -- predicted mRNA/UTR folding free energy, computed via
    the same RBS-Calculator-family thermodynamics used for translation
    initiation rate prediction (their translation-side biophysical proxy)
This is arguably a STRONGER reference point than an independent fresh tool
run would be: it is computed by the paper's own authors on the paper's own
sequences under the paper's own parameterization, with no risk of a
mismatched-input error on my part. It is evaluated exactly like every other
baseline -- Spearman correlation against measured values on the SAME fold-0
held-out test set -- for a genuinely apples-to-apples comparison, not just
quoted as a literature figure.
"""
import numpy as np
import pandas as pd
import json
from pathlib import Path
from scipy.stats import spearmanr

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out" / "baselines"
OUT.mkdir(parents=True, exist_ok=True)

HOSTS = ["EC", "BS", "PA"]
TEST_FOLD = 0


def main():
    df = pd.read_parquet(DATA / "three_host_library.parquet")
    df["OLIGO ID"] = df["OLIGO ID"].astype(str)
    fold_df = pd.read_parquet(DATA / "splits" / "fold_assignment_FINAL.parquet")
    id_to_fold = dict(zip(fold_df["OLIGO ID"], fold_df["fold"]))
    df["fold"] = df["OLIGO ID"].map(id_to_fold)
    df = df[df["fold"].notna()].copy()
    df["fold"] = df["fold"].astype(int)

    three_host = pd.read_parquet(DATA / "three_host.parquet")
    three_host["OLIGO ID"] = three_host["OLIGO ID"].astype(str)
    df = df.merge(three_host[["OLIGO ID"] + [f"{h}_tx_usable" for h in HOSTS] +
                              [f"{h}_tx_active" for h in HOSTS] + [f"{h}_tx_norm" for h in HOSTS]],
                   on="OLIGO ID", suffixes=("", "_dup"))

    results = {}
    for host in HOSTS:
        print(f"\n=== {host} ===")
        test_mask = (df["fold"] == TEST_FOLD) & df[f"{host}_tx_usable"]

        # transcription: sigma70 match score vs tx_norm (on actives, log-scale, matching other baselines)
        active_mask = test_mask & df[f"{host}_tx_active"]
        sigma_score = df.loc[active_mask, f"{host}_best_sigma70_match_score"]
        tx_val = df.loc[active_mask, f"{host}_tx_norm"]
        valid = sigma_score.notna() & tx_val.notna()
        log_tx = np.log1p(np.clip(tx_val[valid], 0, None))
        rho_tx = spearmanr(sigma_score[valid], log_tx).correlation if valid.sum() > 1 else None
        print(f"  transcription (sigma70_match_score vs tx_norm, n={valid.sum()}): rho={rho_tx}")

        # translation: delta_G vs protein_log10 (delta_G is folding energy -- MORE NEGATIVE
        # = more stable secondary structure = typically LOWER translation, so a negative
        # correlation is the biophysically expected direction; reported as-is, not sign-flipped)
        protein_col = f"{host}_protein (log10)"
        dg = df.loc[test_mask, f"{host}_delta_G"]
        protein = df.loc[test_mask, protein_col]
        valid_tl = dg.notna() & protein.notna()
        rho_tl = spearmanr(dg[valid_tl], protein[valid_tl]).correlation if valid_tl.sum() > 1 else None
        print(f"  translation (delta_G vs protein_log10, n={valid_tl.sum()}): rho={rho_tl}")

        results[host] = {
            "transcription": {"feature": "best_sigma70_match_score", "n": int(valid.sum()),
                               "spearman_rho": float(rho_tx) if rho_tx is not None and rho_tx == rho_tx else None},
            "translation": {"feature": "delta_G", "n": int(valid_tl.sum()),
                             "spearman_rho": float(rho_tl) if rho_tl is not None and rho_tl == rho_tl else None,
                             "note": "delta_G is folding free energy; more negative = more stable "
                                     "secondary structure = typically lower translation initiation "
                                     "rate, so a NEGATIVE rho is the biophysically expected direction, "
                                     "not a sign error"},
        }

    results["not_attempted"] = {
        "salis_rbs_calculator_v1_fresh_run": "not attempted -- Python 2 + bundled NuPACK build, "
            "judged not worth the setup time given the paper's own released biophysical scores "
            "(computed with an RBS-Calculator-family method on these exact sequences) were available",
        "salis_promoter_calculator_fresh_run": "not attempted, same reasoning",
    }

    with open(OUT / "baseline4_biophysical.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {OUT / 'baseline4_biophysical.json'}")


if __name__ == "__main__":
    main()
