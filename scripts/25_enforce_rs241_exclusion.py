"""
GATE 2 - Task 3.4: Enforce the charter's non-negotiable rule -- RS241 is never
trained on -- MECHANICALLY, not by discipline.

CRITICAL FINDING: 207 of the 241 RS241 oligo IDs are physically present in the
main 29,249-sequence three-host library (RS241 was "selected" from existing
library members per Johns et al.'s own methods text, Gate 1). Before this
fix, data/splits/fold_assignment_FINAL.parquet assigned these 207 oligos to
ordinary folds like any other library member -- meaning a naive "train on
everything not in the held-out test fold" loop would have trained on ~4/5 of
them for any given fold, directly violating the charter rule.

Fix: remove all 241 RS241 oligo IDs (not just the 207 that happen to overlap
by ID -- all of them, for defense in depth) from the main-library fold
assignment entirely. They are neither train nor test in the main three-host
benchmark; they exist ONLY in data/rs241.parquet's own six-host evaluation,
which is a completely separate evaluation track per the charter.
"""
import pandas as pd
from pathlib import Path
import json

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"
SPLITS = DATA / "splits"


def main():
    rs241 = pd.read_parquet(DATA / "rs241.parquet")
    rs241_ids = set(rs241["id"].astype(str))
    print(f"RS241 total ids: {len(rs241_ids)}")

    fold_df = pd.read_parquet(SPLITS / "fold_assignment_FINAL.parquet")
    before_n = len(fold_df)
    overlap = set(fold_df["OLIGO ID"]) & rs241_ids
    print(f"RS241 ids physically present in main-library fold assignment BEFORE fix: {len(overlap)}")

    fold_df_clean = fold_df[~fold_df["OLIGO ID"].isin(rs241_ids)].copy()
    after_n = len(fold_df_clean)
    print(f"Fold assignment rows: {before_n} -> {after_n} (removed {before_n - after_n})")

    remaining_overlap = set(fold_df_clean["OLIGO ID"]) & rs241_ids
    assert len(remaining_overlap) == 0, f"RS241 exclusion FAILED: {len(remaining_overlap)} still present"
    print("VERIFIED: zero RS241 ids remain in the main-library fold assignment")

    fold_df_clean.to_parquet(SPLITS / "fold_assignment_FINAL.parquet", index=False)

    fold_loads = fold_df_clean["fold"].value_counts().sort_index()
    print(f"\nFinal fold loads after RS241 removal:\n{fold_loads.to_string()}")

    with open(OUT / "rs241_exclusion_fix.json", "w") as f:
        json.dump({
            "rs241_total_ids": len(rs241_ids),
            "rs241_ids_found_in_main_library_before_fix": len(overlap),
            "rows_removed": before_n - after_n,
            "final_fold_loads": fold_loads.to_dict(),
            "verification": "zero RS241 ids remain in fold_assignment_FINAL.parquet, confirmed by set intersection",
        }, f, indent=2, default=str)
    print(f"\nWrote {OUT / 'rs241_exclusion_fix.json'}")


if __name__ == "__main__":
    main()
