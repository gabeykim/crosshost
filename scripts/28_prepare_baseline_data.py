"""
GATE 3 - Shared data prep for all baselines.

Builds, per primary host (EC/BS/PA), a cached array bundle:
  - one-hot encoded sequences (N, 4, 165)
  - transcription: usable mask, active mask (rna_count>0), tx_norm (regression target)
  - translation: usable mask, active mask, protein_log10 (regression target)
  - fold assignment (0-4) from the frozen split

TRANSLATION FLOOR-VALUE FINDING (discovered while building this cache, not
assumed): protein_log10 is heavily pinned at a single repeated value for the
weaker-translating hosts -- 90.5% of B. subtilis's "usable" (non-null) values
round to exactly 2.00, 67.4% of E. coli's round to ~1.87, but only 9.9% of
P. aeruginosa's round to its mode. This tracks the known host ranking
(BS weakest translator < EC < PA strongest, Gate 1) and is almost certainly
the paper's own pseudo-value/floor convention for constructs with
insufficient FACS-seq signal (analogous to the documented transcription
pseudo-value scheme, Gate 1 Methods), not real quantitative measurements.

Fix applied: for each host, the floor value is detected as the most common
rounded (2dp) protein_log10 value IF it accounts for >5% of usable rows
(clearly non-physical for a continuous measurement; a real distribution mode
would not concentrate this much mass at one exact point). Two consequences:
  - tl_usable_for_regression EXCLUDES floor-pinned rows -- training a
    regressor against a target that is 90% identical constant values is not
    meaningful and would produce spuriously "good" but useless MSE.
  - translation ACTIVE is redefined as protein_log10 STRICTLY ABOVE the
    floor value (a principled "detectable above background" cutoff) rather
    than an arbitrary median split, which was the original (weaker) design
    of this script before the floor was found.
"""
import numpy as np
import pandas as pd
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data"
CACHE = DATA / "baseline_cache"
CACHE.mkdir(exist_ok=True)

HOSTS = ["EC", "BS", "PA"]
SEQ_LEN = 165
BASE_TO_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}
TEST_FOLD = 0  # matches every baseline script (30-34) -- the floor value below
                # must be fit on the training pool (fold != TEST_FOLD) ONLY, per
                # standing rule SR3 / audit_leakage.py check 5. An earlier version
                # of this script computed it over ALL folds including the test
                # fold -- a real (if narrow) fold-locality violation, caught and
                # fixed here rather than left as a known issue. See
                # out/GATE3_MEMO.md "WHAT I COULD NOT DO" / normalization check.


def one_hot(seq):
    arr = np.zeros((4, SEQ_LEN), dtype=np.float32)
    for i, b in enumerate(seq):
        idx = BASE_TO_IDX.get(b)
        if idx is not None:
            arr[idx, i] = 1.0
    return arr


def main():
    df = pd.read_parquet(DATA / "three_host.parquet")
    df["OLIGO ID"] = df["OLIGO ID"].astype(str)
    fold_df = pd.read_parquet(DATA / "splits" / "fold_assignment_FINAL.parquet")
    id_to_fold = dict(zip(fold_df["OLIGO ID"], fold_df["fold"]))
    df["fold"] = df["OLIGO ID"].map(id_to_fold)
    df = df[df["fold"].notna()].copy()  # drop RS241-excluded rows (not in fold assignment)
    df["fold"] = df["fold"].astype(int)
    print(f"Rows with a valid fold assignment: {len(df)} (RS241-excluded rows dropped)")

    for host in HOSTS:
        print(f"\n=== {host} ===")
        seqs = df["regulatory_sequence"].values
        onehot = np.stack([one_hot(s) for s in seqs]).astype(np.float32)

        tx_usable = df[f"{host}_tx_usable"].values.astype(bool)
        tx_active = df[f"{host}_tx_active"].values.astype(bool)
        tx_norm = df[f"{host}_tx_norm"].values.astype(np.float32)

        tl_usable = df[f"{host}_tl_usable"].values.astype(bool)
        protein_log10 = df[f"{host}_protein_log10"].values.astype(np.float32)
        fold_arr = df["fold"].values.astype(int)

        # detect floor/pseudo-value -- FOLD-LOCAL: fit using the training pool
        # (fold != TEST_FOLD) only, never the test fold, then apply the SAME
        # fitted threshold to label both train and test rows (correct
        # methodology: the threshold is a fitted statistic, the label
        # application is not).
        train_pool_for_fit = tl_usable & (fold_arr != TEST_FOLD)
        usable_protein_vals = protein_log10[train_pool_for_fit]
        rounded = np.round(usable_protein_vals, 2)
        vals, counts = np.unique(rounded, return_counts=True)
        mode_idx = np.argmax(counts)
        mode_val, mode_frac = vals[mode_idx], counts[mode_idx] / len(usable_protein_vals)
        is_floor_artifact = mode_frac > 0.05
        floor_value = mode_val if is_floor_artifact else None

        if is_floor_artifact:
            at_floor = tl_usable & (np.round(protein_log10, 2) == floor_value)
            tl_usable_for_regression = tl_usable & ~at_floor
            tl_active = tl_usable & (np.round(protein_log10, 2) > floor_value)
        else:
            tl_usable_for_regression = tl_usable.copy()
            tl_active = tl_usable & (protein_log10 > np.median(usable_protein_vals))

        fold = df["fold"].values.astype(int)
        oligo_ids = df["OLIGO ID"].values

        print(f"  tx_usable: {tx_usable.sum()}, tx_active: {tx_active.sum()}")
        print(f"  tl_usable: {tl_usable.sum()}, floor_value={floor_value} (frac={mode_frac:.3f}, "
              f"is_artifact={is_floor_artifact}), tl_usable_for_regression: {tl_usable_for_regression.sum()}, "
              f"tl_active: {tl_active.sum()}")
        print(f"  fold distribution: {np.bincount(fold)}")

        np.savez_compressed(
            CACHE / f"{host}_baseline_data.npz",
            onehot=onehot, oligo_ids=oligo_ids, fold=fold,
            tx_usable=tx_usable, tx_active=tx_active, tx_norm=tx_norm,
            tl_usable=tl_usable, tl_usable_for_regression=tl_usable_for_regression,
            tl_active=tl_active, protein_log10=protein_log10,
            tl_floor_value=floor_value if floor_value is not None else np.nan,
        )
        print(f"  Wrote {CACHE / f'{host}_baseline_data.npz'}")


if __name__ == "__main__":
    main()
