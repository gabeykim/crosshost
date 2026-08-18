"""Frozen fold splits. Genome-blocked + exact-duplicate safety net +
near-duplicate safety net (>=0.85 identity) -- see README.md for the full
methodology. Max train-test sequence identity across all 5 folds: 0.8485.

RS241 IDs are structurally absent from this table (never trained on --
verified by crosshost.audit.check_no_leakage(), which is the same check
this project's own audit_leakage.py runs)."""
from pathlib import Path
import pandas as pd

_DATA_ROOT = Path(__file__).resolve().parent.parent / "data" / "core"

N_FOLDS = 5
PRIMARY_HOSTS = ["EC", "BS", "PA"]
RS241_HOSTS = ["SE", "VN", "CG"]


def load_splits():
    """Returns a DataFrame with columns ['OLIGO ID', 'fold'], fold in 0..4.
    29,042 rows (three_host_library's 29,249 minus 207 RS241-overlapping ids
    removed from the main fold assignment)."""
    p = _DATA_ROOT / "fold_assignment_FINAL.parquet"
    if not p.exists():
        raise FileNotFoundError(f"{p} not found -- see README.md `make fetch-data`.")
    return pd.read_parquet(p)


def loho_train_test_split(library_df, splits_df, held_out_host, test_fold):
    """Leave-one-host-out helper matching this project's own evaluation
    convention exactly: train on the OTHER primary hosts' rows with
    fold != test_fold; test on held_out_host's rows with fold == test_fold.
    Returns (train_df, test_df), both subsets of library_df merged with fold."""
    merged = library_df.merge(splits_df, on="OLIGO ID", how="inner")
    train_hosts = [h for h in PRIMARY_HOSTS if h != held_out_host]
    train_df = merged[merged["fold"] != test_fold]
    test_df = merged[merged["fold"] == test_fold]
    return train_df, test_df, train_hosts
