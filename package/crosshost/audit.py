"""Leakage checks a user can run against their OWN usage of the splits --
e.g. after writing a custom train/test loop, confirm no RS241 leaked into
training and no fold was crossed. This is a lightweight, importable subset
of the project's own scripts/audit_leakage.py (which additionally checks
sequence-identity thresholds and table-hash provenance against the
maintainers' own frozen data -- run that full script, included in this
package's `scripts/` directory, for the complete audit)."""
import pandas as pd
from .splits import load_splits
from .data import load_rs241


def check_no_rs241_leakage(train_oligo_ids):
    """Raises AssertionError if any RS241 id appears in your training set.
    train_oligo_ids: iterable of OLIGO ID values (str or int, compared as str)."""
    rs241 = load_rs241()
    rs241_ids = set(rs241["id"].astype(str))
    train_ids = set(str(x) for x in train_oligo_ids)
    overlap = rs241_ids & train_ids
    assert not overlap, (
        f"{len(overlap)} RS241 id(s) found in your training set -- RS241 must never "
        f"be trained on. First few: {list(overlap)[:5]}"
    )
    return True


def check_no_fold_crossing(oligo_ids, expected_fold):
    """Raises AssertionError if any of oligo_ids belongs to a fold other
    than expected_fold, per the frozen split assignment."""
    splits = load_splits()
    splits = splits.set_index("OLIGO ID")
    ids = [str(x) for x in oligo_ids]
    splits.index = splits.index.astype(str)
    bad = [i for i in ids if i in splits.index and splits.loc[i, "fold"] != expected_fold]
    assert not bad, f"{len(bad)} id(s) belong to a different fold than expected_fold={expected_fold}"
    return True
