"""Package smoke tests -- loaders, evaluate() API, held-out scoring,
leakage checks. Run with `pytest tests/ -v` from the package/ directory."""
import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import crosshost
from crosshost.audit import check_no_rs241_leakage


def test_load_three_host_library():
    df = crosshost.load_three_host_library()
    assert len(df) > 29000
    assert "regulatory_sequence" in df.columns
    assert (df["regulatory_sequence"].str.len() == 165).all()


def test_load_rs241():
    df = crosshost.load_rs241()
    assert len(df) == 241


def test_load_genomic_features():
    df = crosshost.load_genomic_features()
    assert df.shape[0] == 6
    assert set(["EC", "BS", "PA", "SE", "VN", "CG"]).issubset(set(df.index))


def test_load_physiology_features():
    df = crosshost.load_physiology_features()
    assert df.shape[0] == 6


def test_load_splits():
    df = crosshost.load_splits()
    assert set(df["fold"].unique()) == {0, 1, 2, 3, 4}


def test_evaluate_api_runs_and_random_predictions_score_near_zero():
    np.random.seed(0)
    n = 500
    preds = [{"transcription": {"active_prob": np.random.rand(n), "strength_pred": np.random.randn(n)}} for _ in range(5)]
    targs = [{"transcription": {"active_true": np.random.randint(0, 2, n),
                                  "strength_true": np.where(np.random.rand(n) > 0.3, np.random.randn(n), np.nan)}} for _ in range(5)]
    result = crosshost.evaluate(preds, targs, readouts=("transcription",))
    rho = result["transcription"]["spearman_rho"]["mean"]
    assert rho is not None
    assert abs(rho) < 0.15, f"random predictions should score near-zero rho, got {rho}"
    auc = result["transcription"]["auc"]["mean"]
    assert 0.35 < auc < 0.65, f"random predictions should score near-0.5 AUC, got {auc}"


def test_rs241_leakage_check_passes_on_real_split():
    splits = crosshost.load_splits()
    train_ids = splits[splits.fold != 0]["OLIGO ID"].tolist()
    assert check_no_rs241_leakage(train_ids) is True


def test_rs241_leakage_check_catches_injected_leak():
    rs241 = crosshost.load_rs241()
    fake_train = [rs241["id"].iloc[0], rs241["id"].iloc[1], "99999999"]
    with pytest.raises(AssertionError):
        check_no_rs241_leakage(fake_train)


def test_held_out_eval_inputs_have_no_labels():
    import pandas as pd
    inputs = pd.read_parquet(Path(__file__).resolve().parent.parent / "held_out_eval" / "test_inputs.parquet")
    label_prefixes = ("BS_", "EC_", "PA_")
    leaked = [c for c in inputs.columns if c.startswith(label_prefixes) and c not in ("EC", "BS", "PA")]
    assert not leaked, f"held-out test_inputs.parquet must not contain label columns, found: {leaked}"
