"""CROSSHOST: a benchmark for predicting bacterial regulatory DNA activity
across host species. See README.md for the full description and
out/PAPER_FRAMING.md (in the source repository) for the scientific results
this package's baseline scores and splits were used to produce.
"""
from .data import (
    load_three_host_library,
    load_rs241,
    load_genomic_features,
    load_physiology_features,
)
from .splits import load_splits, N_FOLDS, PRIMARY_HOSTS, RS241_HOSTS
from .evaluate import evaluate, bootstrap_ci_90

__version__ = "1.0.0"

__all__ = [
    "load_three_host_library",
    "load_rs241",
    "load_genomic_features",
    "load_physiology_features",
    "load_splits",
    "evaluate",
    "bootstrap_ci_90",
    "N_FOLDS",
    "PRIMARY_HOSTS",
    "RS241_HOSTS",
]
