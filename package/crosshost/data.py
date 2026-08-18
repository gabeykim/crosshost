"""Dataset loaders for the CROSSHOST benchmark.

All tables are derived, with attribution, from Johns et al. 2018
(Nat. Methods 15:323-329, BioProject PRJNA431139). See data/core/CITATION.md.
Raw Springer Nature-copyrighted supplementary tables are NOT redistributed
in this package -- only the derived activity/feature values used by the
benchmark, per the licensing quarantine documented in README.md.
"""
from pathlib import Path
import pandas as pd

_DATA_ROOT = Path(__file__).resolve().parent.parent / "data" / "core"


def _data_path(name):
    p = _DATA_ROOT / name
    if not p.exists():
        raise FileNotFoundError(
            f"{p} not found. If you cloned this repo without Git LFS or the data "
            f"release, run `make fetch-data` first (see README.md)."
        )
    return p


def load_three_host_library():
    """The primary 3-host (E. coli / B. subtilis / P. aeruginosa) regulatory
    sequence library, 29,249 rows before RS241 exclusion. Columns include
    `regulatory_sequence` (165bp), and per-host `{H}_tx_usable`/`{H}_tx_active`/
    `{H}_tx_norm` (transcription) and `{H}_tl_usable`/`{H}_protein_log10`
    (translation) for H in {EC, BS, PA}. See README.md for the exact
    usable/active definitions (dna_count!=0 & (rna+dna)>=15; active =
    usable & rna_count>0; translation floor-value correction applied)."""
    return pd.read_parquet(_data_path("three_host.parquet"))


def load_rs241():
    """The RS241 held-out-host evaluation set (241 sequences spanning all 6
    hosts; 207 have recoverable sequence text). NEVER used for training --
    enforced both by this package's split loader (RS241 ids are absent from
    the fold assignment table) and by `crosshost.audit.check_no_leakage()`.
    """
    return pd.read_parquet(_data_path("rs241.parquet"))


def load_genomic_features():
    """37-dimensional genome-encoded host feature vector, one row per host
    (6 hosts), NOT z-scored (z-score across the full 6-host table yourself
    if replicating this project's model training -- see README.md)."""
    return pd.read_parquet(_data_path("hosts_genomic.parquet")).set_index("host")


def load_physiology_features():
    """6-dimensional depth-matched proteome-derived physiology proxy vector,
    one row per host. See data/core/CITATION.md for per-host provenance
    (PaxDb for 5 hosts, PRIDE PXD027874 for V. natriegens -- heterogeneous
    source types, disclosed per standing rule SR2)."""
    return pd.read_parquet(_data_path("hosts_physiology.parquet")).set_index("host")
