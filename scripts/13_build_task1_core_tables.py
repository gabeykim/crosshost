"""
GATE 2 - Task 1: Build the core frozen data tables for the three-host library
and RS241, from the Gate 1 / Gate 1.5 outputs.

Per amendment C3: translation N = 7,887 (our own computed value), NOT the
paper's stated 8,898. See comment at TRANSLATION_QC below.

Per amendment C1/Gate 1 finding: NO barcode deduplication is applied -- the
released library is already 1:1 OLIGO ID : intgen_id (genomic locus). This is
ASSERTED in code (not just assumed) so a future change to the input data would
fail loudly rather than silently reintroducing a leakage risk.
"""
import pandas as pd
import numpy as np
from pathlib import Path
import hashlib
import json

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"

MANIFEST_PATH = DATA / "MANIFEST.json"


def sha256_of_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def build_three_host():
    df = pd.read_parquet(DATA / "three_host_library.parquet")

    # --- ASSERTION: 1:1 OLIGO ID : intgen_id (genomic locus). This is the Gate 1
    # finding that overturned the charter's barcode-dedup assumption (amendment C1).
    # If this ever fails, STOP -- it means the input data changed and the barcode-
    # replicate structure the charter originally worried about may now be real.
    n_oligo = df["OLIGO ID"].nunique()
    n_locus = df["intgen_id"].nunique()
    assert n_oligo == len(df), f"Duplicate OLIGO IDs found: {n_oligo} unique vs {len(df)} rows"
    assert n_locus == len(df), (
        f"OLIGO ID is no longer 1:1 with intgen_id ({n_locus} unique loci vs {len(df)} rows) -- "
        f"barcode-replicate structure may now be present. Re-run Gate 1's barcode analysis "
        f"(scripts/03_task1_task2_analysis.py Q5) before proceeding; do NOT silently dedupe."
    )
    print(f"ASSERTION PASSED: {n_oligo} rows, all 1:1 OLIGO ID:intgen_id (no barcode dedup needed)")

    # Per-host QC (paper's own rule, verified Gate 1): usable transcription =
    # dna_count!=0 AND (rna+dna)>=15. Usable translation = protein(log10) not null.
    # "Active" (a stricter flag than "usable") = rna_count>0 among usable rows.
    rows = {"OLIGO ID": df["OLIGO ID"], "intgen_id": df["intgen_id"],
            "regulatory_sequence": df["Regulatory Sequence"],
            "genome_id": df["genome_id"], "phylum": df["phylum"], "class": df["class"],
            "family": df["family"], "genus": df["genus"], "order": df["order"],
            "strain_name": df["[strain_name]"], "gc_content": df["Regulatory Sequence"].apply(
                lambda s: (s.count("G") + s.count("C")) / len(s) if isinstance(s, str) and len(s) else np.nan)}

    for host in ["BS", "EC", "PA"]:
        rna = df[f"{host}_rna_count"]
        dna = df[f"{host}_dna_count"]
        tx_usable = (dna != 0) & ((rna + dna) >= 15)
        tx_active = tx_usable & (rna > 0)
        protein = df[f"{host}_protein (log10)"]
        tl_usable = protein.notna()

        rows[f"{host}_tx_raw"] = df[f"{host}_tx_raw"]
        rows[f"{host}_tx_norm"] = df[f"{host}_tx_norm"]
        rows[f"{host}_tx_usable"] = tx_usable
        rows[f"{host}_tx_active"] = tx_active
        rows[f"{host}_protein_log10"] = protein
        rows[f"{host}_tl_usable"] = tl_usable
        rows[f"{host}_rna_count"] = rna
        rows[f"{host}_dna_count"] = dna

    out = pd.DataFrame(rows)

    tx_all3 = out["BS_tx_usable"] & out["EC_tx_usable"] & out["PA_tx_usable"]
    tl_all3 = out["BS_tl_usable"] & out["EC_tl_usable"] & out["PA_tl_usable"]
    n_tx = int(tx_all3.sum())
    n_tl = int(tl_all3.sum())
    print(f"Three-host transcription-usable: {n_tx} (Gate 1 established 11,276)")
    print(f"Three-host translation-usable: {n_tl} (amendment C3: use 7,887, NOT paper's 8,898)")
    assert n_tx == 11276, f"Transcription overlap drifted: {n_tx} != 11276"
    assert n_tl == 7887, f"Translation overlap drifted: {n_tl} != 7887 (amendment C3 value)"

    out["usable_all3_transcription"] = tx_all3
    out["usable_all3_translation"] = tl_all3

    out_path = DATA / "three_host.parquet"
    out.to_parquet(out_path, index=False)
    print(f"Wrote {out_path} ({len(out)} rows, {len(out.columns)} columns)")
    return out_path


def build_rs241():
    rs241 = pd.read_parquet(DATA / "rs241.parquet")
    pairwise = json.load(open(OUT / "task_a_rs241_pairwise_results.json"))

    HOSTS = ["EC", "BS", "PA", "SE", "VN", "CG"]
    out = rs241.copy()
    for h in HOSTS:
        out[f"tx_usable_{h}"] = out[f"tx_log2_{h}"].notna()
        out[f"tl_usable_{h}"] = out[f"tl_log10_{h}"].notna()

    # attach pairwise-usability-derived flags: is this RS usable for each of the
    # Gate 1.5-identified runnable held-out configurations?
    tx_ec_pa = out["tx_usable_EC"] & out["tx_usable_PA"]
    tl_ec_pa = out["tl_usable_EC"] & out["tl_usable_PA"]
    tx_ec_bs_pa = out["tx_usable_EC"] & out["tx_usable_BS"] & out["tx_usable_PA"]
    tl_ec_bs_pa = out["tl_usable_EC"] & out["tl_usable_BS"] & out["tl_usable_PA"]

    for target in ["SE", "VN", "CG"]:
        out[f"usable_heldout_{target}_train_EC_PA_tx"] = out[f"tx_usable_{target}"] & tx_ec_pa
        out[f"usable_heldout_{target}_train_EC_PA_tl"] = out[f"tl_usable_{target}"] & tl_ec_pa
        out[f"usable_heldout_{target}_train_EC_BS_PA_tx"] = out[f"tx_usable_{target}"] & tx_ec_bs_pa
        out[f"usable_heldout_{target}_train_EC_BS_PA_tl"] = out[f"tl_usable_{target}"] & tl_ec_bs_pa

    for target in ["SE", "VN", "CG"]:
        n1 = int(out[f"usable_heldout_{target}_train_EC_PA_tx"].sum())
        n2 = int(out[f"usable_heldout_{target}_train_EC_BS_PA_tx"].sum())
        print(f"  held_out={target}: train=EC+PA N={n1} (tx), train=EC+BS+PA N={n2} (tx)")

    out_path = DATA / "rs241.parquet"  # overwrite with pairwise flags added
    out.to_parquet(out_path, index=False)
    print(f"Wrote {out_path} ({len(out)} rows, {len(out.columns)} columns)")
    return out_path


def build_hosts_physiology():
    phys = pd.read_csv(DATA / "hosts_physiology_final.csv")
    out_path = DATA / "hosts_physiology.parquet"
    phys.to_parquet(out_path, index=False)
    print(f"Wrote {out_path} ({len(phys)} rows, {len(phys.columns)} columns)")
    print(f"  source/proteome_depth columns present: "
          f"{'physiology_source' in phys.columns and 'proteome_depth_used' in phys.columns}")
    return out_path


def update_manifest(paths):
    manifest = {}
    if MANIFEST_PATH.exists():
        manifest = json.load(open(MANIFEST_PATH))
    for p in paths:
        manifest[p.name] = {
            "sha256": sha256_of_file(p),
            "size_bytes": p.stat().st_size,
            "n_rows": len(pd.read_parquet(p)) if p.suffix == ".parquet" else None,
        }
    with open(MANIFEST_PATH, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"\nUpdated {MANIFEST_PATH}")


if __name__ == "__main__":
    print("=== Building data/three_host.parquet ===")
    p1 = build_three_host()
    print("\n=== Building data/rs241.parquet (with pairwise flags) ===")
    p2 = build_rs241()
    print("\n=== Building data/hosts_physiology.parquet ===")
    p3 = build_hosts_physiology()
    update_manifest([p1, p2, p3])
