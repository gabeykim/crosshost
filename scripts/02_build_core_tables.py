"""
Build core Parquet tables from the Johns et al. 2018 Nature Methods supplementary
XLSX files (PMC6065261 / NIHMS945382 supplements 3-7).

Produces:
  data/three_host_library.parquet  -- one row per OLIGO ID (barcode-level), with
                                       metadata + BS/EC/PA transcription & translation
  data/rs241.parquet               -- one row per RS241 id, six-host transcription
                                       (log2) and translation (log10) values
Prints only summary stats (row counts, null counts) -- never dumps raw rows.
"""
import openpyxl
import pandas as pd
from pathlib import Path

RAW = Path(__file__).resolve().parent.parent / "raw"
DATA = Path(__file__).resolve().parent.parent / "data"
DATA.mkdir(exist_ok=True)


def sheet_to_df(path, sheet_name, header=True):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb[sheet_name]
    rows = ws.iter_rows(values_only=True)
    cols = next(rows)
    data = list(rows)
    wb.close()
    df = pd.DataFrame(data, columns=cols)
    return df


def build_three_host_library():
    p3 = RAW / "NIHMS945382-supplement-3.xlsx"
    p4 = RAW / "NIHMS945382-supplement-4.xlsx"

    print("Loading Supplement 3 (Table 1: Oligo Library Info) - Metadata sheet...")
    meta = sheet_to_df(p3, "Metadata")
    print(f"  Metadata shape: {meta.shape}")
    print(f"  Metadata columns: {list(meta.columns)}")

    print("Loading Supplement 3 - Full Constructs sheet...")
    constructs = sheet_to_df(p3, "Full Constructs")
    constructs.columns = ["OLIGO ID", "full_construct_seq"]
    print(f"  Full Constructs shape: {constructs.shape}")

    print("Loading Supplement 4 (Table 2: BS_EC_PA_Express) - BS/EC/PA sheets...")
    host_dfs = {}
    for host in ["BS", "EC", "PA"]:
        hdf = sheet_to_df(p4, host)
        hdf = hdf.rename(columns={c: f"{host}_{c}" for c in hdf.columns if c != "OLIGO ID"})
        host_dfs[host] = hdf
        print(f"  {host} shape: {hdf.shape}, columns: {list(hdf.columns)}")

    # Merge everything on OLIGO ID
    df = meta.merge(constructs, on="OLIGO ID", how="outer", validate="one_to_one")
    for host in ["BS", "EC", "PA"]:
        df = df.merge(host_dfs[host], on="OLIGO ID", how="outer", validate="one_to_one")

    print(f"\nMerged library shape: {df.shape}")

    numeric_cols = []
    for host in ["BS", "EC", "PA"]:
        numeric_cols += [
            f"{host}_rna_count", f"{host}_dna_count", f"{host}_primary TSS",
            f"{host}_fraction_around_primary_TSS", f"{host}_tx_raw", f"{host}_tx_norm",
            f"{host}_protein (log10)", f"{host}_delta_G", f"{host}_best_sigma70_match_score",
        ]
    numeric_cols += ["intgen_lp", "intgen_rp", "gene_id", "core count", "genome_id"]
    for c in numeric_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    out_path = DATA / "three_host_library.parquet"
    df.to_parquet(out_path, index=False)
    print(f"Wrote {out_path} ({out_path.stat().st_size/1e6:.2f} MB)")
    return df


def build_rs241():
    p6 = RAW / "NIHMS945382-supplement-6.xlsx"
    print("\nLoading Supplement 6 (Table 4: RS241) - Log2 Transcription / Log10 Translation...")
    tx = sheet_to_df(p6, "Log2 Transcription")
    tx.columns = ["id"] + [f"tx_log2_{c}" for c in tx.columns[1:]]
    tl = sheet_to_df(p6, "Log10 Translation")
    tl.columns = ["id"] + [f"tl_log10_{c}" for c in tl.columns[1:]]
    print(f"  Log2 Transcription shape: {tx.shape}")
    print(f"  Log10 Translation shape: {tl.shape}")

    df = tx.merge(tl, on="id", how="outer", validate="one_to_one")
    print(f"  Merged RS241 shape: {df.shape}")

    for c in df.columns:
        if c != "id":
            df[c] = pd.to_numeric(df[c], errors="coerce")

    out_path = DATA / "rs241.parquet"
    df.to_parquet(out_path, index=False)
    print(f"Wrote {out_path} ({out_path.stat().st_size/1e6:.2f} MB)")
    return df


if __name__ == "__main__":
    build_three_host_library()
    build_rs241()
