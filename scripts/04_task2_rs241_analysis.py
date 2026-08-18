"""
GATE 1 - Task 2: RS241 six-host subset integrity check.

"Usable value" for RS241 is defined as: the cell in the Log2 Transcription (or
Log10 Translation) sheet is not null/blank for that host. These are already
paper-processed log-scale values (not raw counts), so no additional QC
threshold is reconstructable from the public table -- we report presence/
absence of a value exactly as released.
"""
import pandas as pd
from pathlib import Path
import json

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"

HOSTS = ["EC", "BS", "PA", "SE", "VN", "CG"]  # E. coli, B. subtilis, P. aeruginosa, S. enterica, V. natriegens, C. glutamicum

def main():
    df = pd.read_parquet(DATA / "rs241.parquet")
    print(f"RS241 table shape: {df.shape}")
    print(f"n unique id: {df['id'].nunique()}")
    print(f"Columns: {list(df.columns)}")

    results = {"n_sequences": int(df["id"].nunique())}

    print("\n" + "=" * 70)
    print("PER-HOST USABLE COUNTS -- TRANSCRIPTION (Log2)")
    print("=" * 70)
    tx_cols = {h: f"tx_log2_{h}" for h in HOSTS}
    per_host_tx = {}
    for h, c in tx_cols.items():
        n_usable = df[c].notna().sum()
        per_host_tx[h] = int(n_usable)
        print(f"  {h}: {n_usable} / {len(df)} usable ({100*n_usable/len(df):.1f}%)")

    print("\n" + "=" * 70)
    print("PER-HOST USABLE COUNTS -- TRANSLATION (Log10)")
    print("=" * 70)
    tl_cols = {h: f"tl_log10_{h}" for h in HOSTS}
    per_host_tl = {}
    for h, c in tl_cols.items():
        n_usable = df[c].notna().sum()
        per_host_tl[h] = int(n_usable)
        print(f"  {h}: {n_usable} / {len(df)} usable ({100*n_usable/len(df):.1f}%)")

    print("\n" + "=" * 70)
    print("SIX-HOST OVERLAP (usable in ALL 6 simultaneously)")
    print("=" * 70)
    tx_mask_all6 = pd.concat([df[c].notna() for c in tx_cols.values()], axis=1).all(axis=1)
    tl_mask_all6 = pd.concat([df[c].notna() for c in tl_cols.values()], axis=1).all(axis=1)
    print(f"  Transcription usable in all 6 hosts: {tx_mask_all6.sum()} / {len(df)}")
    print(f"  Translation usable in all 6 hosts:   {tl_mask_all6.sum()} / {len(df)}")

    print("\n" + "=" * 70)
    print("THREE-HOST-ONLY OVERLAP WITHIN RS241 (EC, BS, PA) -- for cross-check vs main library")
    print("=" * 70)
    tx_mask_3 = df[tx_cols["EC"]].notna() & df[tx_cols["BS"]].notna() & df[tx_cols["PA"]].notna()
    tl_mask_3 = df[tl_cols["EC"]].notna() & df[tl_cols["BS"]].notna() & df[tl_cols["PA"]].notna()
    print(f"  Transcription usable in EC+BS+PA: {tx_mask_3.sum()} / {len(df)}")
    print(f"  Translation usable in EC+BS+PA:   {tl_mask_3.sum()} / {len(df)}")

    print("\n" + "=" * 70)
    print("RAGGEDNESS: distribution of 'number of hosts with usable value' per RS")
    print("=" * 70)
    tx_n_hosts = pd.concat([df[c].notna() for c in tx_cols.values()], axis=1).sum(axis=1)
    tl_n_hosts = pd.concat([df[c].notna() for c in tl_cols.values()], axis=1).sum(axis=1)
    print("Transcription: n_hosts_with_value -> count of RSs")
    print(tx_n_hosts.value_counts().sort_index().to_string())
    print("\nTranslation: n_hosts_with_value -> count of RSs")
    print(tl_n_hosts.value_counts().sort_index().to_string())

    results["transcription"] = {
        "per_host_usable": per_host_tx,
        "all_6_hosts_usable": int(tx_mask_all6.sum()),
        "ec_bs_pa_usable": int(tx_mask_3.sum()),
        "n_hosts_distribution": {str(k): int(v) for k, v in tx_n_hosts.value_counts().sort_index().items()},
    }
    results["translation"] = {
        "per_host_usable": per_host_tl,
        "all_6_hosts_usable": int(tl_mask_all6.sum()),
        "ec_bs_pa_usable": int(tl_mask_3.sum()),
        "n_hosts_distribution": {str(k): int(v) for k, v in tl_n_hosts.value_counts().sort_index().items()},
    }

    # Also check whether RS241 ids are a subset of the main 3-host library OLIGO IDs
    main_df = pd.read_parquet(DATA / "three_host_library.parquet")
    main_ids = set(main_df["OLIGO ID"].tolist())
    rs241_ids = set(df["id"].tolist())
    overlap = rs241_ids & main_ids
    print("\n" + "=" * 70)
    print("RS241 id overlap with main 3-host library OLIGO IDs")
    print("=" * 70)
    print(f"  RS241 ids: {len(rs241_ids)}")
    print(f"  RS241 ids also present as OLIGO ID in main library: {len(overlap)}")
    print(f"  RS241 ids NOT in main library: {len(rs241_ids - main_ids)}")
    results["rs241_ids_in_main_library"] = len(overlap)
    results["rs241_ids_not_in_main_library"] = len(rs241_ids - main_ids)

    with open(OUT / "task2_rs241_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nWrote {OUT / 'task2_rs241_results.json'}")


if __name__ == "__main__":
    main()
