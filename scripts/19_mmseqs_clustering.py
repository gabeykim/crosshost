"""
GATE 2 - Task 3: MMseqs2 sequence-identity clustering, the ONLY leakage
defense per amendment C1 (barcode dedup does not apply). Clusters ALL 29,249
regulatory sequences in the library (not just the 11,276 three-host-usable
subset) so that fold assignments are reusable across any downstream task that
uses a different sequence subset (e.g. per-host baselines with different N,
or a translation-only subset) -- every RS's cluster ID is looked up once here
and can be joined against any task-specific filter later.
"""
import subprocess
import pandas as pd
from pathlib import Path
import shutil

DATA = Path(__file__).resolve().parent.parent / "data"
RAW = Path(__file__).resolve().parent.parent / "raw"
SPLITS = DATA / "splits"
SPLITS.mkdir(exist_ok=True)
MMSEQS_DIR = RAW / "mmseqs"
MMSEQS_DIR.mkdir(exist_ok=True)

THRESHOLDS = [0.3, 0.5, 0.7]


def export_fasta():
    df = pd.read_parquet(DATA / "three_host.parquet")
    fasta_path = MMSEQS_DIR / "all_rs.fasta"
    with open(fasta_path, "w") as f:
        for _, row in df.iterrows():
            f.write(f">{row['OLIGO ID']}\n{row['regulatory_sequence']}\n")
    print(f"Exported {len(df)} sequences to {fasta_path}")
    return fasta_path, df


def run_mmseqs_cluster(fasta_path, min_seq_id, coverage=0.8):
    tag = f"id{int(min_seq_id*100)}"
    out_prefix = MMSEQS_DIR / f"cluster_{tag}"
    tmp_dir = MMSEQS_DIR / f"tmp_{tag}"
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir()

    # This MMseqs2 build has no --strand flag for easy-cluster (nucleotide
    # search here implicitly considers both strands, doubling the effective
    # DB size 29,249 -> 55,814, which caused the first attempt to exceed this
    # machine's 8GB RAM: "Cannot fit databases into 7G"). --split-memory-limit
    # forces MMseqs2 to split the prefilter search into RAM-sized chunks
    # instead of loading everything at once.
    # `easy-cluster` (cascaded clustering: linclust redundancy pass + a full
    # sensitive prefilter/search second pass) repeatedly died with "Cannot
    # fit databases into 1G" even after reducing sensitivity, thread count,
    # and split-memory-limit -- this sandboxed environment has 8GB RAM total
    # but a `vm_stat`/`vm.swapusage` check showed swap already ~92% committed
    # (only ~1.5GB genuinely free), and the full search step's k-mer index
    # apparently cannot shrink below that floor regardless of --split.
    # Switched to `easy-linclust` -- MMseqs2's dedicated linear-time,
    # low-memory clustering mode (Steinegger & Soding 2018, Nat Commun,
    # "Clustering huge protein sequence sets in linear time"), which is
    # MMseqs2's own recommended approach for exactly this situation (large
    # sequence sets, constrained memory). Same --min-seq-id/-c semantics.
    cmd = ["mmseqs", "easy-linclust", str(fasta_path), str(out_prefix), str(tmp_dir),
           "--min-seq-id", str(min_seq_id), "-c", str(coverage), "--threads", "2",
           "--cov-mode", "0"]
    print(f"\nRunning: {' '.join(cmd)}")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        print(proc.stdout[-3000:])
        print(proc.stderr[-3000:])
        raise RuntimeError(f"mmseqs easy-cluster failed at min_seq_id={min_seq_id}")

    cluster_tsv = Path(f"{out_prefix}_cluster.tsv")
    clusters = pd.read_csv(cluster_tsv, sep="\t", header=None, names=["representative", "member"])
    clusters["representative"] = clusters["representative"].astype(str)
    clusters["member"] = clusters["member"].astype(str)
    return clusters


def main():
    fasta_path, df = export_fasta()

    summary_rows = []
    for thr in THRESHOLDS:
        clusters = run_mmseqs_cluster(fasta_path, thr)
        cluster_sizes = clusters.groupby("representative").size().sort_values(ascending=False)
        n_clusters = len(cluster_sizes)
        print(f"\n--- min-seq-id={thr} ---")
        print(f"  n_clusters: {n_clusters}")
        print(f"  largest cluster size: {cluster_sizes.iloc[0]}")
        print(f"  cluster size distribution: singletons={int((cluster_sizes==1).sum())}, "
              f"2-10={int(((cluster_sizes>=2)&(cluster_sizes<=10)).sum())}, "
              f"11-100={int(((cluster_sizes>=11)&(cluster_sizes<=100)).sum())}, "
              f">100={int((cluster_sizes>100).sum())}")

        out_path = SPLITS / f"clusters_id{int(thr*100)}.parquet"
        clusters.to_parquet(out_path, index=False)
        print(f"  wrote {out_path}")

        summary_rows.append({
            "min_seq_id": thr, "n_clusters": n_clusters,
            "largest_cluster_size": int(cluster_sizes.iloc[0]),
            "n_singletons": int((cluster_sizes == 1).sum()),
            "mean_cluster_size": float(cluster_sizes.mean()),
            "median_cluster_size": float(cluster_sizes.median()),
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(SPLITS / "clustering_summary.csv", index=False)
    print(f"\n{summary_df.to_string(index=False)}")
    print(f"\nWrote {SPLITS / 'clustering_summary.csv'}")


if __name__ == "__main__":
    main()
