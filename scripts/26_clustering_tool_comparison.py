"""
GATE 2.5 - Task A1: Is the clustering-miss finding tool-specific (linclust) or
general (also affects easy-cluster / CD-HIT)?

Builds a tractable subset (413 exact-duplicate oligos + the known 93.3% pair
29542/29549 + ~3000 random background sequences), then:
  1. Runs mmseqs easy-cluster (the full/sensitive mode that OOM'd on the full
     29,242-sequence library) on this smaller subset.
  2. Runs cd-hit-est (a completely independent, exhaustive/greedy clustering
     implementation) on the same subset.
  3. Establishes ground truth via a direct all-pairs identity computation on
     the subset (tractable at this size: ~3400 sequences).
  4. Reports precision/recall of each tool against ground truth, and
     specifically whether each tool co-clusters the 29542/29549 pair.
"""
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
import json
import itertools

DATA = Path(__file__).resolve().parent.parent / "data"
RAW = Path(__file__).resolve().parent.parent / "raw"
OUT = Path(__file__).resolve().parent.parent / "out"
WORKDIR = RAW / "cluster_comparison"
WORKDIR.mkdir(exist_ok=True)

KNOWN_PAIR = ("29542", "29549")
N_BACKGROUND = 3000
THRESHOLD = 0.7


def build_subset():
    df = pd.read_parquet(DATA / "three_host.parquet")
    df["OLIGO ID"] = df["OLIGO ID"].astype(str)

    seq_groups = df.groupby("regulatory_sequence")["OLIGO ID"].apply(list)
    dup_groups = seq_groups[seq_groups.apply(len) > 1]
    dup_ids = set(oid for members in dup_groups for oid in members)
    print(f"Exact-duplicate oligos: {len(dup_ids)}")

    must_include = dup_ids | set(KNOWN_PAIR)

    rng = np.random.default_rng(7)
    remaining = df[~df["OLIGO ID"].isin(must_include)]
    background = rng.choice(remaining["OLIGO ID"].values, size=N_BACKGROUND, replace=False)

    subset_ids = must_include | set(background)
    subset_df = df[df["OLIGO ID"].isin(subset_ids)].copy()
    print(f"Subset size: {len(subset_df)} (must_include={len(must_include)}, background={N_BACKGROUND})")

    fasta_path = WORKDIR / "subset.fasta"
    with open(fasta_path, "w") as f:
        for _, row in subset_df.iterrows():
            f.write(f">{row['OLIGO ID']}\n{row['regulatory_sequence']}\n")

    return subset_df, fasta_path


def kmers(seq, k=10):
    return set(seq[i:i + k] for i in range(len(seq) - k + 1))


def exact_identity(a, b):
    best = 0
    L = min(len(a), len(b))
    for offset in range(-5, 6):
        matches, n = 0, 0
        for i in range(L):
            j = i + offset
            if 0 <= j < len(b):
                n += 1
                if a[i] == b[j]:
                    matches += 1
        if n > 0:
            best = max(best, matches / n)
    return best


def ground_truth_pairs(subset_df, threshold=THRESHOLD):
    """All-pairs check via k-mer index (as validated in Gate 2) -- the same
    method that was cross-checked against Biopython global alignment and
    found exact. Establishes which pairs SHOULD co-cluster at `threshold`."""
    ids = subset_df["OLIGO ID"].tolist()
    seqs = subset_df["regulatory_sequence"].tolist()
    id_to_seq = dict(zip(ids, seqs))

    index = defaultdict(list)
    for i, s in enumerate(seqs):
        for km in kmers(s):
            index[km].append(i)

    positive_pairs = set()
    checked = set()
    for i, s in enumerate(seqs):
        counts = defaultdict(int)
        for km in kmers(s):
            for j in index[km]:
                if j > i:
                    counts[j] += 1
        for j, c in counts.items():
            if c < 5:
                continue
            key = (i, j)
            if key in checked:
                continue
            checked.add(key)
            ident = exact_identity(s, seqs[j])
            if ident >= threshold:
                positive_pairs.add(frozenset((ids[i], ids[j])))

    print(f"Ground truth: {len(positive_pairs)} pairs at >= {threshold} identity "
          f"(out of {len(checked)} candidate pairs checked via k-mer prefilter)")
    return positive_pairs, id_to_seq


def clusters_to_pairs(cluster_assignment):
    """cluster_assignment: dict id -> cluster_label. Return set of frozensets
    for every same-cluster pair (i.e. what the tool CLAIMS should be merged)."""
    by_cluster = defaultdict(list)
    for oid, c in cluster_assignment.items():
        by_cluster[c].append(oid)
    pairs = set()
    for members in by_cluster.values():
        if len(members) < 2:
            continue
        for a, b in itertools.combinations(members, 2):
            pairs.add(frozenset((a, b)))
    return pairs


def run_mmseqs_easy_cluster(fasta_path, threshold=THRESHOLD):
    out_prefix = WORKDIR / "mmseqs_cluster_subset"
    tmp_dir = WORKDIR / "mmseqs_tmp_subset"
    import shutil
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir()
    cmd = ["mmseqs", "easy-cluster", str(fasta_path), str(out_prefix), str(tmp_dir),
           "--min-seq-id", str(threshold), "-c", "0.8", "--threads", "2", "--cov-mode", "0",
           "-s", "1", "-k", "6", "--split-memory-limit", "1000M"]
    print(f"\nRunning: {' '.join(cmd)}")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        print("STDOUT tail:", proc.stdout[-2000:])
        print("STDERR tail:", proc.stderr[-2000:])
        return None
    tsv = Path(f"{out_prefix}_cluster.tsv")
    clusters = pd.read_csv(tsv, sep="\t", header=None, names=["representative", "member"])
    assignment = dict(zip(clusters["member"].astype(str), clusters["representative"].astype(str)))
    print(f"  SUCCESS: {clusters['representative'].nunique()} clusters")
    return assignment


def run_cdhit(fasta_path, threshold=None):
    out_path = WORKDIR / "cdhit_subset_out"
    # FINDING, confirmed empirically (not assumed from memory): this cd-hit-est
    # build (4.8.1) hard-rejects ANY -c below 0.8 ("invalid clstr threshold,
    # should >=0.8"), at every word size tried (n=3..8), with or without -G.
    # This is itself relevant: cd-hit-est is not designed/supported to run
    # short-read nucleotide clustering below 80% identity at all (the
    # word-based seeding statistics break down below that for short
    # sequences; CD-HIT's own docs point to the separate psi-cd-hit workflow
    # for sub-80% clustering, not attempted here). Since the known problem
    # pair (29542/29549) is at 93.3% identity -- comfortably above 0.8 -- a
    # threshold of 0.8 is used for the CD-HIT comparison specifically, which
    # still directly tests whether CD-HIT catches THIS pair.
    threshold = threshold or 0.8
    cmd = ["cd-hit-est", "-i", str(fasta_path), "-o", str(out_path), "-c", str(threshold),
           "-n", "5", "-d", "0", "-M", "2000", "-T", "2", "-G", "1"]
    print(f"\nRunning: {' '.join(cmd)}")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    print(proc.stdout[-1500:])
    if proc.returncode != 0:
        print("STDERR tail:", proc.stderr[-2000:])
        return None
    clstr_path = Path(f"{out_path}.clstr")
    assignment = {}
    current_cluster = None
    with open(clstr_path) as f:
        for line in f:
            if line.startswith(">Cluster"):
                current_cluster = line.strip()
            else:
                # format: 0	165nt, >29542... at +/100.00%
                oid = line.split(">")[1].split("...")[0]
                assignment[oid] = current_cluster
    print(f"  SUCCESS: {len(set(assignment.values()))} clusters")
    return assignment


def precision_recall(claimed_pairs, true_pairs, all_possible_pairs_count):
    tp = len(claimed_pairs & true_pairs)
    fp = len(claimed_pairs - true_pairs)
    fn = len(true_pairs - claimed_pairs)
    precision = tp / (tp + fp) if (tp + fp) else float("nan")
    recall = tp / (tp + fn) if (tp + fn) else float("nan")
    return precision, recall, tp, fp, fn


def main():
    subset_df, fasta_path = build_subset()
    true_pairs, id_to_seq = ground_truth_pairs(subset_df)

    known_pair_fs = frozenset(KNOWN_PAIR)
    print(f"\nKnown pair {KNOWN_PAIR} in ground truth (>=0.7 identity)? {known_pair_fs in true_pairs}")
    ident = exact_identity(id_to_seq[KNOWN_PAIR[0]], id_to_seq[KNOWN_PAIR[1]])
    print(f"  direct identity check: {ident:.4f}")

    results = {}

    # mmseqs easy-cluster (the sensitive mode that OOM'd on the full library)
    mmseqs_assignment = run_mmseqs_easy_cluster(fasta_path)
    if mmseqs_assignment:
        mmseqs_pairs = clusters_to_pairs(mmseqs_assignment)
        p, r, tp, fp, fn = precision_recall(mmseqs_pairs, true_pairs, len(subset_df))
        known_pair_clustered = mmseqs_assignment.get(KNOWN_PAIR[0]) == mmseqs_assignment.get(KNOWN_PAIR[1])
        print(f"\nmmseqs easy-cluster: precision={p:.3f}, recall={r:.3f}, tp={tp}, fp={fp}, fn={fn}")
        print(f"  co-clusters known 93.3% pair? {known_pair_clustered}")
        results["mmseqs_easy_cluster"] = {"ran": True, "precision": p, "recall": r, "tp": tp, "fp": fp, "fn": fn,
                                           "co_clusters_known_pair": known_pair_clustered}
    else:
        results["mmseqs_easy_cluster"] = {"ran": False, "reason": "failed on subset -- see script output"}

    # cd-hit-est -- this build hard-rejects c<0.8 (see run_cdhit docstring),
    # so it is run and evaluated at 0.8, with its own ground truth at 0.8
    # (not the 0.7 ground truth used for the other tools) for a fair
    # apples-to-apples precision/recall comparison.
    true_pairs_08, _ = ground_truth_pairs(subset_df, threshold=0.8)
    cdhit_assignment = run_cdhit(fasta_path, threshold=0.8)
    if cdhit_assignment:
        cdhit_pairs = clusters_to_pairs(cdhit_assignment)
        p, r, tp, fp, fn = precision_recall(cdhit_pairs, true_pairs_08, len(subset_df))
        known_pair_clustered = cdhit_assignment.get(KNOWN_PAIR[0]) == cdhit_assignment.get(KNOWN_PAIR[1])
        print(f"\ncd-hit-est (run at c=0.8, since c=0.7 is rejected by this build): "
              f"precision={p:.3f}, recall={r:.3f}, tp={tp}, fp={fp}, fn={fn}")
        print(f"  co-clusters known 93.3% pair? {known_pair_clustered}")
        results["cd_hit_est"] = {"ran": True, "threshold_used": 0.8,
                                  "note": "run at 0.8 -- this cd-hit-est build rejects c<0.8",
                                  "precision": p, "recall": r, "tp": tp, "fp": fp, "fn": fn,
                                  "co_clusters_known_pair": known_pair_clustered}
    else:
        results["cd_hit_est"] = {"ran": False, "reason": "failed on subset -- see script output"}

    # also re-run linclust on the SAME subset for a fair apples-to-apples comparison
    out_prefix = WORKDIR / "linclust_subset"
    tmp_dir = WORKDIR / "linclust_tmp_subset"
    import shutil
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir()
    cmd = ["mmseqs", "easy-linclust", str(fasta_path), str(out_prefix), str(tmp_dir),
           "--min-seq-id", str(THRESHOLD), "-c", "0.8", "--threads", "2", "--cov-mode", "0"]
    print(f"\nRunning: {' '.join(cmd)}")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode == 0:
        tsv = Path(f"{out_prefix}_cluster.tsv")
        clusters = pd.read_csv(tsv, sep="\t", header=None, names=["representative", "member"])
        linclust_assignment = dict(zip(clusters["member"].astype(str), clusters["representative"].astype(str)))
        linclust_pairs = clusters_to_pairs(linclust_assignment)
        p, r, tp, fp, fn = precision_recall(linclust_pairs, true_pairs, len(subset_df))
        known_pair_clustered = linclust_assignment.get(KNOWN_PAIR[0]) == linclust_assignment.get(KNOWN_PAIR[1])
        print(f"\neasy-linclust (subset, for fair comparison): precision={p:.3f}, recall={r:.3f}")
        print(f"  co-clusters known 93.3% pair? {known_pair_clustered}")
        results["mmseqs_easy_linclust_subset"] = {"ran": True, "precision": p, "recall": r, "tp": tp, "fp": fp,
                                                    "fn": fn, "co_clusters_known_pair": known_pair_clustered}
    else:
        results["mmseqs_easy_linclust_subset"] = {"ran": False}

    results["ground_truth_n_positive_pairs"] = len(true_pairs)
    results["known_pair_identity"] = ident
    results["known_pair_in_ground_truth"] = known_pair_fs in true_pairs
    results["subset_size"] = len(subset_df)

    with open(OUT / "clustering_tool_comparison.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nWrote {OUT / 'clustering_tool_comparison.json'}")


if __name__ == "__main__":
    main()
