"""
GATE 2 - Task 3.4: FINAL frozen split construction.

Rationale (see Gate 2 memo for full narrative): three distinct leakage
mechanisms were found in this dataset during Task 3.2b/3.3 investigation,
none of which is fully addressed by mmseqs2 clustering or genome-blocking
alone:
  1. Exact duplicates (offset-0 identical strings) -- 181 groups, 134 of them
     spanning DIFFERENT source genomes. mmseqs linclust caught 180/181;
     closed completely by an exact-string-match safety net.
  2. Near-duplicates via substitution only (e.g. 93.3% identity, no shift)
     -- NOT reliably caught by mmseqs linclust at ANY tested threshold
     (0.3/0.5/0.7), confirmed by cross-checking a specific pair (OLIGO
     29542/29549, 93.3% identity per Biopython global alignment, exact
     match to this script's own k-mer+alignment estimate) that remained in
     two separate singleton mmseqs clusters even at min-seq-id=0.7.
  3. "Shifted" near-duplicates via small indels near the 165bp extraction
     window boundary (e.g. OLIGO 20578/32225: not equal as strings, but
     100% identical once a 2bp register shift is applied -- consistent
     with a small indel in the genomic context between two related loci).

Given this, splitting is done in two layers:
  LAYER 1 (primary structural defense): GENOME-BLOCKED assignment. Every
    sequence's entire source genome is assigned to exactly one fold. This is
    chosen over sequence-cluster-blocking as primary because it is the
    stronger, more easily defended guarantee for a reviewer ("no two
    sequences from the same source organism ever appear on both sides of a
    split") and because Task 3.3 showed cluster-based folds still leak
    substantial cross-genome near-duplicate identity (up to 93.3%) that
    genome-blocking has no reason to prevent either -- so cluster-blocking
    bought no leakage-safety advantage over genome-blocking to justify its
    weaker guarantee.
  LAYER 2 (sequence-level safety net, mandatory given Layer 1 alone still
    showed max_identity=1.000 in Task 3.3): an all-pairs near-duplicate
    search via k-mer-index candidate screening + exact (indel-tolerant,
    shift-window) identity computation, covering the FULL 29,249-sequence
    library once. Any two sequences at or above IDENTITY_MERGE_THRESHOLD are
    unioned into the same connected component (via union-find), and every
    component that would otherwise span >1 fold is consolidated into a
    single fold (majority-share rule, same mechanism as the Task 3.2b
    exact-duplicate fix, now generalized).

IDENTITY_MERGE_THRESHOLD = 0.85 was chosen because: it is comfortably above
the highest CROSS-GENOME identity mmseqs clustering successfully used as an
operating threshold (0.7) while being low enough to catch both concretely
observed problem cases (93.3% substitution-only, 100% shifted) with margin.
It is NOT claimed to be a universal safe threshold -- it is validated
empirically below via the final per-fold max-identity report, which is the
number that actually matters.
"""
import pandas as pd
import numpy as np
from pathlib import Path
import json
from collections import defaultdict

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"
SPLITS = DATA / "splits"
N_FOLDS = 5
KMER_K = 10
IDENTITY_MERGE_THRESHOLD = 0.85
MIN_SHARED_KMERS = 5


def kmers(seq, k=KMER_K):
    return set(seq[i:i + k] for i in range(len(seq) - k + 1))


def exact_identity(a, b):
    best = 0
    L = min(len(a), len(b))
    for offset in range(-5, 6):
        matches = 0
        n = 0
        for i in range(L):
            j = i + offset
            if 0 <= j < len(b):
                n += 1
                if a[i] == b[j]:
                    matches += 1
        if n > 0:
            best = max(best, matches / n)
    return best


class UnionFind:
    def __init__(self, items):
        self.parent = {x: x for x in items}

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def all_pairs_near_duplicates(ids, seqs, threshold=IDENTITY_MERGE_THRESHOLD):
    """Full all-pairs near-duplicate search over the WHOLE library via a
    k-mer inverted index (built once) -- not limited to train-vs-test, since
    this defines the merge groups BEFORE fold assignment."""
    n = len(ids)
    kmer_sets = [kmers(s) for s in seqs]
    index = defaultdict(list)
    for i, kset in enumerate(kmer_sets):
        for km in kset:
            index[km].append(i)

    uf = UnionFind(ids)
    n_pairs_checked = 0
    n_merges = 0
    checked_pairs = set()
    for i in range(n):
        candidate_counts = defaultdict(int)
        for km in kmer_sets[i]:
            for j in index[km]:
                if j > i:
                    candidate_counts[j] += 1
        for j, c in candidate_counts.items():
            if c < MIN_SHARED_KMERS:
                continue
            pair_key = (i, j)
            if pair_key in checked_pairs:
                continue
            checked_pairs.add(pair_key)
            n_pairs_checked += 1
            ident = exact_identity(seqs[i], seqs[j])
            if ident >= threshold:
                uf.union(ids[i], ids[j])
                n_merges += 1
    return uf, n_pairs_checked, n_merges


def main():
    df = pd.read_parquet(DATA / "three_host.parquet")
    df["OLIGO ID"] = df["OLIGO ID"].astype(str)
    ids = df["OLIGO ID"].tolist()
    seqs = df["regulatory_sequence"].tolist()
    id_to_genome = dict(zip(df["OLIGO ID"], df["genome_id"]))
    id_to_seq = dict(zip(ids, seqs))

    print(f"Running all-pairs near-duplicate search (k={KMER_K}, min_shared_kmers={MIN_SHARED_KMERS}, "
          f"merge threshold={IDENTITY_MERGE_THRESHOLD}) over {len(ids)} sequences...")
    uf, n_checked, n_merges = all_pairs_near_duplicates(ids, seqs)
    print(f"  {n_checked} candidate pairs checked, {n_merges} pairs merged (>= {IDENTITY_MERGE_THRESHOLD} identity)")

    components = defaultdict(list)
    for oid in ids:
        components[uf.find(oid)].append(oid)
    comp_sizes = sorted((len(v) for v in components.values()), reverse=True)
    print(f"  {len(components)} connected components; largest={comp_sizes[0]}; "
          f"n components with >1 member: {sum(1 for s in comp_sizes if s > 1)}")

    # LAYER 1: genome-blocked assignment
    genome_counts = df["genome_id"].value_counts()
    genomes_sorted = genome_counts.index.tolist()
    fold_loads = np.zeros(N_FOLDS)
    genome_to_fold = {}
    for g in genomes_sorted:
        target = int(np.argmin(fold_loads))
        genome_to_fold[g] = target
        fold_loads[target] += genome_counts[g]
    member_to_fold = {oid: genome_to_fold[id_to_genome[oid]] for oid in ids}

    # LAYER 2: consolidate every near-duplicate connected component into one fold
    n_components_fixed = 0
    n_reassigned = 0
    for comp_root, members in components.items():
        if len(members) < 2:
            continue
        folds_present = set(member_to_fold[m] for m in members)
        if len(folds_present) > 1:
            n_components_fixed += 1
            fold_counts = pd.Series([member_to_fold[m] for m in members]).value_counts()
            target = int(fold_counts.idxmax())
            for m in members:
                if member_to_fold[m] != target:
                    n_reassigned += 1
                member_to_fold[m] = target

    print(f"\nLayer 2 consolidation: {n_components_fixed} near-duplicate components spanned >1 fold, "
          f"{n_reassigned} oligos reassigned")

    fold_loads_final = np.zeros(N_FOLDS)
    for f in member_to_fold.values():
        fold_loads_final[f] += 1
    print(f"Final fold loads: {fold_loads_final.astype(int).tolist()}")

    fold_df = pd.DataFrame({"OLIGO ID": list(member_to_fold.keys()), "fold": list(member_to_fold.values())})
    fold_df.to_parquet(SPLITS / "fold_assignment_FINAL.parquet", index=False)

    # verify: final max train-test identity, full check (not sampled)
    print(f"\nFinal per-fold max train-test identity (full, unsampled check):")

    def compute_max_identity_full(train_ids_l, test_ids_l):
        train_seqs_l = [id_to_seq[i] for i in train_ids_l]
        test_seqs_l = [id_to_seq[i] for i in test_ids_l]
        idx = defaultdict(list)
        for i, s in enumerate(train_seqs_l):
            for km in kmers(s):
                idx[km].append(i)
        best, best_pair, checked = 0.0, None, 0
        for tj, tseq in enumerate(test_seqs_l):
            counts = defaultdict(int)
            for km in kmers(tseq):
                for ti in idx.get(km, ()):
                    counts[ti] += 1
            for ti, c in counts.items():
                if c < MIN_SHARED_KMERS:
                    continue
                checked += 1
                ident = exact_identity(train_seqs_l[ti], tseq)
                if ident > best:
                    best, best_pair = ident, (train_ids_l[ti], test_ids_l[tj])
        return best, best_pair, checked

    final_results = []
    for test_fold in range(N_FOLDS):
        test_ids = fold_df.loc[fold_df.fold == test_fold, "OLIGO ID"].tolist()
        train_ids = fold_df.loc[fold_df.fold != test_fold, "OLIGO ID"].tolist()
        max_id, pair, n_c = compute_max_identity_full(train_ids, test_ids)
        print(f"  fold {test_fold}: test_n={len(test_ids)}, train_n={len(train_ids)}, "
              f"max_identity={max_id:.4f}, pair={pair}")
        final_results.append({"test_fold": test_fold, "test_n": len(test_ids), "train_n": len(train_ids),
                               "max_identity": max_id, "pair": pair})

    with open(OUT / "final_split_results.json", "w") as f:
        json.dump({
            "identity_merge_threshold": IDENTITY_MERGE_THRESHOLD,
            "n_candidate_pairs_checked_allpairs": n_checked,
            "n_pairs_merged": n_merges,
            "n_connected_components": len(components),
            "n_multi_member_components": sum(1 for s in comp_sizes if s > 1),
            "largest_component_size": comp_sizes[0],
            "n_components_spanning_multiple_folds_before_fix": n_components_fixed,
            "n_oligos_reassigned_layer2": n_reassigned,
            "final_fold_loads": fold_loads_final.astype(int).tolist(),
            "final_per_fold_max_identity": final_results,
        }, f, indent=2, default=str)
    print(f"\nWrote {OUT / 'final_split_results.json'}")


if __name__ == "__main__":
    main()
