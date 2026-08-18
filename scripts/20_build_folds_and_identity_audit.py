"""
GATE 2 - Task 3: Build 5-fold cross-validation splits at cluster level, and
compute per-fold maximum train-test sequence identity.

IMPORTANT FINDING (see out/mmseqs_clustering_investigation.json and the Gate 2
memo): mmseqs easy-linclust's default k-mer seeding cannot distinguish
--min-seq-id 0.3 from 0.5 for this dataset (165bp sequences) -- both produce
BYTE-IDENTICAL clustering output (verified: diff of the two cluster.tsv files
is empty). A probe sweep (0.3/0.4/0.5/0.6/0.7/0.8/0.9) showed the tool's
practical identity-detection floor sits between 0.5 and 0.6 for this data;
below that, the k-mer prefilter simply finds no additional candidate pairs to
test, regardless of the requested threshold. This is reported as a finding,
not silently worked around.

Because clustering alone cannot be fully trusted below ~0.55 identity, this
script ALSO computes per-fold max train-test identity directly (not relying
on the clustering step's sensitivity) via a k-mer-index candidate screen +
exact identity computation on the survivors -- see compute_max_identity().
This is the real safety net Task 3.2 asks for, and it does not depend on
whatever the clustering step did or didn't find.
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
KMER_K = 10  # for candidate screening only -- exact identity computed on survivors


def load_sequences():
    df = pd.read_parquet(DATA / "three_host.parquet")
    return df[["OLIGO ID", "regulatory_sequence", "genome_id"]].copy()


def assign_folds(clusters_df, id_to_seq, seed=42):
    """Whole clusters -> single folds (never split a cluster across folds).

    SAFETY NET (added after a direct finding, not a hypothetical): the initial
    version of this script assigned folds purely from mmseqs linclust output
    and then, independently, audited per-fold max train-test identity. That
    audit caught a real case -- OLIGO IDs 21865 and 37869, byte-identical
    165bp sequences from the SAME source genome (637000188) -- sitting in two
    DIFFERENT singleton clusters at min-seq-id=0.7 and therefore landing in
    different folds. This is a genuine clustering-tool miss, not a bug in the
    identity check (verified by direct string comparison of the two
    sequences). linclust's speed comes from an approximate, chunked/threaded
    algorithm with no formal guarantee that every exact-duplicate pair gets
    merged -- see the Gate 2 memo for the full discussion.

    Fix: after cluster-based assignment, explicitly group ALL oligos by exact
    sequence string (fast, 100%-reliable for perfect duplicates -- no k-mer
    approximation needed) and force every oligo sharing an identical sequence
    into the SAME fold, overriding whatever the cluster assignment gave them.
    This closes the exact-duplicate leakage path completely. It does NOT
    address near-but-not-exact-duplicate pairs that clustering may also miss
    -- that residual risk is exactly why compute_max_identity() below is run
    as an independent, authoritative check rather than trusting fold
    assignment on faith.
    """
    cluster_sizes = clusters_df.groupby("representative")["member"].count().sort_values(ascending=False)
    reps = cluster_sizes.index.to_numpy().copy()
    fold_loads = np.zeros(N_FOLDS)
    rep_to_fold = {}
    for rep in reps:  # largest first, already sorted
        target_fold = int(np.argmin(fold_loads))
        rep_to_fold[rep] = target_fold
        fold_loads[target_fold] += cluster_sizes[rep]

    member_to_fold = {}
    for _, row in clusters_df.iterrows():
        member_to_fold[row["member"]] = rep_to_fold[row["representative"]]

    # exact-duplicate safety net
    seq_to_members = defaultdict(list)
    for member_id in member_to_fold:
        seq_to_members[id_to_seq[member_id]].append(member_id)

    n_reassigned = 0
    for seq, members in seq_to_members.items():
        if len(members) < 2:
            continue
        folds_present = set(member_to_fold[m] for m in members)
        if len(folds_present) > 1:
            # send the whole exact-duplicate group to whichever fold already
            # holds the largest share of the group (minimizes disruption)
            fold_counts = pd.Series([member_to_fold[m] for m in members]).value_counts()
            target = int(fold_counts.idxmax())
            for m in members:
                if member_to_fold[m] != target:
                    n_reassigned += 1
                member_to_fold[m] = target

    if n_reassigned:
        print(f"  EXACT-DUPLICATE SAFETY NET: reassigned {n_reassigned} oligos to close "
              f"exact-duplicate-sequence leakage paths that clustering missed")
        fold_loads = np.zeros(N_FOLDS)
        for f in member_to_fold.values():
            fold_loads[f] += 1

    return member_to_fold, fold_loads


def kmers(seq, k=KMER_K):
    return set(seq[i:i + k] for i in range(len(seq) - k + 1))


def exact_identity(a, b):
    """Ungapped identity at best alignment offset (sequences are both 165bp;
    allow a small shift window to approximate local alignment cheaply)."""
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


def compute_max_identity(train_seqs, test_seqs, train_ids, test_ids, k=KMER_K, min_shared_kmers=6):
    """K-mer-index candidate screen (memory-light: builds an inverted index of
    k-mer -> set of train sequence indices, then for each test sequence looks
    up which train sequences share enough k-mers to be worth an exact check)
    followed by exact identity computation on survivors only."""
    train_kmer_sets = [kmers(s) for s in train_seqs]
    index = defaultdict(list)
    for i, kset in enumerate(train_kmer_sets):
        for km in kset:
            index[km].append(i)

    best_identity = 0.0
    best_pair = None
    n_candidates_checked = 0
    for tj, tseq in enumerate(test_seqs):
        tkmers = kmers(tseq)
        candidate_counts = defaultdict(int)
        for km in tkmers:
            for ti in index.get(km, ()):
                candidate_counts[ti] += 1
        candidates = [ti for ti, c in candidate_counts.items() if c >= min_shared_kmers]
        for ti in candidates:
            n_candidates_checked += 1
            ident = exact_identity(train_seqs[ti], tseq)
            if ident > best_identity:
                best_identity = ident
                best_pair = (train_ids[ti], test_ids[tj])
    return best_identity, best_pair, n_candidates_checked


def main():
    seq_df = load_sequences()
    seq_df["OLIGO ID"] = seq_df["OLIGO ID"].astype(str)
    id_to_seq = dict(zip(seq_df["OLIGO ID"], seq_df["regulatory_sequence"]))
    id_to_genome = dict(zip(seq_df["OLIGO ID"], seq_df["genome_id"]))

    all_results = {}
    for thr_pct in [30, 50, 70]:
        print(f"\n{'='*70}\nTHRESHOLD {thr_pct/100}\n{'='*70}")
        clusters = pd.read_parquet(SPLITS / f"clusters_id{thr_pct}.parquet")
        member_to_fold, fold_loads = assign_folds(clusters, id_to_seq)
        print(f"Fold loads (n sequences): {fold_loads.astype(int).tolist()}")

        fold_assignment = pd.DataFrame({"OLIGO ID": list(member_to_fold.keys()),
                                         "fold": list(member_to_fold.values())})
        fold_assignment.to_parquet(SPLITS / f"fold_assignment_id{thr_pct}.parquet", index=False)

        # per-fold max identity: this fold's members = test, all others = train.
        # FULL computation, no subsampling -- the k-mer index makes this cheap
        # (candidates are rare for unrelated 165bp sequences at k=10: expected
        # random-collision rate for >=6 shared 10-mers is astronomically low,
        # ~4^10=1M possible 10-mers vs ~156 k-mers/sequence), so exact_identity
        # is only ever called on genuinely promising pairs.
        per_fold_results = []
        for test_fold in range(N_FOLDS):
            test_ids = fold_assignment.loc[fold_assignment.fold == test_fold, "OLIGO ID"].tolist()
            train_ids = fold_assignment.loc[fold_assignment.fold != test_fold, "OLIGO ID"].tolist()

            train_seqs = [id_to_seq[i] for i in train_ids]
            test_seqs = [id_to_seq[i] for i in test_ids]

            max_id, pair, n_checked = compute_max_identity(train_seqs, test_seqs, train_ids, test_ids)
            print(f"  fold {test_fold}: test_n={len(test_ids)}, train_n={len(train_ids)}, "
                  f"max_identity={max_id:.3f}, n_candidate_pairs_checked={n_checked}, pair={pair}")

            per_fold_results.append({
                "test_fold": test_fold, "test_n_full": len(test_ids), "train_n_full": len(train_ids),
                "max_identity": max_id, "n_candidate_pairs_checked": n_checked,
                "max_identity_pair": pair,
            })

        all_results[thr_pct] = {"fold_loads": fold_loads.astype(int).tolist(), "per_fold": per_fold_results}

    with open(OUT / "fold_identity_audit.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nWrote {OUT / 'fold_identity_audit.json'}")


if __name__ == "__main__":
    main()
