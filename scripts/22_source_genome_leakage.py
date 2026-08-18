"""
GATE 2 - Task 3.3: Source-genome leakage check.

The 29,249 sequences were mined from 184 prokaryotic genomes (Gate 1,
Metadata sheet). Sequence-identity clustering may not fully break relatedness
between sequences from related genomes -- Task 3.2b already found 134/181
exact-duplicate groups span DIFFERENT source genomes, so this is a live risk,
not a hypothetical.
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


def compute_max_identity_cross_genome(train_ids, test_ids, id_to_seq, id_to_genome, min_shared_kmers=6):
    train_seqs = [id_to_seq[i] for i in train_ids]
    test_seqs = [id_to_seq[i] for i in test_ids]
    train_kmer_sets = [kmers(s) for s in train_seqs]
    index = defaultdict(list)
    for i, kset in enumerate(train_kmer_sets):
        for km in kset:
            index[km].append(i)

    best_identity = 0.0
    best_pair = None
    n_checked = 0
    n_cross_genome_checked = 0
    for tj, tseq in enumerate(test_seqs):
        tkmers = kmers(tseq)
        candidate_counts = defaultdict(int)
        for km in tkmers:
            for ti in index.get(km, ()):
                candidate_counts[ti] += 1
        candidates = [ti for ti, c in candidate_counts.items() if c >= min_shared_kmers]
        for ti in candidates:
            n_checked += 1
            if id_to_genome[train_ids[ti]] == id_to_genome[test_ids[tj]]:
                continue  # same-genome pair, not what this check is for
            n_cross_genome_checked += 1
            ident = exact_identity(train_seqs[ti], tseq)
            if ident > best_identity:
                best_identity = ident
                best_pair = (train_ids[ti], test_ids[tj], id_to_genome[train_ids[ti]], id_to_genome[test_ids[tj]])
    return best_identity, best_pair, n_checked, n_cross_genome_checked


def main():
    df = pd.read_parquet(DATA / "three_host.parquet")
    df["OLIGO ID"] = df["OLIGO ID"].astype(str)
    meta = pd.read_parquet(DATA / "three_host_library.parquet")
    meta["OLIGO ID"] = meta["OLIGO ID"].astype(str)
    tax_cols = meta.set_index("OLIGO ID")[["phylum", "class", "family", "genus", "order"]]

    id_to_seq = dict(zip(df["OLIGO ID"], df["regulatory_sequence"]))
    id_to_genome = dict(zip(df["OLIGO ID"], df["genome_id"]))

    # 1. genome distribution
    print("=" * 70)
    print("1. SOURCE GENOME DISTRIBUTION")
    print("=" * 70)
    genome_counts = df["genome_id"].value_counts()
    print(f"n distinct source genomes: {genome_counts.nunique() if False else len(genome_counts)}")
    print(f"total sequences: {len(df)}")
    print(f"\nTop 10 genomes by sequence count:")
    print(genome_counts.head(10).to_string())
    print(f"\nTail shape: median={genome_counts.median():.0f}, "
          f"genomes contributing <10 sequences: {(genome_counts < 10).sum()}, "
          f"genomes contributing 1 sequence: {(genome_counts == 1).sum()}")
    top1_share = genome_counts.iloc[0] / len(df) * 100
    top10_share = genome_counts.head(10).sum() / len(df) * 100
    print(f"top-1 genome share: {top1_share:.2f}% of library, top-10 genomes share: {top10_share:.2f}%")

    # 3. phylogenetic breakdown (using phylum/class from Gate 1 Metadata)
    print("\n" + "=" * 70)
    print("3. PHYLOGENETIC / TAXONOMIC SPREAD")
    print("=" * 70)
    genome_phylum = meta.drop_duplicates("genome_id").set_index("genome_id")["phylum"]
    genome_class = meta.drop_duplicates("genome_id").set_index("genome_id")["class"]
    phylum_counts = genome_phylum.value_counts()
    print(f"n distinct phyla represented among source genomes: {phylum_counts.nunique() if False else len(phylum_counts)}")
    print(f"\nGenomes per phylum:")
    print(phylum_counts.to_string())
    class_counts = genome_class.value_counts()
    print(f"\nGenomes per class (top 15):")
    print(class_counts.head(15).to_string())

    # 2. cross-genome max identity per fold per threshold (cluster-based splits)
    print("\n" + "=" * 70)
    print("2. MAX TRAIN-TEST IDENTITY, DIFFERENT-SOURCE-GENOME PAIRS ONLY (cluster-based folds)")
    print("=" * 70)
    cross_genome_results = {}
    for thr_pct in [30, 50, 70]:
        fold_assignment = pd.read_parquet(SPLITS / f"fold_assignment_id{thr_pct}.parquet")
        print(f"\n--- threshold {thr_pct/100} ---")
        per_fold = []
        for test_fold in range(N_FOLDS):
            test_ids = fold_assignment.loc[fold_assignment.fold == test_fold, "OLIGO ID"].tolist()
            train_ids = fold_assignment.loc[fold_assignment.fold != test_fold, "OLIGO ID"].tolist()
            max_id, pair, n_checked, n_cross = compute_max_identity_cross_genome(
                train_ids, test_ids, id_to_seq, id_to_genome)
            print(f"  fold {test_fold}: max_cross_genome_identity={max_id:.3f}, "
                  f"n_candidates={n_checked}, n_cross_genome_candidates={n_cross}, pair={pair}")
            per_fold.append({"test_fold": test_fold, "max_cross_genome_identity": max_id,
                              "pair": pair, "n_candidates": n_checked, "n_cross_genome_candidates": n_cross})
        cross_genome_results[thr_pct] = per_fold

    # 4. source-genome-blocked splitting alternative
    print("\n" + "=" * 70)
    print("4. SOURCE-GENOME-BLOCKED SPLITTING (alternative scheme)")
    print("=" * 70)
    rng = np.random.default_rng(42)
    genomes_sorted = genome_counts.index.tolist()  # largest first, already sorted by value_counts
    fold_loads = np.zeros(N_FOLDS)
    genome_to_fold = {}
    for g in genomes_sorted:
        target = int(np.argmin(fold_loads))
        genome_to_fold[g] = target
        fold_loads[target] += genome_counts[g]
    print(f"Fold loads (n sequences) under genome-blocked assignment: {fold_loads.astype(int).tolist()}")

    df["genome_fold"] = df["genome_id"].map(genome_to_fold)
    genome_fold_df = df[["OLIGO ID", "genome_fold"]].copy()
    genome_fold_df.to_parquet(SPLITS / "fold_assignment_genome_blocked.parquet", index=False)

    # class balance / active-fraction survival check under genome-blocked folds
    print(f"\nActive-fraction survival check per fold (usable_all3_transcription rate):")
    for f in range(N_FOLDS):
        sub = df[df["genome_fold"] == f]
        rate = sub["usable_all3_transcription"].mean()
        print(f"  fold {f}: n={len(sub)}, usable_all3_transcription rate={rate:.3f}")

    print(f"\nMax train-test identity under genome-blocked folds (all pairs, not just cross-genome "
          f"since folds ARE genome-blocked so every train-test pair is automatically cross-genome):")
    genome_blocked_results = []
    for test_fold in range(N_FOLDS):
        test_ids = df.loc[df.genome_fold == test_fold, "OLIGO ID"].tolist()
        train_ids = df.loc[df.genome_fold != test_fold, "OLIGO ID"].tolist()
        max_id, pair, n_checked, n_cross = compute_max_identity_cross_genome(
            train_ids, test_ids, id_to_seq, id_to_genome)
        print(f"  fold {test_fold}: test_n={len(test_ids)}, train_n={len(train_ids)}, "
              f"max_identity={max_id:.3f}, pair={pair}")
        genome_blocked_results.append({"test_fold": test_fold, "test_n": len(test_ids),
                                        "train_n": len(train_ids), "max_identity": max_id, "pair": pair})

    results = {
        "n_source_genomes": int(len(genome_counts)),
        "top10_genomes": genome_counts.head(10).to_dict(),
        "top1_share_pct": float(top1_share),
        "top10_share_pct": float(top10_share),
        "n_genomes_single_sequence": int((genome_counts == 1).sum()),
        "phylum_counts": phylum_counts.to_dict(),
        "class_counts_top15": class_counts.head(15).to_dict(),
        "cross_genome_max_identity_by_threshold": cross_genome_results,
        "genome_blocked_fold_loads": fold_loads.astype(int).tolist(),
        "genome_blocked_max_identity": genome_blocked_results,
    }
    with open(OUT / "source_genome_leakage_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nWrote {OUT / 'source_genome_leakage_results.json'}")


if __name__ == "__main__":
    main()
