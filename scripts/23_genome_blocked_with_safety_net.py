"""
GATE 2 - Task 3.3/3.4: Genome-blocked fold assignment WITH the exact-duplicate
safety net applied (the naive genome-blocked scheme in script 22 showed
max_identity=1.000 in 4/5 folds -- direct evidence that genome-blocking alone
does not prevent cross-genome exact-duplicate leakage, exactly as the 134
different-genome duplicate groups from Task 3.2b predict). This script
reruns genome-blocked assignment with the same majority-fold consolidation
fix used for cluster-based folds, then re-measures max identity to see
whether the residual leakage ceiling matches the cluster-based scheme's
~0.933, which would indicate a shared underlying cause (a family of
near-but-not-exact duplicate sequences neither scheme's safety net addresses).
"""
import pandas as pd
import numpy as np
from pathlib import Path
import json
from collections import defaultdict
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
leakmod = import_module("22_source_genome_leakage")

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"
SPLITS = DATA / "splits"
N_FOLDS = 5


def main():
    df = pd.read_parquet(DATA / "three_host.parquet")
    df["OLIGO ID"] = df["OLIGO ID"].astype(str)
    id_to_seq = dict(zip(df["OLIGO ID"], df["regulatory_sequence"]))
    id_to_genome = dict(zip(df["OLIGO ID"], df["genome_id"]))

    genome_counts = df["genome_id"].value_counts()
    genomes_sorted = genome_counts.index.tolist()
    fold_loads = np.zeros(N_FOLDS)
    genome_to_fold = {}
    for g in genomes_sorted:
        target = int(np.argmin(fold_loads))
        genome_to_fold[g] = target
        fold_loads[target] += genome_counts[g]

    member_to_fold = {oid: genome_to_fold[gid] for oid, gid in id_to_genome.items()}

    # exact-duplicate safety net (same logic as script 20's assign_folds)
    seq_to_members = defaultdict(list)
    for member_id in member_to_fold:
        seq_to_members[id_to_seq[member_id]].append(member_id)

    n_reassigned = 0
    n_groups_fixed = 0
    for seq, members in seq_to_members.items():
        if len(members) < 2:
            continue
        folds_present = set(member_to_fold[m] for m in members)
        if len(folds_present) > 1:
            n_groups_fixed += 1
            fold_counts = pd.Series([member_to_fold[m] for m in members]).value_counts()
            target = int(fold_counts.idxmax())
            for m in members:
                if member_to_fold[m] != target:
                    n_reassigned += 1
                member_to_fold[m] = target

    print(f"Exact-duplicate safety net on genome-blocked folds: {n_groups_fixed} groups fixed, "
          f"{n_reassigned} oligos reassigned")

    fold_loads_after = np.zeros(N_FOLDS)
    for f in member_to_fold.values():
        fold_loads_after[f] += 1
    print(f"Fold loads after safety net: {fold_loads_after.astype(int).tolist()}")

    fold_df = pd.DataFrame({"OLIGO ID": list(member_to_fold.keys()), "fold": list(member_to_fold.values())})
    fold_df.to_parquet(SPLITS / "fold_assignment_genome_blocked_safe.parquet", index=False)

    print(f"\nMax train-test identity, genome-blocked + safety net:")
    results = []
    for test_fold in range(N_FOLDS):
        test_ids = fold_df.loc[fold_df.fold == test_fold, "OLIGO ID"].tolist()
        train_ids = fold_df.loc[fold_df.fold != test_fold, "OLIGO ID"].tolist()
        max_id, pair, n_checked, n_cross = leakmod.compute_max_identity_cross_genome(
            train_ids, test_ids, id_to_seq, id_to_genome)
        print(f"  fold {test_fold}: test_n={len(test_ids)}, train_n={len(train_ids)}, "
              f"max_identity={max_id:.3f}, pair={pair}")
        results.append({"test_fold": test_fold, "test_n": len(test_ids), "train_n": len(train_ids),
                         "max_identity": max_id, "pair": pair})

    # active-fraction survival check
    df["genome_fold_safe"] = df["OLIGO ID"].map(dict(zip(fold_df["OLIGO ID"], fold_df["fold"])))
    print(f"\nActive-fraction survival check:")
    balance = []
    for f in range(N_FOLDS):
        sub = df[df["genome_fold_safe"] == f]
        rate = sub["usable_all3_transcription"].mean()
        print(f"  fold {f}: n={len(sub)}, usable_all3_transcription rate={rate:.3f}")
        balance.append({"fold": f, "n": len(sub), "usable_all3_transcription_rate": float(rate)})

    with open(OUT / "genome_blocked_safe_results.json", "w") as f:
        json.dump({
            "n_groups_fixed": n_groups_fixed, "n_reassigned": n_reassigned,
            "fold_loads": fold_loads_after.astype(int).tolist(),
            "max_identity_per_fold": results, "balance_check": balance,
        }, f, indent=2, default=str)
    print(f"\nWrote {OUT / 'genome_blocked_safe_results.json'}")


if __name__ == "__main__":
    main()
