"""
GATE 1.5 - Task A: Correct RS241 usability analysis for leave-one-host-out evaluation.

Gate 1 computed the SIX-WAY intersection (usable in all 6 hosts simultaneously),
which over-restricts: leave-one-host-out only requires a usable value in the held-out
host plus usable values in whichever hosts are used for training (not necessarily all
5 others). This script recomputes per-host counts, the full pairwise matrix, and the
actually-runnable held-out-host configurations.

"Usable value" definition: IDENTICAL to Gate 1 Task 1/2 -- the cell in the Log2
Transcription (or Log10 Translation) sheet of Supplementary Data Table 4 (RS241) is
not null. No additional threshold is reconstructable from this table (these are
already paper-processed log-scale values, not raw counts), so "usable" = present,
exactly as in Gate 1.
"""
import pandas as pd
from pathlib import Path
import json
import itertools

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"

HOSTS = ["EC", "BS", "PA", "SE", "VN", "CG"]
HOST_NAMES = {
    "EC": "E. coli", "BS": "B. subtilis", "PA": "P. aeruginosa",
    "SE": "S. enterica", "VN": "V. natriegens", "CG": "C. glutamicum",
}
PRIMARY_HOSTS = ["EC", "BS", "PA"]
ADDED_HOSTS = ["SE", "VN", "CG"]


def analyze(df, cols, readout_name):
    print(f"\n{'='*70}\n{readout_name}\n{'='*70}")

    masks = {h: df[cols[h]].notna() for h in HOSTS}

    # 1. Per-host usable counts
    print("\n-- 1. Per-host usable counts --")
    per_host = {}
    for h in HOSTS:
        n = int(masks[h].sum())
        per_host[h] = n
        print(f"  {h} ({HOST_NAMES[h]}): {n} / {len(df)}")

    # 2. Pairwise matrix
    print("\n-- 2. Pairwise usable-in-both matrix --")
    pairwise = {}
    header = "        " + "  ".join(f"{h:>6s}" for h in HOSTS)
    print(header)
    for h1 in HOSTS:
        row = {}
        row_str = f"{h1:>6s}: "
        for h2 in HOSTS:
            n = int((masks[h1] & masks[h2]).sum())
            row[h2] = n
            row_str += f"{n:>6d}  "
        pairwise[h1] = row
        print(row_str)

    # 3. Held-out-host viability
    print("\n-- 3. Held-out-host viability (target host + >=2 / >=3 OTHER hosts usable) --")
    viability = {}
    for target in HOSTS:
        others = [h for h in HOSTS if h != target]
        n_others_usable = pd.concat([masks[h] for h in others], axis=1).sum(axis=1)
        target_usable = masks[target]

        n_ge2 = int((target_usable & (n_others_usable >= 2)).sum())
        n_ge3 = int((target_usable & (n_others_usable >= 3)).sum())

        # which other hosts contribute most, among rows where target usable
        sub = df.loc[target_usable]
        contrib = {h: int(sub[cols[h]].notna().sum()) for h in others}
        contrib_sorted = dict(sorted(contrib.items(), key=lambda kv: -kv[1]))

        viability[target] = {
            "n_target_usable": int(target_usable.sum()),
            "n_target_and_ge2_others": n_ge2,
            "n_target_and_ge3_others": n_ge3,
            "other_host_contribution_when_target_usable": contrib_sorted,
        }
        tag = "PRIMARY" if target in PRIMARY_HOSTS else "ADDED"
        print(f"  [{tag}] {target} ({HOST_NAMES[target]}): "
              f"n_usable={target_usable.sum()}, "
              f"n_usable_and_>=2_others={n_ge2}, "
              f"n_usable_and_>=3_others={n_ge3}")
        print(f"      other-host coverage when {target} usable: {contrib_sorted}")

    # 4. Best achievable evaluation configurations: for each held-out host, find the
    # single best training-host-set (from among all non-empty subsets of the other 5)
    # that maximizes N = count of RSs usable in target AND in ALL chosen training hosts.
    # We report the full ranked list restricted to "realistic" training sets: all
    # 5 others, and the best subset of size >=2 (since a real LOHO run typically trains
    # on multiple/all remaining hosts, not just 1).
    print("\n-- 4. Best achievable evaluation configurations --")
    configs = []
    for target in HOSTS:
        others = [h for h in HOSTS if h != target]
        target_mask = masks[target]
        # try all non-empty subsets of "others" (5 hosts -> 31 subsets), find N for each
        best_for_target = []
        for r in range(1, len(others) + 1):
            for subset in itertools.combinations(others, r):
                combined = target_mask.copy()
                for h in subset:
                    combined = combined & masks[h]
                n = int(combined.sum())
                best_for_target.append({"held_out": target, "train_hosts": list(subset), "n": n})
        # keep: the all-5-training-hosts config, and the best config using >=2 training hosts,
        # and the best config using all others minus worst 1 (i.e. best 4-of-5)
        all5 = [c for c in best_for_target if len(c["train_hosts"]) == 5][0]
        best_ge2 = max([c for c in best_for_target if len(c["train_hosts"]) >= 2], key=lambda c: c["n"])
        best_4of5 = max([c for c in best_for_target if len(c["train_hosts"]) == 4], key=lambda c: c["n"])
        configs.append({"held_out": target, "config": "train_on_all_5_others", **all5})
        configs.append({"held_out": target, "config": "best_4_of_5_others", **best_4of5})
        configs.append({"held_out": target, "config": "best_subset_ge2_others", **best_ge2})

    configs_sorted = sorted(configs, key=lambda c: -c["n"])
    for c in configs_sorted:
        print(f"  held_out={c['held_out']:>3s} ({HOST_NAMES[c['held_out']]:15s}) "
              f"config={c['config']:25s} train={c['train_hosts']} N={c['n']}")

    return {
        "per_host_usable": per_host,
        "pairwise_matrix": pairwise,
        "held_out_viability": viability,
        "ranked_configs": configs_sorted,
    }


def main():
    df = pd.read_parquet(DATA / "rs241.parquet")
    print(f"RS241 table: {df.shape}")

    tx_cols = {h: f"tx_log2_{h}" for h in HOSTS}
    tl_cols = {h: f"tl_log10_{h}" for h in HOSTS}

    # Sanity check: definition identical to Gate 1 -- both are just .notna() on the
    # same columns used in scripts/04_task2_rs241_analysis.py. Re-verify the known
    # Gate 1 anchor numbers (per-host counts, all-6 intersection) reproduce exactly.
    print("\n" + "#" * 70)
    print("SANITY CHECK vs Gate 1 Task 2 (scripts/04_task2_rs241_analysis.py)")
    print("#" * 70)
    tx_all6 = pd.concat([df[c].notna() for c in tx_cols.values()], axis=1).all(axis=1).sum()
    tl_all6 = pd.concat([df[c].notna() for c in tl_cols.values()], axis=1).all(axis=1).sum()
    print(f"Recomputed all-6-hosts-usable: transcription={tx_all6} (Gate1 reported 111), "
          f"translation={tl_all6} (Gate1 reported 126)")
    assert tx_all6 == 111, f"MISMATCH: expected 111, got {tx_all6}"
    assert tl_all6 == 126, f"MISMATCH: expected 126, got {tl_all6}"
    print("MATCH CONFIRMED -- usable-value definition is identical to Gate 1 Task 1/2 "
          "(non-null cell in the paper-processed log2/log10 sheets, no additional threshold).")

    results = {}
    results["transcription"] = analyze(df, tx_cols, "TRANSCRIPTION (Log2)")
    results["translation"] = analyze(df, tl_cols, "TRANSLATION (Log10)")

    with open(OUT / "task_a_rs241_pairwise_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nWrote {OUT / 'task_a_rs241_pairwise_results.json'}")


if __name__ == "__main__":
    main()
