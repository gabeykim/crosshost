"""
GATE 10 - Task 5: feasibility verdict for a future host-count sweep
(Gate 11's proposed design, NOT run or evaluated here). Descriptive
combinatorics only -- no models, no hypothesis test.
"""
import json
import itertools
import numpy as np
import pandas as pd
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

SPECIES = ["Ec", "Ef", "Se", "Ko", "Pa", "Pp", "Vn", "Bs", "Cg", "Ll"]


def n_usable_in_all(df, host_subset):
    mask = np.ones(len(df), dtype=bool)
    for h in host_subset:
        mask &= df[f"{h}_usable"].values
    return int(mask.sum())


def main():
    df = pd.read_parquet(DATA / "drafts.parquet")

    print("=" * 80)
    print("N sequences usable in ALL k hosts, by host count k -- max 200 subsets sampled per k for k>4 (C(10,5)=252 total, exhaustive up to k=6)")
    print("=" * 80)

    by_k = {}
    for k in range(2, 11):
        combos = list(itertools.combinations(SPECIES, k))
        if len(combos) > 300:
            rng = np.random.default_rng(0)
            idx = rng.choice(len(combos), size=300, replace=False)
            combos = [combos[i] for i in sorted(idx)]
        ns = [n_usable_in_all(df, c) for c in combos]
        by_k[k] = {"n_subsets_checked": len(combos), "n_total_possible_subsets": len(list(itertools.combinations(SPECIES, k))),
                   "min_N": int(min(ns)), "median_N": int(np.median(ns)), "max_N": int(max(ns)),
                   "best_subset": list(combos[int(np.argmax(ns))]), "best_subset_N": int(max(ns)),
                   "worst_subset": list(combos[int(np.argmin(ns))]), "worst_subset_N": int(min(ns))}
        print(f"k={k}: N range [{by_k[k]['min_N']}, {by_k[k]['max_N']}], median={by_k[k]['median_N']}, "
              f"best={by_k[k]['best_subset']} (N={by_k[k]['best_subset_N']}), "
              f"checked {by_k[k]['n_subsets_checked']}/{by_k[k]['n_total_possible_subsets']} subsets")

    # Volume-controlled sweep: total_N = k * per_host_N. Find the largest total_N
    # such that per_host_N = total_N/k is achievable (<= N usable in all k hosts,
    # for SOME k-host subset) across the FULL range k=2..10.
    print("\n" + "=" * 80)
    print("VOLUME-CONTROLLED SWEEP: largest total-N supporting the full k=2..10 range")
    print("=" * 80)
    # max achievable total_N at each k = k * max_N_at_k (using the BEST subset at each k)
    max_total_by_k = {k: by_k[k]["best_subset_N"] * k for k in range(2, 11)}
    for k in range(2, 11):
        print(f"  k={k}: max per-host-N={by_k[k]['best_subset_N']}, max total-N (k*perhost)={max_total_by_k[k]}")
    largest_feasible_total = min(max_total_by_k.values())
    limiting_k = min(max_total_by_k, key=max_total_by_k.get)
    print(f"\nLargest total-N supporting ALL k=2..10 (limited by k={limiting_k}, "
          f"max_N={by_k[limiting_k]['best_subset_N']}): total_N={largest_feasible_total}")
    print("Row counts at each k under this volume-controlled design:")
    volume_controlled_rows = {}
    for k in range(2, 11):
        per_host_n = largest_feasible_total // k
        volume_controlled_rows[k] = {"per_host_N": per_host_n, "total_N": per_host_n * k}
        print(f"  k={k}: {k} hosts x {per_host_n} sequences/host = {per_host_n*k} total")

    # Host-subset sampling: how many distinct k-host subsets are available at each k?
    n_subsets_by_k = {k: len(list(itertools.combinations(SPECIES, k))) for k in range(2, 11)}

    summary = {
        "max_training_hosts_available": 10,
        "n_usable_by_host_count": by_k,
        "volume_controlled_sweep": {
            "feasible": True,
            "largest_total_N_supporting_full_k2_to_k10_range": largest_feasible_total,
            "limiting_host_count": limiting_k,
            "per_k_row_counts": volume_controlled_rows,
            "caveat": "This is the theoretical max given the BEST subset at the limiting k; using a FIXED subset "
                      "across all k values (needed to separate 'how many hosts' from 'which hosts' cleanly) will "
                      "have a smaller usable N than this best-case figure -- see host-subset-sampling section.",
        },
        "n_distinct_subsets_by_host_count": n_subsets_by_k,
    }
    with open(RESULTS / "gate10_feasibility.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\nWrote {RESULTS / 'gate10_feasibility.json'}")


if __name__ == "__main__":
    main()
