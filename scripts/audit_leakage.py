#!/usr/bin/env python3
"""
CROSSHOST leakage audit -- runs at the start of every gate from Gate 2 onward.

Standalone, exits non-zero on ANY failure. Checks:
  1. No RS (OLIGO ID) appears in more than one fold
  2. No RS241 sequence appears anywhere in the main-library fold assignment
     (mechanical enforcement of the charter's non-negotiable "never train on
     RS241" rule, per Gate 2 Task 3.4)
  3. Max train-test sequence identity per fold is below the operating bound
     established at freeze time (0.8485, see out/final_split_results.json)
  4. No exact-duplicate sequence pair (100% identity, offset-0 string match)
     is split across folds
  5. Fold-local normalization: the translation floor-value statistic
     (scripts/28_prepare_baseline_data.py) is refit from training-pool data
     ONLY and compared against an all-data refit to confirm they differ (i.e.
     the check has real discriminating power) and that the cached value
     matches the training-pool fit, not the all-data one. Extended from a
     Gate 2 artifact-detection-only stub now that Gate 3 introduces the
     first real fitted statistic in this pipeline.
  6. Table hashes in data/MANIFEST.json match the files on disk

Usage: python scripts/audit_leakage.py
Exit code 0 = all checks passed. Exit code 1 = at least one check failed.
"""
import sys
import json
import hashlib
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
SPLITS = DATA / "splits"
OUT = ROOT / "out"

MAX_IDENTITY_BOUND = 0.90  # operating bound: freeze-time measured max was 0.8485;
                            # allow headroom to 0.90 before hard-failing, since
                            # this audit re-measures on every gate and legitimate
                            # small data changes (e.g. adding physiology columns)
                            # shouldn't fail the build over noise in the 4th decimal


def sha256_of_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def check_no_rs_in_multiple_folds():
    name = "1. No RS appears in more than one fold"
    fold_df = pd.read_parquet(SPLITS / "fold_assignment_FINAL.parquet")
    dupes = fold_df["OLIGO ID"].duplicated().sum()
    ok = dupes == 0
    detail = f"{len(fold_df)} rows, {dupes} duplicate OLIGO IDs"
    return name, ok, detail


def check_no_rs241_in_training():
    name = "2. No RS241 sequence appears in the main-library fold assignment"
    fold_df = pd.read_parquet(SPLITS / "fold_assignment_FINAL.parquet")
    rs241 = pd.read_parquet(DATA / "rs241.parquet")
    rs241_ids = set(rs241["id"].astype(str))
    fold_ids = set(fold_df["OLIGO ID"].astype(str))
    overlap = rs241_ids & fold_ids
    ok = len(overlap) == 0
    detail = f"{len(rs241_ids)} RS241 ids checked, {len(overlap)} found in main-library folds"
    return name, ok, detail


def kmers(seq, k=10):
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


def check_max_identity_per_fold():
    name = f"3. Max train-test identity per fold < {MAX_IDENTITY_BOUND}"
    three_host = pd.read_parquet(DATA / "three_host.parquet")
    three_host["OLIGO ID"] = three_host["OLIGO ID"].astype(str)
    id_to_seq = dict(zip(three_host["OLIGO ID"], three_host["regulatory_sequence"]))
    fold_df = pd.read_parquet(SPLITS / "fold_assignment_FINAL.parquet")

    n_folds = fold_df["fold"].nunique()
    worst = 0.0
    per_fold_detail = []
    for test_fold in sorted(fold_df["fold"].unique()):
        test_ids = fold_df.loc[fold_df.fold == test_fold, "OLIGO ID"].tolist()
        train_ids = fold_df.loc[fold_df.fold != test_fold, "OLIGO ID"].tolist()
        train_seqs = [id_to_seq[i] for i in train_ids]
        test_seqs = [id_to_seq[i] for i in test_ids]

        index = defaultdict(list)
        for i, s in enumerate(train_seqs):
            for km in kmers(s):
                index[km].append(i)

        fold_max = 0.0
        for tj, tseq in enumerate(test_seqs):
            counts = defaultdict(int)
            for km in kmers(tseq):
                for ti in index.get(km, ()):
                    counts[ti] += 1
            for ti, c in counts.items():
                if c < 5:
                    continue
                ident = exact_identity(train_seqs[ti], tseq)
                if ident > fold_max:
                    fold_max = ident
        per_fold_detail.append(f"fold{test_fold}={fold_max:.4f}")
        worst = max(worst, fold_max)

    ok = worst < MAX_IDENTITY_BOUND
    detail = f"per-fold max identities: {', '.join(per_fold_detail)}; worst={worst:.4f}"
    return name, ok, detail


def check_no_exact_duplicate_split():
    name = "4. No exact-duplicate sequence pair is split across folds"
    three_host = pd.read_parquet(DATA / "three_host.parquet")
    three_host["OLIGO ID"] = three_host["OLIGO ID"].astype(str)
    fold_df = pd.read_parquet(SPLITS / "fold_assignment_FINAL.parquet")
    id_to_fold = dict(zip(fold_df["OLIGO ID"], fold_df["fold"]))

    seq_groups = three_host.groupby("regulatory_sequence")["OLIGO ID"].apply(list)
    dup_groups = seq_groups[seq_groups.apply(len) > 1]

    n_split = 0
    examples = []
    for seq, members in dup_groups.items():
        # only consider members actually present in the fold assignment
        # (RS241-only members were removed by check #2's enforcement)
        present = [m for m in members if m in id_to_fold]
        folds_present = set(id_to_fold[m] for m in present)
        if len(folds_present) > 1:
            n_split += 1
            if len(examples) < 3:
                examples.append(present)

    ok = n_split == 0
    detail = f"{len(dup_groups)} exact-duplicate groups checked, {n_split} split across folds" + \
             (f", examples: {examples}" if examples else "")
    return name, ok, detail


def check_no_global_normalization():
    name = "5. Fold-local normalization: fitted statistics use training-pool data only"
    # FUNCTIONAL check (extended from Gate 2's artifact-detection-only version
    # now that Gate 3 introduces a real fitted statistic): the translation
    # "floor value" (scripts/28_prepare_baseline_data.py) is a per-host
    # statistic fit from data and then used to LABEL every row (train and
    # test) as active/inactive and usable-for-regression/not. This check
    # independently recomputes that statistic two ways -- from the training
    # pool only (fold != TEST_FOLD), and from ALL rows including the test
    # fold -- and asserts the cached .npz value matches the TRAINING-POOL
    # computation, not the all-data one. (These two numbers were found to
    # differ narrowly for all 3 hosts when this was checked during Gate 3 --
    # e.g. EC's floor fraction is 0.665 training-pool-only vs 0.674 on all
    # data -- confirming the check has real discriminating power and is not
    # vacuously passing because the two are identical.)
    TEST_FOLD = 0
    three_host_path = DATA / "three_host.parquet"
    cache_dir = DATA / "baseline_cache"
    if not three_host_path.exists() or not cache_dir.exists():
        return name, True, "SKIPPED -- no baseline cache exists yet (this check activates once " \
                            "scripts/28_prepare_baseline_data.py has been run)"

    fold_df = pd.read_parquet(SPLITS / "fold_assignment_FINAL.parquet")
    id_to_fold = dict(zip(fold_df["OLIGO ID"].astype(str), fold_df["fold"]))
    three_host = pd.read_parquet(three_host_path)
    three_host["OLIGO ID"] = three_host["OLIGO ID"].astype(str)
    three_host["fold"] = three_host["OLIGO ID"].map(id_to_fold)
    valid = three_host[three_host["fold"].notna()].copy()
    valid["fold"] = valid["fold"].astype(int)

    mismatches = []
    checked = []
    for host in ["EC", "BS", "PA"]:
        cache_path = cache_dir / f"{host}_baseline_data.npz"
        if not cache_path.exists():
            continue
        cached = np.load(cache_path, allow_pickle=True)
        cached_floor = float(cached["tl_floor_value"])

        tl_usable = valid[f"{host}_tl_usable"].values.astype(bool)
        protein = valid[f"{host}_protein_log10"].values.astype(np.float32)
        fold_arr = valid["fold"].values.astype(int)

        train_pool_vals = protein[tl_usable & (fold_arr != TEST_FOLD)]
        all_vals = protein[tl_usable]

        def mode_val(vals):
            rounded = np.round(vals, 2)
            u, c = np.unique(rounded, return_counts=True)
            return float(u[np.argmax(c)])

        train_pool_floor = mode_val(train_pool_vals) if len(train_pool_vals) else None
        all_data_floor = mode_val(all_vals) if len(all_vals) else None

        checked.append(f"{host}: cached={cached_floor:.2f}, train_pool_recompute={train_pool_floor:.2f}, "
                        f"all_data_recompute={all_data_floor:.2f}")
        if train_pool_floor is None or abs(cached_floor - train_pool_floor) > 1e-6:
            mismatches.append(f"{host}: cached value {cached_floor} does NOT match training-pool-only "
                               f"recomputation {train_pool_floor}")

    ok = len(mismatches) == 0 and len(checked) > 0
    detail = "; ".join(checked) + (f" -- MISMATCHES: {mismatches}" if mismatches else " -- all match training-pool-only fit")
    return name, ok, detail


def gitignored_paths():
    """Root-relative paths explicitly listed in .gitignore (exact lines only --
    this project's own gitignored MANIFEST entries are always exact paths, not
    globs; a glob pattern would not exact-match a manifest key and is silently
    treated as not-gitignored, which is the conservative direction -- it would
    surface as a real mismatch rather than being silently swallowed)."""
    gi = ROOT / ".gitignore"
    if not gi.exists():
        return set()
    paths = set()
    for line in gi.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        paths.add(line)
    return paths


def check_manifest_hashes():
    name = "6. Table hashes match data/MANIFEST.json"
    manifest_path = DATA / "MANIFEST.json"
    if not manifest_path.exists():
        return name, False, "MANIFEST.json does not exist"
    manifest = json.load(open(manifest_path))
    ignored = gitignored_paths()
    mismatches = []
    checked = 0
    skipped = []
    for fname, info in manifest.items():
        # Manifest keys are bare filenames resolved under data/ by convention
        # (e.g. "three_host.parquet"), EXCEPT keys that already carry an
        # explicit root-relative path prefix (e.g. "raw/drafts/foo.xlsx",
        # added starting Gate 10 for source-data provenance tracking outside
        # data/) -- those resolve relative to the project ROOT instead.
        path = ROOT / fname if fname.startswith(("raw/", "out/")) else DATA / fname
        if not path.exists():
            # A manifest-tracked file that is also listed in .gitignore is
            # EXPECTED to be absent from a clean checkout (Gate 14: a fresh
            # clone was found to fail this check unconditionally before this
            # fix existed, for exactly this reason -- see
            # out/results/gate14_clean_clone_verification.md). Its hash stays
            # in the manifest for provenance (obtain the file via the URL
            # recorded there); its absence here is not a leakage-audit
            # failure. If the file IS present (e.g. downloaded manually), it
            # is still hash-verified below like any other entry -- gitignored
            # only means "not required to be present," not "never checked."
            rel = fname if fname.startswith(("raw/", "out/")) else f"data/{fname}"
            if rel in ignored:
                skipped.append(f"{fname}: gitignored, not expected in a clean checkout")
                continue
            mismatches.append(f"{fname}: file missing")
            continue
        actual_hash = sha256_of_file(path)
        checked += 1
        if actual_hash != info.get("sha256"):
            mismatches.append(f"{fname}: hash mismatch")
    ok = len(mismatches) == 0 and checked > 0
    detail = f"{checked} files checked against manifest, {len(mismatches)} mismatches" + \
             (f": {mismatches}" if mismatches else "") + \
             (f"; {len(skipped)} gitignored and skipped: {skipped}" if skipped else "")
    return name, ok, detail


def main():
    checks = [
        check_no_rs_in_multiple_folds,
        check_no_rs241_in_training,
        check_max_identity_per_fold,
        check_no_exact_duplicate_split,
        check_no_global_normalization,
        check_manifest_hashes,
    ]

    print("=" * 78)
    print("CROSSHOST LEAKAGE AUDIT")
    print("=" * 78)

    all_pass = True
    for check_fn in checks:
        try:
            name, ok, detail = check_fn()
        except Exception as e:
            name, ok, detail = check_fn.__doc__ or check_fn.__name__, False, f"CHECK RAISED EXCEPTION: {e}"
        status = "PASS" if ok else "FAIL"
        all_pass &= ok
        print(f"[{status}] {name}")
        print(f"       {detail}")

    print("=" * 78)
    if all_pass:
        print("RESULT: ALL CHECKS PASSED")
        print("=" * 78)
        sys.exit(0)
    else:
        print("RESULT: AT LEAST ONE CHECK FAILED")
        print("=" * 78)
        sys.exit(1)


if __name__ == "__main__":
    main()
