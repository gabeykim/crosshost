"""Maintainer-run scoring script for held-out submissions (no hosted server
yet -- see README.md in this directory for the submission process).

Usage:
    python score_submission.py path/to/submission.parquet

submission.parquet must have columns: 'OLIGO ID', and for each readout you
are scoring, '{readout}_active_prob' and '{readout}_strength_pred' (e.g.
'transcription_active_prob', 'transcription_strength_pred'). Missing a
readout's columns skips scoring for that readout.

Requires the PRIVATE labels file (out/gate8_held_out_eval_private_labels.parquet
in the source repository -- NOT distributed with this package). Maintainer-only.
"""
import sys
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from crosshost.evaluate import evaluate

PRIVATE_LABELS = Path(__file__).resolve().parent.parent.parent / "out" / "gate8_held_out_eval_private_labels.parquet"

READOUT_TARGET_COLS = {
    "transcription": {"active": "{h}_tx_active", "strength": "{h}_tx_norm", "usable": "{h}_tx_usable"},
    "translation": {"active": None, "strength": "{h}_protein_log10", "usable": "{h}_tl_usable"},
}
HOSTS = ["EC", "BS", "PA"]


def main(submission_path):
    if not PRIVATE_LABELS.exists():
        raise FileNotFoundError(
            f"{PRIVATE_LABELS} not found. This script is maintainer-only and requires the "
            f"private held-out labels file, which is never distributed with this package."
        )
    labels = pd.read_parquet(PRIVATE_LABELS)
    labels["OLIGO ID"] = labels["OLIGO ID"].astype(str)
    sub = pd.read_parquet(submission_path)
    sub["OLIGO ID"] = sub["OLIGO ID"].astype(str)
    merged = labels.merge(sub, on="OLIGO ID", how="inner")
    print(f"Scoring {len(merged)} / {len(labels)} held-out rows with a submitted prediction "
          f"({len(labels) - len(merged)} rows in the held-out set were not covered by this submission).")

    for readout in ["transcription", "translation"]:
        prob_col = f"{readout}_active_prob"
        strength_col = f"{readout}_strength_pred"
        if prob_col not in merged.columns and strength_col not in merged.columns:
            print(f"  {readout}: not scored (no prediction columns found)")
            continue
        # per-host, since the held-out set spans all 3 primary hosts and each host has its own label columns
        for host in HOSTS:
            spec = READOUT_TARGET_COLS[readout]
            usable_col = spec["usable"].format(h=host)
            usable = merged[usable_col].fillna(False).astype(bool)
            if usable.sum() == 0:
                continue
            sub_h = merged[usable]
            preds_per_fold = [{readout: {
                "active_prob": sub_h[prob_col].values if prob_col in sub_h.columns else np.full(len(sub_h), 0.5),
                "strength_pred": sub_h[strength_col].values if strength_col in sub_h.columns else np.zeros(len(sub_h)),
            }}]
            active_true = sub_h[spec["active"].format(h=host)].values if spec["active"] else (sub_h[spec["strength"].format(h=host)].values > 0)
            strength_true = sub_h[spec["strength"].format(h=host)].values
            targets_per_fold = [{readout: {"active_true": active_true, "strength_true": strength_true}}]
            result = evaluate(preds_per_fold, targets_per_fold, readouts=(readout,))
            r = result[readout]["spearman_rho"]
            print(f"  {readout} {host}: n={usable.sum()}, spearman_rho={r['mean']}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python score_submission.py path/to/submission.parquet")
        sys.exit(1)
    main(sys.argv[1])
