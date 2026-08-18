"""
GATE 8 - Task 2b: build the withheld-label held-out evaluation split.

DESIGN DECISION, disclosed: RS241 cannot serve as a genuine withheld test
set for future submitters -- its labels have already been used to compute
and publish results throughout this project (Gates 5, 5.5, 6, 7;
out/PAPER_FRAMING.md cites RS241 zero-shot numbers directly). Withholding
already-published labels would be theater, not a real evaluation barrier.

Instead: fold 4 of the existing frozen 5-fold split (5,813 rows, arbitrary
but fixed choice -- any fold would do) is designated the official withheld
test fold. Its activity-value label columns are stripped from the PUBLIC
package (package/held_out_eval/test_inputs.parquet: sequence + host
metadata only). The true labels are written OUTSIDE the distributable
package directory (out/gate8_held_out_eval_private_labels.parquet) --
explicitly NOT for inclusion in any public release, git push, Zenodo
archive, or HF dataset upload. package/held_out_eval/README.md documents
this and the scoring protocol.

Folds 0-3 remain fully public (in data/core/three_host.parquet) for model
development -- this mirrors the standard practice (ProteinGym, DREAM
challenges) of public dev data + private test labels, adapted to a project
with no submission server: scoring is done by the maintainer running
score_submission.py against the private file on request.
"""
import pandas as pd
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = Path(__file__).resolve().parent.parent / "out"
PKG = Path(__file__).resolve().parent.parent / "package"
HELD_OUT_DIR = PKG / "held_out_eval"
HELD_OUT_DIR.mkdir(parents=True, exist_ok=True)

WITHHELD_FOLD = 4


def main():
    lib = pd.read_parquet(DATA / "three_host.parquet")
    splits = pd.read_parquet(DATA / "splits" / "fold_assignment_FINAL.parquet")
    lib["OLIGO ID"] = lib["OLIGO ID"].astype(str)
    splits["OLIGO ID"] = splits["OLIGO ID"].astype(str)
    merged = lib.merge(splits, on="OLIGO ID", how="inner")

    held_out = merged[merged.fold == WITHHELD_FOLD].copy()
    label_cols = [c for c in held_out.columns if any(c.startswith(p) for p in ["BS_", "EC_", "PA_"])]
    input_cols = [c for c in held_out.columns if c not in label_cols and c != "fold"]

    test_inputs = held_out[input_cols].reset_index(drop=True)
    test_labels = held_out[["OLIGO ID"] + label_cols].reset_index(drop=True)

    test_inputs.to_parquet(HELD_OUT_DIR / "test_inputs.parquet", index=False)
    test_labels.to_parquet(OUT / "gate8_held_out_eval_private_labels.parquet", index=False)

    print(f"Withheld fold: {WITHHELD_FOLD}")
    print(f"Test inputs (PUBLIC): {len(test_inputs)} rows, {len(input_cols)} cols -> {HELD_OUT_DIR / 'test_inputs.parquet'}")
    print(f"Test labels (PRIVATE, do not distribute): {len(test_labels)} rows, {len(label_cols)} cols -> "
          f"{OUT / 'gate8_held_out_eval_private_labels.parquet'}")
    print(f"\nNOTE: data/core/three_host.parquet (public dev data) still contains fold-4 rows WITH labels, "
          f"since it is the existing full-library table used elsewhere in this project. The held-out-eval "
          f"protocol's actual guarantee is that held_out_eval/test_inputs.parquet (the file a submitter "
          f"would use) carries no labels -- submitters are expected to use the held-out protocol's files, "
          f"not re-derive fold 4 from the full public library. This is disclosed as a real, non-cryptographic "
          f"limitation of a solo-maintainer benchmark with no hosted server: a sufficiently motivated user "
          f"COULD look up fold 4's labels in the public dev table. See out/KNOWN_ISSUES.md.")


if __name__ == "__main__":
    main()
