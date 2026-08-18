"""
GATE 7 - Task A: verify the predict_seqonly/predict_film batching fix is
numerically identical to the original unbatched path, and benchmark the
speedup.

NOTE on the benchmark baseline: the UNBATCHED path's cost at 23,232 rows
(8,454s / 2.35hr) was already directly measured during Gate 7's diagnosis
(a real run, not an estimate -- see out/GATE7_MEMO.md). Re-running the
unbatched path at that size here to "confirm" it again would cost another
2+ hours for no new information, so this script reuses that measured value
as the "before" baseline and only benchmarks the NEW chunked path fresh
(fast, since the fix's whole point is that it no longer falls off the
cliff). The 5,810-row unbatched baseline (28.5s) is also reused from the
same diagnostic run.
"""
import sys
import time
import json
import numpy as np
import torch
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
so = import_module("58_sequence_only_model")
m40 = import_module("40_film_cnn_model")

OUT = Path(__file__).resolve().parent.parent / "out"
RESULTS = OUT / "results"
MODELS_DIR = OUT / "models"
DATA = Path(__file__).resolve().parent.parent / "data"

MEASURED_UNBATCHED = {  # from the Gate 7 diagnostic run, out/GATE7_MEMO.md
    "seqonly_5810_rows_sec": 28.474313020706177,
    "seqonly_23232_rows_sec": 8454.029106855392,
}


def unbatched_seqonly_reference(model, seq):
    """Exact copy of the ORIGINAL (pre-fix) predict_seqonly body, kept here
    only as a reference implementation for the identity check -- not used
    anywhere else, and not reintroducing the bug into the shared module."""
    model.eval()
    dev = next(model.parameters()).device
    with torch.no_grad():
        X = torch.tensor(seq, dtype=torch.float32).to(dev)
        preds = model(X)
        return {k: v.cpu().numpy() for k, v in preds.items()}


def unbatched_film_reference(model, seq, host_vec):
    model.eval()
    dev = next(model.parameters()).device
    with torch.no_grad():
        X = torch.tensor(seq, dtype=torch.float32).to(dev)
        H = torch.tensor(host_vec, dtype=torch.float32).to(dev)
        preds = model(X, H)
        return {k: v.cpu().numpy() for k, v in preds.items()}


def verify_identity():
    print("=" * 78)
    print("STEP 1: numerical identity check, chunked vs unbatched, n=2000 (below the cliff)")
    print("=" * 78)
    results = {}

    d = dict(np.load(DATA / "baseline_cache" / "EC_baseline_data.npz", allow_pickle=True))
    seq_small = d["onehot"][:2000]

    model_so = so.SequenceOnlyCNN()
    model_so.load_state_dict(torch.load(MODELS_DIR / "seqonly_loho_EC_fold0.pt", map_location="cpu"))
    ref = unbatched_seqonly_reference(model_so, seq_small)
    chunked = so.predict_seqonly(model_so, seq_small, batch_size=1024)
    max_diffs = {k: float(np.max(np.abs(ref[k] - chunked[k]))) for k in ref}
    print(f"  predict_seqonly max abs diff per head: {max_diffs}")
    results["predict_seqonly"] = max_diffs

    g_z, gcols, p_z, pcols, imputed = m40.load_host_features()
    seq_te, hv_te, targets_te, masks_te, _ = m40.build_pooled_arrays(
        ["BS"], g_z, fold_filter_fn=lambda f: f == 0)
    seq_small2 = seq_te[:2000]
    hv_small2 = hv_te[:2000]
    model_film = m40.FiLMSequenceCNN(host_dim=len(gcols), use_film=True)
    model_film.load_state_dict(torch.load(MODELS_DIR / "loho_BS_genomic_fold0.pt", map_location="cpu"))
    ref2 = unbatched_film_reference(model_film, seq_small2, hv_small2)
    chunked2 = m40.predict_film(model_film, seq_small2, hv_small2, batch_size=1024)
    max_diffs2 = {k: float(np.max(np.abs(ref2[k] - chunked2[k]))) for k in ref2}
    print(f"  predict_film max abs diff per head: {max_diffs2}")
    results["predict_film"] = max_diffs2

    overall_max = max(max(v.values()) for v in results.values())
    results["overall_max_abs_diff"] = overall_max
    results["verdict"] = "NUMERICALLY IDENTICAL (float32 precision)" if overall_max < 1e-4 else "DIFFERS -- INVESTIGATE BEFORE TRUSTING"
    print(f"\n  Overall max abs diff: {overall_max:.2e}")
    print(f"  VERDICT: {results['verdict']}")
    return results


def benchmark_chunked():
    print("\n" + "=" * 78)
    print("STEP 2: benchmark the CHUNKED path (unbatched baseline reused from Gate 7 diagnosis)")
    print("=" * 78)
    d = dict(np.load(DATA / "baseline_cache" / "EC_baseline_data.npz", allow_pickle=True))
    model = so.SequenceOnlyCNN()
    model.load_state_dict(torch.load(MODELS_DIR / "seqonly_loho_EC_fold0.pt", map_location="cpu"))

    results = {}
    for n_rows, label in [(5810, "5810_rows"), (23232, "23232_rows")]:
        seq = d["onehot"][:n_rows]
        t0 = time.time()
        _ = so.predict_seqonly(model, seq, batch_size=1024)
        elapsed = time.time() - t0
        baseline_key = f"seqonly_{n_rows}_rows_sec"
        baseline = MEASURED_UNBATCHED.get(baseline_key)
        speedup = baseline / elapsed if baseline else None
        results[label] = {"n_rows": n_rows, "chunked_sec": elapsed, "unbatched_sec_measured_earlier": baseline,
                           "speedup": speedup}
        print(f"  n={n_rows}: chunked={elapsed:.2f}s, unbatched(measured earlier)={baseline}, speedup={speedup}")

    ratio_check = results["23232_rows"]["chunked_sec"] / results["5810_rows"]["chunked_sec"]
    row_ratio = 23232 / 5810
    print(f"\n  Chunked scaling check: 4x rows -> {ratio_check:.2f}x time (linear would be {row_ratio:.2f}x)")
    results["scaling_check"] = {"row_ratio": row_ratio, "time_ratio": ratio_check,
                                 "roughly_linear": abs(ratio_check - row_ratio) < row_ratio}

    if torch.backends.mps.is_available():
        print("\n  MPS benchmark:")
        model_mps = so.SequenceOnlyCNN().to("mps")
        model_mps.load_state_dict(torch.load(MODELS_DIR / "seqonly_loho_EC_fold0.pt", map_location="mps"))
        for n_rows, label in [(5810, "5810_rows"), (23232, "23232_rows")]:
            seq = d["onehot"][:n_rows]
            t0 = time.time()
            _ = so.predict_seqonly(model_mps, seq, batch_size=1024)
            elapsed = time.time() - t0
            results[f"mps_{label}"] = {"n_rows": n_rows, "chunked_sec": elapsed}
            print(f"    n={n_rows}: chunked={elapsed:.2f}s")
    else:
        print("\n  MPS not available in this environment -- CPU-only benchmark reported.")
        results["mps_available"] = False

    return results


def main():
    identity = verify_identity()
    if identity["overall_max_abs_diff"] >= 1e-4:
        print("\nABORTING benchmark -- identity check failed, fix is not trustworthy as-is.")
        with open(RESULTS / "gate7_batching_fix_verification.json", "w") as f:
            json.dump({"identity": identity}, f, indent=2)
        return
    bench = benchmark_chunked()
    with open(RESULTS / "gate7_batching_fix_verification.json", "w") as f:
        json.dump({"identity": identity, "benchmark": bench}, f, indent=2, default=str)
    print(f"\nWrote {RESULTS / 'gate7_batching_fix_verification.json'}")


if __name__ == "__main__":
    main()
