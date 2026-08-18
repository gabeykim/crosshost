"""
GATE 3.5 - Task 3: measure the ACTUAL cost of the specified CNN architecture
(scripts/29_model.py, 213,058 params) in this environment, at Gate-4-scale
data, and extrapolate to Gate 4 (12-24 fits) and Gate 7 (60-120 ensemble
fits).

MAJOR FINDING this task exists to report: this machine has Apple Silicon
(torch.backends.mps.is_available() == True), which Gate 3 never used or
checked -- Gate 3's "no BLAS/MKL" finding was CPU-specific and did not rule
out MPS. This changes the compute-feasibility picture substantially.

Benchmarks THREE things:
  1. Raw forward+backward step cost, CPU vs MPS, at realistic Gate-4-scale N
     (~45,000 -- the actual size of a two-primary-host pooled training pool,
     not the N<=3000 calibration-curve values Gate 3 benchmarked).
  2. A REAL end-to-end train_model() call (with early stopping, on real BS
     transcription data) at both a small and the full-scale N, to confirm
     the raw-step extrapolation against actual convergence behavior.
  3. Network reachability (for the cloud-GPU feasibility question) -- checked,
     not assumed.
"""
import sys
import time
import subprocess
from pathlib import Path
import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
model_mod = import_module("29_model")

DATA = Path(__file__).resolve().parent.parent / "data"
CACHE = DATA / "baseline_cache"
OUT = Path(__file__).resolve().parent.parent / "out"


def check_devices():
    info = {
        "cuda_available": torch.cuda.is_available(),
        "mps_available": torch.backends.mps.is_available(),
        "mkl_available": torch.backends.mkl.is_available(),
        "mkldnn_available": torch.backends.mkldnn.is_available(),
        "torch_version": torch.__version__,
    }
    return info


def check_network():
    results = {}
    for name, url in [("generic_https", "https://aws.amazon.com"), ("dns_resolution_only", None)]:
        if url is None:
            continue
        try:
            out = subprocess.run(["curl", "-s", "-m", "5", "-o", "/dev/null", "-w", "%{http_code}", url],
                                  capture_output=True, text=True, timeout=10)
            results[name] = {"url": url, "http_code": out.stdout.strip()}
        except Exception as e:
            results[name] = {"url": url, "error": str(e)}
    return results


def raw_step_benchmark(N, batch_sizes=(256, 1024)):
    """Synthetic data, isolates pure forward+backward+optimizer-step cost,
    no data loading / early-stopping overhead -- a clean per-epoch estimate."""
    X = np.random.rand(N, 4, 165).astype(np.float32)
    y_active = np.random.randint(0, 2, N).astype(np.float32)
    y_strength = np.random.rand(N).astype(np.float32)
    active_mask = np.random.randint(0, 2, N).astype(bool)

    results = {}
    devices = ["cpu"] + (["mps"] if torch.backends.mps.is_available() else [])
    for device_name in devices:
        for batch_size in batch_sizes:
            device = torch.device(device_name)
            torch.manual_seed(0)
            model = model_mod.SequenceCNN().to(device)
            opt = torch.optim.Adam(model.parameters(), lr=1e-3)
            bce = torch.nn.BCEWithLogitsLoss()
            mse = torch.nn.MSELoss()
            Xt = torch.tensor(X, dtype=torch.float32).to(device)
            ya = torch.tensor(y_active).to(device)
            ys = torch.tensor(y_strength).to(device)
            am = torch.tensor(active_mask).to(device)
            n = len(Xt)
            perm = torch.randperm(n)
            for i in range(0, min(batch_size * 3, n), batch_size):
                idx = perm[i:i + batch_size]
                opt.zero_grad()
                al, st = model(Xt[idx], None)
                loss = bce(al, ya[idx])
                loss.backward()
                opt.step()
            if device_name == "mps":
                torch.mps.synchronize()
            t0 = time.time()
            n_batches = 0
            for i in range(0, n, batch_size):
                idx = perm[i:i + batch_size]
                opt.zero_grad()
                al, st = model(Xt[idx], None)
                loss = bce(al, ya[idx])
                am_b = am[idx]
                if am_b.sum() > 0:
                    loss = loss + mse(st[am_b], ys[idx][am_b])
                loss.backward()
                opt.step()
                n_batches += 1
            if device_name == "mps":
                torch.mps.synchronize()
            t1 = time.time()
            epoch_time = t1 - t0
            key = f"{device_name}_batch{batch_size}"
            results[key] = {"device": device_name, "batch_size": batch_size, "N": N,
                             "n_batches": n_batches, "epoch_time_sec": epoch_time,
                             "sixty_epoch_fit_min": epoch_time * 60 / 60}
            print(f"  {key}: N={N}, {n_batches} batches, epoch_time={epoch_time:.2f}s, "
                  f"60-epoch fit = {epoch_time:.2f} min")
    return results


def real_end_to_end_benchmark():
    """One REAL train_model() call with early stopping, on real B. subtilis
    transcription data (BS chosen since it's the smallest primary host --
    conservative/fast sanity check), CPU vs MPS, to confirm the raw-step
    extrapolation against actual convergence dynamics (not just synthetic
    per-epoch cost)."""
    d = dict(np.load(CACHE / "BS_baseline_data.npz", allow_pickle=True))
    usable = d["tx_usable"]
    idx = np.where(usable)[0][:2000]  # a moderate, fast-to-fit real slice
    X = d["onehot"][idx]
    y_active = d["tx_active"][idx].astype(np.float32)
    y_strength = d["tx_norm"][idx].astype(np.float32)
    active_mask = d["tx_active"][idx]

    results = {}
    devices = ["cpu"] + (["mps"] if torch.backends.mps.is_available() else [])
    for device_name in devices:
        t0 = time.time()
        model = model_mod.train_model(X, y_active, y_strength, active_mask,
                                       epochs=60, batch_size=64, seed=0, device=device_name)
        t1 = time.time()
        prob, strength = model_mod.predict(model, X)
        results[device_name] = {"N": len(idx), "wall_clock_sec": t1 - t0}
        print(f"  REAL fit, device={device_name}, N={len(idx)}: {t1 - t0:.1f}s wall-clock "
              f"(early-stopped at or before 60 epochs)")
    return results


def extrapolate(per_fit_minutes, label):
    for n_fits in [12, 24, 60, 120]:
        total_min = per_fit_minutes * n_fits
        print(f"    {label}: {n_fits} fits -> {total_min:.0f} min = {total_min/60:.1f} hr = {total_min/60/24:.2f} days")


def main():
    print("=== Device / acceleration availability ===")
    dev_info = check_devices()
    for k, v in dev_info.items():
        print(f"  {k}: {v}")

    print("\n=== Network reachability (for cloud-GPU feasibility question) ===")
    net_info = check_network()
    for k, v in net_info.items():
        print(f"  {k}: {v}")

    GATE4_REALISTIC_N = 45000  # ~ two pooled primary hosts' full training data (LOHO)
    print(f"\n=== Raw per-epoch step cost at Gate-4-realistic scale (N={GATE4_REALISTIC_N}) ===")
    raw_results = raw_step_benchmark(GATE4_REALISTIC_N)

    print("\n=== Real end-to-end train_model() call (real data, early stopping) ===")
    real_results = real_end_to_end_benchmark()

    print("\n=== Extrapolation to Gate 4 / Gate 7 fit counts ===")
    cpu_key = "cpu_batch256"
    mps_key = "mps_batch1024" if "mps_batch1024" in raw_results else None
    print("  CPU (batch=256):")
    extrapolate(raw_results[cpu_key]["sixty_epoch_fit_min"], "CPU")
    if mps_key:
        print("  MPS (batch=1024, best observed):")
        extrapolate(raw_results[mps_key]["sixty_epoch_fit_min"], "MPS")

    output = {
        "device_info": dev_info,
        "network_check": net_info,
        "raw_step_benchmark_N": GATE4_REALISTIC_N,
        "raw_step_benchmark": raw_results,
        "real_end_to_end_benchmark_N": 2000,
        "real_end_to_end_benchmark": real_results,
        "gate4_realistic_pooled_n_per_loho_fit": {
            "held_out_BS_train_EC_PA_tx": 24412 + 21306,
            "held_out_EC_train_BS_PA_tx": 15697 + 21306,
            "held_out_PA_train_EC_BS_tx": 24412 + 15697,
        },
    }
    import json
    with open(OUT / "gate3_5_cnn_compute_benchmark.json", "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nWrote {OUT / 'gate3_5_cnn_compute_benchmark.json'}")


if __name__ == "__main__":
    main()
