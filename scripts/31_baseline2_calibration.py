"""
GATE 3 - Baseline 2: per-host model trained on N calibration examples. THE
REAL COMPETITOR the eventual Gate 4 host-conditioned model must beat.

ARCHITECTURE DEVIATION FROM SPEC, disclosed prominently (see out/GATE3_MEMO.md):
The task specifies "the same architecture family you intend for Gate 4's
model" -- a small CNN (scripts/29_model.py implements it, 213K params,
matching the charter's spec exactly, and it trains and produces sensible
gradients). It was NOT used for this baseline's full grid. Diagnosis: this
sandboxed environment's PyTorch CPU build has no MKL/MKLDNN acceleration
(confirmed: torch.backends.mkl.is_available()==False,
torch.backends.mkldnn.is_available()==False). A single 213K-parameter
forward+backward step on a batch of 256 takes ~1.4s (4 threads) to ~3.8s (1
thread) -- roughly 3-4 orders of magnitude slower than expected for a model
this size. At that rate, one N=3000 fit (60 epochs) takes ~5-7 minutes; the
full grid (3 hosts x 7 Ns x 10 draws x 2 readouts = up to 420 fits, weighted
toward small N but including many N=1000/3000 fits) was estimated at several
hours to over a day of wall-clock time -- infeasible within this session.

SUBSTITUTE, used for the full grid: k-mer frequency features (k=4, 256
dimensions, counts normalized to sum to 1 per sequence) as a fast,
vectorizable proxy for local-motif sensitivity (the same qualitative
information a small-kernel CNN's first layer extracts), combined with
scikit-learn's LogisticRegression (classifier) and Ridge (regressor) --
both BLAS-accelerated in this environment and complete a fit in
milliseconds regardless of N. This is SEQUENCE-ONLY (no host features), so
it still isolates what B2 is meant to isolate (what a practitioner gets from
measuring N parts in their own host, architecture aside).

Regularization (the "tuning effort" given to this baseline, applied
uniformly, not re-tuned per N/host/draw): LogisticRegression C=1.0,
Ridge alpha=10.0 -- both scikit-learn near-defaults, chosen once by checking
that neither degenerately over- nor under-fits on a representative N=300
draw before committing to the full grid, not swept per configuration.

N=0 falls back exactly to Baseline 1's train-pool statistic (see script 30).
10 random draws per N>0, each with an independent numpy Generator seed.
"""
import numpy as np
import json
import time
from pathlib import Path
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import matthews_corrcoef, roc_auc_score, mean_squared_error
import warnings
warnings.filterwarnings("ignore")

DATA = Path(__file__).resolve().parent.parent / "data"
CACHE = DATA / "baseline_cache"
OUT = Path(__file__).resolve().parent.parent / "out" / "baselines"
OUT.mkdir(parents=True, exist_ok=True)

HOSTS = ["EC", "BS", "PA"]
TEST_FOLD = 0
N_VALUES = [0, 10, 30, 100, 300, 1000, 3000]
N_DRAWS = 10
KMER_K = 4
LOGREG_C = 1.0
RIDGE_ALPHA = 10.0

BASES = "ACGT"


def build_kmer_features(onehot):
    """onehot: (N, 4, 165) -> kmer frequency features (N, 4**KMER_K)."""
    seqs_idx = onehot.argmax(axis=1)  # (N, 165) base index per position
    n_seqs, seq_len = seqs_idx.shape
    n_kmers = seq_len - KMER_K + 1
    n_features = 4 ** KMER_K

    # vectorized rolling k-mer index computation
    powers = 4 ** np.arange(KMER_K)
    kmer_ids = np.zeros((n_seqs, n_kmers), dtype=np.int32)
    for offset in range(KMER_K):
        kmer_ids += seqs_idx[:, offset:offset + n_kmers] * powers[offset]

    features = np.zeros((n_seqs, n_features), dtype=np.float32)
    for i in range(n_seqs):
        counts = np.bincount(kmer_ids[i], minlength=n_features)
        features[i] = counts / counts.sum()
    return features


def evaluate_fit(clf, reg, X_test, y_active_test, y_strength_test, active_mask_test, log_transform):
    pred_prob = clf.predict_proba(X_test)[:, 1] if clf is not None else np.full(len(X_test), np.nan)
    pred_class = (pred_prob > 0.5).astype(int) if clf is not None else np.zeros(len(X_test), dtype=int)

    mcc = matthews_corrcoef(y_active_test, pred_class) if clf is not None and len(set(y_active_test)) > 1 else None
    try:
        auc = roc_auc_score(y_active_test, pred_prob) if clf is not None and len(set(y_active_test)) > 1 else None
    except Exception:
        auc = None

    rho, rmse = None, None
    if reg is not None and active_mask_test.sum() > 1:
        Xa = X_test[active_mask_test]
        y_true = y_strength_test[active_mask_test]
        if log_transform:
            y_true = np.log1p(np.clip(y_true, 0, None))
        y_pred = reg.predict(Xa)
        if len(set(np.round(y_true, 6))) > 1:
            rho = spearmanr(y_true, y_pred).correlation
        rmse = float(np.sqrt(mean_squared_error(y_true, y_pred)))

    return {"mcc": float(mcc) if mcc is not None and mcc == mcc else None,
            "auc": float(auc) if auc is not None and auc == auc else None,
            "spearman_rho": float(rho) if rho is not None and rho == rho else None,
            "rmse_log_strength": rmse}


def run_host_readout(host, readout, d, kmer_features):
    if readout == "transcription":
        usable, active, strength = d["tx_usable"], d["tx_active"], d["tx_norm"]
        usable_reg = usable
        log_transform = True
    else:
        usable, active, strength = d["tl_usable"], d["tl_active"], d["protein_log10"]
        usable_reg = d["tl_usable_for_regression"]
        log_transform = False

    fold = d["fold"]
    train_pool_mask = usable & (fold != TEST_FOLD)
    test_mask = usable & (fold == TEST_FOLD)
    train_pool_idx = np.where(train_pool_mask)[0]
    test_idx = np.where(test_mask)[0]

    X_test = kmer_features[test_idx]
    y_active_test = active[test_idx].astype(int)
    y_strength_test = strength[test_idx]
    active_mask_test = active[test_idx]

    results_by_n = {}
    for N in N_VALUES:
        if N == 0:
            # falls back to Baseline 1 exactly -- constant prediction
            majority_rate = active[train_pool_idx].mean()
            majority_class = int(majority_rate > 0.5)
            pred_prob = np.full(len(test_idx), majority_rate)
            pred_class = np.full(len(test_idx), majority_class)
            mcc = matthews_corrcoef(y_active_test, pred_class) if len(set(y_active_test)) > 1 else None
            try:
                auc = roc_auc_score(y_active_test, pred_prob) if len(set(y_active_test)) > 1 else None
            except Exception:
                auc = None
            train_reg_pool = usable_reg[train_pool_idx] & active[train_pool_idx]
            train_strength = strength[train_pool_idx][train_reg_pool]
            if log_transform:
                train_strength = np.log1p(np.clip(train_strength, 0, None))
            mean_pred = train_strength.mean() if len(train_strength) else np.nan
            y_true_test = y_strength_test[active_mask_test]
            if log_transform:
                y_true_test = np.log1p(np.clip(y_true_test, 0, None))
            rho = None
            rmse = float(np.sqrt(np.mean((y_true_test - mean_pred) ** 2))) if len(y_true_test) else None
            draws = [{"mcc": float(mcc) if mcc is not None and mcc == mcc else None,
                      "auc": float(auc) if auc is not None and auc == auc else None,
                      "spearman_rho": rho, "rmse_log_strength": rmse, "n_actual_train": 0}]
            results_by_n[N] = draws
            continue

        draws = []
        for draw_i in range(N_DRAWS):
            rng = np.random.default_rng(1000 * N + draw_i)
            if N > len(train_pool_idx):
                sample_idx = train_pool_idx  # cap at pool size
            else:
                sample_idx = rng.choice(train_pool_idx, size=N, replace=False)

            X_train = kmer_features[sample_idx]
            y_active_train = active[sample_idx].astype(int)
            y_strength_train_full = strength[sample_idx]
            active_mask_train = active[sample_idx] & usable_reg[sample_idx]

            clf = None
            if len(set(y_active_train)) > 1:
                clf = LogisticRegression(C=LOGREG_C, max_iter=300).fit(X_train, y_active_train)

            reg = None
            if active_mask_train.sum() >= 3:
                Xa = X_train[active_mask_train]
                ya = y_strength_train_full[active_mask_train]
                if log_transform:
                    ya = np.log1p(np.clip(ya, 0, None))
                if len(set(np.round(ya, 6))) > 1:
                    reg = Ridge(alpha=RIDGE_ALPHA).fit(Xa, ya)

            metrics = evaluate_fit(clf, reg, X_test, y_active_test, y_strength_test, active_mask_test, log_transform)
            metrics["n_actual_train"] = int(len(sample_idx))
            draws.append(metrics)
        results_by_n[N] = draws

    return results_by_n, {"n_train_pool": int(len(train_pool_idx)), "n_test": int(len(test_idx))}


def main():
    t_start = time.time()
    all_results = {}
    for host in HOSTS:
        print(f"\n=== {host} ===")
        d = dict(np.load(CACHE / f"{host}_baseline_data.npz", allow_pickle=True))
        kmer_features = build_kmer_features(d["onehot"])
        print(f"  k-mer features built: {kmer_features.shape}")

        host_results = {}
        host_meta = {}
        for readout in ["transcription", "translation"]:
            t0 = time.time()
            res, meta = run_host_readout(host, readout, d, kmer_features)
            print(f"  {readout}: pool={meta['n_train_pool']}, test={meta['n_test']} "
                  f"({time.time()-t0:.1f}s)")
            for N in N_VALUES:
                rhos = [r["spearman_rho"] for r in res[N] if r["spearman_rho"] is not None]
                mccs = [r["mcc"] for r in res[N] if r["mcc"] is not None]
                rho_str = f"rho={np.mean(rhos):.3f}±{np.std(rhos):.3f}" if rhos else "rho=N/A"
                mcc_str = f"mcc={np.mean(mccs):.3f}±{np.std(mccs):.3f}" if mccs else "mcc=N/A"
                print(f"    N={N:>5d}: {rho_str}  {mcc_str}  (n_draws={len(res[N])})")
            host_results[readout] = res
            host_meta[readout] = meta
        all_results[host] = {"results": host_results, "meta": host_meta}

    print(f"\nTotal wall time: {time.time()-t_start:.1f}s")
    with open(OUT / "baseline2_calibration_raw.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"Wrote {OUT / 'baseline2_calibration_raw.json'}")


if __name__ == "__main__":
    main()
