"""
GATE 1.5 - Task B1/B2: Download PaxDb v6.0 protein abundance data for 5 of the 6
crosshost benchmark hosts (V. natriegens is absent from PaxDb -- see script 08).

API discovered by reverse-engineering the PaxDb website's Vue.js bundle (app.js),
which references https://api.pax-db.org. The OpenAPI spec is served at
https://api.pax-db.org/try-api/openapi.json. No registered API key was obtained;
all requests use the documented public fallback `x-api-key: test` scope (shared,
rate-limited pool, explicitly sanctioned for testing by the API's own error
message). Base path is /v6.

Datasets selected (see raw/paxdb_datasets_*.json for the full per-species listing
this was chosen from):
  - PRIMARY pick per host: the highest-scoring "Integrated" (weighted-average-
    across-all-studies) dataset PaxDb provides for that species -- this is PaxDb's
    own best single estimate of "typical" abundance.
  - CROSS-CHECK pick for B. subtilis, P. aeruginosa, S. enterica: a single-study
    2025 dataset (Abele et al., Mol Cell Proteomics, MassIVE MSV000096603 -- a
    303-species DSMZ-strain-collection bacterial proteome atlas) that happens to
    cover exactly these 3 of our 6 hosts under one consistent lab/method (though
    NOT one consistent growth medium -- see memo comparability section). Used to
    sanity-check whether the derived metrics are sensitive to source choice.
    E. coli and C. glutamicum were checked and are NOT in this atlas under PaxDb's
    current per-species dataset listing (confirmed by grep, see terminal history).
"""
import requests
import json
import time
from pathlib import Path

RAW = Path(__file__).resolve().parent.parent / "raw"
API = "https://api.pax-db.org/v6"
HEADERS = {"Accept": "application/json", "User-Agent": "Mozilla/5.0"}

DATASETS = {
    # host_tag: [(dataset_id, species_id, label, num_abundances), ...]
    "EC": [(3740039012, 511145, "integrated", 3747)],
    "BS": [(4187424622, 224308, "integrated", 4052), (945251378, 224308, "abele2025", 3278)],
    "PA": [(1789377356, 208964, "integrated", 5034), (938070019, 208964, "abele2025", 4509)],
    "SE": [(1793422220, 99287, "integrated", 2620), (1071551919, 99287, "abele2025", 2492)],
    "CG": [(3268883846, 196627, "integrated", 1227)],
}


def fetch_dataset_abundances(dataset_id, num_abundances, page_size=1000):
    all_rows = []
    start = 0
    while start < num_abundances:
        end = min(start + page_size, num_abundances)
        url = f"{API}/abundances/dataset/{dataset_id}?start={start}&end={end}&sort=-abundance"
        resp = requests.get(url, headers=HEADERS, timeout=30)
        resp.raise_for_status()
        payload = resp.json()
        rows = payload.get("data", [])
        if not rows:
            break
        all_rows.extend(rows)
        start = end
        time.sleep(0.15)  # be polite to the shared test-key rate-limited pool
    return all_rows


def main():
    for host, dsets in DATASETS.items():
        for dataset_id, species_id, label, n in dsets:
            out_path = RAW / f"paxdb_abundances_{host}_{label}_{dataset_id}.json"
            if out_path.exists():
                print(f"SKIP (already downloaded): {out_path.name}")
                continue
            print(f"Downloading {host} / {label} / dataset {dataset_id} (expect ~{n} proteins)...")
            rows = fetch_dataset_abundances(dataset_id, n)
            print(f"  got {len(rows)} rows")
            with open(out_path, "w") as f:
                json.dump({"host": host, "dataset_id": dataset_id, "species_id": species_id,
                           "label": label, "n_expected": n, "data": rows}, f)
            print(f"  wrote {out_path}")


if __name__ == "__main__":
    main()
