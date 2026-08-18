"""
GATE 2 - Task 2: Download RefSeq reference genomes for all 6 hosts.

Assembly accessions identified via NCBI Datasets REST API v2
(api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{taxid}/dataset_report),
restricted to reference/representative complete genomes. All 6 are RefSeq
"Complete Genome" assemblies -- exact accessions below, [VERIFIED] via the API
response saved to raw/genomes/*_dataset_report.json.
"""
import requests
import json
import subprocess
from pathlib import Path

RAW = Path(__file__).resolve().parent.parent / "raw" / "genomes"
RAW.mkdir(parents=True, exist_ok=True)

HOSTS = {
    "EC": {"accession": "GCF_000005845.2", "taxid": 511145, "name": "Escherichia coli K-12 MG1655"},
    "BS": {"accession": "GCF_000009045.1", "taxid": 224308, "name": "Bacillus subtilis 168"},
    "PA": {"accession": "GCF_000006765.1", "taxid": 208964, "name": "Pseudomonas aeruginosa PAO1"},
    "SE": {"accession": "GCF_000006945.2", "taxid": 99287, "name": "Salmonella enterica Typhimurium LT2"},
    "CG": {"accession": "GCF_000011325.1", "taxid": 196627, "name": "Corynebacterium glutamicum ATCC 13032"},
    "VN": {"accession": "GCF_001456255.1", "taxid": 1219067, "name": "Vibrio natriegens ATCC 14048"},
}


def ftp_path_for_accession(acc):
    # RefSeq FTP layout: GCF/000/005/845/GCF_000005845.2_ASM584v2/
    prefix = acc.split("_")[1].split(".")[0]  # 000005845
    triplets = [prefix[0:3], prefix[3:6], prefix[6:9]]
    return f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/{triplets[0]}/{triplets[1]}/{triplets[2]}"


def find_full_dirname(base_url, acc):
    # directory listing to find the full assembly-name-suffixed folder
    resp = requests.get(base_url + "/", timeout=30)
    resp.raise_for_status()
    import re
    matches = re.findall(rf'href="({re.escape(acc)}[^"/]*)/"', resp.text)
    if not matches:
        raise RuntimeError(f"Could not find directory for {acc} at {base_url}")
    return matches[0]


def download_file(url, out_path):
    if out_path.exists() and out_path.stat().st_size > 0:
        print(f"  SKIP (exists): {out_path.name}")
        return
    resp = requests.get(url, timeout=60, stream=True)
    resp.raise_for_status()
    with open(out_path, "wb") as f:
        for chunk in resp.iter_content(chunk_size=1 << 16):
            f.write(chunk)
    print(f"  downloaded {out_path.name} ({out_path.stat().st_size/1e6:.2f} MB)")


def main():
    manifest = {}
    for host, info in HOSTS.items():
        acc = info["accession"]
        print(f"\n=== {host}: {info['name']} ({acc}) ===")
        base = ftp_path_for_accession(acc)
        dirname = find_full_dirname(base, acc)
        full_url = f"{base}/{dirname}"
        print(f"  directory: {full_url}")

        host_dir = RAW / host
        host_dir.mkdir(exist_ok=True)

        files = {
            "genomic.fna.gz": f"{dirname}_genomic.fna.gz",
            "protein.faa.gz": f"{dirname}_protein.faa.gz",
            "genomic.gff.gz": f"{dirname}_genomic.gff.gz",
            "cds_from_genomic.fna.gz": f"{dirname}_cds_from_genomic.fna.gz",
            "rna_from_genomic.fna.gz": f"{dirname}_rna_from_genomic.fna.gz",
        }
        for label, fname in files.items():
            url = f"{full_url}/{fname}"
            out_path = host_dir / label
            try:
                download_file(url, out_path)
            except Exception as e:
                print(f"  FAILED {label}: {e}")

        manifest[host] = {**info, "directory_name": dirname, "source_url": full_url}

    with open(RAW / "genome_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"\nWrote {RAW / 'genome_manifest.json'}")


if __name__ == "__main__":
    main()
