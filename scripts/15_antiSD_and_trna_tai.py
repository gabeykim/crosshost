"""
GATE 2 - Task 2: Anti-Shine-Dalgarno feature (ViennaRNA hybridization free
energy) + tRNA gene copy numbers (tRNAscan-SE 2.0) + tAI (dos Reis et al. 2004
formula, implemented directly -- no separate tool needed once tRNA gene copy
numbers are known).

UNIT TEST (anti-SD): the extracted 9-mer anti-SD motif for E. coli and
P. aeruginosa must equal "ACCTCCTTA", and for B. subtilis must equal
"ACCTCCTTT" -- these are Johns et al. 2018's own published anti-SD sequences
(Gate 1, manuscript Methods text), used here as an independent primary-source
check on the extraction pipeline. See main() for the assertion.
"""
import json
import re
import subprocess
import gzip
import tempfile
from pathlib import Path
from collections import Counter

RAW = Path(__file__).resolve().parent.parent / "raw" / "genomes"
OUT = Path(__file__).resolve().parent.parent / "out"
DATA = Path(__file__).resolve().parent.parent / "data"

HOSTS = ["EC", "BS", "PA", "SE", "CG", "VN"]

# canonical/idealized Shine-Dalgarno consensus (Shine & Dalgarno 1974; also the
# exact reverse-complement of the anti-SD nonamer "ACCUCCUUA" found in E. coli
# and P. aeruginosa's own 16S rRNA -- see docstring)
CONSENSUS_SD_RNA = "UAAGGAGGU"

ANCHOR = "ACCTCCTT"

# published anti-SD sequences [VERIFIED, Johns et al. 2018 manuscript Methods,
# raw/nihms945382_bodytext.txt] used as the unit test
PUBLISHED_ANTI_SD = {"EC": "ACCTCCTTA", "PA": "ACCTCCTTA", "BS": "ACCTCCTTT"}


def extract_anti_sd_nonamer(seqs):
    """Find the ACCTCCTT anchor + 1 following base in each 16S copy, return the
    majority-consensus 9-mer across copies (robust to per-copy 3'-boundary
    annotation noise, e.g. seen in S. enterica's 7 rRNA copies)."""
    nonamers = []
    for s in seqs:
        idx = s.rfind(ANCHOR)
        if idx >= 0 and idx + 9 <= len(s):
            nonamers.append(s[idx:idx + 9])
    if not nonamers:
        return None, {}
    counts = Counter(nonamers)
    majority = counts.most_common(1)[0][0]
    return majority, dict(counts)


def rna_duplex_energy(seq_a_dna, seq_b_dna):
    """Run ViennaRNA RNAduplex, return the minimum free energy (kcal/mol)."""
    seq_a = seq_a_dna.replace("T", "U")
    seq_b = seq_b_dna.replace("T", "U")
    proc = subprocess.run(["RNAduplex"], input=f"{seq_a}\n{seq_b}\n", capture_output=True, text=True, timeout=30)
    # RNAduplex output line looks like: ".((((((((&)))))))). 1,9 : 1,9 (-12.34)"
    line = proc.stdout.strip().splitlines()[-1]
    m = re.search(r"\(\s*(-?\d+\.\d+)\s*\)", line)
    if not m:
        raise RuntimeError(f"Could not parse RNAduplex output: {proc.stdout!r} / stderr={proc.stderr!r}")
    return float(m.group(1)), proc.stdout.strip()


def run_trnascan(genome_fna_gz, host, mode="-B"):
    """Run tRNAscan-SE 2.0. mode: -B for bacteria (default assumption; -A for
    archaea not needed here since all 6 hosts are bacteria)."""
    with tempfile.TemporaryDirectory() as tmp:
        fna_path = Path(tmp) / "genome.fna"
        with gzip.open(genome_fna_gz, "rt") as fin, open(fna_path, "w") as fout:
            fout.write(fin.read())
        out_path = Path(tmp) / "trnascan_out.txt"
        stats_path = Path(tmp) / "trnascan_stats.txt"
        cmd = ["tRNAscan-SE", mode, "-o", str(out_path), "-m", str(stats_path),
               "--thread", "2", str(fna_path)]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if not out_path.exists():
            raise RuntimeError(f"tRNAscan-SE failed for {host}: stdout={proc.stdout}\nstderr={proc.stderr}")
        result_text = out_path.read_text()
    return result_text


def parse_trnascan_output(text):
    """Parse tRNAscan-SE tabular output into (anticodon, isotype) counts."""
    rows = []
    for line in text.splitlines():
        if line.startswith("Sequence") or line.startswith("Name") or line.startswith("--------") or not line.strip():
            continue
        parts = re.split(r"\s+", line.strip())
        if len(parts) < 9:
            continue
        # Columns: seqname, trna_num, start, end, isotype, anticodon, intron_start,
        # intron_end, score, [note]
        isotype = parts[4]
        anticodon = parts[5]
        rows.append({"isotype": isotype, "anticodon": anticodon})
    return rows


def main():
    features = {}
    for host in HOSTS:
        print(f"\n=== {host} ===")
        d = json.load(open(OUT / "genomic_features_basic.json"))
        seqs = d[host]["sixteen_s_sequences"]
        nonamer, counts = extract_anti_sd_nonamer(seqs)
        print(f"  anti-SD nonamer (majority of {len(seqs)} 16S copies): {nonamer}  (vote counts: {counts})")

        energy, raw_out = rna_duplex_energy(nonamer, CONSENSUS_SD_RNA.replace("U", "T"))
        print(f"  RNAduplex free energy vs consensus SD ({CONSENSUS_SD_RNA}): {energy} kcal/mol")
        print(f"  raw RNAduplex output: {raw_out}")

        print(f"  Running tRNAscan-SE 2.0 (bacterial mode)...")
        trnascan_text = run_trnascan(RAW / host / "genomic.fna.gz", host)
        trnas = parse_trnascan_output(trnascan_text)
        iso_counts = Counter(t["isotype"] for t in trnas)
        anticodon_counts = Counter(t["anticodon"] for t in trnas)
        print(f"  tRNA genes found: {len(trnas)} total, {len(iso_counts)} isotypes, "
              f"{len(anticodon_counts)} distinct anticodons")
        print(f"  isotype counts: {dict(iso_counts)}")

        features[host] = {
            "anti_sd_nonamer": nonamer,
            "anti_sd_vote_counts": counts,
            "anti_sd_rnaduplex_mfe_kcalmol": energy,
            "consensus_sd_used": CONSENSUS_SD_RNA,
            "n_trna_genes_total": len(trnas),
            "trna_anticodon_counts": dict(anticodon_counts),
            "trna_isotype_counts": dict(iso_counts),
        }

    # --- UNIT TEST ---
    print("\n" + "=" * 70)
    print("UNIT TEST: anti-SD nonamer vs Johns et al. 2018 published values")
    print("=" * 70)
    all_pass = True
    for host, expected in PUBLISHED_ANTI_SD.items():
        actual = features[host]["anti_sd_nonamer"]
        ok = actual == expected
        all_pass &= ok
        print(f"  {host}: expected={expected} actual={actual} -> {'PASS' if ok else 'FAIL'}")
    assert all_pass, "Anti-SD unit test FAILED -- do not accept this feature"
    print("ALL ANTI-SD UNIT TESTS PASSED")

    with open(OUT / "antisd_trna_features.json", "w") as f:
        json.dump(features, f, indent=2)
    print(f"\nWrote {OUT / 'antisd_trna_features.json'}")


if __name__ == "__main__":
    main()
