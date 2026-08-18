"""
GATE 2 - Task 2: tRNA Adaptation Index (tAI), dos Reis et al. 2003/2004.

Implementation ported directly (line-for-line logic) from the CANONICAL R
source, fetched live from the author's own package repository:
  https://raw.githubusercontent.com/mariodosreis/tai/master/R/tAI.R
  (get.ws function; s-values and codon-grouping formula)
This was necessary because the exact wobble s-value constants are not
something to reconstruct from memory -- they are specific fitted parameters
from dos Reis, Savva & Wernisch (2004) NAR 32:5036-44, and using the wrong
values would silently produce plausible-looking but wrong tAI numbers.

UNIT TEST: this script's output is checked element-by-element against R's own
get.ws() output, computed by actually running the fetched R source in this
same environment (R was installed via conda specifically for this
cross-check) on the package's bundled reference E. coli K-12 tRNA gene vector
(87 tRNA genes, inst/extdata/ecolik12.trna). See main() -- the assertion
requires an exact (float-tolerance) match on all 60 values, not just
"close enough."
"""
import json
import re
import numpy as np
import pandas as pd
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "out"
DATA = Path(__file__).resolve().parent.parent / "data"

# TCAG-nested codon order, EXACTLY matching the R package's ecolik12.trna
# vector ordering (verified against codonR's codonM Perl script and the
# get.ws() index arithmetic: stop codons removed at positions 11,12,15 (1-idx)
# = TAA,TAG,TGA; Met at position 36 (1-idx) = ATG; bacterial ATA-Ile special
# case at position 35 (1-idx) = ATA -- all confirmed to land exactly where
# this ordering places them).
BASES = ["T", "C", "A", "G"]
CODONS_64 = [b1 + b2 + b3 for b1 in BASES for b2 in BASES for b3 in BASES]
assert len(CODONS_64) == 64
assert CODONS_64[10] == "TAA" and CODONS_64[11] == "TAG" and CODONS_64[14] == "TGA"  # 0-indexed = R's 11,12,15
assert CODONS_64[35] == "ATG"  # 0-indexed = R's 36
assert CODONS_64[34] == "ATA"  # 0-indexed = R's 35

# optimised s-values, dos Reis et al. 2004 (fetched from get.ws() in the
# canonical R source -- see docstring)
S_VALUES = [0.0, 0.0, 0.0, 0.0, 0.41, 0.28, 0.9999, 0.68, 0.89]
P_VALUES = [1 - s for s in S_VALUES]

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def reverse_complement(seq):
    return "".join(COMPLEMENT[b] for b in reversed(seq))


def build_trna_vector_from_anticodons(anticodon_counts):
    """dos Reis's tRNA[i] = count of tRNA genes whose anticodon is the exact
    Watson-Crick complement of codon i (i.e. anticodon = reverse_complement(codon)).
    anticodon_counts: dict of {anticodon_seq (DNA, T not U): count} from tRNAscan-SE."""
    vec = np.zeros(64)
    for i, codon in enumerate(CODONS_64):
        wc_anticodon = reverse_complement(codon)
        vec[i] = anticodon_counts.get(wc_anticodon, 0)
    return vec


def get_ws(trna, sking=1, s=None):
    """Direct port of get.ws() from tAI.R. trna: length-64 array in CODONS_64
    order. sking: 1=prokaryote (enables the bacterial ATA-Ile special case)."""
    if s is None:
        s = S_VALUES
    p = [1 - si for si in s]

    W = []
    for i in range(0, 64, 4):  # 0-indexed groups of 4: NNT,NNC,NNA,NNG
        W.append(p[0] * trna[i] + p[4] * trna[i + 1])       # NNT
        W.append(p[1] * trna[i + 1] + p[5] * trna[i])       # NNC
        W.append(p[2] * trna[i + 2] + p[6] * trna[i])       # NNA
        W.append(p[3] * trna[i + 3] + p[7] * trna[i + 2])   # NNG
    W = np.array(W)

    W[35] = p[3] * trna[35]  # index 36 in R (1-idx) = Met, 0-idx 35

    if sking == 1:
        W[34] = p[8]  # index 35 in R (1-idx) = bacterial ATA-Ile, 0-idx 34

    # remove stop codons (0-idx 10,11,14) and Met (0-idx 35)
    keep = [i for i in range(64) if i not in (10, 11, 14, 35)]
    W = W[keep]
    codons_60 = [CODONS_64[i] for i in keep]

    w = W / W.max()
    zero_mask = (w == 0)
    if zero_mask.sum() > 0:
        nonzero = w[~zero_mask]
        gm = np.exp(np.sum(np.log(nonzero)) / len(nonzero))
        w[zero_mask] = gm

    return w, codons_60


def validate_against_r_reference():
    """UNIT TEST: run this Python implementation on the R package's own
    bundled reference E. coli tRNA vector (87 genes) and require an exact
    match against R's own get.ws() output (computed live via Rscript, see
    scripts/16 docstring for the exact command run to produce the reference
    CSV)."""
    ref_trna_path = Path("/tmp/tai_r/ecolik12.trna")
    ref_ws_path = Path("/tmp/tai_r/ecoli_ws_reference.csv")
    if not ref_trna_path.exists() or not ref_ws_path.exists():
        raise RuntimeError(
            "Reference files not found. Run: mamba install -n crosshost -c conda-forge r-base, "
            "then fetch tAI.R and ecolik12.trna from mariodosreis/tai and run get.ws() in R "
            "to produce /tmp/tai_r/ecoli_ws_reference.csv (see terminal history for exact commands)."
        )
    ref_trna = np.array([float(x) for x in ref_trna_path.read_text().split()])
    ref_ws = pd.read_csv(ref_ws_path)["w"].values

    my_ws, my_codons = get_ws(ref_trna, sking=1)

    assert len(my_ws) == len(ref_ws) == 60, f"length mismatch: {len(my_ws)} vs {len(ref_ws)}"
    max_diff = np.max(np.abs(my_ws - ref_ws))
    print(f"Max absolute difference vs R reference output: {max_diff:.2e}")
    assert max_diff < 1e-8, f"Python get_ws() does NOT match R's get.ws() -- max diff {max_diff}"
    print("UNIT TEST PASSED: Python get_ws() exactly matches R get.ws() on reference E. coli data "
          f"(87 tRNA genes, {len(my_ws)} codon weights, max diff {max_diff:.2e})")


def get_tai(codon_freqs, w, codons_60):
    """Geometric mean of w across a sequence's codon usage (dos Reis get.tai())."""
    log_w = np.log(w)
    total = 0.0
    n = 0
    for codon, weight in zip(codons_60, log_w):
        c = codon_freqs.get(codon, 0)
        total += c * weight
        n += c
    return np.exp(total / n) if n > 0 else np.nan


def main():
    print("=" * 70)
    print("UNIT TEST: validating Python tAI implementation against R reference")
    print("=" * 70)
    validate_against_r_reference()

    print("\n" + "=" * 70)
    print("Computing tAI weights and genome-wide tAI per host")
    print("=" * 70)

    antisd = json.load(open(OUT / "antisd_trna_features.json"))
    basic = json.load(open(OUT / "genomic_features_basic.json"))

    results = {}
    for host in ["EC", "BS", "PA", "SE", "CG", "VN"]:
        anticodon_counts_dna = {ac.replace("U", "T"): n
                                 for ac, n in antisd[host]["trna_anticodon_counts"].items()}
        trna_vec = build_trna_vector_from_anticodons(anticodon_counts_dna)
        w, codons_60 = get_ws(trna_vec, sking=1)

        codon_freqs = basic[host]["codon_counts"]
        genome_tai = get_tai(codon_freqs, w, codons_60)

        ranked = sorted(zip(codons_60, w), key=lambda x: -x[1])
        print(f"\n{host}: total tRNA genes (WC-matchable)={int(trna_vec.sum())}, "
              f"genome-wide mean tAI={genome_tai:.4f}")
        print(f"  top 5 codons by tAI weight: {ranked[:5]}")
        print(f"  bottom 5 codons by tAI weight: {ranked[-5:]}")

        results[host] = {
            "trna_vector_64": trna_vec.tolist(),
            "codon_weights_w": dict(zip(codons_60, w.tolist())),
            "genome_wide_tai": float(genome_tai),
            "top5_codons": ranked[:5],
        }

    # secondary unit test: NOT a claim from memory (an earlier draft of this
    # script asserted CTG=1.0 as a "well-known optimal codon" fact -- this was
    # WRONG, caught by checking the actual R reference output directly rather
    # than trusting recollection; the real top codon in the reference 87-tRNA
    # E. coli K-12 vector is AAA (Lys), w=1.0, computed directly from R's own
    # get.ws()). This test instead re-derives the same fact independently: it
    # asserts our full 6-host pipeline's own EC weights (computed from our
    # tRNAscan-SE run, not the reference vector) rank AAA at or near the top,
    # since E. coli's tRNA gene complement is what it is regardless of source.
    ec_w = results["EC"]["codon_weights_w"]
    print(f"\nSecondary unit test: E. coli AAA (Lys) tAI weight = {ec_w['AAA']:.4f} "
          f"(expect 1.0 -- verified against R's own get.ws() output on the reference "
          f"87-tRNA E. coli vector, NOT a memorized claim, see /tmp/tai_r/ecoli_ws_reference.csv)")
    assert abs(ec_w["AAA"] - 1.0) < 1e-9, "AAA should have the maximum tAI weight (1.0) in E. coli"
    print("PASSED")

    with open(OUT / "tai_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {OUT / 'tai_results.json'}")


if __name__ == "__main__":
    main()
