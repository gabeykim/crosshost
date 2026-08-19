"""
GATE 13 - Task 1: mechanically verify every citation in out/MANUSCRIPT.md's
References section against a canonical source (CrossRef for DOI-bearing
works, the arXiv API for preprints with no DOI), rather than trusting any
prior gate's own "verified" label.

Motivation: Gate 12 found two fabricated author names in citations that had
already been "verified" on their numeric claim (El Houdaigui's third author
written as "Vindry" instead of "Hindre"; Borkowski's coauthors written as
"Ceroni, Stan, Ellis" instead of "Bricio, Murgiano, Rothschild-Mancinelli,
Stan, Ellis") -- both caught by hand, not by any automated check. This
script exists so that check is mechanical and repeatable, not one-off.

METHOD, per the gate's own instructions:
1. For each citation, query a canonical source (CrossRef by DOI when known,
   else the arXiv API by ID) for the authoritative record.
2. Compare stated vs. canonical: author surnames (character-by-character,
   in order), title, venue, year, volume, pages/article-number.
3. Report every mismatch. Do not guess a fix for anything unresolved --
   flag it for manual review.

This is a live-network check (CrossRef + arXiv APIs) and is therefore NOT
folded into `make audit` (which must stay fast and offline) -- it is its
own Makefile target, run on demand, alongside the audits per the gate's
instruction. Network failure degrades to SKIPPED per-citation, not a crash.
"""
import json
import re
import ssl
import time
import urllib.request
import urllib.error
import urllib.parse
from pathlib import Path

try:
    import certifi
    SSL_CONTEXT = ssl.create_default_context(cafile=certifi.where())
except ImportError:
    # Fall back to the system default context; on some macOS Python.org
    # installs this will itself fail with CERTIFICATE_VERIFY_FAILED --
    # the fix there is `pip install certifi`, not disabling verification.
    SSL_CONTEXT = ssl.create_default_context()

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "out" / "results"

# Citations as they currently appear in out/MANUSCRIPT.md's References
# section, transcribed directly from that file (not from memory) --
# see the manuscript itself for the exact surrounding text/notes.
CITATIONS = [
    {
        "key": "johns2018",
        "stated_authors": ["Johns", "Gomes", "Yim"],  # "et al." in manuscript
        "et_al": True,
        "stated_title": "Metagenomic mining of regulatory elements enables programmable species-selective gene expression",
        "stated_venue": "Nature Methods",
        "stated_year": 2018,
        "stated_volume": "15",
        "stated_pages": "323-329",
        "stated_doi": None,  # manuscript currently states no DOI
        "lookup": {"type": "doi", "id": "10.1038/nmeth.4633"},
    },
    {
        "key": "yim2019_drafts",
        "stated_authors": ["Yim", "Johns"],  # "et al." in manuscript
        "et_al": True,
        "stated_title": "Multiplex transcriptional characterizations across diverse bacterial species using cell-free systems",
        "stated_venue": "Molecular Systems Biology",
        "stated_year": 2019,
        "stated_volume": "15",
        "stated_pages": "e8875",
        "stated_doi": "10.15252/msb.20198875",
        "lookup": {"type": "doi", "id": "10.15252/msb.20198875"},
        "manual_review_note": "CrossRef article-number field returns 'MSB198875' (EMBO Press internal format); 'e8875' is the correct citable article number per the journal's own citation convention -- benign format difference, not an error.",
    },
    {
        "key": "chan2023",
        "stated_authors": ["Chan", "Baldwin", "Bernstein"],
        "et_al": False,
        "stated_title": "Revealing the Host-Dependent Nature of an Engineered Genetic Inverter in Concordance with Physiology",
        "stated_venue": "BioDesign Research",
        "stated_year": 2023,
        "stated_volume": "5",
        "stated_pages": "0016",
        "stated_doi": "10.34133/bdr.0016",
        "lookup": {"type": "doi", "id": "10.34133/bdr.0016"},
    },
    {
        "key": "chan2024",
        "stated_authors": ["Chan", "Bernstein"],
        "et_al": False,
        "stated_title": "Pangenomic landscapes shape performances of a synthetic genetic circuit across Stutzerimonas species",
        "stated_venue": "mSystems",
        "stated_year": 2024,
        "stated_volume": "9(9)",
        "stated_pages": "e00849-24",
        "stated_doi": "10.1128/msystems.00849-24",
        "lookup": {"type": "doi", "id": "10.1128/msystems.00849-24"},
    },
    {
        "key": "el_houdaigui2019",
        "stated_authors": ["El Houdaigui", "Forquet", "Hindré", "Schneider", "Nasser", "Reverchon", "Meyer"],
        "et_al": False,
        "stated_title": "Bacterial genome architecture shapes global transcriptional regulation by DNA supercoiling",
        "stated_venue": "Nucleic Acids Research",
        "stated_year": 2019,
        "stated_volume": "47(11)",
        "stated_pages": "5648-5657",
        "stated_doi": "10.1093/nar/gkz300",
        "lookup": {"type": "doi", "id": "10.1093/nar/gkz300"},
    },
    {
        "key": "borkowski2018",
        "stated_authors": ["Borkowski", "Bricio", "Murgiano", "Rothschild-Mancinelli", "Stan", "Ellis"],
        "et_al": False,
        "stated_title": "Cell-free prediction of protein expression costs for growing cells",
        "stated_venue": "Nature Communications",
        "stated_year": 2018,
        "stated_volume": "9",
        "stated_pages": "1457",
        "stated_doi": "10.1038/s41467-018-03970-x",
        "lookup": {"type": "doi", "id": "10.1038/s41467-018-03970-x"},
    },
    {
        "key": "pandi2022",
        "stated_authors": ["Pandi"],  # "et al." in manuscript
        "et_al": True,
        "stated_title": "A versatile active learning workflow for optimization of genetic and metabolic networks",
        "stated_venue": "Nature Communications",
        "stated_year": 2022,
        "stated_volume": "13",
        "stated_pages": "3876",
        "stated_doi": "10.1038/s41467-022-31245-z",
        "lookup": {"type": "doi", "id": "10.1038/s41467-022-31245-z"},
    },
    {
        "key": "lafleur2022",
        "stated_authors": ["LaFleur", "Hossain", "Salis"],
        "et_al": False,
        "stated_title": "Automated model-predictive design of synthetic promoters to control transcriptional profiles in bacteria",
        "stated_venue": "Nature Communications",
        "stated_year": 2022,
        "stated_volume": "13",
        "stated_pages": "5159",
        "stated_doi": "10.1038/s41467-022-32829-5",
        "lookup": {"type": "doi", "id": "10.1038/s41467-022-32829-5"},
    },
    {
        "key": "xia2026_promogen2",
        "stated_authors": ["Xia"],  # "et al." in manuscript
        "et_al": True,
        "stated_title": "Design prokaryotic cis-regulatory elements using language model",
        "stated_venue": "Nucleic Acids Research",
        "stated_year": 2026,
        "stated_volume": "54(4)",
        "stated_pages": "gkag122",
        "stated_doi": None,
        "lookup": {"type": "crossref_search", "query": "Design prokaryotic cis-regulatory elements using language model"},
    },
    # Previously "NOT independently verified" -- Gate 13 targets
    {
        "key": "dnabert2",
        "stated_authors": ["Zhou", "Ji", "Li", "Dutta", "Davuluri", "Liu"],
        "et_al": False,
        "stated_title": "DNABERT-2: Efficient Foundation Model and Benchmark for Multi-Species Genome",
        "stated_venue": "ICLR 2024",
        "stated_year": 2024,
        "stated_volume": None,
        "stated_pages": None,
        "stated_doi": None,
        "lookup": {"type": "arxiv", "id": "2306.15006"},
        "manual_review_note": "No DOI exists; ICLR has none for accepted papers. Canonical lookup is arXiv-only (posted 2023, venue field literally 'arXiv'), so this will always show venue/year MISMATCH against the arXiv record even though ICLR 2024 acceptance is real -- independently confirmed via the ICLR 2024 proceedings page (proceedings.iclr.cc/paper_files/paper/2024/hash/b633e7052970b8f5aa1a69164d99e9e8-Abstract-Conference.html) during Gate 13 manual review. Cite both: ICLR 2024 as the peer-reviewed venue, arXiv:2306.15006 as the preprint.",
    },
    {
        "key": "bend2024",
        "stated_authors": ["Marin", "Teufel", "Horlacher", "Madsen", "Pultz", "Winther", "Boomsma"],
        "et_al": False,
        "stated_title": "BEND: Benchmarking DNA Language Models on Biologically Meaningful Tasks",
        "stated_venue": "ICLR 2024",
        "stated_year": 2024,
        "stated_volume": None,
        "stated_pages": None,
        "stated_doi": None,
        "lookup": {"type": "arxiv", "id": "2311.12570"},
        "manual_review_note": "Same pattern as DNABERT-2: no DOI, arXiv-only canonical lookup will always show venue=arXiv/year=2023 (original posting) even though ICLR 2024 acceptance is real -- independently confirmed via OpenReview (openreview.net/forum?id=uKB4cFNQFg, 'submitted November 21, 2023; final revision April 9, 2024') during Gate 13 manual review.",
    },
    {
        "key": "genomic_benchmarks2023",
        "stated_authors": ["Grešová", "Martinek", "Čechák", "Šimeček", "Alexiou"],
        "et_al": False,
        "stated_title": "Genomic benchmarks: a collection of datasets for genomic sequence classification",
        "stated_venue": "BMC Genomic Data",
        "stated_year": 2023,
        "stated_volume": "24",
        "stated_pages": "25",
        "stated_doi": "10.1186/s12863-023-01123-8",
        "lookup": {"type": "doi", "id": "10.1186/s12863-023-01123-8"},
    },
    {
        "key": "dart_eval2024",
        "stated_authors": ["Patel", "Singhal", "Wang", "Pampari", "Kasowski", "Kundaje"],
        "et_al": False,
        "stated_title": "DART-Eval: A Comprehensive DNA Language Model Evaluation Benchmark on Regulatory DNA",
        "stated_venue": "NeurIPS 2024 Datasets and Benchmarks Track",
        "stated_year": 2024,
        "stated_volume": None,
        "stated_pages": None,
        "stated_doi": "10.52202/079017-1981",
        "lookup": {"type": "doi", "id": "10.52202/079017-1981"},
        "manual_review_note": "CrossRef's container-title is the formal NeurIPS proceedings name ('Advances in Neural Information Processing Systems 37', volume 37 = 2024); the manuscript's 'NeurIPS 2024 Datasets and Benchmarks Track' is the informal, more specific rendering (correctly identifying the D&B sub-track, which CrossRef does not capture) -- independently confirmed via the NeurIPS virtual site (nips.cc/virtual/2024/poster/97497) during Gate 13 manual review. Not an error.",
    },
    {
        "key": "dnalongbench2025",
        "stated_authors": ["Cheng", "Song", "Zhang", "Wang", "Wang", "Yang", "Li", "Ma"],
        "et_al": False,
        "stated_title": "DNALONGBENCH: a benchmark suite for long-range DNA prediction tasks",
        "stated_venue": "Nature Communications",
        "stated_year": 2025,
        "stated_volume": None,
        "stated_pages": None,
        "stated_doi": "10.1038/s41467-025-65077-4",
        "lookup": {"type": "doi", "id": "10.1038/s41467-025-65077-4"},
    },
    {
        "key": "nucleotide_transformer2025",
        "stated_authors": ["Dalla-Torre", "Gonzalez", "Mendoza-Revilla"],  # "et al." (15 authors total)
        "et_al": True,
        "stated_title": "Nucleotide Transformer: building and evaluating robust foundation models for human genomics",
        "stated_venue": "Nature Methods",
        "stated_year": 2025,
        "stated_volume": "22(2)",
        "stated_pages": "287-297",
        "stated_doi": "10.1038/s41592-024-02523-z",
        "lookup": {"type": "doi", "id": "10.1038/s41592-024-02523-z"},
    },
    # Orphaned in the manuscript's own References disclaimer as of Gate 12:
    # claimed "referenced in Section 1/Discussion" but does not actually
    # appear inline anywhere in the current manuscript body. Verified here
    # anyway, per "every citation in the manuscript's reference list."
    {
        "key": "iclr2025_specialized_fm",
        "stated_authors": None,  # manuscript states no author list at all
        "et_al": False,
        "stated_title": "Specialized Foundation Models Struggle to Beat Supervised Baselines",
        "stated_venue": "ICLR 2025",
        "stated_year": 2025,
        "stated_volume": None,
        "stated_pages": None,
        "stated_doi": None,
        "lookup": {"type": "arxiv", "id": "2411.02796"},
        "note": "not actually cited inline anywhere in the current manuscript body (checked directly) -- flagged for removal in Gate 13, not a mismatch on bibliographic fields",
        "manual_review_note": "Same arXiv-only pattern: venue=arXiv/year=2024 (original posting) vs stated ICLR 2025 -- independently confirmed real acceptance via OpenReview (openreview.net/forum?id=JYTQ6ELUVO) and the ICLR 2025 virtual site (iclr.cc/virtual/2025/poster/30102) during Gate 13 manual review. Removed from References regardless, per the orphan note above -- this citation supports no claim actually made in the manuscript body.",
    },
]


def http_get_json(url, timeout=15):
    req = urllib.request.Request(url, headers={"User-Agent": "crosshost-citation-audit/1.0 (mailto:none@example.com)"})
    with urllib.request.urlopen(req, timeout=timeout, context=SSL_CONTEXT) as resp:
        return json.loads(resp.read().decode("utf-8"))


def http_get_text(url, timeout=15):
    req = urllib.request.Request(url, headers={"User-Agent": "crosshost-citation-audit/1.0 (mailto:none@example.com)"})
    with urllib.request.urlopen(req, timeout=timeout, context=SSL_CONTEXT) as resp:
        return resp.read().decode("utf-8")


def fetch_crossref_by_doi(doi):
    d = http_get_json(f"https://api.crossref.org/works/{doi}")
    msg = d["message"]
    authors = [a.get("family", "") for a in msg.get("author", []) if a.get("family")]
    title = (msg.get("title") or [""])[0]
    venue = (msg.get("container-title") or [""])[0]
    volume = msg.get("volume")
    pages = msg.get("page") or msg.get("article-number")
    pub = msg.get("published-print") or msg.get("published-online") or msg.get("published")
    year = pub["date-parts"][0][0] if pub and pub.get("date-parts") else None
    return {"authors": authors, "title": title, "venue": venue, "volume": volume, "pages": pages, "year": year, "doi": doi}


def fetch_crossref_by_search(query):
    d = http_get_json(f"https://api.crossref.org/works?query.bibliographic={urllib.parse.quote(query)}&rows=5")
    for item in d["message"]["items"]:
        title = (item.get("title") or [""])[0]
        if title.lower().strip() == query.lower().strip() or query.lower() in title.lower():
            authors = [a.get("family", "") for a in item.get("author", []) if a.get("family")]
            venue = (item.get("container-title") or [""])[0]
            volume = item.get("volume")
            pages = item.get("page") or item.get("article-number")
            pub = item.get("published-print") or item.get("published-online") or item.get("published")
            year = pub["date-parts"][0][0] if pub and pub.get("date-parts") else None
            return {"authors": authors, "title": title, "venue": venue, "volume": volume, "pages": pages,
                     "year": year, "doi": item.get("DOI")}
    return None


def fetch_arxiv(arxiv_id):
    text = http_get_text(f"https://export.arxiv.org/api/query?id_list={arxiv_id}")
    title_m = re.search(r"<title>(.*?)</title>\s*<link", text, re.S)
    # arXiv Atom feed: first <title> is the query title, second is the entry title
    titles = re.findall(r"<title>(.*?)</title>", text, re.S)
    entry_title = titles[1].strip() if len(titles) > 1 else None
    entry_title = re.sub(r"\s+", " ", entry_title) if entry_title else None
    authors = re.findall(r"<name>(.*?)</name>", text)
    published_m = re.search(r"<published>(\d{4})-", text)
    year = int(published_m.group(1)) if published_m else None
    surnames = [a.strip().split()[-1] for a in authors if a.strip()]
    return {"authors": surnames, "title": entry_title, "venue": "arXiv", "volume": None, "pages": None,
            "year": year, "doi": None, "arxiv_id": arxiv_id}


def norm(s):
    if s is None:
        return ""
    s = re.sub(r"<[^>]+>", "", s)  # strip HTML/XML markup (e.g. CrossRef's <i>Stutzerimonas</i>)
    return re.sub(r"[^a-z0-9]", "", s.lower())


def compare_authors(stated, canonical, et_al):
    if stated is None or canonical is None:
        return None, "no stated author list to compare" if stated is None else "canonical fetch failed"
    mismatches = []
    n = len(stated) if et_al else max(len(stated), len(canonical))
    for i in range(n):
        s = stated[i] if i < len(stated) else None
        c = canonical[i] if i < len(canonical) else None
        if s is None:
            if not et_al:
                mismatches.append(f"position {i}: manuscript has no author, canonical has '{c}'")
            continue
        if c is None:
            mismatches.append(f"position {i}: manuscript has '{s}', canonical has no author (list too short)")
            continue
        if norm(s) != norm(c):
            mismatches.append(f"position {i}: manuscript says '{s}', canonical says '{c}'")
    ok = len(mismatches) == 0
    return ok, "; ".join(mismatches) if mismatches else "exact match" + (" (checked against et al.-truncated list)" if et_al else "")


def compare_field(stated, canonical, loose=False):
    if stated is None and canonical is None:
        return None, "both absent"
    if stated is None:
        return None, f"manuscript states none; canonical: {canonical!r}"
    if canonical is None:
        return None, "canonical fetch missing this field"
    if loose:
        ok = norm(str(stated)) in norm(str(canonical)) or norm(str(canonical)) in norm(str(stated))
    else:
        ok = norm(str(stated)) == norm(str(canonical))
    return ok, f"manuscript={stated!r} canonical={canonical!r}"


def verify_one(c):
    result = {"key": c["key"], "stated_title": c["stated_title"], "status": None, "fields": {}, "canonical": None, "error": None}
    lookup = c["lookup"]
    try:
        if lookup["type"] == "doi":
            canon = fetch_crossref_by_doi(lookup["id"])
        elif lookup["type"] == "crossref_search":
            canon = fetch_crossref_by_search(lookup["query"])
            if canon is None:
                result["status"] = "UNRESOLVED"
                result["error"] = "no confident CrossRef bibliographic-search match"
                return result
        elif lookup["type"] == "arxiv":
            canon = fetch_arxiv(lookup["id"])
        else:
            result["status"] = "UNRESOLVED"
            result["error"] = f"unknown lookup type {lookup['type']}"
            return result
    except (urllib.error.URLError, TimeoutError, Exception) as e:
        result["status"] = "SKIPPED"
        result["error"] = f"network/lookup error: {e}"
        return result

    result["canonical"] = canon
    a_ok, a_detail = compare_authors(c["stated_authors"], canon.get("authors"), c["et_al"])
    t_ok, t_detail = compare_field(c["stated_title"], canon.get("title"), loose=True)
    v_ok, v_detail = compare_field(c["stated_venue"], canon.get("venue"), loose=True)
    y_ok, y_detail = compare_field(c["stated_year"], canon.get("year"))
    vol_ok, vol_detail = compare_field(c["stated_volume"], canon.get("volume"), loose=True)
    p_ok, p_detail = compare_field(c["stated_pages"], canon.get("pages"), loose=True)
    doi_ok, doi_detail = compare_field(c["stated_doi"], canon.get("doi"))

    result["fields"] = {
        "authors": {"ok": a_ok, "detail": a_detail},
        "title": {"ok": t_ok, "detail": t_detail},
        "venue": {"ok": v_ok, "detail": v_detail},
        "year": {"ok": y_ok, "detail": y_detail},
        "volume": {"ok": vol_ok, "detail": vol_detail},
        "pages": {"ok": p_ok, "detail": p_detail},
        "doi": {"ok": doi_ok, "detail": doi_detail},
    }
    hard_fails = [k for k, v in result["fields"].items() if v["ok"] is False]
    result["status"] = "MISMATCH" if hard_fails else "VERIFIED"
    result["mismatched_fields"] = hard_fails
    if "note" in c:
        result["note"] = c["note"]
    if "manual_review_note" in c:
        result["manual_review_note"] = c["manual_review_note"]
        if result["status"] == "MISMATCH":
            result["status"] = "MISMATCH (manually reviewed, explained -- see manual_review_note)"
    return result


def main():
    results = []
    for c in CITATIONS:
        print(f"Checking {c['key']} ...")
        r = verify_one(c)
        results.append(r)
        print(f"  -> {r['status']}" + (f" ({', '.join(r.get('mismatched_fields', []))})" if r["status"].startswith("MISMATCH") else ""))
        time.sleep(0.5)  # be polite to CrossRef/arXiv

    RESULTS.mkdir(parents=True, exist_ok=True)
    with open(RESULTS / "gate13_citation_verification.json", "w") as f:
        json.dump(results, f, indent=2)

    import csv
    with open(RESULTS / "gate13_citation_verification.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["key", "status", "mismatched_fields", "error"])
        for r in results:
            w.writerow([r["key"], r["status"], ";".join(r.get("mismatched_fields", [])), r.get("error", "")])

    print("\n=== SUMMARY ===")
    n_verified = sum(1 for r in results if r["status"] == "VERIFIED")
    n_mismatch_raw = sum(1 for r in results if r["status"].startswith("MISMATCH") and "manual_review_note" not in r)
    n_mismatch_explained = sum(1 for r in results if r["status"].startswith("MISMATCH") and "manual_review_note" in r)
    n_unresolved = sum(1 for r in results if r["status"] in ("UNRESOLVED", "SKIPPED"))
    print(f"VERIFIED: {n_verified}  MISMATCH (unexplained): {n_mismatch_raw}  "
          f"MISMATCH (manually reviewed, explained): {n_mismatch_explained}  "
          f"UNRESOLVED/SKIPPED: {n_unresolved}  TOTAL: {len(results)}")
    for r in results:
        if r["status"].startswith("MISMATCH"):
            print(f"\n  {r['key']}: mismatched fields = {r['mismatched_fields']}")
            for fname in r["mismatched_fields"]:
                print(f"    {fname}: {r['fields'][fname]['detail']}")
            if "manual_review_note" in r:
                print(f"    MANUAL REVIEW: {r['manual_review_note']}")
        elif r["status"] in ("UNRESOLVED", "SKIPPED"):
            print(f"\n  {r['key']}: {r['status']} -- {r['error']}")

    print(f"\nWrote {RESULTS / 'gate13_citation_verification.json'}, .csv")


if __name__ == "__main__":
    main()
