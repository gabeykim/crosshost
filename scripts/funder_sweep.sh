#!/usr/bin/env bash
# Funding-landscape sweep for the cross-host / chassis-effect bioinformatics benchmark project.
# Hits public funder APIs directly for 8 query terms across 5 sources, saves raw JSON to ../raw/,
# and prints a compact summary to stdout.
#
# Usage: bash funder_sweep.sh
# Run from anywhere; paths below are relative to this script's location.

set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"
RAW_DIR="$ROOT_DIR/raw"
mkdir -p "$RAW_DIR"
cd "$ROOT_DIR" || exit 1

TERMS=(
  "chassis effect"
  "host context gene expression prediction"
  "cross-species promoter prediction"
  "broad host range synthetic biology"
  "genetic part portability"
  "regulatory element prediction bacteria"
  "multi-host expression prediction"
  "machine learning promoter bacteria"
)

slugify() { echo "$1" | tr ' ' '_' | tr -cd 'a-zA-Z0-9_'; }
urlencode() { python3 -c "import urllib.parse,sys; print(urllib.parse.quote(sys.argv[1]))" "$1"; }

# ---------------------------------------------------------------------------
# 1. NIH RePORTER  (POST, JSON body, advanced_text_search across title+terms+abstract)
#    Schema check performed first: a trivial {"criteria":{"fiscal_years":[2024]}} query
#    confirmed the endpoint/shape; advanced_text_search with search_field
#    "projecttitle,terms,abstracttext" is the correct field for full-text title+abstract search.
#    This is the FULL method (not the simpler criteria.terms fallback) -- it worked on first try.
# ---------------------------------------------------------------------------
echo "=== NIH RePORTER ==="
# schema-confirmation query (trivial)
curl -s -X POST https://api.reporter.nih.gov/v2/projects/search \
  -H "Content-Type: application/json" \
  -d '{"criteria":{"fiscal_years":[2024]},"include_fields":["ProjectTitle","ProjectNum"],"offset":0,"limit":3}' \
  -o "$RAW_DIR/nih_test_trivial.json"

for term in "${TERMS[@]}"; do
  slug=$(slugify "$term")
  curl -s -X POST https://api.reporter.nih.gov/v2/projects/search \
    -H "Content-Type: application/json" \
    -d "{\"criteria\":{\"advanced_text_search\":{\"operator\":\"and\",\"search_field\":\"projecttitle,terms,abstracttext\",\"search_text\":\"$term\"},\"fiscal_years\":[2024,2025,2026]},\"include_fields\":[\"ProjectTitle\",\"AbstractText\",\"FiscalYear\",\"OrgName\",\"ProjectNum\",\"AwardAmount\",\"ProjectStartDate\",\"ProjectEndDate\",\"PiName\"],\"offset\":0,\"limit\":25}" \
    -o "$RAW_DIR/nih_reporter_${slug}.json"
  echo "  [$term] -> nih_reporter_${slug}.json"
done

# ---------------------------------------------------------------------------
# 2. NSF Awards API (GET, keyword=). NOTE: api.nsf.gov redirects http->https (301); use https
#    directly or curl -L. KNOWN LIMITATION (confirmed empirically in this sweep): the `keyword`
#    param only reliably filters for SINGLE distinctive words (e.g. "bacteriophage", "chassis").
#    Multi-word phrases (2+ words) degrade to an apparently unfiltered, date-sorted result set --
#    confirmed by comparing against a no-keyword control query which returned the same
#    most-recent-award ordering. This is documented in the output report as a fallback/limitation.
# ---------------------------------------------------------------------------
echo "=== NSF Awards API ==="
for term in "${TERMS[@]}"; do
  slug=$(slugify "$term")
  encoded=$(urlencode "$term")
  curl -s -L "https://api.nsf.gov/services/v1/awards.json?keyword=${encoded}&printFields=id,title,abstractText,startDate,expDate,piFirstName,piLastName,awardeeName" \
    -o "$RAW_DIR/nsf_${slug}.json"
  echo "  [$term] -> nsf_${slug}.json"
done

# ---------------------------------------------------------------------------
# 3. DOE / Agile BioFoundry. No dedicated award-search API was found for DOE Office of Science
#    or the Agile BioFoundry consortium specifically. Fallback used per task instructions:
#    OSTI.gov public API (v1), which indexes DOE-funded publications/technical reports (NOT an
#    awards database -- records reference contract/award numbers but this is publication metadata).
#    Exact-phrase search via quoted q= param; publication_date_start/end restrict to FY2024+.
# ---------------------------------------------------------------------------
echo "=== DOE / OSTI.gov (publications fallback) ==="
for term in "${TERMS[@]}"; do
  slug=$(slugify "$term")
  encoded=$(urlencode "\"${term}\"")
  curl -s "https://www.osti.gov/api/v1/records?q=${encoded}&publication_date_start=01/01/2024&publication_date_end=12/31/2026&rows=25" \
    -o "$RAW_DIR/osti_${slug}.json"
  echo "  [$term] -> osti_${slug}.json"
done
# Supplementary: Agile BioFoundry sponsor-acknowledgment search (not one of the 8 core terms,
# but used to directly enumerate FY2024+ Agile BioFoundry outputs for manual review).
curl -s "https://www.osti.gov/api/v1/records?q=%22Agile+BioFoundry%22&publication_date_start=01/01/2024&publication_date_end=12/31/2026&rows=25" \
  -o "$RAW_DIR/osti_agile_biofoundry_supplementary.json"

# ---------------------------------------------------------------------------
# 4. CORDIS (EU). Public search API at cordis.europa.eu/search/en?q=...&format=json works.
#    Restrict to contenttype='project' to avoid news/results-in-brief articles being mixed in.
# ---------------------------------------------------------------------------
echo "=== CORDIS ==="
for term in "${TERMS[@]}"; do
  slug=$(slugify "$term")
  q="contenttype='project' AND '${term}'"
  encoded=$(urlencode "$q")
  curl -s "https://cordis.europa.eu/search/en?q=${encoded}&format=json&num=25" \
    -o "$RAW_DIR/cordis_${slug}.json"
  echo "  [$term] -> cordis_${slug}.json"
done

# ---------------------------------------------------------------------------
# 5. UKRI Gateway to Research. Endpoint: gtr.ukri.org/gtr/api/projects?q=...
#    KNOWN LIMITATION (confirmed empirically, same pattern as NSF): the `q` param filters well
#    for single distinctive words but multi-word phrases return huge totalSize counts that look
#    like relevance-ranked-but-noisy OR matching rather than strict AND/phrase matching. Unlike
#    NSF, GTR's default sort does appear to be relevance-based (top results stayed thematically
#    on-topic in spot checks), so results were still manually reviewed for genuine hits.
#    Funding dates come from the FUND link's start/end epoch-ms timestamps, not top-level fields.
# ---------------------------------------------------------------------------
echo "=== UKRI Gateway to Research ==="
for term in "${TERMS[@]}"; do
  slug=$(slugify "$term")
  encoded=$(urlencode "$term")
  curl -s -H "Accept: application/json" "https://gtr.ukri.org/gtr/api/projects?q=${encoded}&fetchSize=25" \
    -o "$RAW_DIR/gtr_${slug}.json"
  echo "  [$term] -> gtr_${slug}.json"
done

echo
echo "Done. Raw JSON saved under $RAW_DIR"
