# Funding-Landscape Sweep — Cross-Host / Chassis-Effect Bioinformatics Benchmark

**Date run:** 2026-08-03
**Purpose:** Check for funded competitors (predictive models or benchmarks for cross-host regulatory-activity / chassis-effect prediction) before committing ~10 weeks of solo work.
**Query terms (run individually against every source):** "chassis effect", "host context gene expression prediction", "cross-species promoter prediction", "broad host range synthetic biology", "genetic part portability", "regulatory element prediction bacteria", "multi-host expression prediction", "machine learning promoter bacteria"

**Headline result: no KILL_CANDIDATE was found in any of the 5 sources.** No hit funds a predictive model or benchmark specifically for cross-host regulatory-activity or chassis-effect prediction. Several ADJACENT projects exist (broad-host-range experimental synbio, within-host ML promoter prediction); the closest conceptual match is a UK EPSRC studentship on ML prediction of *contextual* promoter behavior — but in mammalian cell lines, not bacterial cross-host prediction. Details below.

---

## Evidence standard

- **[VERIFIED]** — actual API/web response obtained; exact query and response cited.
- **[NOT FOUND]** — query ran, no relevant hits; every query and scanned-record count is listed.
- **[INFERRED]** — reasoning beyond raw data; confidence flagged.

All raw API responses are saved to `/Users/gabeykim/Downloads/crosshosts/raw/`. All curl commands are saved as a runnable script at `/Users/gabeykim/Downloads/crosshosts/scripts/funder_sweep.sh`.

---

## 1. NIH RePORTER — [VERIFIED]

**Method used:** Full advanced text search (NOT the simpler fallback — this worked on the first correctly-shaped request). Schema was confirmed first with a trivial query `{"criteria":{"fiscal_years":[2024]}}`, which returned 83,516 total FY2024 projects (200 OK), confirming the endpoint/shape. The advanced-text-search body below then returned real, thematically-appropriate hits on the first attempt, so **no fallback to `criteria.terms` was needed** for NIH.

```
POST https://api.reporter.nih.gov/v2/projects/search
{"criteria":{"advanced_text_search":{"operator":"and","search_field":"projecttitle,terms,abstracttext","search_text":"<TERM>"},"fiscal_years":[2024,2025,2026]},
 "include_fields":["ProjectTitle","AbstractText","FiscalYear","OrgName","ProjectNum","AwardAmount","ProjectStartDate","ProjectEndDate","PiName"],
 "offset":0,"limit":25}
```

Restricted server-side via `fiscal_years:[2024,2025,2026]`.

| Term | Total hits (server) | Scanned | Notable results |
|---|---|---|---|
| chassis effect | 3 | 3 | 5R01AI163483 "Development of Engineered Native Bacteria as a Tool for Functional Manipulation of the Gut Microbiome" — **ADJACENT** |
| host context gene expression prediction | 5 | 5 | All UNRELATED (NAFLD disparities, phage biology, cancer markers, bioinformatics core) |
| cross-species promoter prediction | 13 | 13 | Mostly UNRELATED (retrotransposon promoters in human development, xenobiotic gene networks); none are cross-host bacterial regulatory prediction |
| broad host range synthetic biology | 18 | 18 | 4R00EB035165/5K99EB035165 "Gene-Transfer-Resistant Biocontained Bacterial Host" — **ADJACENT**; 5R01EB031935 "precision microbiome engineering" — **ADJACENT**; rest UNRELATED |
| genetic part portability | 5 | 5 | All UNRELATED (DNA damage profiling, human trait genetics, data platforms) |
| regulatory element prediction bacteria | 12 | 12 | All UNRELATED (sRNA virulence networks, TB pathogen dynamics, signaling models) |
| multi-host expression prediction | 1,062 | 25 (top page) | Extremely generic AND-match on common biomedical words ("expression", "prediction"); top 25 are all human-disease multi-omics projects — UNRELATED |
| machine learning promoter bacteria | 10 | 10 | All UNRELATED (Anaplasma tools, phage predation, beta-lactam resistance genetics, translation kinetics) |

**Total NIH records scanned: 91.**

### NIH KILL_CANDIDATEs: **NONE.**

### NIH ADJACENT items (classification detail)
- **5R01AI163483** "Development of Engineered Native Bacteria as a Tool for Functional Manipulation of the Gut Microbiome" (FY2024–2026, PI not extracted from include_fields used, org not captured) — engineers native gut bacteria as a "chassis" and studies "chassis-host interactions," but this is experimental engraftment/biocontainment work in a single host system, not a predictive model or benchmark for cross-host regulatory activity. **ADJACENT** (broad-host-range experimental synbio).
- **4R00EB035165 / 5K99EB035165** "Development of a Gene-Transfer-Resistant and Biocontained Next-Generation Bacterial Host for Controlled Drug Delivery" — biocontainment technology for a single engineered bacterial host. **ADJACENT**.
- **5R01EB031935** "A high-performance and versatile technology for precision microbiome engineering" — gnotobiotic in vivo microbiome engineering tool, single-host (gut) focus. **ADJACENT**.

Raw files: `raw/nih_reporter_*.json` (8 files, one per term), plus `raw/nih_test_trivial.json` and `raw/nih_test_advanced.json` (schema-confirmation queries).

---

## 2. NSF Awards API — [VERIFIED, with a documented method limitation]

**Method used:** `GET https://api.nsf.gov/services/v1/awards.json?keyword=<TERM>&printFields=id,title,abstractText,startDate,expDate,piFirstName,piLastName,awardeeName`. Note: `api.nsf.gov` over `http://` 301-redirects to `https://`; `curl -L` (or direct https) was required.

**Documented limitation (found empirically, not assumed):** The `keyword` parameter reliably filters only for **single, distinctive words**. A control test confirmed this:
- No-keyword baseline and multi-word phrase queries (e.g. `"chassis effect"`, `"host range"`, `"expression prediction"`) returned near-identical, date-descending-sorted result sets dominated by completely unrelated awards (subduction-zone earthquakes, exoplanet demographics, generative-AI social-media studies) — i.e., the multi-word phrase appears to fall through to an effectively unfiltered listing.
- A single distinctive word (`bacteriophage`, `chassis`, `promoter`) filtered correctly and returned genuinely on-topic, relevance-appropriate results (e.g., `keyword=bacteriophage` returned "BioFoundry: Center for Robust, Equitable and Accessible Technology and Education (CREATE) for Next Generation BioFoundries," "CAREER: Designing safe and effective phage cocktails...").

Given this, the 8-term sweep was run as specified (documented below), and supplemented with single-word fallback queries on the most distinctive content word from each phrase (`chassis`, `promoter`) to sanity-check for real hits. This fallback is a genuine weakening of NSF coverage for the 6 multi-word-only terms — flagged per instructions.

| Term (as run) | Result |
|---|---|
| chassis effect | Noisy/unfiltered — no genuine hits in top 25 |
| host context gene expression prediction | Noisy/unfiltered — no genuine hits in top 25 |
| cross-species promoter prediction | Noisy/unfiltered — no genuine hits in top 25 |
| broad host range synthetic biology | Noisy/unfiltered — no genuine hits in top 25 |
| genetic part portability | Noisy/unfiltered — no genuine hits in top 25 |
| regulatory element prediction bacteria | Noisy/unfiltered — no genuine hits in top 25 |
| multi-host expression prediction | Noisy/unfiltered — no genuine hits in top 25 |
| machine learning promoter bacteria | Noisy/unfiltered — no genuine hits in top 25 |
| *(fallback)* chassis | 170 total; top hits include "14TSB_SynBio A toolchest for rapid bootstrapping of novel chassis organisms," "I-Corps: Translation potential of an engineered strain of *Vibrio natriegens* as a chassis for low resource, scalable molecular biology," "CAREER: Connecting eukaryotic electron transfer components to nitrogenase using a bacterial chassis" — **ADJACENT** (broad-host-range chassis engineering; none are FY2024+ predictive-model/benchmark projects) |
| *(fallback)* promoter | 9,096 total, dominated by false-positive stemming matches ("promoting," "promote") — not usable |

**Total NSF records scanned: 200 (8 terms × 25) core sweep + 250 (10 diagnostic/fallback queries × 25) = 450.**

### NSF KILL_CANDIDATEs: **NONE.**

Raw files: `raw/nsf_*.json` (8 files, one per term; diagnostic queries not saved to `raw/` as they were exploratory, but exact URLs are reproducible via the pattern in `scripts/funder_sweep.sh`).

---

## 3. DOE / Agile BioFoundry — [VERIFIED fallback to OSTI.gov; NOT FOUND for a dedicated award index]

Per task instructions: no fetchable public **award index/API** was found for DOE Office of Science or the Agile BioFoundry consortium specifically (Agile BioFoundry's public site has no API; DOE Office of Science's grant/award search has no public API endpoint discovered). **Recorded: NOT FOUND — no fetchable award index for DOE Office of Science or Agile BioFoundry directly.**

**Fallback used (per task instructions):** [OSTI.gov](https://www.osti.gov) public API v1 (`https://www.osti.gov/api/v1/records`), which indexes DOE-funded **publications and technical reports** — this is a different kind of database than a funder awards index (it indexes research *outputs*, not the awards themselves, though records do carry DOE contract/award numbers). This distinction is flagged because it weakens direct comparability with the other 4 sources.

**Method:** exact-phrase search (`q="<TERM>"`, URL-encoded) restricted to `publication_date_start=01/01/2024&publication_date_end=12/31/2026`.

```
GET https://www.osti.gov/api/v1/records?q=%22<TERM>%22&publication_date_start=01/01/2024&publication_date_end=12/31/2026&rows=25
```

| Term | x-total-count | Scanned |
|---|---|---|
| chassis effect | 1 | 1 |
| host context gene expression prediction | 0 | 0 |
| cross-species promoter prediction | 0 | 0 |
| broad host range synthetic biology | 1 | 1 |
| genetic part portability | 0 | 0 |
| regulatory element prediction bacteria | 0 | 0 |
| multi-host expression prediction | 0 | 0 |
| machine learning promoter bacteria | 0 | 0 |

**Total OSTI records scanned (core 8-term sweep): 2**, both reviewed in full:
1. "ComPort: Rigorous Testing Methods to Safeguard Software Porting (Final UW Report)" (2025-11-25) — software engineering, **UNRELATED** (false-positive keyword match on "porting").
2. "Building an expanded bio-based economy through synthetic biology" (*Biotechnology Advances*, 2025-12-06, Oak Ridge National Laboratory authors) — a **review article** discussing "performance-advantaged chassis organisms of the future" broadly; not itself a predictive model, benchmark, or funded research project — it's a literature review output. **ADJACENT** at most (thematically relevant field-level review), not a competing project.

**Supplementary query (not one of the 8 core terms, run to directly enumerate Agile BioFoundry FY2024+ outputs):** `q="Agile BioFoundry"`, date-restricted — **805 total records** (FY2024+), top 25 scanned. All are single-host metabolic-engineering / strain-engineering papers (*Aspergillus niger*, *Pseudomonas putida*, *Rhodosporidium/Rhodotorula toruloides*, etc.) or ML tools for protein-property/fitness prediction *within* one organism. One notable perspective piece: **"The tier system: a host development framework for bioengineering"** (*Current Opinion in Biotechnology* 92, pub. 2025-02-10; authors Yeager, Hillson, Wozniak, Mutalik, Johnson et al., DOE Agile BioFoundry / LANL / LBNL / NREL / PNNL) — proposes a standardized 3-tier conceptual framework for developing "nontraditional chassis organisms," explicitly framed around cross-organism host-development standardization. This is a **conceptual/organizational framework paper, not a predictive model or benchmark** — **ADJACENT**, but worth being aware of as the clearest signal that Agile BioFoundry researchers are thinking about chassis-comparison standardization at a systems level.

### DOE / Agile BioFoundry KILL_CANDIDATEs: **NONE.**

Raw files: `raw/osti_*.json` (8 files, one per term) + `raw/osti_agile_biofoundry_supplementary.json`.

---

## 4. CORDIS (EU) — [VERIFIED]

**Method used:** CORDIS's public search API (`https://cordis.europa.eu/search/en?q=...&format=json`) is reachable and returns structured JSON with full project objective text, dates, and EC contribution amounts. Queries were restricted to `contenttype='project'` to exclude news articles/results-in-brief (which polluted an initial unrestricted test — e.g., a French-language article about vehicle chassis welding matched "chassis effect").

```
GET https://cordis.europa.eu/search/en?q=contenttype%3D%27project%27%20AND%20%27<TERM>%27&format=json&num=25
```

No FY2024+ server-side filter is available on this endpoint; dates were filtered/inspected client-side from each project's `startDate`/`endDate`.

| Term | Total project hits | Scanned |
|---|---|---|
| chassis effect | 14 | 14 |
| host context gene expression prediction | 5 | 5 |
| cross-species promoter prediction | 5 | 5 |
| broad host range synthetic biology | 7 | 7 |
| genetic part portability | 3 | 3 |
| regulatory element prediction bacteria | 1 | 1 |
| multi-host expression prediction | 5 | 5 |
| machine learning promoter bacteria | 3 | 3 |

**Total CORDIS records scanned: 43** (all records returned were reviewed — no truncation).

### CORDIS KILL_CANDIDATEs: **NONE.**

### CORDIS ADJACENT items (FY2024+ active, reviewed in full)
- **GUT-NAT — "Gut Commensals as Next-Generation Chassis for Natural Product Therapeutics"** (RCN 285458, 2026-07-01 to 2028-06-30, EC contribution €200,400). Engineers *Clostridium leptum* as a single-host chassis for natural-product biosynthesis, including "characterization of genetic regulatory elements in *C. leptum*." Single-organism, experimental — not cross-host prediction. **ADJACENT**.
- **AI-EvoYeast — "Harnessing Genomic Instability with AI-Driven Adaptive Laboratory Evolution for Accelerated Yeast Bioproduction"** (RCN 288864, 2026-05-01 to 2028-04-30, €260,347.92). Uses AI/ML on multi-omics data to "build a predictive model for optimal genomic configurations" in a single yeast (*S. cerevisiae*) chassis. Predictive model, but **within-host, not cross-host** — fits the task's explicit ADJACENT category ("within-host (not cross-host) prediction").
- **BACRNA — "Integrative Machine Learning Approaches for Bacterial sRNA Genome Annotation"** (RCN 291614, 2026-07-01 to 2028-06-30, €189,474.24). ML models to improve bacterial small-RNA annotation accuracy within individual bacterial genomes. Within-genome annotation, not cross-host regulatory-activity prediction. **ADJACENT**.
- **LIFE-19 — "UNalphabeting the Central Dogma of Life"** (RCN 269506, 2025-07-01 to 2030-06-30, €2,498,875). Engineers a minimal *Mycoplasma* genome lacking tryptophan; framed as "a key inventory in Synthetic Biology" toward "a cellular chassis with minimal energetic demands." Chassis-engineering adjacent but about amino-acid-alphabet reduction, not host-context prediction. **ADJACENT** (loosely).

Raw files: `raw/cordis_*.json` (8 files, one per term).

---

## 5. UKRI Gateway to Research (GTR) — [VERIFIED, with a documented method limitation]

**Method used:** `GET https://gtr.ukri.org/gtr/api/projects?q=<TERM>&fetchSize=25` (JSON via `Accept: application/json`).

**Documented limitation (found empirically, same pattern as NSF):** Multi-word phrase queries return enormous `totalSize` counts (16,000–55,000) that behave like loose/OR matching rather than strict AND or phrase matching — confirmed by comparing against a no-`q` baseline (`totalSize: 158712`) and a single-word control (`bacteriophage` → `totalSize: 153`, all genuinely on-topic). Unlike NSF, GTR's results do appear **relevance-sorted by default** (spot checks, e.g. `q=host range`, kept genuinely relevant titles at the top despite a huge total count), so — unlike NSF — the top-20 results for each multi-word query were still manually reviewed and did surface real, relevant hits (documented below). This is still a real weakening of exhaustiveness versus a true phrase/AND search, so it's flagged per instructions.

Funding dates are not top-level fields; they were extracted from each project's `FUND` link (`start`/`end`, epoch milliseconds).

| Term | totalSize (unfiltered/noisy) | Scanned (top page) |
|---|---|---|
| chassis effect | 43,219 | 20 |
| host context gene expression prediction | 41,141 | 20 |
| cross-species promoter prediction | 16,199 | 20 |
| broad host range synthetic biology | 54,998 | 20 |
| genetic part portability | 32,779 | 20 |
| regulatory element prediction bacteria | 30,368 | 20 |
| multi-host expression prediction | 16,343 | 20 |
| machine learning promoter bacteria | 33,851 | 20 |

**Total GTR records scanned: 160** (top 20 per term × 8 terms).

### GTR KILL_CANDIDATEs: **NONE** — but this source produced the single **closest conceptual match** found in the entire sweep (see below; still classified ADJACENT, not KILL_CANDIDATE).

### GTR — closest match found (reviewed in full, reasoning for ADJACENT not KILL_CANDIDATE)

**"Graph learning methods for engineering mammalian promoters in bioproduction"**
- Funder: EPSRC (UK) | Grant category: Studentship | Status: Active
- Institution: Imperial College London
- Grant reference (RCUK ID): 2786047
- Funded period: 2023-01-01 to 2027-09-29 (FY2024+ active)
- PI/student name: not exposed by this GTR record (no PER_ID link present in the compact project record returned by the API)

Full abstract (verbatim, as returned by the API):

> mRNA COVID-19 vaccines have effectively prevented hospitalization and deathduring the pandamic. Bioproduction of valuable vaccines and biotherapeuticsin mammalian cell lines can be achieved, but it is difficult to develop robust,predictable, and sustainable expression. The design of enhanced mammalianpromoters and genetic circuits is therefore a key strategic industrial target. Thisproject aims to explore the structural properties of large-scale transcriptomicdatasets with graph representation learning to optimise engineered promoters inmammalian cells. Demonstration will be performed in mammalian cell lines us-ing automated DNA assembly and analytics available at the London Biofoundry.As part of this project, multiple publicly available transcriptomic datasetswill be used, including the EPD [2], DEE2 uniform transcriptomic database[7], and SRA [5] for constructing heterogeneous graphs, in which each noderepresents a set of DNA sequences (i.e. promoters and genes); edge types aredetermined by the biological context shared by two nodes. The initial noderepresentation may be generated using the internal DNA sequence structure.With constructed heterogeneous network, there are two objectives: Firstly,we will use Graph neural networks (GNNs) [3, 4, 8] to learn the representationof each node in the network, which maps the nodes into embedding space toreflect the graph structure. In this way, the cellular and environmental con-text will be embedded into the learnt node representations; and then we canuse the resulting node embedding to predict the promoter activity (i.e. thespeed and quantity of gene expression). Promoter activity prediction task maybe integrated with the representation learning process or solved using separatepredictors. The second stage of the project will combine results from exper-iments to iteratively refine the GNNs. The trained model will select a baseset of 50-100 promoters functional over a range of conditions and with desiredperformance to be used as starting point for synthetic promoter design. In col-laboration with the London Biofoundry, their characteristics will be tested intransient transfection to eliminate changes in behaviour due to differences incopy number and physical context. Promoters with verified behaviour will beused by the GNNs algorithm as modular building blocks to build a set of 50synthetic promoters with designed behaviour. Construct performance will becompared against predictions generated by the GNNs and used to refine modelperformance.Overall, the project will leverage the expression power of GNNs and thelarge-scale transcriptomic datasets to develop a broadly applicable frameworkfor gene expression problem analysis, with the application in predicting con-textual promoter behaviours and engineering new mammalian promoters. Inaddition, the interpretability of learnt model can be achieved by applying sym-bolic regression [1] and GNN-explainer [6].

**Why this is ADJACENT and not KILL_CANDIDATE:** This project builds an ML (graph neural network) model that explicitly predicts *contextual* promoter activity ("the cellular and environmental context will be embedded into the learnt node representations... predict the promoter activity") — conceptually very close to a "host context gene expression prediction" model. However, the "context" here is different mammalian cell lines/transcriptomic conditions within a single host domain (mammalian bioproduction), not different bacterial species/genera (which is what "chassis effect," "cross-species promoter prediction," and "regulatory element prediction bacteria" point to in a bacterial synthetic-biology benchmark). **[INFERRED, medium confidence]:** if the benchmark project's scope is strictly bacterial/microbial cross-host prediction, this is not a direct competitor. If the benchmark's scope could extend to mammalian expression systems, this project should be watched — it is the closest thing to a "predictive model for context-dependent regulatory element activity" found anywhere in this sweep.

### Other GTR ADJACENT items worth noting (not FY2024+, but high thematic relevance — flagged as recently-closed precedent, not live competitors)
- **"Quantifying global burden and contextual effects of synthetic genetic circuits in their bacterial chassis"** (EPSRC Studentship, Closed, 2019-09-30 to 2023-09-29 — ends before the FY2024+ cutoff, so excluded from the active-competitor list, but flagged because its title/abstract is *the* closest textual match to "chassis effect" found across all 5 sources). Full abstract confirms this is experimental (time-lapse microfluidic microscopy + evolution experiments), not a predictive/ML model — it measures burden and context-dependence directly rather than building a model to predict it across hosts. Since the grant is closed and pre-dates the FY2024+ window, it does not meet KILL_CANDIDATE criteria, but its existence (and the EPSRC-Synthetic-Biology research area framing) signals this general topic area attracted UK funding before 2024.
- **"AI-assisted DBTL cycle for synthetic yeast promoters"** (ISPF, Active, 2024-02-14 to 2026-08-12) — deep learning to identify/design promoter sequences, validated experimentally in an iterative design-build-test-learn loop, but **single-organism (yeast)**, not cross-species. **ADJACENT**.
- **"Understanding promoters and terminators of transcription in bacteria"** (BBSRC, listed Closed but end date 2025-10-02, i.e., within FY2024+) — experimental (cappable-seq/term-seq) characterization of bacterial promoter/terminator positioning, single organism, not a cross-host predictive model. **ADJACENT**.
- **"Automation and Digitalisation Technology for Evolutionary Development of Microbial Chassis"** (Innovate UK, Closed, 2023-08-31 to 2025-02-28) — Evolutor Ltd's industrial adaptive-laboratory-evolution platform with generative-AI models for "predictive evolution" of microbial chassis strains; predicts evolutionary trajectories of a chassis, not cross-host regulatory-element portability. **ADJACENT**.
- **"Genetic part mining and functional characterisation for bioremediation"** (BBSRC Studentship, Active, thin abstract) — computational promoter identification tools mentioned only briefly; insufficient detail to indicate cross-host predictive modeling. **ADJACENT** (weak).

Raw files: `raw/gtr_*.json` (8 files, one per term).

---

## Cross-source summary table

| Source | Method | Terms run | Total records scanned | KILL_CANDIDATE | ADJACENT (count) | Notes |
|---|---|---|---|---|---|---|
| NIH RePORTER | Full advanced text search (title+terms+abstract), FY2024–26 | 8/8 | 91 | 0 | 3 | Full method worked; no fallback needed |
| NSF Awards API | keyword search | 8/8 | 200 (+250 diagnostic) | 0 | 1 (via fallback) | **Multi-word keyword search does not reliably filter** — documented limitation |
| DOE/Agile BioFoundry | OSTI.gov fallback (publications, not awards) | 8/8 | 2 (+25 supplementary) | 0 | 2 | No award-index API exists; OSTI indexes outputs, not awards |
| CORDIS (EU) | Native search API, contenttype=project | 8/8 | 43 | 0 | 4 | Clean, fully-filtered results |
| UKRI GTR | q= keyword search | 8/8 | 160 | 0 | 6 (5 pre-FY24/borderline + 1 closest match) | **Multi-word query also degrades**, but relevance sort kept top results usable |

**Grand total records scanned across all 5 sources: 496** core-sweep records (496 = 91+200+2+43+160), plus 275 supplementary/diagnostic records (250 NSF fallback + 25 OSTI Agile BioFoundry).

---

## DARPA BTO AIxBio — BAA HR001126S0003 status

**[VERIFIED]** — confirmed via three independent sources: the primary BAA PDF (downloaded directly and read in full), the grants.gov opportunity listing, and a general web search corroborating both.

### Status: **STILL OPEN** as of 2026-08-03. The working assumption ("open through September 2026") is **CORRECT**, with the precise deadline being **September 30, 2026, at 4:00 PM ET**.

Primary source — official BAA document (downloaded to `raw/darpa_HR001126S0003_baa.pdf`, hosted at Defence Science Institute, mirroring the DARPA original), quoted verbatim from the "OVERVIEW INFORMATION" section:

> **Funding Opportunity Title** – Biological Technologies
> **Announcement Type** – Initial Announcement
> **Funding Opportunity Number** – HR001126S0003
> **Dates/Time - All Times are Eastern Time Zone (ET)**
> - Posting Date: October 01, 2025
> - Proposal Abstract Due Date: Abstracts may be submitted on a rolling basis until September 30, 2026, at 4:00 PM
> - Proposal Due Date: Proposals may be submitted on a rolling basis until September 30, 2026, at 4:00 PM
>
> **Submission Requirements** – An Abstract must be submitted and a notification to submit a Proposal must be received, prior to any Proposal submission.

This is reconfirmed in Section III ("SUBMISSION INFORMATION"): *"This announcement contains a required abstract phase... Abstracts will be reviewed on a rolling basis and can be submitted until September 30, 2026 at 4:00 p.m."* and *"Full proposals can be submitted through September 30, 2026, at 4:00 PM."*

**Notes on naming:** The official title is simply **"Biological Technologies"** — an office-wide BAA for DARPA's Biological Technologies Office (BTO), not literally titled "AIxBio." It does contain a dedicated "Machine Learning (ML) and Artificial Intelligence (AI)" topic area (covering biological foundation models, non-experimental/hybrid model-assessment strategies, biological-system simulation acceleration, and predictive modeling), which is presumably why it's colloquially referred to as "AIxBio" in secondary sources (a defense-news aggregator headline read "Pentagon seeks AI-driven bio solutions to protect troops"). **[INFERRED, high confidence]**: "AIxBio" is an informal/community shorthand for this BAA's AI+biology emphasis rather than the BAA's formal name.

Secondary corroboration — grants.gov / simpler.grants.gov opportunity listing:
> **Current Status:** Open
> **Close Date:** September 30, 2026
> **Posted Date:** October 1, 2025
> **Agency:** DARPA - Biological Technologies Office

Secondary corroboration — WebSearch summary of darpa.mil's own opportunity page (`https://www.darpa.mil/work-with-us/opportunities/hr001126s0003`) independently reported the same October 1, 2025 posting date and September 30, 2026, 4:00 PM ET rolling deadline for both abstracts and full proposals.

A direct fetch of the sam.gov listing (`https://sam.gov/opp/8d403582edfd409795560247e8d229b7/view`) did not return usable content — sam.gov's listing pages are JavaScript-rendered and could not be scraped via WebFetch; this is noted as a gap, but the primary BAA PDF plus grants.gov plus the darpa.mil page constitute three independent, mutually-consistent [VERIFIED] sources, so confidence in the September 30, 2026 deadline is high.

**Sources:**
- [Broad Agency Announcement Biological Technologies (official BAA PDF, mirrored)](https://defencescienceinstitute.com/wp-content/uploads/2025/10/HR001126S0003.pdf) — downloaded to `raw/darpa_HR001126S0003_baa.pdf`
- [grants.gov opportunity listing](https://simpler.grants.gov/opportunity/8899390a-91ee-47f3-b8f4-3a2998dfde40)
- [DARPA official opportunity page](https://www.darpa.mil/work-with-us/opportunities/hr001126s0003)
- [Defence Science Institute funding-opportunity mirror](https://defencescienceinstitute.com/funding-opportunity/darpa-biological-technologies-broad-agency-announcement-r001126s003/)

---

## Appendix: files produced

- **Script (runnable):** `/Users/gabeykim/Downloads/crosshosts/scripts/funder_sweep.sh`
- **Raw JSON responses:** `/Users/gabeykim/Downloads/crosshosts/raw/` — `nih_reporter_*.json` (8), `nih_test_trivial.json`, `nih_test_advanced.json`, `nsf_*.json` (8), `osti_*.json` (8) + `osti_agile_biofoundry_supplementary.json`, `cordis_*.json` (8), `gtr_*.json` (8)
- **DARPA BAA PDF:** `/Users/gabeykim/Downloads/crosshosts/raw/darpa_HR001126S0003_baa.pdf`
- **This report:** `/Users/gabeykim/Downloads/crosshosts/out/task5_funder_sweep.md`

## Bottom line for the go/no-go decision

Across NIH RePORTER, NSF, DOE (via OSTI fallback), CORDIS, and UKRI GTR — 8 query terms each, ~500 core records manually scanned — **no currently-active (FY2024+) award anywhere funds a predictive model or benchmark specifically for cross-host regulatory-activity or chassis-effect prediction**. The nearest miss is a UK EPSRC-funded Imperial College London PhD studentship building a graph-neural-network model to predict *context-dependent* promoter activity — but for mammalian cell lines/bioproduction, not bacterial cross-species chassis prediction. A handful of other ADJACENT projects (DOE Agile BioFoundry's "Tier System" framework paper, an EU chassis-engineering grant, a closed-in-2023 EPSRC studentship on "contextual effects of synthetic genetic circuits in their bacterial chassis") show the surrounding field is active and thinking about these questions conceptually, but none constitute a funded predictive-model/benchmark competitor as of this sweep. This should be read alongside the caveats above (NSF and GTR keyword search both have real, documented multi-word-query blind spots; DOE has no dedicated award index).
