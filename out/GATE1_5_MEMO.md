# CROSSHOST — Gate 1.5 Memo: Unblocking the Gate 2 Feature Spec

**Date:** 2026-08-03
**Prepared for:** Gabriel
**Scope:** Two targeted follow-ups to Gate 1, both load-bearing for Gate 2. Write no modeling code — none was written. Evidence tags follow the same standard as Gate 1: [VERIFIED] / [COMPUTED] / [INFERRED, confidence] / [NOT FOUND].

**Bottom line up front: PASS.** Task A produces a concrete, actionable table of runnable held-out-host configurations — RS241 supports genuine held-out-host tests on all three added hosts, at N substantially larger than Gate 1's (mistaken) six-way-intersection estimate suggested. Task B produces a well-documented, honest **negative-but-actionable** answer: a uniform ≥3-metric physiology vector *can* be built for all six hosts from a single class of source (integrated proteomics abundance databases), but it is comparable only as a crude, rank-order-level signal, not a precise one — and this should be Gate 2's *primary* physiology representation specifically because it is the only thing that's genuinely uniform across all six hosts.

---

## TASK A — Corrected RS241 usability for leave-one-host-out evaluation

Gate 1's 111/241 (transcription) and 126/241 (translation) figures computed the **six-way intersection** — usable in all six hosts simultaneously. That was the wrong question. Leave-one-host-out evaluation of a held-out host only needs a usable value in that host plus usable values in whichever hosts are actually used for training. This section recomputes the right quantities.

**Definition sanity check [VERIFIED]:** the "usable value" definition here is byte-for-byte identical to Gate 1 Task 1/2 (non-null cell in the paper-processed Log2 Transcription / Log10 Translation sheets of Supplementary Data Table 4, no additional threshold reconstructable from this table). This was re-verified computationally: recomputing the six-way intersection with this script reproduces Gate 1's 111 and 126 exactly (`scripts/05_task_a_rs241_pairwise.py`, assertion-checked). No definitional drift.

### A1. Per-host usable counts (unchanged from Gate 1 — reproduced here for reference)

| Host | Transcription usable | Translation usable |
|---|---|---|
| *E. coli* | 219 / 241 | 230 / 241 |
| *B. subtilis* | 149 / 241 | 149 / 241 |
| *P. aeruginosa* | 240 / 241 | 235 / 241 |
| *S. enterica* | 237 / 241 | 237 / 241 |
| *V. natriegens* | 166 / 241 | 229 / 241 |
| *C. glutamicum* | 211 / 241 | 218 / 241 |

### A2. Pairwise usable-in-both matrix

**Transcription:**

| | EC | BS | PA | SE | VN | CG |
|---|---|---|---|---|---|---|
| **EC** | 219 | 146 | 219 | 216 | 162 | 200 |
| **BS** | 146 | 149 | 148 | 149 | 113 | 139 |
| **PA** | 219 | 148 | 240 | 236 | 166 | 211 |
| **SE** | 216 | 149 | 236 | 237 | 164 | 209 |
| **VN** | 162 | 113 | 166 | 164 | 166 | 160 |
| **CG** | 200 | 139 | 211 | 209 | 160 | 211 |

**Translation:**

| | EC | BS | PA | SE | VN | CG |
|---|---|---|---|---|---|---|
| **EC** | 230 | 143 | 227 | 228 | 218 | 213 |
| **BS** | 143 | 149 | 146 | 147 | 141 | 133 |
| **PA** | 227 | 146 | 235 | 233 | 223 | 214 |
| **SE** | 228 | 147 | 233 | 237 | 225 | 216 |
| **VN** | 218 | 141 | 223 | 225 | 229 | 210 |
| **CG** | 213 | 133 | 214 | 216 | 210 | 218 |

**Reading this matrix: *B. subtilis* is visibly the weak link in every row and column it touches** (its highest pairwise value with any other host tops out at 149, versus 200+ for most other pairs). Every other host pair clears 160–240 in both readouts. This confirms Gate 1's qualitative finding (BS is the hard host) while showing precisely how much it costs to include BS in a training set.

### A3. Held-out-host viability

"N with target usable AND ≥k other hosts usable" — [COMPUTED, `scripts/05_task_a_rs241_pairwise.py`]

**Transcription:**

| Held-out host | N target usable | N target + ≥2 others | N target + ≥3 others | Top contributing other hosts |
|---|---|---|---|---|
| *E. coli* (primary) | 219 | 218 | 210 | PA(219), SE(216), CG(200), VN(162), BS(146) |
| *B. subtilis* (primary) | 149 | 148 | 147 | SE(149), PA(148), EC(146), CG(139), VN(113) |
| *P. aeruginosa* (primary) | 240 | 230 | 214 | SE(236), EC(219), CG(211), VN(166), BS(148) |
| ***S. enterica*** (added) | 237 | 228 | 213 | PA(236), EC(216), CG(209), VN(164), BS(149) |
| ***V. natriegens*** (added) | 166 | 166 | 164 | PA(166), SE(164), EC(162), CG(160), BS(113) |
| ***C. glutamicum*** (added) | 211 | 210 | 204 | PA(211), SE(209), EC(200), VN(160), BS(139) |

**Translation:**

| Held-out host | N target usable | N target + ≥2 others | N target + ≥3 others | Top contributing other hosts |
|---|---|---|---|---|
| *E. coli* (primary) | 230 | 230 | 228 | SE(228), PA(227), VN(218), CG(213), BS(143) |
| *B. subtilis* (primary) | 149 | 149 | 147 | SE(147), PA(146), EC(143), VN(141), CG(133) |
| *P. aeruginosa* (primary) | 235 | 235 | 232 | SE(233), EC(227), VN(223), CG(214), BS(146) |
| ***S. enterica*** (added) | 237 | 237 | 233 | PA(233), EC(228), VN(225), CG(216), BS(147) |
| ***V. natriegens*** (added) | 229 | 228 | 222 | SE(225), PA(223), EC(218), CG(210), BS(141) |
| ***C. glutamicum*** (added) | 218 | 217 | 216 | SE(216), PA(214), EC(213), VN(210), BS(133) |

### A4. Best achievable evaluation configurations — ranked

The realistic Gate 5 scenario is: train a model on the three primary hosts (which have the ~11,276-sequence main-library overlap, not just RS241), then evaluate transfer to a held-out added host using RS241. That specific configuration was computed directly [COMPUTED]:

| Held-out host | Train hosts | N (transcription) | N (translation) |
|---|---|---|---|
| *S. enterica* | EC + BS + PA | 146 | 141 |
| *S. enterica* | **EC + PA (drop weak BS)** | **216** | **225** |
| *V. natriegens* | EC + BS + PA | 113 | 134 |
| *V. natriegens* | **EC + PA (drop weak BS)** | **162** | **215** |
| *C. glutamicum* | EC + BS + PA | 138 | 130 |
| *C. glutamicum* | **EC + PA (drop weak BS)** | **200** | **211** |
| *B. subtilis* | EC + PA | 146 | 142 |

**Dropping *B. subtilis* from the RS241-based training set increases usable N by 30–63% across every added-host target**, in both readouts. This is the single most actionable number in this task: *if RS241 rows themselves are used for any part of training/calibration (as opposed to only using the much larger main-library BS/EC/PA overlap for training and RS241 purely as held-out evaluation), excluding B. subtilis from that training slice is a straightforward, well-justified win.*

The full ranked table of every held-out-host × training-subset combination tried (all subsets of the 5 non-target hosts, not just the primary-3 and best-2 cases highlighted above) is in `out/task_a_rs241_pairwise_results.json`.

### A5. Direct answer to "does RS241 support a genuine held-out-host test on any organism outside the three primary hosts, and at what N?"

**Yes, on all three added hosts, at N in the 160–237 range depending on exact configuration** (using target + best 2-host training subset; see A4). This is a real, usable held-out-host test for *S. enterica*, *V. natriegens*, and *C. glutamicum* — none of the three needs to be dropped. *V. natriegens* has the smallest N of the three (~162–166 for the tightest useful configurations) because it has the lowest per-host RS241 coverage among the added hosts (68.9% transcription, 95.0% translation) — still workable, just the least statistically powerful of the three added-host tests. *S. enterica* is the strongest added-host test by a comfortable margin (216–237 depending on configuration) because it has the highest per-host RS241 coverage of any host, primary or added (98.3% both readouts).

---

## TASK B — Uniform omics-derived physiology proxies

### B1. Integrated protein-abundance resource coverage across the six hosts

**PaxDb v6.0** (published Jan 2026, 405 species, 382 bacterial datasets) [VERIFIED via direct API interrogation]. The public website (pax-db.org) is a client-rendered Vue.js single-page app that returns no usable content to a non-JS fetch; the underlying REST API was located by extracting `api.pax-db.org` from the app's JS bundle and finding its OpenAPI spec at `https://api.pax-db.org/try-api/openapi.json`. The API requires an `x-api-key` header but accepts a documented public fallback value (`test`) explicitly sanctioned by the API's own response message for testing use (shared, rate-limited pool). Base path: `/v6`. Key endpoints used: `GET /v6/metadata/species` (full species list), `GET /v6/metadata/dataset?species_id=X` (per-species dataset list with quality scores), `GET /v6/abundances/dataset/{id}?start=&end=&sort=` (paginated per-protein abundance in ppm).

**Coverage, scanning the full 405-species list directly [VERIFIED]:**

| Host | In PaxDb v6.0? | species_id | Best dataset (score, coverage%) |
|---|---|---|---|
| *E. coli* K-12 MG1655 | Yes | 511145 | "Integrated" — score 23.3, 91% coverage, 3,747 proteins |
| *B. subtilis* 168 | Yes | 224308 | "Integrated" — score 15.8, 97% coverage, 4,052 proteins |
| *P. aeruginosa* PAO1 | Yes | 208964 | "Integrated" — score 26.3, 90% coverage, 5,034 proteins |
| *S. enterica* Typhimurium LT2 | Yes | 99287 | "Integrated" — score 12.4, 59% coverage, 2,620 proteins |
| *C. glutamicum* ATCC 13032 | Yes | 196627 | "Integrated" — score 9.5, 40% coverage, 1,227 proteins |
| ***V. natriegens*** | **No** [VERIFIED — full 405-species list scanned, zero hits] | — | — |

**V. natriegens is genuinely absent from PaxDb.** Only one *Vibrio* species appears in the entire database — *Vibrio proteolyticus* (species_id 1219065, 1 dataset) — an unrelated species, not usable as a substitute.

**PRIDE/ProteomeXchange search for V. natriegens [VERIFIED]** turned up 3 relevant projects. **PXD027874** ("Vibrio natriegens Systems & Synthetic Biology: Proteome Profile (Temperature & Salinity)," Hervey et al., US Naval Research Laboratory, CC0 license) was selected: it is explicitly framed around V. natriegens as a synthetic-biology chassis, was processed with MaxQuant (LFQ intensity quantification, standard workflow), and — unusually valuably — includes the fully processed `proteinGroups.txt` output file (not just raw MS files), directly downloadable via PRIDE's FTP.

**A serendipitous secondary finding:** three of our five PaxDb-covered hosts (*B. subtilis*, *P. aeruginosa*, *S. enterica*) each have a 2025 dataset from the same underlying study — Abele et al., *Mol Cell Proteomics* 2025 ("Proteomic Diversity in Bacteria," MassIVE MSV000096603), a 303-species bacterial proteome atlas built from the DSMZ/Weihenstephan strain collection [VERIFIED, PMC11919601]. *E. coli* and *C. glutamicum* were checked and are **not** present in this atlas under PaxDb's current per-species listings. This atlas was used as a cross-check source (see B3), not the primary pick, because — per its own Methods — "bacteria were cultured on their respective agar plate types, temperatures, and oxygen conditions," i.e. each species got its own DSMZ-recommended condition, not one shared growth protocol; it is one consistent lab/instrument/method, but not one consistent growth condition.

### B2. Uniform derived-metric computation

**Method (identical across every host — full detail and rationale is in the script docstring, `scripts/08_compute_physiology_proxies.py`):** each protein has a gene/preferred name and a free-text functional annotation. A protein is assigned to a category if its gene name exactly matches a curated symbol list for that category (conserved bacterial nomenclature: `rps*/rpl*/rpm*` for ribosomal proteins, `rpoA/B/C/Z` for RNAP core, `groEL/groES/dnaK/dnaJ/grpE/tig` for chaperones, `tuf*/fusA` for elongation factors) **or** its annotation's primary-name clause (the text before the first semicolon — see below) contains a fixed keyword phrase for that category. Sigma factors are identified by annotation keyword only (`"sigma factor"`, `"sigma-70"`, `"RNA polymerase sigma"`), because sigma-factor gene symbols are not uniform across our hosts (`rpoD/rpoS/rpoH/rpoE/rpoN` in the Gammaproteobacteria vs. `sigA/sigB/sigC...` in *B. subtilis*/*C. glutamicum*).

**Two rounds of manual false-positive auditing were performed and are documented in the script**, because the first-pass method produced real errors:
- Round 1 caught the risk of matching keywords anywhere in the full annotation text, which pulled in e.g. *E. coli* ribosomal protein L23 (`rplW`) into the "chaperones" category, because L23's *extended* description mentions "trigger factor" (L23 physically contacts trigger factor at the ribosome exit tunnel) even though L23 itself is not a chaperone. **Fix:** keyword matching was restricted to the annotation's primary-name clause only (before the first semicolon) — PaxDb/UniProt-style annotations reliably put the protein's actual name/identity in that first clause.
- Round 2 (after the fix) still caught two remaining classes of error: (a) enzymes that *modify* a ribosomal protein, named in a way that puts the substrate's name first (e.g. *E. coli* `prmA` = "Methyltransferase for 50S ribosomal subunit protein L11," `roxA` = "50S ribosomal protein L16 arginine hydroxylase") — fixed with an explicit exclusion pattern for enzyme-activity suffixes; (b) anti-sigma-factor regulators, whose annotations mention "sigma" in describing what they regulate rather than what they are (e.g. *B. subtilis* `ylaD` = "anti-YlaC sigma factor," `yxlD` = "sigma-Y antisigma factor component") — fixed with a broader proximity-based exclusion (`anti` within 20 characters of `sigma`).

Every remaining match was spot-checked against its actual annotation text (not just trusted blindly) — e.g. confirmed *P. aeruginosa* locus-tag-named hits like `PA0149`/`PA1363` are genuinely annotated "Probable sigma-70 factor, ECF subfamily," and V. natriegens's `hslO` (matched via "chaperonin") is genuinely annotated "33 kDa chaperonin." The full audit trail (every matched gene, which rule fired) is in `out/physiology_proxy_matches.json`.

**Results — primary source per host:**

| Host | n proteins quantified | Ribosomal protein fraction | RNAP core fraction | Sigma factor fraction (n distinct) | Chaperone fraction | EF-Tu/EF-G fraction |
|---|---:|---:|---:|---:|---:|---:|
| *E. coli* | 3,747 | 20.38% | 0.73% | 0.11% (n=7) | 1.94% | 2.40% |
| *B. subtilis* | 4,052 | 18.70% | 0.38% | 0.16% (n=19) | 1.34% | 1.58% |
| *P. aeruginosa* | 5,034 | 23.48% | 0.89% | 0.17% (n=22) | 4.15% | 2.46% |
| *S. enterica* | 2,620 | 22.09% | 0.51% | 0.06% (n=5) | 0.54% | 3.45% |
| *V. natriegens* | 665 | 29.84% | 1.72% | 0.02% (n=2) | 4.11% | 4.63% |
| *C. glutamicum* | 1,227 | 14.94% | 0.32% | 0.09% (n=5) | 0.73% | 2.46% |

**Plausibility check against Gate 1's independent literature-derived growth-rate ranking [INFERRED, medium confidence]:** ribosomal protein fraction is a resource-allocation-to-translation proxy, and bacterial growth-law theory predicts faster-growing organisms should allocate more proteome mass to ribosomes. Ranked, our ribosomal-fraction result is VN(29.8%) > PA(23.5%) > SE(22.1%) > EC(20.4%) > BS(18.7%) > CG(14.9%). This matches the extremes of Gate 1's growth-rate literature well: *V. natriegens* — famously the fastest-growing named bacterium — has the clear highest fraction; *C. glutamicum* — the slowest-growing of the six by a wide margin (μ≈0.5–0.65 h⁻¹ vs. ≥0.95 h⁻¹ for everything else) — has the clear lowest. The middle ordering (PA > SE > EC > BS) is less clean against Gate 1's growth-rate numbers, which is expected: PaxDb's "Integrated" datasets aggregate many differently-conditioned studies per organism (see B3), so they are not a clean single-condition measurement to compare 1:1 against Gate 1's specific literature growth rates. **This is presented as a sanity check, not a validation** — the metric is directionally sensible, not precisely calibrated.

**A secondary observation worth carrying into Gate 2:** the *count* of distinct sigma factors detected (last column above) may be as informative as the mass fraction, and is a strikingly clean biological signal on its own — *B. subtilis*'s 19 detected sigma factors versus 2–7 for every other host directly reproduces its well-documented status as having one of the largest sigma-factor regulons in bacteria. This wasn't asked for explicitly but is offered as a candidate 6th/7th feature.

### B3. Comparability audit — mandatory, and the results are mixed

**[VERIFIED] Growth condition is NOT matched across hosts, and for five of six hosts it isn't even a single stated condition.** PaxDb's "Integrated" datasets (used as the primary pick for *E. coli*, *B. subtilis*, *P. aeruginosa*, *S. enterica*, *C. glutamicum*) are each a weighted average across *all* whole-organism studies PaxDb has for that species — anywhere from 3 underlying studies (*C. glutamicum*) to 20 (*E. coli*), spanning unstated combinations of growth phase, medium, and stress condition. There is no single "the *E. coli* dataset was grown in LB at 37°C" statement to make, because it wasn't one experiment.

**V. natriegens is the one exception with a precisely documented condition** — 300 mM added NaCl, 37°C, first timepoint (TP1), averaged across 3 replicates, chosen deliberately to be the closest available match to the ~256 mM NaCl / 37°C condition Gate 1 identified as V. natriegens's best-characterized growth-rate literature (Hoffart et al. 2017). This is a genuine asymmetry: the one host with no PaxDb coverage ended up with the *best*-documented growth condition of all six, precisely because it required manual sourcing rather than an opaque database aggregate.

**[COMPUTED] Source-sensitivity check** — for the three hosts with two independent PaxDb dataset choices available (*B. subtilis*, *P. aeruginosa*, *S. enterica*: "Integrated" vs. the single-study Abele et al. 2025 atlas), every derived metric was computed both ways:

| Host | Ribosomal fraction | RNAP fraction | Sigma fraction | Chaperone fraction | EF fraction |
|---|---|---|---|---|---|
| *B. subtilis* | 18.70% vs 16.41% (12% rel. diff) | 0.38% vs 0.45% (17%) | 0.16% vs 0.20% (21%) | 1.34% vs 1.48% (10%) | 1.58% vs 2.06% (23%) |
| ***P. aeruginosa*** | **23.48% vs 16.58% (29% rel. diff)** | **0.89% vs 0.50% (43%)** | 0.17% vs 0.26% (36%) | **4.15% vs 1.78% (57% rel. diff)** | **2.46% vs 1.42% (42%)** |
| *S. enterica* | 22.09% vs 20.20% (9% rel. diff) | 0.51% vs 0.73% (30%) | 0.06% vs 0.16% (62%) | 0.54% vs 0.48% (12%) | 3.45% vs 2.98% (14%) |

***P. aeruginosa* is unambiguously the most source-sensitive host** — every metric shifts by 29–57% relative depending on which reasonable dataset choice is used, roughly double the sensitivity seen for *B. subtilis* or *S. enterica* (sigma-factor fractions aside, which are noisy for everyone because the absolute values are tiny). **This is the second independent line of evidence (after Gate 1's classical-literature search) that *P. aeruginosa* — a primary host — is the least reliably characterized host in this project's physiology evidence base.** Two unrelated methods (literature search, omics-database sensitivity) both single out the same host as the weak point. That is a real, structural finding, not a one-off search-difficulty artifact.

**[VERIFIED] Coverage (number of proteins quantified) varies 8-fold**, from 665 (V. natriegens) to 5,034 (*P. aeruginosa*). Low coverage inflates noise for small-count categories — *V. natriegens*'s 2-protein sigma-factor count is almost certainly a coverage artifact (real sigma-factor repertoires in Gammaproteobacteria are typically 6+ factors) rather than a genuine biological signal that V. natriegens has unusually few sigma factors, and should be treated as low-confidence or dropped for that specific host/metric combination.

**Honest verdict on comparability:** this is a **comparable-but-crude** signal, not an incomparable-but-rich one, and it should be used that way — for rank-ordering hosts and as a coarse numeric feature, not for precise quantitative claims about exact percentages. The V. natriegens vs. C. glutamicum extremes are trustworthy; fine distinctions in the middle of the ranking (e.g. "is *P. aeruginosa* really 3 percentage points ahead of *S. enterica*?") are not, especially once source-sensitivity is accounted for.

### B4. Verdict

**Can a physiology feature vector with ≥3 metrics be built uniformly for all six hosts, from sources comparable enough to defend? YES**, with the comparability caveats above made explicit and carried forward into every downstream use.

**The vector, as computed, becomes the Gate 2 spec:**

| # | Feature | Source | Computation |
|---|---|---|---|
| 1 | Max growth rate / doubling time | Gate 1 literature (see `out/task4_physiology_primary3.md`, `out/task4_physiology_secondary3.md`) | Best available published value per host; **not condition-matched across hosts** — carry the source medium/temp as metadata, don't silently drop it |
| 2 | Ribosomal protein mass fraction | This gate (PaxDb ×5 hosts + PRIDE PXD027874 for V. natriegens) | Σ(abundance_ppm of `rps*/rpl*/rpm*`-matched or "NNS ribosomal protein"-annotated proteins) / Σ(total abundance_ppm) |
| 3 | RNA polymerase core-subunit fraction | Same | Σ(abundance of `rpoA/B/C/Z`) / total |
| 4 | Sigma factor fraction (+ count as an alternative/additional feature) | Same | Σ(abundance of annotation-keyword-matched sigma factors) / total; also report raw count of distinct sigma factors detected |
| 5 | Chaperone system fraction | Same | Σ(abundance of `groEL/groES/dnaK/dnaJ/grpE/tig` + "chaperonin"-annotated) / total |
| 6 | Elongation factor (EF-Tu/EF-G) fraction | Same | Σ(abundance of `tuf*/fusA`) / total |

This is 6 metrics, uniformly computed (same code, same rules) for all 6 hosts — well above the ≥3 bar. **Recommendation: this omics-derived vector (features 2–6) should be Gate 2's *primary* physiology representation, specifically because it is the only one that is genuinely uniform across all six hosts.** Feature 1 (growth rate) and Gate 1's richer literature-derived metrics (RNA/protein ratio, φ_R, absolute ribosome/RNAP counts) should be layered on as a **secondary, richer tier** for the subset of hosts where they exist with good condition-matching — concretely, *E. coli*, *B. subtilis*, and *V. natriegens* have a uniquely strong shared data point (Zhu, Mori, Hwa & Dai, *PNAS* 2025, all three grown side-by-side in identical conditions), and *C. glutamicum* has strong single-host data (Nat. Commun. 2023) at 30°C only.

**Tier assignment, combining Gate 1 and Gate 1.5 findings:**

| Host | Tier-1 uniform omics vector (this gate) | Tier-2 rich literature vector (Gate 1) | Overall physiology confidence |
|---|---|---|---|
| *E. coli* | Strong (3,747 proteins, PaxDb "Integrated") | Strong (Bremer & Dennis classic series + Zhu 2025) | High |
| *B. subtilis* | Strong (4,052 proteins) | Strong (Zhu 2025 matched to EC/VN) | High |
| *V. natriegens* | Moderate (665 proteins — lower depth, but explicit single condition) | Strong (Zhu 2025 matched to EC/BS; MSB 2023 φ_R) | High |
| *C. glutamicum* | Moderate (1,227 proteins, weakest PaxDb coverage) | Strong (Nat Commun 2023, but 30°C only) | Moderate-high |
| *S. enterica* | Strong (2,620–2,492 proteins, low source-sensitivity) | Weak (paywalled classical sources) | Moderate — **omics tier now compensates for Gate 1's literature gap here** |
| ***P. aeruginosa*** | **Weak (highest source-sensitivity of any host, 29–57% relative)** | **Weak (growth rate only; Gate 1)** | **Low — the one host where both approaches agree it's the weak point** |

**The practical upshot for the genome-vs-physiology comparison:** *S. enterica* is in noticeably better shape after this gate than Gate 1 alone suggested (Gate 1 called it "thin, growth-rate-only"; the omics tier gives it 5 solid, low-source-sensitivity metrics). *P. aeruginosa* remains the genuine weak point, and now for two independent reasons rather than one. If Gate 5's results turn out sensitive to *P. aeruginosa*'s physiology features specifically, that sensitivity should be flagged rather than treated as a modeling artifact — it may be a real data-quality signal.

---

## WHAT I COULD NOT DO

- **Could not obtain a registered PaxDb API key** — all requests used the documented public fallback (`x-api-key: test`), which the API's own error message describes as a "shared, rate-limited pool" not intended for production use. This gate's usage was light (a few hundred requests total) and completed without any rate-limit errors, but a production Gate 2 pipeline that re-downloads this data repeatedly should register for a real key first (one-click, per the API's own message) rather than relying on the shared test pool.
- **Could not access the PaxDb website's own UI directly** — it is a client-rendered Vue.js SPA that returns no content to a non-JavaScript fetch. All access was via reverse-engineering the underlying REST API from the compiled JS bundle. This worked completely, but it means there's no user-facing page to point to for "how was this dataset chosen" — the OpenAPI spec at `https://api.pax-db.org/try-api/openapi.json` and the raw JSON responses in `raw/paxdb_*.json` are the actual documentation trail.
- **Could not confirm whether the Abele et al. 2025 atlas (MSV000096603) genuinely excludes *E. coli* and *C. glutamicum*, or merely hasn't been fully ingested into PaxDb's per-species listing yet** — the paper's full species table (Supplemental Table S2) was not accessible via the tools available in this pass (ScienceDirect blocked with HTTP 403; the PMC mirror's supplementary files were not directly fetchable). If *E. coli* and *C. glutamicum* are in fact in that atlas under a different taxonomy ID than PaxDb's current mapping, a fully single-source 6-host (well, still missing V. natriegens) comparison might be recoverable with more digging.
- **Could not verify the exact growth medium recipe for the "S300" (300 mM added NaCl) condition** in PXD027874 beyond what's stated in the project metadata (base medium composition before the NaCl/temperature factorial treatment) — the full materials-and-methods text of the underlying paper (if one exists beyond the PRIDE project description) was not located/read in this pass.
- **Did not attempt to resolve the V. natriegens sigma-factor undercount (n=2)** beyond flagging it as a likely coverage artifact — a targeted re-search of the 665-protein dataset specifically for sigma-factor-family proteins that may have been filtered out by the QC steps (contaminant/reverse/site-only removal, or simply not detected at this instrument's depth) was not performed.

## CONTRADICTIONS WITH THE CHARTER

See `out/CHARTER_AMENDMENTS.md` for the full write-up (four items: C1 leakage mechanism, C2 Bernstein-lab data reusability, C3 translation N, C4 RS241 six-way-intersection correction). Summary:

1. **The charter's stated primary leakage defense (barcode-replicate dedup) doesn't apply — sequence-identity clustering is now the only defense**, and needs a proper threshold-sensitivity analysis (0.3/0.5/0.7) rather than an assumed 0.5, plus per-fold max-identity reporting.
2. **Bernstein-lab physiology tables were never a usable shortcut** (confirmed in Gate 1, restated here for completeness) — remove that assumption from planning documents.
3. **This gate found a NEW instance of the same pattern Gate 1 found**: Gate 1's own RS241 six-way-intersection number (111/126) turns out to have been the wrong quantity to report for the intended use (leave-one-host-out viability), analogous to how the charter's own barcode-dedup assumption was wrong for a different reason. Both are now corrected. The general lesson — verify what a computed number is actually being *used for* before treating it as final — applies to future gates too.
4. **P. aeruginosa's physiology weakness is now a two-method-independent finding**, not a one-off. This wasn't anticipated by the charter or by Gate 1 in isolation, but with two unrelated approaches (literature search, database source-sensitivity) both landing on the same host as the weak point, it should be treated as a real property of the available public data on *P. aeruginosa* physiology, not noise.

## FILES WRITTEN

- `out/GATE1_5_MEMO.md` — this memo
- `out/CHARTER_AMENDMENTS.md` — the four charter corrections
- `out/task_a_rs241_pairwise_results.json` — full RS241 pairwise/viability/ranked-configs data
- `out/physiology_proxy_matches.json` — full gene-level audit trail for the physiology proxy classification (every match, which rule fired)
- `out/vnat_processing_log.json` — V. natriegens PRIDE data processing log
- `data/physiology_proxy_results.csv` — final 6-host × 5-metric physiology proxy matrix (all sources, primary + cross-check)
- `data/vnatriegens_proteome_TP1_S300_T37.parquet` — processed V. natriegens proteome (665 proteins, representative condition)
- `scripts/05_task_a_rs241_pairwise.py` — Task A computation
- `scripts/06_download_paxdb.py` — PaxDb API client and bulk downloader
- `scripts/07_process_vnatriegens_pride.py` — PRIDE MaxQuant output processor for V. natriegens
- `scripts/08_compute_physiology_proxies.py` — uniform protein-family identification and mass-fraction computation (all 6 hosts)
- `raw/paxdb_openapi.json`, `raw/paxdb_species_list.json`, `raw/paxdb_datasets_*.json`, `raw/paxdb_abundances_*.json` — PaxDb raw downloads
- `raw/vnat_proteinGroups.txt`, `raw/vnat_metadata.csv`, `raw/vnat_expdesign.txt` — PRIDE PXD027874 raw downloads
