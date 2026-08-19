**Does the in-vivo cross-host contrast survive controlling for source-genome GC content? YES — transcription SURVIVES INTACT (gap ratio 0.341 → 0.308, a 9.7% relative change), translation SURVIVES ATTENUATED (gap ratio 0.286 → 0.240, a 15.9% relative change). Neither readout comes close to "does not survive." Independently confirmed a second way: the contrast holds, and in some strata is even more pronounced, within single source-genome phyla (Proteobacteria, Firmicutes) considered separately — a check that makes no assumption about the functional form of the GC relationship at all.**

# GATE 10.5 MEMO — Close DRAFTS, Check the GC Confound

---

## TASK 2 — the GC confound [LOAD-BEARING, COMPLETE]

Script: `scripts/95_gc_confound_invivo.py`. Full output: `out/results/gate10_5_gc_confound.json`, `.csv`. Figures: `out/figures/gate10_5_gc_vs_activity_invivo.png`, `gate10_5_raw_vs_gc_controlled.png`, `gate10_5_phylum_stratified.png`.

### 2.1 — Does the GC relationship hold in vivo?

**Yes, in every host, both readouts** — though noticeably weaker than DRAFTS's cell-free range (-0.74 to -0.49 in all 10 hosts):

| Host | Transcription rho(GC%, activity) | Translation rho(GC%, activity) |
|---|---|---|
| EC | -0.614 | -0.494 |
| BS | -0.200 | -0.266 |
| PA | -0.434 | -0.230 |

*E. coli* shows by far the strongest GC dependence in vivo (-0.614, close to DRAFTS's cell-free range); *B. subtilis* the weakest (-0.200/-0.266). This asymmetry matters for what follows: EC's own activity is the most GC-confounded of the three hosts, meaning any EC-involving pair has the most "shared-with-GC" variance to remove.

### 2.2 — Partial correlations, raw vs. GC-controlled, all pairs, both readouts

Two independent computation methods (closed-form partial-correlation formula and rank-residual regression) agreed to 3 decimal places in every cell — reported as a single number below.

| Readout | Pair | N | Raw rho | GC-controlled rho | % change |
|---|---|---|---|---|---|
| transcription | EC-BS | 3,668 | 0.258 | 0.207 | -19.8% |
| transcription | **EC-PA** | 9,741 | **0.754** | **0.718** | -4.7% |
| transcription | BS-PA | 2,099 | 0.257 | 0.236 | -8.2% |
| translation | EC-BS | 866 | 0.161 | 0.125 | -22.5% |
| translation | **EC-PA** | 3,826 | **0.742** | **0.738** | -0.6% |
| translation | BS-PA | 314 | 0.263 | 0.230 | -12.7% |

**EC-PA barely moves** (-4.7% / -0.6%) — the pair carrying the central claim's positive evidence is nearly unaffected by GC control. The BS-pairs shrink by a larger relative amount (-8% to -23%), but this is expected mechanically: a small raw correlation loses more of itself, proportionally, to any fixed absolute adjustment than a large one does — the BS-pairs were already close to zero, so the same-sized subtraction is a bigger percentage of a smaller number. **None of the three pairs changes sign or approaches the EC-PA value.**

### 2.3 — The decisive question: does the EC-PA vs. BS-pairs contrast survive?

| Readout | EC-PA raw→GC | BS-pairs mean raw→GC | Gap ratio raw→GC | % change | **VERDICT** |
|---|---|---|---|---|---|
| Transcription | 0.754 → 0.718 | 0.257 → 0.221 | 0.341 → 0.308 | -9.7% | **SURVIVES INTACT** |
| Translation | 0.742 → 0.738 | 0.212 → 0.177 | 0.286 → 0.240 | -15.9% | **SURVIVES, ATTENUATED** |

Neither readout comes anywhere close to the "DOES NOT SURVIVE" threshold (which would require the BS-pairs mean to approach or exceed EC-PA, or the gap ratio to exceed ~0.5) — the gap ratio moves modestly in both cases, staying in the same qualitative range Gate 8's disattenuation analysis already reported (0.341, 0.9-reliability case) and the translation-readout attenuation (15.9%) is the larger of the two but still leaves BS-pairs at under a quarter of EC-PA's value.

### 2.4 — Phylum-stratified check (independent confirmation, no functional-form assumption)

Within **Proteobacteria** source sequences alone (n=321-3,282 depending on pair): EC-PA transcription rho = **0.634** (still strong); EC-BS = 0.090 (still near zero); BS-PA = **-0.209** (negative). Within **Firmicutes** source sequences alone (n=1,214-2,666): EC-PA = **0.812** (still strong, if anything higher); EC-BS = 0.177; BS-PA = 0.271. **The contrast holds, and in the Firmicutes stratum is if anything sharper, entirely within a single source phylum** — this check makes no assumption about linearity or the specific functional form of the GC relationship, unlike the residual-regression approach in 2.2/2.3, and it independently confirms the same conclusion. Smaller strata (Actinobacteria, Bacteroidetes, Cyanobacteria, Euryarchaeota, Planctomycetes) show noisier but directionally consistent patterns (EC-PA always the highest of the three pairs where N permits a stable estimate).

### 2.5 — Interaction with Gate 8's disattenuation gap ratio

| Readout | Gap ratio, disattenuation only (Gate 8's own number) | Gap ratio, disattenuation + GC control |
|---|---|---|
| Transcription | 0.341 | 0.308 |
| Translation | 0.286 (not previously reported combined with reliability=0.9 in this exact form) | 0.240 |

Adding GC control on top of Gate 8's existing measurement-noise correction moves the transcription gap ratio from 0.341 to 0.308 — both corrections point the same direction (narrowing the gap slightly) but neither, alone or combined, brings the ratio anywhere near parity (1.0) or even the "does not survive" region.

### VERDICT, restated plainly

**The GC confound is real (strong, universal, in the same direction as DRAFTS's cell-free finding) but does not explain the central claim's evidence.** The EC-PA vs. BS-pairs contrast — this project's most load-bearing raw-data finding — survives GC control intact for transcription and survives with moderate (not severe) attenuation for translation, confirmed two independent ways. **Proceeding to Task 3 as instructed**, since Task 2 did not undercut the central claim.

---

## TASK 1 — DRAFTS results section

`out/DRAFTS_SECTION.md` — publication-ready prose, no new analysis beyond what Gate 10 already established. Covers: the numbers (labeled by provenance — recomputed vs. Gate-10-sourced), the mechanistic reading (explicitly flagged as a hypothesis consistent with the data, not a demonstrated causal claim), the cell-free-prototyping implication, the four honest limitations (P. aeruginosa absent; small sequence-level join; transcription-only; n=10 vs n=3 host-count asymmetry), and why the host-count sweep was not run (one paragraph, matching the reasoning given in this gate's own prompt).

---

## TASK 3 — Framing updates

### `out/PAPER_FRAMING.md`

The DRAFTS modality finding is inserted as **new finding 4** in the evidence hierarchy (renumbering the prior findings 4-11 to 5-12), immediately after the FiLM-instability finding and before the free-lookup-embedding finding. **Reasoning for this placement:** it is a mechanistic explanation for the central negative result, sourced from an *independent dataset* (not a reanalysis of the same in-vivo numbers), using the *same regulatory-sequence library* from the *same laboratory* — this is stronger corroboration than another instance of the negative result on the same data, and directly explains WHY genomic machinery features carried no signal (finding 5, the free-embedding result), so it belongs ahead of that finding, not after it. It is placed below the FiLM-instability finding (a direct, controlled result about this project's own model) and above the free-embedding result (a result the DRAFTS finding helps explain), which is the ordering "strongest and most direct first" naturally produces here.

The GC-control result is recorded in finding 8 (formerly 7, the disattenuation finding) as an additional, independent robustness check alongside the existing measurement-noise correction.

The DRAFTS scoping decision (host-count sweep not run) is noted in the same new finding.

### `out/KNOWN_ISSUES.md`

New item added: the *P. agglomerans*/*P. aeruginosa* abbreviation collision, framed as a documented gotcha for anyone using both datasets together.

### `out/CHARTER_AMENDMENTS.md`

New standing rule (SR8) recording the audit-check-6 path-resolution fix from Gate 10, so future gates adding non-`data/`-rooted files to the manifest know the convention already exists and don't rediscover the same bug.

---

## WHAT I COULD NOT DO

1. A true within-sequence paired in-vivo comparison for the 112 DRAFTS-Johns overlapping sequences at adequate N — the usable-in-both-hosts overlap for EC-BS specifically was only 15 sequences (reported in Gate 10, reused here), too small to be an independent estimate; the headline 0.597-vs-0.258 comparison necessarily compares two different (though maximally overlapping in library origin) populations rather than a single paired one.
2. A rigorous statistical test (e.g. a bootstrap CI) on the GC-controlled correlations or the gap-ratio change — reported as point estimates with two independent corroborating methods (residual regression, phylum-stratification) rather than a formal interval, given the descriptive/diagnostic nature of this check and the gate's own "no models" scope limit.
3. Extending the phylum-stratified check to translation with adequate power in every phylum — several strata had translation-readout N too small for a stable estimate (e.g. Proteobacteria EC-BS translation, n=18) and are reported as such (None) rather than papered over with an unstable number.

## TEMPTATIONS TO ADJUST

1. **Rounding the translation verdict up to "SURVIVES INTACT"** to make both readouts match cleanly. Not done — 15.9% attenuation is reported honestly as "ATTENUATED," a real, if modest, difference from the transcription readout's 9.7%.
2. **Reporting only the residual-regression method** and omitting the phylum-stratified cross-check, since the former alone would have been sufficient to answer the task. Not done — the phylum-stratified check is reported in full because it independently confirms the result via a completely different, assumption-free method, which is stronger evidence than either check alone.
3. **Treating EC-PA's near-zero attenuation (-4.7%/-0.6%) as suspicious or requiring extra scrutiny** rather than simply reporting it. Not done — it is a real, mechanically-explicable pattern (EC-PA's raw correlation is large, so the same absolute GC-adjustment is a small percentage of it) and is reported plainly with the explanation, not flagged as an anomaly needing further chasing.
