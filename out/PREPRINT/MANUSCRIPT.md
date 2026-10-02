---
header-includes: |
  \usepackage{float}
  \floatplacement{figure}{H}
  \usepackage{caption}
  \captionsetup{labelformat=empty,font=small,justification=raggedright,singlelinecheck=false}
---

# Cell-free transcription does not reproduce the bacterial chassis effect: evidence that host specificity resides in cellular context rather than transcription machinery

**Gabriel Kim**

Stanford University, Stanford, CA, USA\
Correspondence: gabeykim@stanford.edu

---

## Abstract

A regulatory DNA sequence characterized in one bacterial species often behaves differently in another — the chassis effect. We compare cell-free (DRAFTS, Yim et al. 2019) against in-vivo (Johns et al. 2018) measurements of the same 165 bp library. Pooled, the two modalities are indistinguishable (**0.616 cell-free, 0.655 in vivo**). Restricted to sequences detectably active in both hosts, the in-vivo correlation falls to **0.258** while the cell-free correlation rises to **0.677**, a 2.6-fold gap that survives partial-correlation control for source-genome GC composition (2.63× raw, 2.69× controlled). In living cells, cross-host agreement is carried substantially by agreement on which sequences are silent: **26.0% of these sequences are co-active in vivo against 93.6% in lysate**. Screening constructs for one fixed chassis via lysate is therefore reasonable; ranking a part across candidate hosts is not. **This comparison rests on a single species pair**, the only one present in both datasets. We interpret it as evidence, not proof, that host specificity resides substantially in cellular context rather than in the transcriptional machinery lysates retain. That interpretation predicts host-descriptor features encoding machinery composition carry no signal, and we tested it: a pre-registered hypothesis failed on all 8 primary comparisons; an arbitrary host-identity tag matched or beat both a 37-feature genomic and a 6-feature physiology vector in all 12 cells; and sequence-only prediction was not distinguishably beaten in 16 of 18 tests across three conditioning mechanisms. **n_hosts ≤ 6** throughout. We release the benchmark, splits, evaluation code, and baselines, including the negative result and a retraction.

---

## 1. Introduction

Synthetic biology relies on characterized genetic parts, most of them characterized in one or a handful of model organisms — overwhelmingly *E. coli* — with behavior in other hosts assumed rather than measured. The literature calls the resulting mismatch the **chassis effect**, and describes it as a barrier that hinders prediction of function from genetic composition alone and drives costly repetitions of the design-build-test cycle (Chan, Baldwin & Bernstein, 2023); identifying the biological properties that underpin it remains an open challenge (Chan & Bernstein, 2024).

Part of the field routes around the question. Orthogonal-machinery approaches — T7 RNA-polymerase-based expression modules and similar systems — partially neutralize host variation for users willing to adopt an orthogonal transcription/translation system. That helps users who can redesign a circuit around orthogonal machinery, not users who need to predict how an existing, non-orthogonal part will behave in a new host.

One direct, falsifiable hypothesis about the chassis effect's cause comes from work on engineered circuits across related *Pseudomonas* and *Halopseudomonas* hosts and across six closely related *Stutzerimonas* strains. That work reports that "hosts exhibiting more similar metrics of growth and molecular physiology also exhibit more similar performance of the genetic inverter" (Chan, Baldwin & Bernstein, 2023) — that host physiology, measured directly, is what explains chassis effects.

We test a sharper version of that hypothesis using a comparison across measurement modalities. Johns et al. (2018) measured transcription and translation activity for tens of thousands of regulatory sequences mined from 184 bacterial genomes across three primary recipient hosts and, for a 241-sequence subset, three more — to our knowledge the only public dataset measuring the same regulatory sequences across multiple bacterial hosts at scale. DRAFTS (Yim et al., 2019), from the same laboratory, measured cell-free (TXTL) transcription of the same library across ten species.

Comparing the two (Section 3.2) motivates a testable prediction: if host specificity is substantially contextual rather than mechanistic, host-descriptor features built from genome- and proteome-encoded machinery composition should carry no signal for cross-host prediction. We pre-registered a test of this before training any model. It held on every primary comparison. We then made four independent attempts to find a case where it did not. Two initially appeared to be counterexamples. One survives narrowly, in a single cell; the other did not survive its control (Section 3.8).

**n_hosts ≤ 6.** Every claim here is bounded by that number; the unit of generalization is the host, not the sequence.

Five major benchmark suites — BEND (Marin et al., 2024), Genomic Benchmarks (Grešová et al., 2023), DART-Eval (Patel et al., 2024), DNALONGBENCH (Cheng et al., 2025), and the Nucleotide Transformer suite (Dalla-Torre et al., 2025) — define tasks on eukaryotic genomes only; we found no bacterial cross-host regulatory activity task among them. We built CROSSHOST, a frozen versioned benchmark on the Johns data with a small demonstration model, to carry out this test and leave a reusable artifact.

---

## 2. Methods

### 2.1 Data

**Johns et al. (2018)**, *Nature Methods* 15:323–329 (BioProject PRJNA431139): 29,249 regulatory-sequence (RS) records, 165 bp each, mined from 184 bacterial genomes and measured in *E. coli*, *B. subtilis*, and *P. aeruginosa* for transcription (RNA-seq-derived) and translation (FACS-seq GFP). Per-host usable transcription counts: EC 24,613, BS 15,848, PA 21,473. A 241-sequence extension panel (RS241) adds *S. enterica*, *V. natriegens*, and *C. glutamicum*; 34 RS241 IDs have no recoverable sequence text, leaving 207 usable.

Active fractions depend on the denominator, and the two in use are not interchangeable. Over each host's own usable population — the denominator we use, since baselines and folds are host-specific — they are EC 61.1%, BS 28.0%, PA 83.3%. Over the three-host union set (n = 11,265–11,276), the denominator Johns et al. use, they are EC 52.0%, BS 18.9%, PA 83.8%; our recomputation on that population gives 17.9% for BS, a 1.0-point difference we have not traced. The ordering is identical either way; the percentages are not comparable across denominators.

**Translation floor artifact.** `protein_log10` is pinned at a per-host floor value for 66.5% (EC), 89.9% (BS), and 10.1% (PA) of nominally usable rows, consistent with a below-detection convention rather than measurement. We exclude floor-pinned rows from the regression target; corrected usable N is 9,146 / 1,101 / 17,630.

**DRAFTS** (Yim, Johns et al., 2019, *Molecular Systems Biology* 15:e8875), from the same laboratory, measured cell-free transcription of the same library across ten species (1,047 sequences at full coverage). Files came from the PMC Open Access bulk package and were hash-verified. DRAFTS is used descriptively — never for model training, never in contact with the frozen splits.

### 2.2 Three dataset collisions, and the restriction regimes

The two datasets share a laboratory and a source library but not their conventions, and three collisions between them change results materially.

**"Usable" means different things.** DRAFTS marks template-present, zero-RNA rows as `no_RNA_counts` and excludes them from its usable set; Johns retains the equivalent rows with `tx_norm` exactly 0. **DRAFTS-usable is therefore the analogue of Johns-active**, and a comparison built on the shared column name `n_shared_usable` would silently compare a co-active cell-free correlation against a pooled in-vivo one. We accordingly report every correlation under one of two labeled regimes:

- **Pooled** — all sequences measured in both hosts of a pair, including those measured as inactive. For DRAFTS this requires assigning activity 0 to `no_RNA_counts` rows, which carry a sentinel rather than a numeric zero; that is an imputation by definition, not a measurement.
- **Co-active** — sequences detectably active in both hosts.

Section 3.1's primary figures are co-active. A full regime audit of every reported correlation ships with the benchmark.

**The ID spaces are disjoint.** DRAFTS Oligo IDs run 13097–14477, ours 14478–43726, despite both libraries deriving from the same 184-genome mining effort. Joins must use sequence text, which recovers 112 of DRAFTS's 1,047 sequences as also present in the three-host library.

Separately, DRAFTS's two-letter code `Pa` denotes *Pantoea agglomerans*, not *Pseudomonas aeruginosa* — which is absent from DRAFTS entirely. A script matching on the bare code will silently merge two unrelated organisms.

### 2.3 Host features

**Genomic (37-D):** sigma-factor complement, anti-Shine-Dalgarno hybridization free energy, tRNA copy number and tAI, codon usage, GC content, and chaperone/heme/RNAP subunit counts. Each feature was unit-tested against a published reference value before use — for example, the *E. coli* anti-Shine-Dalgarno sequence and published *E. coli* tAI ranks.

**Physiology (6-D):** ribosomal-protein, RNAP core, sigma-factor, chaperone, and elongation-factor proteome fractions plus growth rate, from public reference proteomes (primarily PaxDb), depth-matched to correct a 7-fold variation in native coverage depth. This is a coarse, database-derived composition proxy, not a measurement of dynamic physiological state — closer in kind to the genomic vector's machinery annotations than to the membrane potential, supercoiling, resource competition, or growth-phase transitions implicated by the mechanistic hypothesis in Section 3.4.

Both vectors are z-scored across all six hosts rather than fold-locally, since host covariate information is known at deployment time for a new host; sequence-level label statistics are always fit fold-locally.

### 2.4 Architecture

A four-layer 1D convolutional trunk (kernel widths 15/9/5/3, channels 128/128/64/64, matching the 165 bp part length), BatchNorm, a residual connection on the final conv layer, global average pooling, and two two-stage output heads (transcription, translation; each a classifier plus a strength regressor on actives) — 213,956 parameters with no conditioning pathway (`SequenceOnlyCNN`, the primary baseline).

The pre-registered conditioning mechanism is FiLM: a small MLP maps the host vector to per-channel (γ, β) applied multiplicatively at conv layers 2 and 3. Three alternatives were added later to test whether FiLM's own weakness rather than conditioning per se explained the negative result: concatenation, per-host output heads combined by prediction-space averaging, and per-host heads using the nearest training host. Those three across six (host, readout) cells give the 18 comparisons in Section 3.5; FiLM is the comparator, not one of the three.

### 2.5 Splitting

Sequences are deduplicated to unique text before splitting. The primary leakage defense is genome-blocked sequence-homology clustering (MMseqs2, `--min-seq-id 0.5 -c 0.8`) with an exact- and near-duplicate safety net — a k-mer all-pairs check force-merging any connected component spanning folds — added after clustering alone was found to miss byte-identical pairs at every threshold tested.

**Maximum train–test sequence identity across all five folds is 0.8485.** We report this as a measured number rather than as a claim that clustering was applied, since approximate clustering carries no completeness guarantee and the 184 source genomes include phylogenetically close pairs.

RS241 is never trained on, enforced by an assertion in a standing leakage audit: 0 of 241 RS241 IDs appear in any main-library fold. Within-host quantile normalization is fit on training folds only, inside the fold loop, verified by the same audit.

### 2.6 Pre-registration and decision rule

`PREREGISTRATION.md`, committed 2026-08-04, before any host-conditioned training. The commit timestamp is independently verifiable in the public repository's git history.

**H-MAIN** (the kill gate, dual-arm per readout): the cross-host model at N=100 host-specific calibration examples matches or beats a per-host-only baseline at N=3,000 (transcription) or its measured saturation ceiling, N≈300 (translation), on *B. subtilis* held out.

**H-SCIENCE:** genomic and physiology host features differ measurably in leave-one-host-out performance; direction not predicted, either outcome a result.

**H-DIAGNOSTIC:** the host-biology representation outperforms a free per-host lookup embedding under a pre-specified fallback (mean of training-host embeddings).

One amendment, dated 2026-08-05 and made before H-MAIN was evaluated, widened the decision interval from an 80% t-interval to a 90% percentile bootstrap (10,000 resamples, hierarchical: fold identities resampled, then within-fold draws), after fold-to-fold variance was found to exceed within-fold draw variance by a factor inconsistent with the t-interval's assumptions. The change makes MET harder to declare, and was not revisited after results were seen.

**Decision rule.** H-MAIN is MET only if the model's 90% CI lower bound exceeds the baseline's 90% CI upper bound. Overlapping intervals are NOT MET regardless of which mean is higher.

### 2.7 AI assistance

Data processing, model implementation, evaluation code, audit scripts, and manuscript drafting were carried out with the assistance of an AI system (Claude, Anthropic), operating under task specifications and editorial direction from the author. The author designed the study, specified the pre-registration, made all interpretive decisions, verified every reported number against committed outputs, and is responsible for all claims.

---

## 3. Results

### 3.1 Cross-host measurement structure (Figures 1–2)

Correlations are reported under the two regimes defined in Section 2.2, and the choice changes the answer.

**Co-active:** EC–PA transcription ρ = 0.754 (n = 9,741); EC–BS 0.258 (n = 3,668); BS–PA 0.257 (n = 2,099). Translation shows the same ordering (0.742 / 0.161 / 0.263).

**Pooled:** EC–PA 0.621, EC–BS 0.655, BS–PA 0.508 for transcription — a different ordering, since pooled correlation is inflated by agreement on which sequences are inactive.

We use co-active as primary because it isolates strength from silence. Section 3.2 shows that the difference between the two is where the host effect lives.

The qualitative co-active pattern is Johns et al.'s own finding — their Supplementary Figs. S13 and S15 report it, and their main text states that "between recipients, only *E. coli* and *P. aeruginosa* showed significant correlations." Our contribution is a larger-scale quantification, explicit regime labeling, and two robustness checks not previously applied to this comparison.

We summarize the contrast as a **gap ratio**: the mean of the two *B. subtilis*-pair correlations divided by the EC–PA correlation, within a readout. Under the co-active regime the raw transcription gap ratio is 0.341. Gap ratios and percentage changes are computed from unrounded correlation values; figure labels are rounded to two decimals, so small changes may not be visible in the figures.

**The gap survives measurement noise.** Anchored by a measured *E. coli* transcription reliability of 0.912 (five growth-condition replicates, cross-checked to 0.929 by an independent method — the only host×readout combination with usable replicate structure), disattenuation leaves the gap ratio unchanged at 0.341 when both sides are corrected at reliability 0.9, and moves it to 0.515 only at the extreme of the sensitivity grid, where *B. subtilis* is assigned a pessimistic 0.5 while the others retain 0.9. Growth-condition replicates capture biological variation alongside measurement error, so this reliability estimate is conservative.

![**Figure 1 — Raw cross-host measurement correlation, with disattenuation correction.** Spearman ρ for all three pairs, both readouts, co-active regime, raw and disattenuation-corrected. Co-active N per pair: transcription EC–BS 3,668, EC–PA 9,741, BS–PA 2,099; translation 866, 3,826, 314. No resampled interval — these are population correlations, and the correction is a deterministic formula shown at three reliability levels (0.5/0.7/0.9) as a sensitivity grid. The two EC–PA bars at reliability 0.5 and 0.7 are clipped at ρ = 1 (the analysis stores min(corrected, 1.0)); their uncapped values — 1.51 and 1.08 transcription, 1.48 and 1.06 translation — are printed above the ceiling line. A disattenuated value above 1 means the assumed reliability is too low for that pair.](figures/Figure1_disattenuation.png)

**The gap survives a source-composition confound.** DRAFTS revealed a universal source-genome GC–activity relationship in cell-free data (ρ −0.49 to −0.74 across ten hosts); the same relationship holds more weakly in vivo (ρ −0.20 to −0.61). Partial-correlation control with GC held constant — confirmed to three decimals by a closed-form formula and by rank-residual regression independently — leaves the transcription gap ratio at 0.308 (from 0.341, a 9.7% change) and translation at 0.240 (from 0.286, 15.9%). EC–PA itself moves only −4.7% (0.754 → 0.718) for transcription and under 1% for translation; the *B. subtilis* pairs shrink more in relative terms from a smaller base (EC–BS 0.258 → 0.207), which is mechanically expected. Phylum stratification agrees independently. In all five source-genome phyla where every pair clears the n > 20 reporting threshold, EC–PA transcription stays above both *B. subtilis* pairs: Proteobacteria 0.634 (n = 3,282) against EC–BS 0.090 (n = 321) and BS–PA -0.209 (n = 327); Firmicutes 0.812 (n = 2,666) against EC–BS 0.177 (n = 2,532) and BS–PA 0.271 (n = 1,214); Actinobacteria 0.693 (n = 824) against EC–BS 0.092 (n = 97) and BS–PA -0.174 (n = 119); Bacteroidetes 0.536 (n = 842) against EC–BS -0.010 (n = 146) and BS–PA 0.018 (n = 81); Cyanobacteria 0.543 (n = 455) against EC–BS 0.332 (n = 62) and BS–PA -0.022 (n = 44). The *B. subtilis* pairs range -0.21 to 0.33 across strata rather than sitting near zero in every one, and the five EC–PA values (0.536–0.812) bracket the unstratified co-active 0.754, so the stratification reveals no anomaly. The margin is narrowest in Cyanobacteria, the smallest stratum, where EC–BS reaches 0.332 on n = 62.

![**Figure 2a — GC-composition control.** Raw vs. GC-controlled Spearman ρ, all three pairs and both readouts, partial-correlation control confirmed independently by closed-form formula and rank-residual regression.](figures/Figure2a_gc_control.png)

![**Figure 2b — Phylum stratification.** The same contrast computed within each source-genome phylum separately, as an assumption-free check on the GC control. All five strata in which every host pair clears the n > 20 threshold are shown — Proteobacteria, Firmicutes, Actinobacteria, Bacteroidetes and Cyanobacteria — with per-stratum N on every bar.](figures/Figure2b_phylum_stratified.png)

### 3.2 The modality contrast (Figure 3)

DRAFTS characterized transcription from the same library using cell-free lysates, which retain RNA polymerase, sigma factors, and ribonucleotides but lack an intact membrane, supercoiling homeostasis, macromolecular resource competition, and growth-phase physiology.

**Pooled, the two modalities agree.** For *E. coli*–*B. subtilis*, pooled transcription correlation is **0.616 cell-free** (n = 862) versus **0.655 in vivo** (n = 14,088).

**Conditioning on co-activity separates them.** Restricted to sequences active in both hosts, the in-vivo correlation falls to **0.258** (n = 3,668), a 61% drop, while the cell-free correlation *rises* to **0.677** (n = 807), a 10% increase. That is a **2.6-fold** modality gap under matched restriction, against no gap when pooled. The two movements are not symmetric and should not be read as such: the in-vivo collapse discards 10,420 of 14,088 sequences and is the robust half of the contrast, whereas the cell-free rise follows from dropping 55 of 862 and is small enough to be fragile. All four figures are population Spearman correlations over the full stated N, not resampled estimates, so no interval is attached to them.

![**Figure 3 — The DRAFTS modality comparison.** Cross-host Spearman ρ for *E. coli*–*B. subtilis*, cell-free vs. in vivo, under the pooled and co-active regimes. Pooled n = 862 (cell-free) and 14,088 (in vivo); co-active n = 807 and 3,668.](figures/Figure3_drafts_modality_comparison.png)

**The comparison rests on a single species pair.** Only *E. coli* and *B. subtilis* appear in both datasets; *P. aeruginosa* is absent from DRAFTS, and the three RS241 hosts DRAFTS does cover have at most 207 in-vivo sequences, of which the co-active subset would be smaller. We did not compute cross-modality correlations for them, and note that doing so is the most direct route to widening this comparison. The claim should be read as one well-characterized instance, not a survey.

The mechanism is visible in how much each modality loses to the restriction. Of the EC–BS sequences measured in both hosts, **93.6% (807 of 862) are co-active in lysate against 26.0% (3,668 of 14,088) in vivo** (transcription) — the same operation on the same pair. Cell-free activity is also continuous rather than bimodal, with no exact zeros among usable rows and no floor pile-up (the largest single rounded value accounts for 0.60% of rows); across the ten DRAFTS species the below-detection class is under 5% in nine, with *S. enterica* at 45.9% the exception. In vivo, cross-host agreement is carried substantially by agreement about what does not work. In lysate there is little silence to agree about.

**Implications for cell-free prototyping.** Cell-free measurement is a reasonable proxy for behavior *within* a host — within-species in-vitro/in-vivo agreement runs 0.69–0.90 across seven species (*E. coli* 0.901, *B. subtilis* 0.693), recomputed from DRAFTS's released source data and close to but not identical with their stated 0.71–0.90. What it does not reproduce is host-to-host difference, because the dominant in-vivo difference is which sequences are silent. A practitioner screening constructs for one fixed chassis via lysate should expect results broadly consistent with in vivo. A practitioner ranking a part across candidate hosts should not.

### 3.3 Magnitude, robustness, and the published objections

**DRAFTS's authors do not make this comparison.** We searched DRAFTS's full text, all six Appendix tables, all fourteen Appendix figures, all five source-data spreadsheets' structure, and the peer-review correspondence for any comparison of cross-species correlation between the two modalities. None exists. DRAFTS reports within-species in-vitro/in-vivo fidelity (their Appendix Fig. S5; 94.2% of *E. coli* activities within 1-log variation, their Fig. 1D) and within-cell-free cross-species structure (their Fig. 3) as separate findings, adjacent in their Discussion but never combined.

**The magnitude is regime-dependent.** Across six restriction levels on the 862-sequence shared set, cell-free EC–BS ranges 0.386–0.677. Each row below is compared against the in-vivo figure computed under the matching restriction: pooled against pooled, every restricted regime against the in-vivo co-active 0.258.

| Regime | Cell-free ρ (n) | In-vivo ρ (n) | Ratio |
|---|---|---|---|
| Pooled | 0.616 (862) | 0.655 (14,088) | 0.94 |
| Co-active | 0.677 (807) | 0.258 (3,668) | 2.6× |
| Both above 25th percentile | 0.647 (499) | — | — |
| Both above median | 0.592 (287) | — | — |
| Active-fraction matched (per-host) | 0.509 (218) | 0.258 | 2.0× |
| Top 26% of both (uniform quantile) | 0.386 (150) | 0.258 | 1.5× |

The per-host matched variant restricts each cell-free host to the top fraction matching that host's own in-vivo active fraction (*E. coli* 61.1%, *B. subtilis* 28.0%), retaining 25.3% of shared sequences — close to the in-vivo co-active retention of 26.04%. The uniform variant applies 26.04% to each host independently, which compounds to 17.4% retention and is therefore the stricter of the two. Both take quantiles over the DNA-adequate population rather than the co-detected one. The intermediate quantile cuts have no in-vivo counterpart, since the in-vivo restriction is a detection threshold rather than a quantile. Under the harshest cut the contrast falls to 1.5×, and 4 of the 45 DRAFTS species pairs fall below 0.258 — all four involving *L. lactis*. Those four should be read with caution for a second reason: Yim et al. excluded ScrFI/CCNGG-containing sequences from their *L. lactis* analysis as restriction artifacts, so a *L. lactis* pair may carry residual artifact that has nothing to do with the modality contrast. The 45-pair band is likewise regime-specific: 0.623–0.911 co-active, 0.364–0.892 pooled. EC–BS (0.677) is the third-lowest of the 45 under the regime we use, so the reported pair is conservative rather than favorable.

**The contrast survives a source-composition control.** The GC partial-correlation control applied to the in-vivo pairs in Section 3.1 applies equally to this comparison, where it is the check that speaks most directly to a selection artifact, so we apply it here. Controlling for source-genome GC within each modality, the cell-free co-active correlation falls from 0.677 to 0.556 (a 17.9% reduction) and the in-vivo co-active correlation from 0.258 to 0.207 (19.8%); the ratio moves from 2.63× to **2.69×**. Pooled, the cell-free correlation falls from 0.616 to 0.561 (8.9%) and the in-vivo from 0.655 to 0.517 (21.1%), moving the ratio from 0.94× to 1.09× — still no modality difference. Under co-active restriction the two modalities attenuate by nearly the same proportion, so the co-active contrast — the one this paper's claim rests on — is not an artifact of source-genome GC composition. The pooled attenuations differ considerably, but the pooled comparison shows no modality difference before or after the control, so nothing in the argument turns on them. The closed-form partial correlation and an independent rank-residual regression agree (0.556 against 0.548 cell-free; 0.207 against 0.205 in vivo). One scope limit: the two modalities are measured on different sequence populations, so GC is controlled within each separately — this is not a paired control.

Prior work has noted that cell-free and in-vivo measurements can diverge. Pandi et al. (2022) report a cell-free/in-vivo yield correlation of 0.41 for 20 optimized constructs, concluding that "the optimal candidates are not necessarily directly transferable in vivo." That concerns a different quantity — absolute yield within a single host — but supports the general caution this section sharpens by identifying which axis of fidelity holds and which does not.

**The strongest published argument the other way.** Yim et al.'s own Fig. 3 reports that cell-free cross-species correlation clusters into distinct gram-negative and gram-positive groups (their Fig. 3C), and that more phylogenetically related species share more similar transcription profiles, at r = 0.73 against 16S rRNA similarity (their Fig. 3D). Read as a claim about host specificity, that is the most direct published challenge to the framing here. The two results concern different quantities: their gradient describes how cross-host agreement *varies* with phylogenetic similarity within the cell-free system, while the claim here is about its *level* under matched restriction. Both can hold at once, and here they do — the 45-pair cell-free band runs 0.623–0.911 co-active and never approaches the in-vivo co-active 0.258 for the one pair measured in both systems, so a phylogenetic gradient inside that band does not account for the modality difference. The direct test is to recompute their gradient under the restriction regimes used here; the released data supports it, and **we have not run it**. We flag it as the open check rather than a settled point, because the gradient could as easily strengthen under co-active restriction as weaken. It is one of two checks this paper leaves open; the other, cross-modality correlations for the three RS241 hosts DRAFTS covers, is stated in Limitations item 2.

### 3.4 The mechanistic hypothesis

Lysates retain the enzymatic core of transcription but remove the membrane, chromosomal architecture, growth-phase physiology, and resource-competition context of a living cell. The measured consequence is that the large in-vivo silent class largely disappears, and with it the component of cross-host agreement that co-activity conditioning removes.

This is consistent with established links between each removed element and transcriptional output. DNA supercoiling acts as a global transcriptional regulator whose effect on a promoter depends quantitatively on its discriminator GC content (El Houdaigui et al., 2019). DRAFTS supplied its library as midiprepped plasmid, so those templates are supercoiled; what a lysate lacks is the gyrase/topoisomerase homeostasis and chromosomal context that regulate supercoiling in a living cell, leaving it fixed at whatever the preparation produced rather than dynamically maintained. That is a weaker claim than absence, and it is the accurate one. Resource competition between a gene of interest and the host's own machinery is measurable from lysate capacity assays and predicts in-vivo burden well (in vitro/in vivo capacity R² = 0.74; Borkowski et al., 2018). Growth-phase physiology is, by construction, a property only a dividing cell has.

**This is a mechanistic hypothesis consistent with the data, not a demonstrated causal claim.** The individual mechanisms are established; the link from them to the modeling failure in Section 3.5 is an inference. If host specificity resides substantially in context-dependent silencing rather than in fixed machinery differences, then host-descriptor features built from machinery and composition annotations — genomic or proteomic — should carry no signal for cross-host prediction, however well constructed. Confirming this directly would require titrating physiological context back into a lysate, or measuring the same library in vivo across a matched growth-condition series in multiple hosts.

**The strongest competing account comes from the dataset's own authors.** Yim et al. observe that their DNA library distribution is far more uniform in vitro than in vivo, and attribute that to the decoupling of gene expression from cell fitness in a lysate. Their statement is about library representation rather than about silent sequences, so extending it to activity is our reading, not their claim — but extended, it predicts what we observe: in a living cell, expression that is costly can be selected against or shut down, producing a large silent class, while a lysate has no fitness for expression to couple to and silences almost nothing. That account and the contextual one above make the same prediction for every quantity available here. Separating them requires a per-sequence fitness or burden measurement in each host, which neither dataset provides; GC, cross-host activity and the silent fraction are consistent with both. We name it as a competing explanation of equal standing rather than one this paper rules out. Both accounts nonetheless place the difference in cellular context rather than in transcriptional machinery — expression coupled to growth is a property of a living cell, not of its polymerase. What the data available here cannot settle is which contextual mechanism dominates, not whether the difference is contextual.

### 3.5 Testing the prediction (Figures 4–5)

Four independent lines of evidence test the prediction above, and each comes out as Section 3.4 says it should.

![**Figure 4 — The pre-registered kill gate (H-MAIN).** Model vs. per-host baseline across all 8 primary comparisons, with 90% percentile bootstrap intervals (10,000 resamples, hierarchical over folds then draws). None met the pre-registered bar.](figures/Figure4_hmain_kill_gate.png)

**An arbitrary host tag matches or beats both feature vectors in every cell tested (pre-registered, H-DIAGNOSTIC).** A free per-host lookup embedding carrying zero biological content matched or beat both the 37-D genomic and the 6-D physiology vectors at zero-shot in all **12 of 12** tested (host, readout, variant) cells — zero clean wins for the biological features, and 10 cells where the free embedding won outright on point estimate. Neither vector carries detectable content beyond an arbitrary host label at this n_hosts.

**A feature-group ablation finds no group carrying signal.** A five-group ablation of the genomic vector (sigma-factor, anti-Shine-Dalgarno, tAI/codon usage, RNAP subunit count, chaperone/heme) found no group's removal changing rank correlation beyond fold-to-fold noise (max |Δ| = 0.089). This is scoped to the genomic vector; the physiology vector was not ablated.

**Sequence-only prediction is not distinguishably beaten under three conditioning mechanisms, in 16 of 18 tests.** Across 3 mechanisms × 6 (host, readout) cells, the two exceptions are both at *E. coli* transcription, where concatenation reaches 0.555 and per-host-heads-averaged 0.615 against sequence-only 0.367 and FiLM 0.371 — the one combination with the largest usable N and the only measured reliability estimate. Every *B. subtilis* cell, both readouts, is statistically indistinguishable across sequence-only and all three mechanisms. No mechanism recovers signal for the host the central question concerns.

![**Figure 5 — Zero-shot performance of every system.** Zero-shot Spearman ρ for all nine systems — sequence-only CNN, a free per-host embedding, genomic + FiLM, physiology + FiLM, the three conditioning mechanisms (concatenation, per-host heads averaged and nearest), DNABERT-2 and PromoGen2 — across all 6 (host, readout) cells, 90% bootstrap intervals, on one shared y-scale. Fold-to-fold variance is Figure 7.](figures/Figure5_conditioning_mechanisms.png)

**The pre-registered kill gate failed on all 8 primary comparisons.** Neither readout, neither host-feature variant, on either primary held-out host met the bar. The closest was *B. subtilis* transcription with genomic features: model ρ = 0.213 [0.154, 0.260] against baseline ρ = 0.218 [0.160, 0.251] — NOT MET, point estimate favoring the baseline, and not revisited under any alternative mechanism. Several cells failed by wide margins (*P. aeruginosa* genomic transcription: 0.104 against 0.447).

### 3.6 Calibration (Figure 6)

The sequence-only model's zero-shot classifier is well calibrated on *E. coli* (ECE 0.053 transcription, 0.099 translation) and severely miscalibrated on *B. subtilis* (0.429, 0.344) and *P. aeruginosa* (0.474, 0.431) — roughly a 3.5- to 9-fold degradation depending on the pairing. Split-conformal interval coverage for the strength regressor stayed near nominal regardless (77–94% empirical against 80/90% targets), since conformal intervals guarantee marginal coverage independent of the underlying model's calibration; the ECE failure is a point-probability problem that coverage numbers alone would not reveal.

![**Figure 6 — Off-distribution calibration collapse.** Reliability diagrams for the zero-shot active/inactive classifier, all three primary hosts and both readouts, with Expected Calibration Error per panel.](figures/Figure6_calibration_collapse.png)

Published cross-dataset predictors show comparable off-distribution behavior: LaFleur, Hossain & Salis (2022) report R² = 0.80 on their own data but 0.45–0.65 across three independent external in-vivo datasets. The practical implication is that a model winning on rank correlation is not automatically trustworthy in absolute probability terms cross-host, and any deployment on a meaningfully different host should recalibrate first.

### 3.7 Secondary findings (Figure 7)

**Co-activity conditioning reverses the in-vivo ordering of host pairs.** Under the pooled regime EC–BS ranks highest and EC–PA second (Section 3.1); restricting to co-active sequences halves both *B. subtilis* pairs (0.655 → 0.258; 0.508 → 0.257) while raising EC–PA (0.621 → 0.754), moving it from the middle of the ordering to the top. Section 3.2 shows the same operation raises the cell-free correlation rather than lowering it. Range restriction, the classical mechanical cause of lower observed correlation under subsetting, does not explain the in-vivo drop: the restricted subset's interquartile range is 5×–53× *wider*, not narrower.

**FiLM is measurably unstable, independent of whether conditioning helps.** Its fold-to-fold standard deviation of zero-shot Spearman ρ exceeds every alternative's — sequence-only, concatenation, per-host heads — in all six (host, readout) cells, 4.0–7.8× higher at *E. coli* transcription (FiLM 0.200 against 0.031, 0.050, 0.026). A γ/β generator fit from only two training-host vectors is a poor default in this few-domain regime.

![**Figure 7 — FiLM fold-to-fold instability.** Fold-to-fold standard deviation of zero-shot ρ across all 6 (host, readout) cells and all four systems. No interval: this figure reports variance itself, computed across N = 5 folds per cell.](figures/Figure7_film_instability.png)

### 3.8 A retracted metric and a retracted finding (Figure 8)

**The ceiling metric.** A "percent of cross-host measurement-correlation ceiling" metric was proposed and then retired, because the sequence-only model genuinely exceeds the correctly computed ceiling for *B. subtilis* on both readouts (102–111%) — an explainable consequence of disattenuation, since measurement noise depresses a raw pairwise correlation in a way a model trained on thousands of examples is not. It appears nowhere in this paper's results.

**Shift-prediction.** A further attempt reframed the target: instead of absolute activity, predict the *shift* from a reference host's measured value to a target host, given sequence and the reference value as inputs. This initially appeared to succeed for *B. subtilis* transcription — 7 of 8 configurations for that host and readout distinguishably beat a mean-shift constant baseline, which no absolute-level framing had achieved.

It did not survive its control. Spearman correlation between the shift and the reference value itself is strongly negative for most host pairs (as low as −0.812, *P. aeruginosa*→*B. subtilis* transcription): a sequence with a high reference value mechanically has more room to fall than to rise. Across all 11 cells originally reported as wins, a reference-value-only baseline — ordinary least squares, one feature, no sequence input — matched or beat the original model in 10. A second control, retraining the identical architecture with sequences randomly permuted against their targets while keeping reference values intact, reproduced most of the original performance despite the sequence carrying no information (*E. coli*→*B. subtilis* transcription: original ρ = 0.432, shuffled-sequence 0.557 — *higher* with no sequence signal). **We retract this finding.**

The one cell surviving both controls (*P. aeruginosa*→*E. coli* translation with conditioning: 0.341 against reference-only 0.142 and shuffled −0.015) shares no host or readout with the retracted claim and is reported as a separate minor result. Any shift-prediction extension of this dataset should run both controls first.

![**Figure 8 — The shift-prediction retraction.** Original model vs. reference-value-only baseline vs. shuffled-sequence control, Spearman ρ against true shift, for the originally-reported winning cells. N ranges 294–9,609 per cell (for example, EC→BS transcription 3,589; EC→PA 9,609; BS→PA translation 294). Single point estimates on the held-out fold; the retraction rests on the comparisons, not on interval non-overlap.](figures/Figure8_shift_prediction_retraction.png)

### 3.9 Appendix: foundation models

DNABERT-2 (117M parameters; Zhou et al., 2024) and PromoGen2 (148M; Xia et al., 2026), evaluated uniformly via frozen embeddings plus a shallow head, lost to the 213,956-parameter sequence-only model in 35 of 37 statistically distinguishable comparisons across the primary hosts and RS241. We make no capacity claim from this: PromoGen2's own published native zero-shot protocol on the identical dataset reports 0.68/0.52/0.30 (EC/PA/BS transcription), higher than both our embedding-based numbers (0.495/0.490/0.243) and sequence-only itself at every host. The protocol, not either model's capability, explains the result.

---

## 4. Discussion

**Robustness, not proof.** The co-active contrast survived three checks — disattenuation for measurement error; partial-correlation control for source-genome GC composition on the in-vivo pairs, confirmed independently by phylum stratification (Section 3.1); and that same GC control applied to the modality contrast itself (Section 3.3). None was pre-registered; all three were motivated after the fact, by a subsequent review, by a finding in the cell-free dataset, and by a selection-artifact objection respectively, and surviving three does not rule out a fourth. What can be said is that the two confounds a careful reader proposes first have each been tested — GC twice, on the in-vivo pairs and on the modality contrast — and neither survives contact with the data.

**What this says about the physiology hypothesis: nothing decisive, in either direction.** Our physiology vector is a coarse reference-proteome composition proxy, not the direct wet-lab measurement (growth curves, transcriptomes) the Bernstein-lab papers used, and the genomic-feature ablation never touched it. What can be stated precisely is that the regime where Chan & Bernstein (2024) found physiology predictive — six closely related *Stutzerimonas* strains — is the kind of phylogenetically close comparison our own data shows has the highest cross-host measurement agreement, and the only regime where any conditioning mechanism here showed a benefit. These are compatible claims tested in non-overlapping regimes of host similarity, not competing claims tested in the same one.

**For host-conditioned modeling generally**, the instability finding in Section 3.7 — FiLM's fold-to-fold variance 4.0–7.8× worse than the simpler alternatives in a few-training-domain regime, independent of whether conditioning helps at all — is worth carrying into any host-conditioned model built on a handful of domains.

**What would settle the open questions**, both requiring wet-lab capability this study does not have: a direct measurement-reliability estimate for *B. subtilis* and *P. aeruginosa*, and for translation in any host; and new multi-host in-vivo data, ideally on hosts phylogenetically intermediate between *E. coli*/*P. aeruginosa* and *B. subtilis*, which would let a future study separate "conditioning helps only in the Gammaproteobacteria-close regime" from "conditioning helps only where N and reliability are highest" — two explanations three primary hosts cannot distinguish.

---

## 5. Limitations and unresolved gaps

1. **n_hosts ≤ 6.** The generalization unit is the host, not the sequence: three primary hosts with dense coverage, three more at reduced coverage. Claims about "genome-encoded functions" or "bacterial regulatory prediction" as a general class are not supported by six data points. Every interval bootstraps over folds and hosts, not sequences.

2. **The modality contrast rests on one species pair and is regime-dependent.** *E. coli*–*B. subtilis* is the only pair present in both datasets. Cell-free EC–BS ranges 0.386–0.677 across six restriction levels against an in-vivo co-active 0.258 — a ratio between 1.5× and 2.6×. Under the harshest cut, 4 of 45 DRAFTS pairs fall below 0.258, all involving *L. lactis* — a species for which Yim et al. excluded ScrFI/CCNGG-containing sequences as restriction artifacts, so those four pairs may be artifact-affected. The pooled cell-free figures require imputing activity 0 for `no_RNA_counts` rows, a definitional assignment rather than a measurement. We did not compute cross-modality correlations for the three RS241 hosts DRAFTS covers, which is the most direct route to widening the comparison. The mechanism behind the contrast is also unidentified: Yim et al.'s fitness-decoupling account predicts the same observation, and separating it from the cellular-context account requires per-sequence fitness or burden data that neither dataset provides.

3. **No direct reliability estimate exists for *B. subtilis* or *P. aeruginosa*, or for translation in any host.** No replicate or condition-series data exists for them in the released tables. The one measured reliability (0.912, EC transcription) is for the host and readout the central claim depends on least; *B. subtilis*'s own reliability is bounded by a sensitivity grid, which is a different strength of evidence than a measurement.

4. **The translation floor artifact.** `protein_log10` is pinned at a per-host floor for 66.5% / 89.9% / 10.1% of nominally usable rows (EC/BS/PA). Corrected usable-for-regression N: 9,146 / 1,101 / 17,630.

5. **No within-sequence paired in-vivo/cell-free comparison at adequate N.** The DRAFTS–Johns overlap usable in both *E. coli* and *B. subtilis* simultaneously is 15 sequences (ρ = 0.147, reported for completeness and not used as evidence). Section 3.2 compares two populations measured on the same library by the same laboratory. DRAFTS measures transcription only, so those conclusions are scoped to transcriptional host-specificity.

6. **Fold variance exceeds draw variance and grows with N** — in 15 of 23 baseline combinations checked, reaching a 91× ratio at *E. coli* transcription N = 3,000. A single fixed test fold can look precise while being unrepresentative; use the full five-fold rotation when extending this benchmark.

7. **Cross-host calibration failure.** ECE 3.5–9× worse off the training host (Section 3.6). Recalibrate before trusting raw predicted probabilities cross-host.

8. **Regression to the mean is a serious confound for any shift-prediction extension of this dataset.** Run a reference-value-only baseline and a shuffled-sequence control before trusting such a result (Section 3.8).

9. **Three dataset collisions, each of which silently corrupts a naive join (Section 2.2).** "Usable" differs between the two datasets, so joining on `n_shared_usable` compares mismatched restriction regimes. The ID spaces are disjoint, so any join must use sequence text. And DRAFTS's `Pa` is *Pantoea agglomerans*, not *P. aeruginosa* — which is absent from DRAFTS entirely — so a script matching on the bare two-letter code merges two unrelated organisms.

10. **Evo 2 was not evaluated, and the two foundation models that were did not run under their native protocols.** Evo 2's official package requires CUDA on a Hopper GPU; this work ran on Apple Silicon. For the models evaluated instead, a uniform frozen-embedding protocol was chosen for comparability but measurably cost PromoGen2 relative to its native numbers (Section 3.9). The capacity objection is weakened, not closed.

11. **We could not separate two explanations for the transcription/translation asymmetry** — that translation is fundamentally less cross-host-conserved, versus that the FACS-seq readout is too noisy to support this analysis. Our working view is a mix weighted toward noise for *B. subtilis*, stated as inference.

12. **The held-out evaluation split is not cryptographically enforced.** The public development table contains the withheld fold's labels; the protocol relies on convention, not a technical barrier.

---

## Data and Code Availability

All data are public: Johns et al. (2018) via BioProject PRJNA431139 and the paper's supplementary tables; DRAFTS via the PMC6692573 open-access package; reference proteomes via PaxDb.

The CROSSHOST benchmark — frozen splits, evaluation code, nine baseline systems, the leakage, provenance, and restriction-regime audits, and the citation-verification script — is available at https://github.com/gabeykim/crosshost and archived at https://doi.org/10.5281/zenodo.23051834. Derived activity values from Johns et al. are redistributed with attribution; PromoGen2 and DNABERT-2 embeddings were computed from publicly released model weights and are not redistributed.

`make reproduce` (the fast path), `make audit`, and `make verify-citations` were verified from a fresh clone on macOS, Apple Silicon, Python 3.14.2. `make reproduce-full`, the complete from-raw-data pipeline including all model training, was not re-executed end to end, as it would cost the 40+ cumulative hours already expended; every script in the dependency graph has run and produced its output at least once.

The provenance audit matches output filenames against script text, which is a substring check rather than a dependency check. During development it reported no orphans for the nine preprint figures although no committed script produced them, because the figure-embedding script names every file. Producing scripts were subsequently added and the audit now covers them genuinely.

## Acknowledgements

I thank Harris Wang for helpful correspondence about the Johns et al. and DRAFTS datasets, and for connecting me with Sung Sun Yim and Nathan Johns.

## Funding

No funding was received for this work.

## Competing Interests

The author declares no competing interests.

## Author Contributions

G.K. conceived the study, performed all analyses, and wrote the manuscript.

---

## References

Every citation was verified against a canonical source — CrossRef for DOI-bearing works, the arXiv API for preprints, direct venue search where no DOI exists.

- Borkowski, O., Bricio, C., Murgiano, M., Rothschild-Mancinelli, B., Stan, G.-B. & Ellis, T. (2018). Cell-free prediction of protein expression costs for growing cells. *Nature Communications* 9, 1457. DOI 10.1038/s41467-018-03970-x.
- Chan, K., Baldwin, G.S. & Bernstein, H.C. (2023). Revealing the Host-Dependent Nature of an Engineered Genetic Inverter in Concordance with Physiology. *BioDesign Research* 5, 0016. DOI 10.34133/bdr.0016.
- Chan, K. & Bernstein, H.C. (2024). Pangenomic landscapes shape performances of a synthetic genetic circuit across *Stutzerimonas* species. *mSystems* 9(9), e00849-24. DOI 10.1128/msystems.00849-24.
- Cheng, W., Song, Z., Zhang, Y., Wang, S., Wang, D., Yang, M., Li, L. & Ma, J. (2025). DNALONGBENCH: a benchmark suite for long-range DNA prediction tasks. *Nature Communications*. DOI 10.1038/s41467-025-65077-4.
- Dalla-Torre, H., Gonzalez, L., Mendoza-Revilla, J. et al. (2025). Nucleotide Transformer: building and evaluating robust foundation models for human genomics. *Nature Methods* 22(2), 287–297. DOI 10.1038/s41592-024-02523-z.
- El Houdaigui, B., Forquet, R., Hindré, T., Schneider, D., Nasser, W., Reverchon, S. & Meyer, S. (2019). Bacterial genome architecture shapes global transcriptional regulation by DNA supercoiling. *Nucleic Acids Research* 47(11), 5648–5657. DOI 10.1093/nar/gkz300.
- Grešová, K., Martinek, V., Čechák, D., Šimeček, P. & Alexiou, P. (2023). Genomic benchmarks: a collection of datasets for genomic sequence classification. *BMC Genomic Data* 24, 25. DOI 10.1186/s12863-023-01123-8.
- Johns, N.I., Gomes, A.L.C., Yim, S.S., et al. (2018). Metagenomic mining of regulatory elements enables programmable species-selective gene expression. *Nature Methods* 15, 323–329. DOI 10.1038/nmeth.4633.
- LaFleur, T.L., Hossain, A. & Salis, H.M. (2022). Automated model-predictive design of synthetic promoters to control transcriptional profiles in bacteria. *Nature Communications* 13, 5159. DOI 10.1038/s41467-022-32829-5.
- Marin, F.I., Teufel, F., Horlacher, M., Madsen, D., Pultz, D., Winther, O. & Boomsma, W. (2024). BEND: Benchmarking DNA Language Models on Biologically Meaningful Tasks. *ICLR 2024*. arXiv:2311.12570.
- Pandi, A. et al. (2022). A versatile active learning workflow for optimization of genetic and metabolic networks. *Nature Communications* 13, 3876. DOI 10.1038/s41467-022-31245-z.
- Patel, A., Singhal, A., Wang, A., Pampari, A., Kasowski, M. & Kundaje, A. (2024). DART-Eval: A Comprehensive DNA Language Model Evaluation Benchmark on Regulatory DNA. *Advances in Neural Information Processing Systems* 37 (NeurIPS 2024 Datasets and Benchmarks Track). DOI 10.52202/079017-1981. arXiv:2412.05430.
- Xia, Y. et al. (2026). Design prokaryotic cis-regulatory elements using language model. *Nucleic Acids Research* 54(4), gkag122. [PromoGen2; PMC12907563.]
- Yim, S.S., Johns, N.I., et al. (2019). Multiplex transcriptional characterizations across diverse bacterial species using cell-free systems. *Molecular Systems Biology* 15, e8875. DOI 10.15252/msb.20198875. [DRAFTS; PMC6692573.]
- Zhou, Z., Ji, Y., Li, W., Dutta, P., Davuluri, R.V. & Liu, H. (2024). DNABERT-2: Efficient Foundation Model and Benchmark for Multi-Species Genome. *ICLR 2024*. arXiv:2306.15006. [Model weights: `zhihan1996/DNABERT-2-117M`.]
