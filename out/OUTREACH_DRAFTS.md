# CROSSHOST — Outreach Drafts

**SEND NOTHING FROM THIS FILE WITHOUT EXPLICIT GO-AHEAD.** These are
drafts for review, not queued messages. Contact names/emails are left as
placeholders where this project has not independently verified a current
address — filling those in and actually sending is a decision for Gabriel,
not this session. **Rewritten in full this gate** (Gate 11) to reflect the
final manuscript state (DRAFTS section, GC-control robustness) and a
reordered, expanded recipient list — the four drafts from Gate 9 are
superseded, not appended to.

**A naming trap, flagged explicitly because two of the five recipients
below could be confused with each other:** "Harris Wang lab, Columbia"
(recipient 1 — created the Johns et al. 2018 and DRAFTS libraries this
entire project is built on) and "Wang lab, Tsinghua" (recipient 5 —
DeepCROSS) are **two different, unrelated people who happen to share a
common surname.** Neither draft below has independently verified the
Tsinghua PI's given name — check before sending, and do not let the shared
surname cause the two contacts, asks, or drafts to bleed into each other.

---

## 1. Harris Wang lab, Columbia — two specific, answerable questions

**Highest-value contact by a wide margin.** This lab produced both Johns
et al. (2018), the in-vivo dataset this project's central claim rests on,
and DRAFTS (Yim et al., 2019), the cell-free dataset acquired in Gate 10.
Two concrete, narrow asks, both already identified and quantified by this
project's own work — not a general "can you help" request. Either answer
changes what this project can state with confidence, not just what it
would be nice to know.

**NOTE:** contact information not independently verified this session —
confirm the current corresponding author and address before sending.

> **Subject: Two specific questions about the 2018/2019 regulatory-sequence datasets (replicate reliability, RS241 gaps)**
>
> Hi [Harris Wang / corresponding author — verify current contact],
>
> I've built a benchmark (CROSSHOST) on your 2018 *Nature Methods*
> regulatory-sequence library and, more recently, brought in your 2019
> DRAFTS cell-free dataset (*Molecular Systems Biology*) as an independent
> check — full preprint/code attached/linked. Two specific things came up
> that you may have direct answers for, faster than I could reconstruct
> them from the released tables:
>
> **1. Measurement reliability for *B. subtilis* and *P. aeruginosa*, in
> vivo.** I estimated *E. coli* transcription reliability directly from
> five growth-condition measurements in your 2018 Supplementary Data Table
> 3 (~0.91–0.93, cross-checked two ways), and separately found a
> biological-replicate reliability figure for *E. coli* in your 2019
> DRAFTS Appendix (Table S6: 0.89 in vivo, 0.93 cell-free) — a useful
> second confirmation, but still *E. coli* only. I could not find any
> equivalent replicate or condition-series structure for *B. subtilis* or
> *P. aeruginosa* specifically, in either dataset's released tables. Does
> any such data exist — even unpublished or informal — for those two
> hosts? This is the single gap most likely to change how confidently I
> can state my paper's central finding (a cross-host measurement-agreement
> gap that is largest for *B. subtilis*, the host this reliability
> question matters most for).
>
> **2. The 34 missing RS241 sequences.** 34 of the 241 nominal RS241 oligo
> IDs have no recoverable regulatory-sequence text in the released
> supplementary tables (leaving 207 usable). Is there a released or
> recoverable source for these, or were they excluded from public release
> for a specific reason?
>
> Either answer is useful on its own — if the reliability data doesn't
> exist for those hosts, that's worth stating plainly in my paper too,
> rather than left ambiguous, which is currently what I do.
>
> [name]

---

## 2. Sung Sun Yim (DRAFTS co-first author) — your dataset, a finding you may not have seen

**Second-highest value: a narrower, more technical ask than #1, to the
person most likely to have already thought about this exact comparison.**
Per the task context this project was given: DRAFTS co-first author, was
at Stanford Bioengineering at the time of publication (2019) — **not
independently verified this session; current affiliation may have
changed in the years since, confirm before sending.**

> **Subject: A cross-modality finding in the DRAFTS dataset you may find interesting**
>
> Hi Dr. Yim [verify current affiliation and preferred contact — address
> may have changed since the 2019 DRAFTS paper],
>
> I've built a benchmark (CROSSHOST) testing cross-host bacterial
> regulatory-activity prediction on the Johns et al. 2018 in-vivo dataset,
> and brought in DRAFTS (your 2019 cell-free dataset) as an independent
> check on one specific question: is the cross-host "chassis effect" that
> dominates in-vivo measurements the same phenomenon in a cell-free
> system?
>
> Short version of what I found: no. For the one pair testable in both
> datasets (*E. coli*–*B. subtilis*), the cell-free cross-host correlation
> is more than double the in-vivo figure (0.60 vs. 0.26), and *B.
> subtilis* — a severe outlier in vivo — is unremarkable across your full
> ten-species cell-free panel (mean correlation 0.78, solidly mid-pack).
> Within-species agreement between the two modalities is separately good
> (0.69–0.90 across your RS234 seven-species in-vitro/in-vivo design), so
> this isn't a general cell-free measurement problem — it specifically
> doesn't reproduce cross-host *differences*.
>
> I read this as evidence that host specificity in this system resides
> substantially in cellular context (membrane, supercoiling, growth
> physiology) rather than in the transcriptional machinery a TXTL lysate
> retains — full writeup attached/linked, clearly flagged as a hypothesis
> consistent with the data, not a demonstrated mechanism.
>
> Has anyone in your group looked at this specific in-vivo-vs-cell-free
> cross-host comparison before? I'd be glad to know if this replicates or
> conflicts with anything you've already seen in the data, and would
> welcome any read on whether the mechanistic interpretation above is
> reasonable.
>
> [name]

---

## 3. de Boer lab (UBC) — GAME benchmark contribution

**Highest-leverage single action for adoption.** GAME (the de Boer lab's
benchmark framework) has a species-mapping layer and, per this project's
own charter-stage research, no microbial cross-host task. This converts
the largest benchmark-side entity with overlapping infrastructure from a
competitor into a distribution channel, and borrows standing this
project's solo, unaffiliated author does not otherwise have.

**NOTE:** This project has not independently re-verified GAME's current
task list or confirmed it still lacks a microbial module as of today —
that claim is carried from planning-stage research and should be
re-checked (visit the current GAME repository/paper) before sending.

> **Subject: Microbial cross-host module for GAME?**
>
> Hi [de Boer lab contact — verify current corresponding author before
> sending],
>
> I've built CROSSHOST, a benchmark for predicting bacterial regulatory
> DNA activity across host species — three primary hosts (*E. coli*, *B.
> subtilis*, *P. aeruginosa*) with dense coverage plus three more at
> reduced N, built on the Johns et al. 2018 metagenomic regulatory-sequence
> library, with frozen genome-blocked splits, a pre-registered hypothesis
> test, a full baseline suite, and a pip-installable package with a
> held-out evaluation split.
>
> The central finding is a negative result, tested unusually thoroughly:
> host-conditioned prediction does not beat a sequence-only model under
> three distinct conditioning mechanisms, two host-feature representations,
> and two target reframings — and an independent cell-free dataset on the
> same sequence library offers a candidate mechanistic explanation. Full
> writeup attached/linked. I think that makes it a useful addition to a
> benchmark ecosystem specifically because it's a clean, well-controlled
> negative result in a domain (bacteria) your framework doesn't currently
> cover, alongside a species-mapping layer that (as I understand it,
> though I'd want to confirm against your current codebase) is
> architecturally similar to what a microbial module would need.
>
> Would you be open to a microbial cross-host module contribution, or
> pointing me to how GAME handles new-species/new-domain submissions? Happy
> to share the full preprint, code, and data — everything is already
> packaged and tested from a clean install.
>
> [name]

---

## 4. Bernstein lab — precise, scoped, not a refutation

**Two papers, two different findings, must be addressed distinctly** (per
this project's own verification, `out/GATE8_5_MEMO.md` Task 2c): the 2023
BioDesign Research paper (6 hosts, 3 genera, direct physiology measurement,
no transcriptomes) and the 2024 mSystems paper (6 *Stutzerimonas* strains,
direct physiology + transcriptomes). The outreach names both and is
explicit that this project tested a coarser, database-derived
operationalization of "physiology," in a different (more phylogenetically
distant) host regime, and is not claiming to have tested their hypothesis
on their own terms.

> **Subject: A benchmark testing the genome-vs-physiology chassis-effect question — not a refutation, a scoping note**
>
> Hi Dr. Bernstein and Kevin [Chan — verify current contact/affiliation],
>
> Your 2023 BioDesign Research paper and 2024 mSystems paper on chassis
> effects and physiology motivated a benchmark I've built (CROSSHOST) to
> test the genome-vs-physiology question directly using the Johns et al.
> 2018 cross-host regulatory-activity dataset (*E. coli*, *B. subtilis*,
> *P. aeruginosa*, plus three more hosts at reduced N).
>
> I want to be precise about what this does and doesn't say relative to
> your work, because I think an imprecise framing would do your findings a
> disservice: my physiology feature vector is a coarse, 6-dimension,
> reference-proteome-database-derived proxy (PaxDb-style aggregate
> abundance fractions), not the direct wet-lab measurement (growth curves,
> plasmid copy number, transcriptomes) either of your papers used. And my
> host panel spans two phyla (Firmicutes vs. Gammaproteobacteria) at
> substantial phylogenetic distance, while your 2024 panel is six closely
> related *Stutzerimonas* strains. Neither genomic nor physiology-proxy
> conditioning beat a sequence-only model in my setup — but I don't think
> this tests your hypothesis on its own terms, and I say so explicitly in
> the paper.
>
> What I think IS a genuinely interesting scoping observation: my own data
> independently shows cross-host measurement agreement is highest between
> phylogenetically close hosts (exactly the *E. coli*/*P. aeruginosa*
> Gammaproteobacteria pair in my panel) and drops sharply for the more
> distant *B. subtilis* comparisons — surviving two robustness checks aimed
> directly at it (a measurement-noise correction, and independently, a
> real source-genome GC-composition confound I found via an unrelated
> cell-free dataset from the Wang/Johns lab). That pattern would predict
> your 2024 close-relative regime is close to the best case for
> physiology-based prediction to work, and mine is closer to the worst
> case. That's consistent with, not contradictory to, what you found.
>
> I'd value your read on whether this framing is fair, and whether a
> direct (not proxy) physiology measurement on a more phylogenetically
> distant panel is something either of our groups might be positioned to
> test. Full preprint/code attached/linked.
>
> [name]

---

## 5. Wang lab (Tsinghua) — DeepCROSS, generous framing, benchmark inclusion

**Context:** DeepCROSS (Wang lab, *Nat Commun* 2025) does cross-species
regulatory-sequence design across 2 bacteria, has its own MPRA capability,
and already uses the Johns et al. data — the group most likely to have
standing to extend this benchmark with genuinely new multi-host data, which
this project cannot generate (no wet lab, by charter). Generous framing:
offer inclusion, not competition. **See the naming-collision note at the
top of this file — this is a different Wang lab from recipient 1.**

> **Subject: CROSSHOST — a cross-host bacterial regulatory-activity benchmark, DeepCROSS-relevant**
>
> Hi [Wang lab contact — verify current corresponding author, and confirm
> given name; do not conflate with the Columbia Wang lab this project also
> contacts separately],
>
> I've built CROSSHOST, a benchmark for cross-host bacterial regulatory
> DNA activity prediction, built on the Johns et al. 2018 dataset your
> DeepCROSS work also draws on. The central result is a negative one:
> genome- or physiology-conditioned models don't beat a sequence-only
> baseline for held-out hosts in my setup, tested under several independent
> conditioning mechanisms and two foundation-model comparisons, with a
> cell-free dataset (DRAFTS, from a different lab that also worked with the
> Johns library) offering a candidate mechanistic explanation.
>
> I think DeepCROSS is positioned to do something I can't: you have your
> own MPRA capability and can generate new multi-host regulatory data,
> where I'm limited to existing published datasets (by design — no new
> experimental data is in scope for this project). If useful, I'd like to
> offer inclusion of DeepCROSS as a baseline/comparison system in the
> benchmark's public leaderboard-equivalent (a held-out evaluation split
> with withheld labels), and would welcome any interest in using CROSSHOST
> as an evaluation target for future DeepCROSS-generated data. No
> expectation either way — sharing in case it's useful. Full
> preprint/code/package attached/linked.
>
> [name]

---

## Sending checklist (for Gabriel, not to be acted on by this session)

- [ ] Verify each recipient's current affiliation/contact — all five
      drafts above contain a bracketed verification note; none has been
      checked this session.
- [ ] Double-check the Harris Wang (Columbia) / Wang lab (Tsinghua) name
      collision before sending either — see the note at the top of this
      file.
- [ ] Attach or link the actual preprint once posted to bioRxiv (see
      `out/SUBMISSION_CHECKLIST.md`, `out/SHIP_CHECKLIST.md`) — none of
      these drafts should be sent referencing a preprint that doesn't yet
      have a public URL.
- [ ] Decide sending order — Harris Wang lab first is the strongest default
      (their answers could materially affect the manuscript, e.g. the
      reliability disclosure in Limitations item 3, before it's finalized
      elsewhere), then Yim (a related, narrower ask to the same broad
      group — consider sending together or close in time), then de Boer
      lab, then Bernstein lab and Wang lab (Tsinghua) in either order.
- [ ] None of these have been sent. This file is drafts only.
