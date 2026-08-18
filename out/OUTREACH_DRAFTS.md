# CROSSHOST — Outreach Drafts

**SEND NOTHING FROM THIS FILE WITHOUT EXPLICIT GO-AHEAD.** These are
drafts for review, not queued messages. Contact names/emails are left as
placeholders where this project has not independently verified a current
address — filling those in and actually sending is a decision for Gabriel,
not this session.

Each draft states what was built, what is being asked, and nothing beyond
what this project's own results support — no overclaiming, per instruction.

---

## 1. de Boer lab (UBC) — GAME benchmark contribution

**Priority: highest-leverage single action available.** GAME (the de Boer
lab's benchmark framework) has a species-mapping layer and, per this
project's own charter-stage research, no microbial cross-host task. This
converts the largest benchmark-side entity with overlapping infrastructure
from a competitor into a distribution channel, and borrows standing this
project's solo, unaffiliated author does not otherwise have.

**NOTE:** This project has not independently re-verified GAME's current
task list or confirmed it still lacks a microbial module as of today — that
claim is carried from planning-stage research and should be re-checked
(visit the current GAME repository/paper) before sending, since benchmark
suites add tasks over time.

> **Subject: Microbial cross-host module for GAME?**
>
> Hi [de Boer lab contact — verify current corresponding author before
> sending],
>
> I've built CROSSHOST, a benchmark for predicting bacterial regulatory DNA
> activity across host species — three primary hosts (*E. coli*, *B.
> subtilis*, *P. aeruginosa*) with dense coverage plus three more at
> reduced N, built on the Johns et al. 2018 metagenomic regulatory-sequence
> library, with frozen genome-blocked splits, a pre-registered hypothesis
> test, a full baseline suite, and a pip-installable package with a
> held-out evaluation split.
>
> The central finding is a negative result, tested unusually thoroughly:
> host-conditioned prediction does not beat a sequence-only model under
> three distinct conditioning mechanisms, two host-feature representations,
> and two target reframings — full writeup attached/linked. I think that
> makes it a useful addition to a benchmark ecosystem specifically because
> it's a clean, well-controlled negative result in a domain (bacteria)
> your framework doesn't currently cover, alongside a species-mapping
> layer that (as I understand it, though I'd want to confirm against your
> current codebase) is architecturally similar to what a microbial module
> would need.
>
> Would you be open to a microbial cross-host module contribution, or
> pointing me to how GAME handles new-species/new-domain submissions? Happy
> to share the full preprint, code, and data — everything is already
> packaged and tested from a clean install.
>
> [name]

---

## 2. Bernstein lab — precise, scoped, not a refutation

**Two papers, two different findings, must be addressed distinctly** (per
this project's own verification, `out/GATE8_5_MEMO.md` Task 2c): the 2023
BioDesign Research paper (6 hosts, 3 genera, direct physiology measurement,
no transcriptomes) and the 2024 mSystems paper (6 *Stutzerimonas* strains,
direct physiology + transcriptomes). The outreach should name both and be
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
> distant *B. subtilis* comparisons — which would predict your 2024
> close-relative regime is close to the best case for physiology-based
> prediction to work, and mine is closer to the worst case. That's
> consistent with, not contradictory to, what you found.
>
> I'd value your read on whether this framing is fair, and whether a
> direct (not proxy) physiology measurement on a more phylogenetically
> distant panel is something either of our groups might be positioned to
> test. Full preprint/code attached/linked.
>
> [name]

---

## 3. Wang lab (Tsinghua) — DeepCROSS, generous framing, benchmark inclusion

**Context:** DeepCROSS (Wang lab, *Nat Commun* 2025) does cross-species
regulatory-sequence design across 2 bacteria, has its own MPRA capability,
and already uses the Johns et al. data — the group most likely to have
standing to extend this benchmark with genuinely new multi-host data, which
this project cannot generate (no wet lab, by charter). Generous framing:
offer inclusion, not competition.

> **Subject: CROSSHOST — a cross-host bacterial regulatory-activity benchmark, DeepCROSS-relevant**
>
> Hi [Wang lab contact — verify current corresponding author],
>
> I've built CROSSHOST, a benchmark for cross-host bacterial regulatory
> DNA activity prediction, built on the Johns et al. 2018 dataset your
> DeepCROSS work also draws on. The central result is a negative one:
> genome- or physiology-conditioned models don't beat a sequence-only
> baseline for held-out hosts in my setup, tested under several independent
> conditioning mechanisms and two foundation-model comparisons.
>
> I think DeepCROSS is positioned to do something I can't: you have your
> own MPRA capability and can generate new multi-host regulatory data,
> where I'm limited to the existing 2018 dataset (by design — no new
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

## 4. Johns et al. (2018) authors — two specific, answerable questions

**Highest expected value per the task's own framing: "the single request
most likely to strengthen the paper."** Two concrete, narrow asks the
original authors are uniquely positioned to answer, both already
identified and quantified by this project's own work — not a general
"can you help" request.

> **Subject: Two specific questions about the 2018 regulatory-sequence dataset (RS241 gaps, reliability estimation)**
>
> Hi [corresponding author(s) — Nathan Johns / Harris Wang, verify current
> contact],
>
> I've built a benchmark (CROSSHOST) on your 2018 Nature Methods
> regulatory-sequence library, testing whether host-conditioned models can
> predict cross-host transcription/translation activity — full
> preprint/code attached/linked. Two specific things came up in the process
> that you may have direct answers for, faster than I could reconstruct
> them from the released tables:
>
> **1. RS241 sequence gap.** 34 of the 241 nominal RS241 oligo IDs have no
> recoverable regulatory-sequence text in the released supplementary
> tables (leaving 207 usable). Is there a released or recoverable source
> for these 34 sequences we're missing, or were they excluded from public
> release for a specific reason?
>
> **2. Measurement reliability for *B. subtilis* and translation.** I found
> usable replicate-like structure (five *E. coli* growth conditions,
> Supplementary Data Table 3) that let me estimate *E. coli* transcription
> measurement reliability directly (~0.91–0.93, cross-checked two ways).
> I could not find any equivalent replicate/condition-series structure for
> *B. subtilis* or *P. aeruginosa* specifically, or for translation in any
> host, in the released tables. Does any such data exist — even
> unpublished or informal — that would let a direct reliability estimate be
> computed for those cells? This is the single gap most likely to change
> how confidently I can state my paper's central finding (a cross-host
> measurement-agreement gap that's largest for *B. subtilis*).
>
> Either answer is useful — if the data doesn't exist, that's worth stating
> plainly in my paper too, rather than left ambiguous.
>
> [name]

---

## Sending checklist (for Gabriel, not to be acted on by this session)

- [ ] Verify each recipient's current affiliation/contact — all four drafts
      above contain a bracketed verification note; none has been checked
      this session.
- [ ] Attach or link the actual preprint once posted to bioRxiv (see
      `out/SUBMISSION_CHECKLIST.md`) — none of these drafts should be sent
      referencing a preprint that doesn't yet have a public URL.
- [ ] Decide sending order — Johns et al. first is a reasonable default
      (their answers could materially affect the manuscript before it's
      finalized elsewhere), then de Boer lab, then Bernstein lab and Wang
      lab in either order.
- [ ] None of these have been sent. This file is drafts only.
