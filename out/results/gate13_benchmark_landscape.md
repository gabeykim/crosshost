# Gate 13, Task 2 — Benchmark-landscape coverage, per suite

**Purpose.** The manuscript's Introduction previously claimed, as an unverified working understanding, that BEND, Genomic Benchmarks, DART-Eval, DNALONGBENCH, and the Nucleotide Transformer suite are "predominantly built on human, animal, or plant genomes." That is an exhaustive negative (no bacterial cross-host task exists among them) that was never checked. This gate checks it directly, retrieving each suite's own paper/repo and reporting what organisms and task types it actually covers.

**Method.** For each suite: identified the primary paper via WebSearch, fetched it (or its abstract/full text) via WebFetch, and searched specifically for any bacterial-genome task, any cross-species/cross-host comparison task, and the full organism list. Bibliographic details for each paper cross-checked mechanically against CrossRef/arXiv via `scripts/97_verify_citations.py` (Task 1) — see `out/results/gate13_citation_verification.json` for the raw comparison.

**Result, stated up front: no suite examined contains a bacterial task of any kind, let alone a bacterial cross-host regulatory-activity task. The stop condition ("if any suite turns out to contain something closer to a bacterial cross-host task than expected, stop and report") is not triggered.**

---

## BEND (Benchmarking DNA Language Models on Biologically Meaningful Tasks)

**Citation:** Marin, F.I., Teufel, F., Horlacher, M., Madsen, D., Pultz, D., Winther, O. & Boomsma, W. (2024). BEND: Benchmarking DNA Language Models on Biologically Meaningful Tasks. *ICLR 2024*. arXiv:2311.12570.

**Source checked:** arXiv abstract page (arxiv.org/abs/2311.12570) and OpenReview (openreview.net/forum?id=uKB4cFNQFg).

**Coverage:** Tasks defined **exclusively on the human genome** — the paper's own abstract states its tasks are "defined on the human genome." Seven tasks total: two long-range (gene finding, enhancer annotation), three genome-scale (histone modification, CpG methylation, chromatin accessibility), two zero-shot (noncoding variant effect for expression and disease).

**Bacterial or cross-host content: none.**

---

## Genomic Benchmarks

**Citation:** Grešová, K., Martinek, V., Čechák, D., Šimeček, P. & Alexiou, P. (2023). Genomic benchmarks: a collection of datasets for genomic sequence classification. *BMC Genomic Data* 24, 25. DOI 10.1186/s12863-023-01123-8.

**Source checked:** PMC full text (PMC10150520).

**Coverage:** Eight/nine datasets (the paper's own count varies slightly by counting convention) spanning **four organisms: human (primary, 6 datasets), mouse, roundworm (*C. elegans*), and fruit fly (*D. melanogaster*)** — the paper's own words: *"our datasets include primarily human data, but also mouse, roundworm and fruit fly."* Task types: exon–intron detection, splice-site detection, gene-expression prediction, MPRA cis-regulatory activity, enhancer classification, promoter prediction, MPRA SNP perturbation, eQTL mapping.

**Bacterial or cross-host content: none.** All four organisms are eukaryotes; no task compares activity of the same sequence across organisms — each dataset is single-organism.

---

## DART-Eval (A Comprehensive DNA Language Model Evaluation Benchmark on Regulatory DNA)

**Citation:** Patel, A., Singhal, A., Wang, A., Pampari, A., Kasowski, M. & Kundaje, A. (2024). DART-Eval: A Comprehensive DNA Language Model Evaluation Benchmark on Regulatory DNA. *Advances in Neural Information Processing Systems* 37 (NeurIPS 2024 Datasets and Benchmarks Track). DOI 10.52202/079017-1981. arXiv:2412.05430.

**Source checked:** arXiv full HTML (arxiv.org/html/2412.05430v2).

**Coverage: entirely human.** Positive-element set derived from 2.3 million ENCODE candidate cis-regulatory elements (cCREs); five human cell lines used throughout (GM12878, H1ESC, HEPG2, IMR90, K562); a caQTL task uses variants from African lymphoblastoid cell lines (still human). The paper's own stated limitation: *"our current evaluations are limited to tasks involving short, local sequence contexts."*

**Bacterial or cross-host content: none.** This is the suite named specifically for regulatory DNA, and it is the most directly comparable of the five to this project's own task framing (regulatory-sequence activity prediction) — its complete absence of any cross-species or bacterial component sharpens, rather than weakens, this project's novelty claim.

---

## DNALONGBENCH (A Benchmark Suite for Long-Range DNA Prediction Tasks)

**Citation:** Cheng, W., Song, Z., Zhang, Y., Wang, S., Wang, D., Yang, M., Li, L. & Ma, J. (2025). DNALONGBENCH: a benchmark suite for long-range DNA prediction tasks. *Nature Communications*. DOI 10.1038/s41467-025-65077-4.

**Source checked:** PMC full text (PMC12627797).

**Coverage:** Five tasks, all with dependencies up to 1 Mb: enhancer–target gene interaction (human, K562 cells), eQTL prediction (9 human tissues), 3D genome/contact-map prediction (5 human cell lines), regulatory sequence activity prediction (**human and mouse**), transcription initiation signal prediction (human). **Organisms: human and mouse only.**

**A "cross-species" phrase appears in the paper but refers to something unrelated** — a reference to "sequential regulatory activity prediction across chromosomes" within a single species (checked directly), not a cross-organism comparison. The paper's own framing scopes every task as operating "within a cell," i.e., eukaryotic cellular context, not bacterial.

**Bacterial or cross-host content: none.**

---

## The Nucleotide Transformer suite

**Citation:** Dalla-Torre, H., Gonzalez, L., Mendoza-Revilla, J. et al. (2025). Nucleotide Transformer: building and evaluating robust foundation models for human genomics. *Nature Methods* 22(2), 287–297. DOI 10.1038/s41592-024-02523-z. (Published online Nov 2024; print issue Feb 2025 — both dates real, print date used as canonical per standard bibliographic convention.)

**Source checked:** PubMed abstract, CrossRef record.

**Coverage:** The *pretraining* corpus spans human (individual + reference) and 850 other species' genomes — broad, but pretraining breadth is not the same as task breadth. The actual **downstream benchmark is 18 curated tasks across four categories: promoter (human/mouse), enhancer (human), splice site (human/multi-species), and histone modification (yeast).** This is the broadest organism range of the five suites examined — it includes yeast, a eukaryotic microbe — but still no bacterial task and no task comparing the same regulatory sequence's activity across multiple recipient organisms (its "multi-species" splice-site task evaluates splice-site *recognition* across species' own genomes independently, not cross-host transfer of the same sequence).

**Bacterial or cross-host content: none.**

---

## Summary table

| Suite | Organisms (tasks, not pretraining data) | Bacterial task? | Cross-host task (same sequence, multiple recipient species)? |
|---|---|---|---|
| BEND | Human only | No | No |
| Genomic Benchmarks | Human, mouse, *C. elegans*, *D. melanogaster* | No | No |
| DART-Eval | Human only | No | No |
| DNALONGBENCH | Human, mouse | No | No |
| Nucleotide Transformer | Human, mouse, yeast | No | No |

**Conclusion, and the exact replacement sentence used in `out/MANUSCRIPT.md` / `out/PREPRINT/MANUSCRIPT.md`:**

> "The benchmark suites we examined — BEND, Genomic Benchmarks, DART-Eval, DNALONGBENCH, and the Nucleotide Transformer suite — define tasks on human, mouse, *C. elegans*, *D. melanogaster*, and yeast genomes; we found no bacterial cross-host regulatory activity task among them."

This is a bounded, checkable claim (five named suites, each individually verified above) rather than the prior exhaustive negative ("no existing benchmark suite poses a cross-host prediction task" — a claim about *every* benchmark suite that exists, which no finite check could support). It is defensible because every suite named was actually read, not assumed.
