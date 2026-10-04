**Were the artifacts source defects or rendering defects, and are any others of the same class present? Neither — all three are artifacts of how the PDF's text was extracted or copied. The Markdown source is correct in all three cases and so is the rendered page. Two of the three do not exist in the PDF's text layer at all. The third is one of seventeen identical instances, which is itself the proof that it is systematic extraction behaviour rather than a defect. Nothing was edited.**

# GATE 37 MEMO — Three reported typographic artifacts

Both audits clean at start and end. **No manuscript change, no PDF change, no number changed.** md5 `6e5bce0c3ea0dd27b99c0a68e1b4aac6`, unchanged since Gate 33.

*(Numbering note: `out/state.json` goes from `gate_35` to `gate_37`. No Gate 36 ran in this session — flagging in case one was lost rather than skipped.)*

---

## TASK 1 — Diagnosis

### Artifact 1 — "Stutzerimonasstrains"

**Source is correct.** Line 30 reads `…across six closely related *Stutzerimonas* strains.` — the space is there.

**It does not appear in the PDF either.** Plain `pdftotext` extraction line 42 reads `across six closely related Stutzerimonas strains`, with the space. Across the whole document: `"Stutzerimonas strains"` occurs **twice**, `"Stutzerimonasstrains"` **zero** times.

### Artifact 2 — "B. subtilissimultaneously"

**Source is correct.** Limitations item 5 reads `…in both *E. coli* and *B. subtilis* simultaneously is 15 sequences…`.

**It does not appear in the PDF either.** Extraction line 597 reads `B. subtilis simultaneously`. `"subtilis simultaneously"` occurs **once**, `"subtilissimultaneously"` **zero** times.

### Artifact 3 — "sequenceonly"

**Source is correct** (`sequence-only`, hyphen present), and **the rendered page is correct**.

This one *does* appear in plain `pdftotext` output, once, in the Abstract. The cause is visible with `-layout`:

```
a 37-feature genomic and a 6-feature physiology vector in all 12 cells; and sequence-
only prediction was not distinguishably beaten in 16 of 18 tests across three condition-
```

LaTeX breaks `sequence-only` at its own hyphen — correct typography. Plain `pdftotext` joins the two lines and drops the trailing hyphen, producing `sequenceonly`. **Confirmed visually**: I rendered page 1 at 300 dpi and cropped to that line. The page shows `…and sequence-` at the line end, **hyphen clearly present**, then `only prediction…` on the next line.

---

## TASK 2 — The right fix is none

All three are correct in the source and correct on the page. Applying any of the proposed remedies — a non-breaking space after an italic close, a zero-width non-joiner, `\mbox{}` around `sequence-only` — would change correct typesetting to work around a text-extraction behaviour. `\mbox{sequence-only}` in particular would **forbid a legitimate line break** and make the abstract's justification worse to fix a string that only exists outside the document.

**No edit made.** This falls under the gate's third option: "Something else → report it."

---

## TASK 3 — The whole class, scanned

### Class A — hyphenated terms broken at their own hyphen: **17 instances**

Every one produces a run-together word in plain `pdftotext` output, and every one is correct on the page:

`sequence-only`, `Johns-active`, `host-specific`, `larger-scale`, `Partial-correlation`, `closed-form`, `cell-free`, `1-log`, `862-sequence`, `gram-positive`, `per-sequence`, `cross-host`, `shuffled-sequence` (×2), `Gammaproteobacteria-close`, `fitness-decoupling`, `Host-Dependent`.

**That the count is 17 and not 1 is the finding.** A rendering defect would not hit every hyphenated term in the document uniformly; a line-joining rule in the extractor does exactly that. The reported `sequenceonly` is not special — it is the one someone happened to notice.

### Class B — italic boundaries: **27 checked, 0 affected**

Every `*italic* word` boundary in the source was checked against the extracted text. All 27 keep their space. The two reported instances are in that set and both survive.

### Figure captions

All 9 embedded caption lines scanned. Nothing of either class.

**Where the reports came from:** not from this PDF. Class A is reproducible from plain `pdftotext`; Class B is not reproducible at all here, and is the signature of a PDF viewer's copy-to-clipboard or an extractor that omits spaces at font changes. Gate 27 established that a previous external review was reading a stale PDF; this one looks like a tooling difference rather than a stale file, because Class B does not depend on which version of the PDF is read.

---

## TASK 4 — Verification

The committed `manuscript.pdf` was **not rebuilt**: its source did not change, and regenerating it would have rewritten 1.2 MB of binary with a new timestamp and no content difference. The build was instead verified by running the documented command to a scratch target.

| check | result |
|---|---|
| build command (to scratch target) | **exit 0**, 0 resource warnings, 9 images |
| committed PDF newer than both copies | **YES** (16:39:52 > 16:39:38) |
| both copies byte-identical | **YES**, md5 `6e5bce0c3ea0dd27b99c0a68e1b4aac6` |
| pages / images | **22** / **9** |
| captions | Figures 1, 2a, 2b, 3, 4, 5, 6, 7, 8 — ascending |
| glyphs | ρ 17, ≤ 3, γ 2, β 2, ≈ 1 |
| `Section X.Y` references | 10, all resolve |
| `Limitations item N` references | 1, resolves (12 items) |
| abstract / body | **250** / **6,950** words |
| `audit_leakage.py` | ALL 6 CHECKS PASSED |
| `audit_provenance.py` | 0 orphans, 177 files |
| `make test` | **9 passed** |
| `make verify-citations` | **11 / 0 / 5 / 0 / 16 — full baseline** |

---

## Verdict

Three reports, three false alarms, zero edits. Two of the three cannot be reproduced from this PDF at all; the third is one of seventeen instances of normal hyphenation being flattened by a text extractor. The source and the page are both correct, and the remedy for all three is to read the page rather than its extraction.
