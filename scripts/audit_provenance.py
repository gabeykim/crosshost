"""
GATE 7 - Task 1.5 (standing check from here on, run alongside audit_leakage.py):
flags any file in out/ or out/results/ that no committed script in scripts/
appears to produce. This is a direct, general-purpose response to the
gate5_5_per_host_ceiling.json failure mode -- a result file that shaped two
written documents (GATE5_5_MEMO.md, PAPER_FRAMING.md) and a Gate 6 prediction,
produced by an unsaved ad-hoc computation with no way to inspect or re-verify
it.

METHOD: for every file under out/ (recursively, excluding out/models -- large
binary checkpoints, not "results" in the sense this check cares about) and
every .py file under scripts/, check whether the target file's basename
appears as a string literal anywhere in any script. This is a plain textual
match, not an AST/write-call analysis -- it will not catch a script that
constructs a filename dynamically (e.g. f"gate6_{model_tag}_loho_results.json"
with model_tag as a variable), so a companion PATTERN_ALLOWLIST below
whitelists known dynamic-filename patterns by regex, checked against scripts'
literal f-string/format templates. A file is an "orphan" only if neither the
exact basename nor any allowlisted dynamic pattern matches any script.

This is a heuristic, disclosed as such -- it will have false positives (a
script that only ever READS a file for further processing, never re-produces
it, e.g. a manually-curated input) and is not a substitute for actually
opening flagged files and checking. Report the list; investigate each one by
hand before trusting or retracting it.
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "out"
SCRIPTS = ROOT / "scripts"

EXCLUDE_DIRS = {"models", "PREPRINT"}  # models: binary checkpoints, not "results" files.
# PREPRINT (Gate 11, SR9): out/PREPRINT/figures/*.png are curated, RENAMED copies
# of already-provenance-verified out/figures/*.png files (Figure1_*.png etc., for
# submission numbering) -- their provenance is documented by filename mapping in
# out/PREPRINT/FIGURE_AUDIT.md, not by direct basename match against scripts/, so
# the textual-match heuristic below would always flag them. Excluding the whole
# directory is deliberate, not a gap: every file under it is a copy or hand-written
# document, never a fresh, only-here result.
EXCLUDE_SUFFIXES = {".md"}  # memos are hand-written documents, not script outputs by design

# OS/editor detritus that is gitignored and is not project content. Gate 15:
# a stray out/.DS_Store tripped this audit as a one-file "orphan" -- the same
# gitignore-unawareness Gate 14 fixed in audit_leakage.py's check 6. Skipping
# these by name is deliberate and narrow: only files that .gitignore already
# excludes AND that no script could ever produce.
EXCLUDE_NAMES = {".DS_Store", "Thumbs.db"}

# Known dynamic-filename templates (model_tag / host / variant substitution)
# seen in scripts/61, 63, 64 etc. -- basename won't literal-match, so allow
# the *_{tag}_* shape explicitly.
DYNAMIC_PATTERNS = [
    re.compile(r"^gate6_(dnabert2|promogen2|evo2)_(loho_results|calibration_curves|rs241_results)\.json$"),
    re.compile(r"^(dnabert2|promogen2|evo2)_(library|rs241)\.npz$"),
    re.compile(r"^fm_(dnabert2|promogen2|evo2)_loho_.+\.pt$"),
]


def all_script_text():
    texts = {}
    for p in sorted(SCRIPTS.glob("*.py")):
        texts[p.name] = p.read_text(errors="ignore")
    return texts


def find_producer(basename, script_texts):
    hits = []
    for script_name, text in script_texts.items():
        if basename in text:
            hits.append(script_name)
    if hits:
        return hits
    for pat in DYNAMIC_PATTERNS:
        if pat.match(basename):
            # confirm at least one script contains the static portion of the pattern
            static_frag = basename.split("_")[0]  # coarse, just for a sanity cross-check
            for script_name, text in script_texts.items():
                if static_frag in text and ("model_tag" in text or "f\"" in text or "f'" in text):
                    hits.append(script_name)
            if hits:
                return hits
            return ["<dynamic-pattern-matched, no confirming script found>"]
    return []


def main():
    script_texts = all_script_text()
    all_files = []
    for p in OUT.rglob("*"):
        if p.is_dir():
            continue
        if any(part in EXCLUDE_DIRS for part in p.relative_to(OUT).parts):
            continue
        if p.suffix in EXCLUDE_SUFFIXES:
            continue
        if p.name in EXCLUDE_NAMES:
            continue
        all_files.append(p)

    orphans = []
    covered = []
    for p in sorted(all_files):
        producers = find_producer(p.name, script_texts)
        rel = p.relative_to(ROOT)
        if producers:
            covered.append((str(rel), producers))
        else:
            orphans.append(str(rel))

    print("=" * 78)
    print("PROVENANCE AUDIT")
    print("=" * 78)
    print(f"{len(covered)} files with an apparent producing script, {len(orphans)} orphans "
          f"(no script references the filename), {len(all_files)} total checked "
          f"(excluding out/models/*.pt and out/*.md)")
    print()
    if orphans:
        print("ORPHAN FILES (no committed script appears to produce these):")
        for o in orphans:
            print(f"  {o}")
    else:
        print("No orphans found.")
    print("=" * 78)
    return orphans


if __name__ == "__main__":
    main()
