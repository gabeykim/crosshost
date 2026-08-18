# This package is a curated export. The authoritative `make reproduce` /
# `make reproduce-full` / `make audit` targets live in the SOURCE
# REPOSITORY's root Makefile (one level up), which operates on scripts/ and
# data/ -- the directories that actually regenerate anything. This
# package/ directory is what those targets PRODUCE for distribution, not a
# self-contained pipeline in its own right (see README.md's Reproducibility
# section for why).
#
# What genuinely IS self-contained here: the crosshost/ package's own
# tests and the held-out scoring script.

.PHONY: test audit-package help

help:
	@echo "This is the DISTRIBUTABLE PACKAGE. For full reproduction, use the"
	@echo "Makefile at the source repository root (../Makefile: audit, test,"
	@echo "reproduce, reproduce-full, package)."
	@echo ""
	@echo "make test           - run this package's own tests (loaders, evaluate API)"
	@echo "make audit-package  - sanity-check this package's shipped data files parse"

test:
	python3 -m pytest tests/ -v

audit-package:
	python3 -c "\
import crosshost; \
lib = crosshost.load_three_host_library(); \
rs = crosshost.load_rs241(); \
g = crosshost.load_genomic_features(); \
p = crosshost.load_physiology_features(); \
s = crosshost.load_splits(); \
print(f'library={len(lib)} rs241={len(rs)} genomic={g.shape} physiology={p.shape} splits={len(s)}'); \
print('All package data files load correctly.')"
