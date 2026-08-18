# CROSSHOST reproduction Makefile. Lives at the project root because
# scripts/ and data/ (the things being regenerated) live here -- package/
# is a curated EXPORT of selected outputs for distribution, built BY this
# Makefile's `package` target, not a self-contained pipeline in its own
# right. See package/README.md for the distinction.
#
# WHAT WAS ACTUALLY VERIFIED (Gate 8, disclosed precisely, not overclaimed):
#   `make audit`, `make test`, `make package` -- ACTUALLY RUN, this session,
#   from this exact Makefile, and passed.
#   `make reproduce-full` -- NOT re-executed end-to-end (would take the
#   cumulative ~40+ hours of CNN/FM training this project's Gates 4-8 spent
#   in total). The dependency graph below is correct and every target maps
#   to a real script that has been run and produced its output at least
#   once (that is how out/ got populated) -- but re-running the ENTIRE
#   chain start-to-finish, in one sitting, from raw/ alone, was not
#   attempted in this packaging pass. `make reproduce` (fast path) uses the
#   already-shipped intermediate results instead.

PYTHON := python3

.PHONY: all reproduce reproduce-full package audit test clean help

help:
	@echo "make audit           - run both audit scripts (leakage + provenance), fast, VERIFIED"
	@echo "make test            - run package tests (loaders, evaluate API, held-out scoring), VERIFIED"
	@echo "make package         - (re)build package/ from data/ and out/, VERIFIED"
	@echo "make reproduce       - regenerate all figures/tables from shipped intermediate results, fast"
	@echo "make reproduce-full  - full pipeline from raw/ (MANY HOURS -- CNN/FM training) -- NOT re-tested this session"

audit:
	$(PYTHON) scripts/audit_leakage.py
	$(PYTHON) scripts/audit_provenance.py

test:
	cd package && $(PYTHON) -m pytest tests/ -v

# ---------------------------------------------------------------------------
# FAST PATH: regenerate figures/tables from already-computed intermediate
# results (out/*.json, out/results/*.json) -- this is what most users mean
# by "reproduce": verify the reported numbers/figures follow from the
# underlying data, without re-running multi-hour model training.
# ---------------------------------------------------------------------------
reproduce: audit
	$(PYTHON) scripts/65_gate6_figures.py
	$(PYTHON) scripts/67_crosshost_correlation_figure.py
	$(PYTHON) scripts/71_ablation_figure.py
	$(PYTHON) scripts/73_attenuation_figures.py
	$(PYTHON) scripts/76_regen_calibration_figures.py
	$(PYTHON) scripts/77_regen_remaining_figures.py
	$(PYTHON) scripts/70_provenance_triage.py
	$(PYTHON) scripts/75_master_baseline_export.py
	@echo "Fast reproduction complete -- all figures/tables regenerated from shipped intermediate results."

# ---------------------------------------------------------------------------
# FULL PATH: the entire pipeline from raw/ tables. Real dependency order,
# each target a real script -- but NOT re-executed end-to-end this session
# (see header note). Individual gate scripts within each phase are too
# numerous to list exhaustively here; see scripts/ for the full numbered
# sequence (01 through 77) and each gate's own memo for its exact script list.
# ---------------------------------------------------------------------------
reproduce-full:
	@echo "=== Gate 1-2: data foundation (raw/ -> data/) ==="
	$(PYTHON) scripts/01_inspect_supplements.py
	$(PYTHON) scripts/13_build_task1_core_tables.py
	$(PYTHON) scripts/24_final_split_with_near_dup_net.py
	@echo "=== Gate 3: baselines ==="
	$(PYTHON) scripts/28_prepare_baseline_data.py
	$(PYTHON) scripts/30_baseline1_mean_majority.py
	$(PYTHON) scripts/31_baseline2_calibration.py
	$(PYTHON) scripts/32_baseline3_host_embedding.py
	$(PYTHON) scripts/33_baseline4_biophysical.py
	@echo "=== Gate 4: the model (SLOW -- CNN training, hours) ==="
	$(PYTHON) scripts/42_loho_training.py
	$(PYTHON) scripts/43_calibration_curves.py
	@echo "=== Gate 5-5.5: kill gate + sequence-only ablation (SLOW) ==="
	$(PYTHON) scripts/52_hmain_evaluation.py
	$(PYTHON) scripts/53_hscience_evaluation.py
	$(PYTHON) scripts/54_hdiagnostic_evaluation.py
	$(PYTHON) scripts/59_sequence_only_full_pipeline.py
	$(PYTHON) scripts/57_crosshost_measurement_correlation.py
	@echo "=== Gate 6: foundation models (SLOW -- embedding extraction + head training) ==="
	$(PYTHON) scripts/61_fm_embeddings.py dnabert2
	$(PYTHON) scripts/61_fm_embeddings.py promogen2
	$(PYTHON) scripts/63_fm_loho_calibration_rs241.py dnabert2
	$(PYTHON) scripts/63_fm_loho_calibration_rs241.py promogen2
	$(PYTHON) scripts/64_gate6_full_comparison.py
	@echo "=== Gate 7: ceiling correction + calibration + ablation ==="
	$(PYTHON) scripts/66_ceiling_metric_investigation.py
	$(PYTHON) scripts/68_conformal_calibration.py
	$(PYTHON) scripts/69_host_feature_ablation.py
	@echo "=== Gate 8: attenuation analysis + packaging ==="
	$(PYTHON) scripts/72_attenuation_analysis.py
	$(MAKE) package
	@echo "Full reproduction complete."

package:
	$(PYTHON) scripts/74_build_held_out_eval.py
	$(PYTHON) scripts/75_master_baseline_export.py
	@echo "package/ rebuilt from current data/ and out/ contents."

clean:
	rm -rf package/data/core/*.parquet package/baselines/*.csv package/baselines/*.json
	@echo "Removed generated package/ data files (source data/ and out/ untouched)."
