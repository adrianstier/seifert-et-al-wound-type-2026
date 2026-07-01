# Coral wound-type analysis — reproducible pipeline
# Usage: make all | make exp1 | make exp2 | make verify | make clean | make clean-derived

.PHONY: all exp1 exp2 verify clean clean-derived

all:
	./scripts/RUN_COMPLETE_ANALYSIS.sh

exp1:
	Rscript scripts/airbrush_dremel_10_15_2025.R
	Rscript scripts/exp1_algae_glmm.R
	Rscript scripts/exp1_na_backtrack_audit.R
	Rscript scripts/exp1_sensitivity_6b.R
	Rscript scripts/exp1_sensitivity_extended.R
	Rscript scripts/exp1_ordinal_clmm.R
	Rscript scripts/exp1_multiplicity.R
	Rscript scripts/exp1_dharma_diagnostics.R

exp2:
	Rscript scripts/exp2_01_fit_glmm_models.R
	Rscript scripts/exp2_02_create_figures.R
	Rscript scripts/exp2_03_create_tables.R
	Rscript scripts/exp2_04_run_dharma_diagnostics.R
	Rscript scripts/exp2_05_create_early_timepoint_figure.R
	Rscript scripts/exp2_06_coenosarc_polyp_lag.R
	Rscript scripts/exp2_07_spline_sensitivity_and_icc.R

# Reproducibility / integrity self-check. Run AFTER `make all`.
verify:
	Rscript scripts/verify_pipeline.R

# Remove regenerated Exp 1 outputs only (raw data + archive untouched).
clean-derived:
	rm -rf output/tables output/figures output/text

# Comprehensive clean: removes ALL regenerated outputs (Exp 1 + Exp 2 +
# combined diagnostics) and the stray repo-root Rplots.pdf.
# Raw data (data/) and provenance (archive/, scripts/archive/) are NEVER touched.
clean: clean-derived
	rm -rf output/exp2_diagnostics output/exp2_figures_main \
	       output/exp2_figures_supplement output/exp2_models output/exp2_tables
	rm -rf output/diagnostics_combined
	rm -f Rplots.pdf
