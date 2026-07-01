#!/bin/bash
###############################################################################
# Complete analysis pipeline — coral wound-type manuscript
#
#   Experiment 1 (2022): authoritative pipeline (lme4 + Firth/brglm2, ordinal
#                        CLMM, sensitivity, endpoint robustness, figures)
#   Experiment 2 (2025): glmmTMB time-spline models, figures, tables, DHARMa
#
# No results are hardcoded. Every reported number is read from the regenerated
# outputs. Canonical summaries:
#   Exp 1: output/text/paper_results_summary.md
#   Exp 2: output/exp2_tables/PUBLICATION_STATISTICS_TABLE.csv
#
# Usage:  ./scripts/RUN_COMPLETE_ANALYSIS.sh   (run from repo root)
###############################################################################
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
START=$(date +%s)
fail=0

run() {  # run <label> <script> <fatal:1|0>
  echo "----------------------------------------------------------------------"
  echo ">>> $1"
  echo "----------------------------------------------------------------------"
  if Rscript "$2"; then
    echo "    OK: $2"
  else
    if [ "${3:-1}" = "1" ]; then
      echo "    FATAL: $2 failed"; exit 1
    else
      echo "    WARN (non-fatal, known issue): $2 failed"; fail=1
    fi
  fi
}

echo "=== EXPERIMENT 1 (2022) — authoritative pipeline ==="
run "Exp 1: GLMM + Firth + CLMM + sensitivity + figures" scripts/airbrush_dremel_10_15_2025.R 1
run "Exp 1: PRIMARY fragment-level + paired analysis"    scripts/exp1_fragment_level_primary.R 1
run "Exp 1: algae/debris colonization GLMM"              scripts/exp1_algae_glmm.R 1
run "Exp 1: NA-backtrack + separation audit"             scripts/exp1_na_backtrack_audit.R 1
run "Exp 1: coral-6b leave-one-out sensitivity"          scripts/exp1_sensitivity_6b.R 1
run "Exp 1: extended sensitivity + day-28 Fisher"        scripts/exp1_sensitivity_extended.R 1
run "Exp 1: ordinal CLMM (3-level healing)"              scripts/exp1_ordinal_clmm.R 1
run "Exp 1: multiplicity (BH/Bonferroni)"                scripts/exp1_multiplicity.R 1
run "Exp 1: DHARMa + temporal + separation diagnostics"  scripts/exp1_dharma_diagnostics.R 1

echo "=== EXPERIMENT 2 (2025) ==="
run "Exp 2: fit GLMMs (glmmTMB, treatment x ns(day,3))"   scripts/exp2_01_fit_glmm_models.R 1
run "Exp 2: figures (Porites 7-outcome panel)"             scripts/exp2_02_create_figures.R 1
run "Exp 2: statistical tables"                            scripts/exp2_03_create_tables.R 1
run "Exp 2: DHARMa diagnostics"                            scripts/exp2_04_run_dharma_diagnostics.R 0
run "Exp 2: early-timepoint figure"                        scripts/exp2_05_create_early_timepoint_figure.R 0
run "Exp 2: coenosarc->polyp lag"                          scripts/exp2_06_coenosarc_polyp_lag.R 0
run "Exp 2: spline-df sensitivity + ICC (all outcomes)"    scripts/exp2_07_spline_sensitivity_and_icc.R 0
END=$(date +%s); MIN=$(( (END-START)/60 )); SEC=$(( (END-START)%60 ))
echo "======================================================================"
echo " PIPELINE DONE in ${MIN}m ${SEC}s"
[ "$fail" = "1" ] && echo " (some non-fatal diagnostic steps warned — see log)"
echo ""
echo " Canonical results (not hardcoded — read these files):"
echo "   Exp 1: output/text/paper_results_summary.md"
echo "   Exp 2: output/exp2_tables/PUBLICATION_STATISTICS_TABLE.csv"
echo "   Figures: output/figures/ (Exp 1) ; output/exp2_figures_main/ (Exp 2)"
echo "======================================================================"
