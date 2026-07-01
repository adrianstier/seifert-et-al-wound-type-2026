# Scripts

Run from the repo root via `./scripts/RUN_COMPLETE_ANALYSIS.sh` or `make all`.

## Experiment 1 (2022) — authoritative

- **`airbrush_dremel_10_15_2025.R`** — the single authoritative Exp 1 pipeline.
  Reads `data/airbrush_dremel.csv` (2026-06-15 independent NA-backtrack re-score). Fits binomial
  GLMMs (`lme4`) for `regenerated` and `healed` (healed separation resolved by the
  re-score → lme4 headline OR ~26), Firth penalized refit (`brglm2`) as the
  regeneration headline and a healed sensitivity check,
  ordinal CLMM on 3-level `healed`, debris model, 6-scenario missing-data
  sensitivity + coral-7b exclusion, day-28 endpoint robustness, `emmeans`
  contrasts, `DHARMa` diagnostics, composition + regeneration figures.
  Writes to `output/{tables,figures,text}/` and prints a narrative summary.
  Migrated from the standalone 2022 wound-type repository.

## Experiment 2 (2025)

- `exp2_01_fit_glmm_models.R` — `glmmTMB` `outcome ~ treatment * ns(day_c,3) + (1|id)`
  for Porites outcomes; standard binomial GLMM for outcomes without separation,
  penalized binomial GLMM (weakly-informative Normal prior: Normal(0,3) slopes,
  Normal(0,10) intercept; Gelman et al. 2008) for outcomes with complete separation
  — fitted on the real 0/1 response. Observation-level random effects (OLRE) are a
  non-selected last-resort fallback; they were not used in the final models.
- `exp2_02_create_figures.R` — manuscript Fig 3D (`figure2_porites_all_outcomes`,
  7-outcome patchwork).
- `exp2_03_create_tables.R` — `output/exp2_tables/PUBLICATION_STATISTICS_TABLE.csv`
  (canonical Exp 2 source of truth).
- `exp2_04_run_dharma_diagnostics.R` — DHARMa for Exp 2 models.
- `exp2_05_create_early_timepoint_figure.R` — early-dynamics supplement.
- `create_combined_dharma_summary.R` — **currently broken (exit 1)**; tracked as
  a known issue, non-fatal in the pipeline.

## Notes

- Superseded Exp 1 scripts (`exp1_01_run_analysis.R`, `exp1_02..05_*`) were moved
  to `archive/scripts_old/` on 2026-05-15 — they used stale pre-backtrack data and a
  weaker model. Do not use.
- Exp 2 includes a multiple-testing sensitivity analysis (Bonferroni and BH-FDR)
  in `output/exp2_tables/MULTIPLE_TESTING_CORRECTION.csv`; unadjusted p is primary
  under the a-priori directed-hypothesis design. See README "Known open items".
