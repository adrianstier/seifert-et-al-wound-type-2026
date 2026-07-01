# =============================================================================
# MAIN ANALYSIS SCRIPT: COMPREHENSIVE GLMM ANALYSIS
#
# Coral Wound Healing Project - Experiment 2 (Porites spp.)
# Species: Porites spp.
# Treatments: Airbrush (tissue removal) vs Dremel (skeletal damage) vs Combo
#
# PURPOSE:
#   Fit GLMMs for ALL outcomes using a principled fallback ladder:
#   1. Standard binomial GLMM (outcomes without complete separation)
#   2. Penalized binomial GLMM with a weakly-informative Normal prior on the
#      fixed effects (Gelman et al. 2008, Stat Med 27:2865) for outcomes with
#      complete separation. This keeps a PROPER binomial likelihood on the
#      real 0/1 response -- it does NOT squeeze/transform the response.
#   3. Observation-level random effect (last-resort overdispersion fallback)
#
# APPROACH:
#   - Try standard GLMM first
#   - If separation detected (non-finite/extreme SEs), refit the SAME model
#     with priors = Normal(0,3) on slopes, Normal(0,10) on the intercept
#   - Extract p-values and effect sizes for all outcomes
#
# OUTPUT:
#   - output/exp2_models/*_summary.txt (full model summaries)
#   - output/exp2_models/ALL_OUTCOMES_SUMMARY.csv (which method worked)
#   - output/exp2_models/*_tidy.csv (tidy fixed effects; consumed by exp2_03)
#
# USAGE:
#   Rscript scripts/1_run_glmm_analysis.R
#
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(here)
  library(glmmTMB)
  library(splines)
  library(emmeans)
  library(broom.mixed)
})

# Output directories
OUT_DIR <- here("output", "exp2_models")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Utility functions
msg <- function(text, level = 1) {
  char <- ifelse(level == 1, "=", "-")
  cat("\n", strrep(char, 80), "\n ", text, "\n", strrep(char, 80), "\n\n", sep = "")
}

# Shared setup: recode_outcome(), prepare_outcome_data(), fit_best_glmm(),
# SEPARATION_PRIORS, model_is_healthy(), DHARMA_SEED, theme_pub(), palette.
# Single source of truth shared with exp2_02..exp2_05.
source(here("scripts", "00_setup.R"))

# Load data
dat_porites <- read_excel(here("data", "Porites microscope characterization - complete.xlsx"), sheet = "data") %>%
  clean_names() %>%
  mutate(
    species = "Porites spp.",
    treatment = tolower(trimws(as.character(treatment))),
    treatment_label = factor(case_when(
      treatment == "airbrush" ~ "Airbrush",
      treatment %in% c("drem", "dremel") ~ "Dremel",
      treatment == "air_drem" ~ "Airbrush-Dremel"
    ), levels = c("Airbrush", "Dremel", "Airbrush-Dremel")),
    day = as.numeric(day),
    id = factor(id)
  ) %>%
  filter(!is.na(treatment_label), !is.na(day))

msg("COMPREHENSIVE GLMM ANALYSIS FOR ALL OUTCOMES", level = 1)

cat("Data loaded:\n")
cat("  Porites: N =", nrow(dat_porites), "obs,", n_distinct(dat_porites$id), "corals\n\n")

# =============================================================================
# Analysis function - tries multiple approaches
# =============================================================================

analyze_with_best_method <- function(data, outcome_col, outcome_name, species_name) {

  msg(paste(species_name, "-", outcome_name), level = 2)

  # Prepare data + fit the separation-robust ladder. fit_best_glmm() lives in
  # 00_setup.R and is the SAME function exp2_04 uses, so the model summarised
  # here for inference is byte-identical to the model DHARMa diagnoses there.
  dat_outcome <- prepare_outcome_data(data, outcome_col)

  if (nrow(dat_outcome) == 0) {
    cat("  [!] No data\n\n")
    return(NULL)
  }

  n_obs <- nrow(dat_outcome)
  n_corals <- n_distinct(dat_outcome$id)

  cat("  N obs:", n_obs, "| N corals:", n_corals, "\n")

  best_model <- fit_best_glmm(dat_outcome)

  # =========================================================================
  # Extract results from best model
  # =========================================================================

  if (!best_model$success) {
    cat("  WARNING: No method succeeded\n\n")
    return(tibble(
      Species = species_name,
      Outcome = outcome_name,
      N_Obs = n_obs,
      Method = "Failed",
      Status = "All methods failed"
    ))
  }

  cat("  -> Using:", best_model$method, "\n")

  # Get model summary
  model <- best_model$model
  model_summary <- summary(model)
  coef_table <- model_summary$coefficients$cond

  # ---------------------------------------------------------------------------
  # PRIMARY interaction result: JOINT test of the full treatment x time block.
  # (The per-coefficient Wald min-p below is retained only as a SECONDARY,
  # parameterization-dependent descriptor.)
  # ---------------------------------------------------------------------------
  joint <- joint_interaction_test(dat_outcome, best_model)
  cat("  -> Joint treatment x time test:", joint$joint_method,
      "| p =", signif(joint$joint_p, 4), "\n")

  # Extract treatment × time interactions
  interaction_rows <- grep("treatment_label.*ns\\(day_c", rownames(coef_table))

  sig_count <- 0L  # initialised so the return tibble is well-defined even
                   # when no interaction terms are present
  secondary_min_p <- NA_real_  # secondary, parameterization-dependent only
  if (length(interaction_rows) > 0) {
    interactions <- coef_table[interaction_rows, , drop = FALSE]
    sig_count <- sum(interactions[, "Pr(>|z|)"] < 0.05, na.rm = TRUE)
    secondary_min_p <- suppressWarnings(min(interactions[, "Pr(>|z|)"], na.rm = TRUE))

    cat("  -> Significant interactions:", sig_count, "of", length(interaction_rows), "\n")

    if (sig_count > 0) {
      sig_rows <- which(interactions[, "Pr(>|z|)"] < 0.05)
      for (i in sig_rows) {
        cat("     *", rownames(interactions)[i], ": p =",
            sprintf("%.4f", interactions[i, "Pr(>|z|)"]), "\n")
      }
    }
  }

  # Save full summary
  safe_species <- gsub(" ", "_", tolower(species_name))
  safe_species <- gsub("\\.", "", safe_species)  # Remove periods
  safe_outcome <- gsub(" ", "_", tolower(outcome_name))
  summary_file <- file.path(OUT_DIR, paste0(safe_species, "_", safe_outcome, ".txt"))

  sink(summary_file)
  cat("Species:", species_name, "\n")
  cat("Outcome:", outcome_name, "\n")
  cat("Method:", best_model$method, "\n")
  cat("N observations:", n_obs, "\n")
  cat("N corals:", n_corals, "\n\n")
  print(model_summary)
  sink()

  # Tidy fixed-effect coefficients: authoritative machine-readable source for
  # exp2_03 tables (replaces fragile text-scraping of the printed summary).
  tidy_fe <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
    mutate(Species = species_name, Outcome = outcome_name,
           Method = best_model$method, .before = 1)
  tidy_file <- file.path(OUT_DIR, paste0(safe_species, "_", safe_outcome, "_tidy.csv"))
  readr::write_csv(tidy_fe, tidy_file)

  cat("  [OK] Saved:", basename(summary_file), "+", basename(tidy_file), "\n\n")

  # Return summary
  tibble(
    Species = species_name,
    Outcome = outcome_name,
    N_Obs = n_obs,
    N_Corals = n_corals,
    Method = best_model$method,
    N_Interactions = length(interaction_rows),
    N_Sig_Interactions = sig_count,
    # PRIMARY interaction result
    Joint_Interaction_P = joint$joint_p,
    Joint_Method = joint$joint_method,
    Joint_Stat = joint$joint_stat,
    Joint_DF = joint$joint_df,
    # SECONDARY (parameterization-dependent) min Wald p over the block
    Secondary_Min_Wald_P = secondary_min_p,
    Status = "Success"
  )
}

# =============================================================================
# Run for all outcomes
# =============================================================================

msg("ANALYZING PORITES SPP.", level = 1)

PORITES_OUTCOMES <- list(
  # Composite "healed" (= coenosarc coverage + no algal plug + polyps in center)
  # dropped from the m = 6 analysis on 2026-05-25; it is redundant with the
  # three constituent outcomes plotted in Figure 3 and inflated the
  # multiple-testing denominator. The figure and the manuscript prose both
  # use the six constituent outcomes below.
  "Coenosarc Coverage" = "coenosarc_coverage",
  "Polyps in Center of Wound" = "polyp_in_center_of_wound",
  "Algal Plug" = "algal_plug",
  "Pink" = "pink",
  "RFP" = "rfp",
  "Yellow Aggregations" = "yellow_aggregations"
)

porites_results <- map_dfr(names(PORITES_OUTCOMES), function(outcome_name) {
  outcome_col <- PORITES_OUTCOMES[[outcome_name]]
  if (outcome_col %in% names(dat_porites)) {
    analyze_with_best_method(dat_porites, outcome_col, outcome_name, "Porites spp.")
  }
})

# =============================================================================
# Compile final summary
# =============================================================================

msg("FINAL SUMMARY", level = 1)

all_results <- porites_results

write_csv(all_results, file.path(OUT_DIR, "ALL_OUTCOMES_SUMMARY.csv"))

print(all_results, n = Inf, width = Inf)

cat("\n")
cat("OVERALL SUCCESS RATE:\n")
cat("  Total outcomes:", nrow(all_results), "\n")
cat("  Successful models:", sum(all_results$Status == "Success", na.rm = TRUE), "\n")
cat("  With significant interactions:", sum(all_results$N_Sig_Interactions > 0, na.rm = TRUE), "\n\n")

cat("METHODS USED:\n")
method_counts <- all_results %>%
  filter(Status == "Success") %>%
  count(Method)
print(method_counts)

cat("\n[OK] All model summaries + tidy CSVs saved to:", OUT_DIR, "\n")
cat("[OK] Summary table saved to: output/exp2_models/ALL_OUTCOMES_SUMMARY.csv\n\n")

msg("ANALYSIS COMPLETE - ALL OUTCOMES MODELED!", level = 1)
