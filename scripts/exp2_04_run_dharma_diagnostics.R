################################################################################
# EXPERIMENT 2: DHARMa MODEL DIAGNOSTICS
#
# Purpose: Run comprehensive DHARMa diagnostics for all Experiment 2 GLMMs
#
# DHARMa checks:
#   - Residual patterns (QQ plot, residuals vs predicted)
#   - Overdispersion
#   - Zero-inflation
#   - Outlier detection
#
# Output:
#   - Diagnostic figures: output/exp2_diagnostics/
#   - Summary table: output/exp2_diagnostics/dharma_summary.csv
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(here)
  library(glmmTMB)
  library(splines)
  library(DHARMa)
  library(patchwork)
})

# Clean ggplot theme for the diagnostic panels (matches exp1_dharma_diagnostics.R).
# Replaces DHARMa's base-graphics plot(), whose stacked titles ("Within-group
# deviation from uniformity", "Levene Test ...") overlap when several models are
# composed in one device.
theme_dharma <- function(base_size = 9) ggplot2::theme_bw(base_size = base_size) +
  ggplot2::theme(panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 plot.title = element_text(size = rel(0.95), face = "bold"),
                 plot.margin = margin(4, 6, 4, 6))

# Build a clean QQ-uniformity + residual-vs-predicted pair for one fitted model.
dharma_panel_pair <- function(sim_res, title) {
  r <- sim_res$scaledResiduals
  n <- length(r)
  qq <- tibble(expected = qunif(ppoints(n)), observed = sort(r))
  rp <- tibble(pred = sim_res$fittedPredictedResponse, res = r)
  p1 <- ggplot(qq, aes(expected, observed)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = 2) +
    geom_point(size = 0.5, alpha = 0.5, colour = "#0072B2") +
    labs(x = "Expected (uniform)", y = "Observed residual",
         title = paste0(title, ": QQ uniformity")) +
    theme_dharma()
  p2 <- ggplot(rp, aes(pred, res)) +
    geom_hline(yintercept = c(0.25, 0.5, 0.75), colour = "grey85", linetype = 2) +
    geom_point(size = 0.5, alpha = 0.4, colour = "#0072B2") +
    geom_smooth(method = "loess", se = FALSE, colour = "#D55E00",
                linewidth = 0.6, na.rm = TRUE) +
    labs(x = "Predicted probability", y = "Scaled residual",
         title = paste0(title, ": residual vs. predicted")) +
    theme_dharma()
  p1 | p2
}

# Create output directory
OUT_DIR <- here("output", "exp2_diagnostics")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("\n", rep("=", 80), "\n")
cat(" EXPERIMENT 2: DHARMa MODEL DIAGNOSTICS\n")
cat(rep("=", 80), "\n\n")

# Shared setup: recode_outcome(), prepare_outcome_data(), fit_best_glmm(),
# DHARMA_SEED. Sourcing this guarantees DHARMa diagnoses the SAME model
# (standard / penalized-prior / OLRE) that exp2_01 used for inference.
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

# This pipeline is Porites-only.

# Run diagnostics for one outcome
run_diagnostics <- function(data, outcome_col, outcome_name, species_name) {

  cat("\n", rep("-", 80), "\n")
  cat(" ", species_name, "-", outcome_name, "\n")
  cat(rep("-", 80), "\n\n")

  # Prepare data with the shared helper — identical recoding/centering to
  # exp2_01 (no inert weight column).
  dat_outcome <- prepare_outcome_data(data, outcome_col)

  if (nrow(dat_outcome) < 20) {
    cat("  [!] Insufficient data (n =", nrow(dat_outcome), ")\n")
    return(NULL)
  }

  # Fit the SAME separation-robust model exp2_01 used for inference, so the
  # residual diagnostics describe the model the manuscript actually reports
  # (the old code always diagnosed the naive standard GLMM instead).
  best <- fit_best_glmm(dat_outcome, verbose = FALSE)

  if (!isTRUE(best$success)) {
    cat("  [x] Model failed to fit (all approaches)\n")
    return(list(summary = tibble(
      Species = species_name, Outcome = outcome_name, N = nrow(dat_outcome),
      Method = "Failed", Status = "Model failed",
      Uniformity_p = NA_real_, Dispersion_p = NA_real_,
      Outliers_p = NA_real_, Pass = FALSE
    ), panel = NULL, outcome = outcome_name))
  }

  m <- best$model
  cat("  -> Diagnosing:", best$method, "\n")

  # DHARMa simulation. Seed fixed in 00_setup.R so the Monte-Carlo p-values
  # are reproducible run-to-run — without this they drifted (the documented
  # Exp 2 p-value instability that made docs disagree, e.g. Pink 0.021/0.062).
  set.seed(DHARMA_SEED)
  sim_res <- simulateResiduals(m, n = 1000, plot = FALSE)

  uniformity_test <- testUniformity(sim_res, plot = FALSE)
  dispersion_test <- testDispersion(sim_res, plot = FALSE)
  outlier_test    <- testOutliers(sim_res, plot = FALSE)

  # Clean ggplot diagnostic panel (replaces the raw base-graphics plot whose
  # stacked titles overlapped). Per-outcome PNG plus a contribution to the
  # combined figure assembled after the loop.
  # Short outcome label keeps the right-hand panel title from truncating; the
  # fitted model/method is given in the figure caption rather than every title.
  fig_label <- dplyr::case_match(outcome_name,
                                 "Polyps in Center of Wound" ~ "Polyp regeneration",
                                 "Yellow Aggregations" ~ "Yellow aggregations",
                                 "Coenosarc Coverage" ~ "Coenosarc coverage",
                                 "Algal Plug" ~ "Algal plug",
                                 "Pink" ~ "Pink aggregations",
                                 .default = outcome_name)
  panel <- dharma_panel_pair(sim_res, fig_label)

  safe_species <- gsub(" ", "_", tolower(species_name))
  safe_species <- gsub("\\.", "", safe_species)
  safe_outcome <- gsub(" ", "_", tolower(outcome_name))
  plot_file <- file.path(OUT_DIR, paste0(safe_species, "_", safe_outcome, ".pdf"))
  ggsave(plot_file, panel, width = 180, height = 70, units = "mm",
         device = cairo_pdf, bg = "white")

  # NOTE: these tests are reported DESCRIPTIVELY. A non-significant test is
  # NOT proof of adequate fit (low power at small N). `no_flag` means simply
  # that no DHARMa test flagged a deviation; it is NOT an inferential gate and
  # excludes no outcome from the manuscript.
  no_flag <- uniformity_test$p.value > 0.05 &&
             dispersion_test$p.value > 0.05 &&
             outlier_test$p.value > 0.05

  cat(sprintf("  -> Uniformity p=%.4f | Dispersion p=%.4f | Outlier p=%.4f%s\n",
              uniformity_test$p.value, dispersion_test$p.value,
              outlier_test$p.value,
              if (!no_flag) "  [flag: inspect residual plot]" else ""))
  cat("  -> Diagnostic plot:", basename(plot_file), "\n")

  list(
    summary = tibble(
      Species = species_name,
      Outcome = outcome_name,
      N = nrow(dat_outcome),
      Method = best$method,
      Status = "Success",
      Uniformity_p = uniformity_test$p.value,
      Dispersion_p = dispersion_test$p.value,
      Outliers_p = outlier_test$p.value,
      Pass = no_flag   # descriptive flag (see NOTE), NOT a pass/fail gate
    ),
    panel = panel,
    outcome = outcome_name
  )
}

# Run for all Porites outcomes
cat("\n", rep("=", 80), "\n")
cat(" PORITES SPP.\n")
cat(rep("=", 80), "\n")

PORITES_OUTCOMES <- list(
  # Six reported Exp 2 outcomes (matches exp2_01 inference and Table S10); the
  # composite "healed" is excluded — it is redundant with coenosarc coverage and
  # is not one of the six outcomes the manuscript analyses.
  "Coenosarc Coverage" = "coenosarc_coverage",
  "Polyps in Center of Wound" = "polyp_in_center_of_wound",
  "Algal Plug" = "algal_plug",
  "Pink" = "pink",
  "RFP" = "rfp",
  "Yellow Aggregations" = "yellow_aggregations"
)

porites_runs <- map(names(PORITES_OUTCOMES), function(outcome_name) {
  outcome_col <- PORITES_OUTCOMES[[outcome_name]]
  if (outcome_col %in% names(dat_porites)) {
    run_diagnostics(dat_porites, outcome_col, outcome_name, "Porites spp.")
  } else NULL
}) %>% compact()

# Combine and save results
all_results <- map_dfr(porites_runs, "summary")
write_csv(all_results, file.path(OUT_DIR, "dharma_summary.csv"))

# Combined clean ggplot diagnostic figure: one [QQ | residual-vs-predicted] row
# per successfully fitted outcome. Replaces the stacked raw base-graphics device
# whose per-model titles overlapped.
panels <- map(porites_runs, "panel") %>% compact()
if (length(panels) > 0) {
  combined <- wrap_plots(panels, ncol = 1) +
    plot_annotation(
      title = "Experiment 2 (Porites) GLMMs: DHARMa simulated residuals (seeded, n = 1000)",
      theme = theme(plot.title = element_text(size = 11, face = "bold")))
  fig_h <- 18 + 32 * length(panels)   # mm; ~32 mm per outcome row + title
  ggsave(file.path(OUT_DIR, "exp2_dharma_residuals.png"), combined,
         width = 180, height = fig_h, units = "mm", dpi = 300, bg = "white", limitsize = FALSE)
  ggsave(file.path(OUT_DIR, "exp2_dharma_residuals.pdf"), combined,
         width = 180, height = fig_h, units = "mm", device = cairo_pdf, bg = "white", limitsize = FALSE)
  cat("OK Combined clean diagnostic figure: exp2_dharma_residuals.{png,pdf}\n")
}

# Print summary
cat("\n", rep("=", 80), "\n")
cat(" DIAGNOSTIC SUMMARY\n")
cat(rep("=", 80), "\n\n")

print(all_results %>% select(Species, Outcome, N, Pass, Uniformity_p, Dispersion_p), n = Inf)

cat("\n")
cat("Overall Results:\n")
cat("  Total models tested:", nrow(all_results), "\n")
cat("  Models passing diagnostics:", sum(all_results$Pass, na.rm = TRUE), "\n")
cat("  Pass rate:", sprintf("%.1f%%", 100 * mean(all_results$Pass, na.rm = TRUE)), "\n\n")

cat("OK Diagnostic plots saved to:", OUT_DIR, "\n")
cat("OK Summary table saved to:", file.path(OUT_DIR, "dharma_summary.csv"), "\n\n")

cat(rep("=", 80), "\n")
cat(" DHARMa DIAGNOSTICS COMPLETE\n")
cat(rep("=", 80), "\n\n")
