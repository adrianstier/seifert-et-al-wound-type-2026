# =============================================================================
# Experiment 2 — Algae spline sensitivity (B9) + per-outcome ICC (B10)
# -----------------------------------------------------------------------------
# Added 2026-05-26.
#
# B9: the headline algae model uses ns(day_c, 3) but the trajectories are
# non-monotonic (algae appear, persist, partially retreat). Re-fit with
# ns(day_c, 4) and ns(day_c, 5) and report whether the treatment × time
# interaction p-value remains non-significant.
#
# B10: the manuscript reports treatment × ns(day_c, 3) + (1|id) for the six
# Porites outcomes but does not report the latent-scale ICC of the colony
# random intercept. This script computes ICC_latent = var(random) / (var(random)
# + π²/3) for each outcome's penalized binomial GLMM (the headline model under
# separation, per scripts/00_setup.R::fit_best_glmm).
#
# Outputs:
#   output/exp2_tables/algae_spline_df_sensitivity.csv
#   output/exp2_tables/porites_icc_latent.csv
#
#   Rscript scripts/exp2_07_spline_sensitivity_and_icc.R
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(readxl); library(dplyr); library(readr)
  library(glmmTMB); library(splines); library(broom.mixed); library(purrr); library(tibble)
})

source(here("scripts", "00_setup.R"))

xlsx_path <- here("data", "Porites microscope characterization - complete.xlsx")
df <- read_excel(xlsx_path, sheet = "data") %>%
  mutate(
    id    = as.character(id),
    day   = as.integer(day),
    day_c = as.numeric(day),
    treatment_label = dplyr::case_match(treatment,
                             "airbrush" ~ "Airbrush",
                             "dremel" ~ "Dremel",
                             "air_drem" ~ "Airbrush-Dremel",
                             .default = treatment),
    treatment_label = factor(treatment_label,
                             levels = c("Airbrush", "Dremel", "Airbrush-Dremel"))
  ) %>%
  filter(!is.na(treatment_label), !is.na(day))

# Map of outcome label → data column
PORITES_OUTCOMES <- list(
  "Coenosarc Coverage"        = "coenosarc_coverage",
  "Polyps in Center of Wound" = "polyp_in_center_of_wound",
  "Algal Plug"                = "algal_plug",
  "Yellow Aggregations"       = "yellow_aggregations",
  "Pink"                      = "pink",
  "RFP"                       = "rfp"
)

# -----------------------------------------------------------------------------
# B9: algae spline sensitivity
# -----------------------------------------------------------------------------
cat("=== B9: Algae spline sensitivity ===\n")
algae_data <- df %>%
  filter(!is.na(algal_plug)) %>%
  mutate(value = if_else(tolower(trimws(algal_plug)) == "yes", 1L, 0L))

cat("N obs:", nrow(algae_data), "| N corals:", n_distinct(algae_data$id), "\n")

spline_p <- function(df_use) {
  out <- list()
  for (k in c(3, 4, 5)) {
    f <- as.formula(paste0("value ~ treatment_label * ns(day_c, ", k, ") + (1 | id)"))
    m <- tryCatch(
      glmmTMB(f, data = df_use, family = binomial(), priors = SEPARATION_PRIORS),
      error = function(e) { cat("  df =", k, ": FAILED:", conditionMessage(e), "\n"); NULL }
    )
    if (is.null(m)) { out[[as.character(k)]] <- NA_real_; next }
    co <- summary(m)$coefficients$cond
    inter_rows <- grep("treatment_label.*ns\\(day_c", rownames(co))
    if (length(inter_rows) == 0) { out[[as.character(k)]] <- NA_real_; next }
    # SECONDARY descriptor only: min of the per-coefficient interaction Wald p
    # across spline-df settings (a robustness check on the headline df = 3
    # choice; the PRIMARY joint treatment x time test lives in exp2_01/03).
    interaction_ps <- co[inter_rows, "Pr(>|z|)"]
    min_p <- min(interaction_ps, na.rm = TRUE)
    out[[as.character(k)]] <- list(min_interaction_p = min_p,
                                    n_inter_terms = length(inter_rows),
                                    n_significant = sum(interaction_ps < 0.05, na.rm = TRUE))
    cat("  ns(day_c, ", k, "): min interaction p = ", sprintf("%.4f", min_p),
        " | n_sig = ", out[[as.character(k)]]$n_significant,
        " of ", length(inter_rows), "\n", sep = "")
  }
  out
}

algae_sens <- spline_p(algae_data)
spline_df_table <- tibble(
  outcome = "Algal Plug",
  spline_df = c(3, 4, 5),
  min_interaction_p = sapply(spline_df, function(k) algae_sens[[as.character(k)]]$min_interaction_p),
  n_interaction_terms = sapply(spline_df, function(k) algae_sens[[as.character(k)]]$n_inter_terms),
  n_significant_p05 = sapply(spline_df, function(k) algae_sens[[as.character(k)]]$n_significant)
)

dir.create(here("output", "exp2_tables"), recursive = TRUE, showWarnings = FALSE)
write_csv(spline_df_table, here("output", "exp2_tables", "algae_spline_df_sensitivity.csv"))
cat("\nWrote: output/exp2_tables/algae_spline_df_sensitivity.csv\n")
print(spline_df_table)

# -----------------------------------------------------------------------------
# B9b: spline-df sensitivity for ALL six Porites outcomes (added 2026-06-19).
# Reviewers may ask whether the headline ns(day_c, 3) choice drives the time
# effects; refit each outcome at df = 3, 4, 5 and report the min interaction p.
# -----------------------------------------------------------------------------
cat("\n=== B9b: spline-df sensitivity, all six outcomes ===\n")
all_spline <- imap_dfr(PORITES_OUTCOMES, function(col, nm) {
  dat <- df %>% filter(!is.na(.data[[col]])) %>%
    mutate(value = if_else(tolower(trimws(.data[[col]])) == "yes", 1L, 0L))
  s <- spline_p(dat)
  map_dfr(c(3, 4, 5), function(k) {
    r <- s[[as.character(k)]]
    if (is.list(r)) tibble(outcome = nm, spline_df = k, min_interaction_p = r$min_interaction_p,
                           n_interaction_terms = r$n_inter_terms, n_significant_p05 = r$n_significant)
    else tibble(outcome = nm, spline_df = k, min_interaction_p = NA_real_,
                n_interaction_terms = NA_integer_, n_significant_p05 = NA_integer_)
  })
})
write_csv(all_spline, here("output", "exp2_tables", "spline_df_sensitivity_all_outcomes.csv"))
cat("Wrote: output/exp2_tables/spline_df_sensitivity_all_outcomes.csv\n")
print(all_spline, n = Inf)

# -----------------------------------------------------------------------------
# B10: per-outcome latent ICC for the colony random intercept
# -----------------------------------------------------------------------------
cat("\n=== B10: Per-outcome ICC (latent scale) ===\n")
icc_rows <- list()

for (outcome_name in names(PORITES_OUTCOMES)) {
  outcome_col <- PORITES_OUTCOMES[[outcome_name]]
  dat <- df %>%
    filter(!is.na(.data[[outcome_col]])) %>%
    mutate(value = if_else(tolower(trimws(.data[[outcome_col]])) == "yes", 1L, 0L))
  if (nrow(dat) == 0) next
  fitres <- fit_best_glmm(dat, verbose = FALSE)
  if (!isTRUE(fitres$success)) {
    cat("  ", outcome_name, ": model failed; skipping ICC\n", sep = "")
    next
  }
  vc <- VarCorr(fitres$model)$cond
  # Extract variance of the 'id' random intercept
  var_id <- NA_real_
  for (nm in names(vc)) {
    re_mat <- attr(vc[[nm]], "stddev")
    if (!is.null(re_mat)) {
      sd_id <- as.numeric(re_mat[1])
      if (nm == "id" || nm == "id:obs_id") var_id <- sd_id^2
    }
  }
  # Latent ICC for logit GLMM: var_id / (var_id + pi^2/3)
  icc <- if (!is.na(var_id)) var_id / (var_id + pi^2 / 3) else NA_real_
  icc_rows[[outcome_name]] <- tibble(
    Outcome = outcome_name,
    Method = fitres$method,
    N_obs = nrow(dat),
    N_colonies = n_distinct(dat$id),
    Var_id_latent = var_id,
    ICC_latent = icc
  )
  cat("  ", outcome_name, ": var(id) = ",
      sprintf("%.4f", var_id), ", ICC_latent = ", sprintf("%.4f", icc),
      " (method: ", fitres$method, ")\n", sep = "")
}

icc_tbl <- bind_rows(icc_rows)
write_csv(icc_tbl, here("output", "exp2_tables", "porites_icc_latent.csv"))
cat("\nWrote: output/exp2_tables/porites_icc_latent.csv\n")
print(icc_tbl)
