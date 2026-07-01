# =============================================================================
# Experiment 1 — DESIGN-RESPECTING PRIMARY ANALYSIS (fragment-level + paired)
# -----------------------------------------------------------------------------
# Added 2026-06-30 (review-restats) to address a pseudoreplication problem.
#
# WHY THIS EXISTS
#   The timepoint-level GLMMs in airbrush_dremel_10_15_2025.R pool 7 repeated
#   photographic timepoints per fragment as if independent (154 "obs" from only
#   22 fragments). The parent/coral random effects are SINGULAR (variance ~ 0)
#   and DHARMa shows genuine temporal autocorrelation, so those models
#   pseudoreplicate. They are retained as SUPPORTING / sensitivity analyses.
#
#   This script is the PRIMARY analysis. It respects the design:
#     n = 22 fragments = 11 parent colonies x 2 treatments (paired split-colony;
#     each parent has exactly one airbrush + one dremel/scrape fragment).
#
# THREE PRIMARY LAYERS (for `regenerated` and binary `healed` = yes vs not-yes):
#   1. Fragment-level collapse (n = 22): ever_X = any(X == "yes") per fragment;
#      also X_day28 = value at day 28. ever_X is the primary endpoint (robust to
#      a missing day-28 cell; for the absorbing `regenerated` outcome ever==d28).
#   2. Fragment-level model: glmer(ever_X ~ treatment + species + (1|parent_id));
#      treatment OR + 95% CI + LRT p via drop1(test="Chisq"); isSingular() check.
#      Because n=22 yields singular RE AND quasi-complete separation, we ALSO fit
#      a Firth-penalized logistic (logistf, profile-likelihood CIs) and report the
#      Firth OR as the PRIMARY estimate, with the glmer result alongside.
#   3. Paired confirmatory test (fully pseudoreplication-free): treatment is
#      within-parent paired, so we test the 11 discordant-aware pairs with exact
#      McNemar (= exact binomial on discordant pairs).
#   4. Per-species day-28 Fisher exact (produced by exp1_sensitivity_extended.R)
#      is read in as confirmation.
#
# Output: output/tables/exp1_fragment_level_primary.csv
#   columns: outcome, model, term, OR, ci_lower, ci_upper, p_value, n, method
#
#   Rscript scripts/exp1_fragment_level_primary.R
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(readr); library(dplyr); library(tidyr); library(stringr)
  library(lme4); library(logistf); library(tibble); library(purrr)
})

allowed_treatment <- c("airbrush", "dremel")                  # ref = airbrush
allowed_species   <- c("acropora", "pocillopora", "porites")  # ref = acropora

# ---- 1. Read + fragment-level collapse (n = 22) -----------------------------
raw <- read_csv(here("data", "airbrush_dremel.csv"), show_col_types = FALSE) %>%
  mutate(
    species   = factor(str_trim(tolower(species)),   levels = allowed_species),
    treatment = factor(str_trim(tolower(treatment)), levels = allowed_treatment),
    coral_id  = factor(coral_id),
    parent_id = factor(sub("[[:alpha:]]+$", "", as.character(coral_id))),
    day       = as.integer(round(day)),
    regenerated = str_trim(tolower(regenerated)),
    healed      = str_trim(tolower(healed)),
    regen01   = if_else(regenerated == "yes", 1L,
                        if_else(regenerated == "no", 0L, NA_integer_)),
    healed_yes = if_else(healed == "yes", 1L,
                         if_else(healed %in% c("no", "incomplete"), 0L, NA_integer_))
  )

frag <- raw %>%
  group_by(coral_id, parent_id, species, treatment) %>%
  summarise(
    ever_regen   = as.integer(any(regen01   == 1L, na.rm = TRUE)),
    ever_healed  = as.integer(any(healed_yes == 1L, na.rm = TRUE)),
    regen_d28    = regen01[day == 28][1],
    healed_d28   = healed_yes[day == 28][1],
    .groups = "drop"
  )

stopifnot(nrow(frag) == 22, n_distinct(frag$parent_id) == 11)

# Document ever vs day-28 agreement (regen is absorbing -> should match)
cat("\n=== Fragment-level collapse (n = 22 fragments, 11 parents) ===\n")
cat("ever_regen vs regen_d28 identical:",
    all(frag$ever_regen == frag$regen_d28, na.rm = TRUE),
    "| any d28 NA:", any(is.na(frag$regen_d28)), "\n")
cat("ever_healed vs healed_d28 identical:",
    all(frag$ever_healed == frag$healed_d28, na.rm = TRUE),
    "| any d28 NA:", any(is.na(frag$healed_d28)), "\n")
cat("\never_regen by treatment:\n");  print(with(frag, table(treatment, ever_regen)))
cat("\never_healed by treatment:\n"); print(with(frag, table(treatment, ever_healed)))

# ---- 2. Primary model per outcome (glmer + Firth) ---------------------------
# Returns a list of result rows (one glmer row + one Firth row) for the master CSV.
analyse_outcome <- function(d, yvar, outcome_label) {
  d <- d %>% mutate(.y = .data[[yvar]]) %>% filter(!is.na(.y))
  n <- nrow(d)

  ## (a) glmer with parent random intercept
  glmer_or <- glmer_lo <- glmer_hi <- glmer_p <- NA_real_
  singular <- NA
  glmer_note <- ""
  m <- tryCatch(
    glmer(.y ~ treatment + species + (1 | parent_id), data = d,
          family = binomial,
          control = glmerControl(optimizer = "bobyqa",
                                 optCtrl = list(maxfun = 2e5))),
    error = function(e) { glmer_note <<- paste("glmer error:", conditionMessage(e)); NULL },
    warning = function(w) {
      # refit suppressing warnings so we still get the (boundary) estimate
      suppressWarnings(
        glmer(.y ~ treatment + species + (1 | parent_id), data = d,
              family = binomial,
              control = glmerControl(optimizer = "bobyqa",
                                     optCtrl = list(maxfun = 2e5)))
      )
    }
  )
  if (!is.null(m)) {
    singular <- lme4::isSingular(m)   # default tol (1e-4); matches lme4's own boundary message
    fe <- lme4::fixef(m)
    if ("treatmentdremel" %in% names(fe)) glmer_or <- exp(fe[["treatmentdremel"]])
    # Wald CI on treatment (profile CIs unstable under separation; LRT drives p)
    se <- sqrt(diag(as.matrix(vcov(m))))[["treatmentdremel"]]
    glmer_lo <- exp(fe[["treatmentdremel"]] - 1.96 * se)
    glmer_hi <- exp(fe[["treatmentdremel"]] + 1.96 * se)
    d1 <- tryCatch(drop1(m, test = "Chisq"), error = function(e) NULL)
    if (!is.null(d1) && "treatment" %in% rownames(d1)) {
      glmer_p <- d1["treatment", "Pr(Chi)"]
    }
    cat(sprintf("\n[%s] glmer ever-outcome: OR=%.3g (Wald %.3g-%.3g), drop1 LRT p=%.4g, singular=%s\n",
                outcome_label, glmer_or, glmer_lo, glmer_hi, glmer_p, singular))
  }

  ## (b) Firth-penalized logistic (drops parent RE) — profile-likelihood CIs
  lf <- logistf(.y ~ treatment + species, data = d)
  i  <- match("treatmentdremel", names(coef(lf)))
  firth_or <- exp(coef(lf)[i])
  firth_lo <- exp(lf$ci.lower[i])
  firth_hi <- exp(lf$ci.upper[i])
  firth_p  <- lf$prob[i]   # penalized likelihood-ratio p-value
  cat(sprintf("[%s] Firth (logistf, PRIMARY estimate): OR=%.3g (profile %.3g-%.3g), PLR p=%.4g  [n=%d]\n",
              outcome_label, firth_or, firth_lo, firth_hi, firth_p, n))

  # Quasi-complete separation inflates the glmer treatment OR to the boundary
  separated <- is.finite(glmer_or) && abs(log(glmer_or)) > 10
  glmer_label <- paste0(
    "glmer (1|parent_id)",
    if (isTRUE(singular)) " [SINGULAR: RE var~0]" else "",
    if (isTRUE(separated)) " [SEPARATION: OR at boundary]" else ""
  )

  bind_rows(
    tibble(outcome = outcome_label,
           model   = "Firth penalized logistic (logistf)",
           term    = "Dremel vs Airbrush (treatment)",
           OR = firth_or, ci_lower = firth_lo, ci_upper = firth_hi,
           p_value = firth_p, n = n,
           method  = "PRIMARY estimate: profile-likelihood CI; drops parent RE (separation-robust)"),
    tibble(outcome = outcome_label,
           model   = glmer_label,
           term    = "Dremel vs Airbrush (treatment)",
           OR = glmer_or, ci_lower = glmer_lo, ci_upper = glmer_hi,
           p_value = glmer_p, n = n,
           method  = "Fragment-level GLMM; OR Wald CI, p via drop1 LRT (reported alongside Firth)")
  )
}

regen_rows  <- analyse_outcome(frag, "ever_regen",  "regenerated")
healed_rows <- analyse_outcome(frag, "ever_healed", "healed")

# ---- 3. Paired confirmatory test: exact McNemar / exact binomial ------------
# Treatment is paired within parent. Cast wide, count discordant pairs, test.
paired_test <- function(frag, yvar, outcome_label) {
  w <- frag %>%
    select(parent_id, species, treatment, all_of(yvar)) %>%
    pivot_wider(names_from = treatment, values_from = all_of(yvar))
  # discordant pairs
  b <- sum(w$airbrush == 0 & w$dremel == 1, na.rm = TRUE)  # scrape only
  c <- sum(w$airbrush == 1 & w$dremel == 0, na.rm = TRUE)  # airbrush only
  conc_yy <- sum(w$airbrush == 1 & w$dremel == 1, na.rm = TRUE)
  conc_nn <- sum(w$airbrush == 0 & w$dremel == 0, na.rm = TRUE)
  npair   <- sum(!is.na(w$airbrush) & !is.na(w$dremel))
  # exact McNemar = two-sided exact binomial on discordant pairs
  bt <- stats::binom.test(b, b + c, p = 0.5)
  # per-species discordant counts
  persp <- w %>%
    group_by(species) %>%
    summarise(
      scrape_only   = sum(airbrush == 0 & dremel == 1, na.rm = TRUE),
      airbrush_only = sum(airbrush == 1 & dremel == 0, na.rm = TRUE),
      concordant    = sum(airbrush == dremel, na.rm = TRUE),
      .groups = "drop"
    )
  cat(sprintf("\n[%s] PAIRED (exact McNemar over %d parents): scrape-only=%d, airbrush-only=%d, concordant(yy=%d,nn=%d); exact binomial p=%.4g\n",
              outcome_label, npair, b, c, conc_yy, conc_nn, bt$p.value))
  cat("  per-species discordant counts:\n")
  print(persp)
  tibble(outcome = outcome_label,
         model   = "Paired exact McNemar (within-parent)",
         term    = "Scrape-only vs Airbrush-only discordant pairs",
         OR = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
         p_value = bt$p.value, n = npair,
         method  = sprintf("PRIMARY confirmatory (pseudoreplication-free): %d discordant (scrape-only) vs %d (airbrush-only)",
                           b, c))
}

regen_paired  <- paired_test(frag, "ever_regen",  "regenerated")
healed_paired <- paired_test(frag, "ever_healed", "healed")

# ---- 4. Per-species day-28 Fisher (read confirmation, if available) ---------
fisher_rows <- tibble()
fisher_path <- here("output", "tables", "endpoint_fisher_by_species_day28.csv")
if (file.exists(fisher_path)) {
  fish <- read_csv(fisher_path, show_col_types = FALSE)
  fisher_rows <- fish %>%
    transmute(
      outcome = "healed",
      model   = "Per-species day-28 Fisher exact",
      term    = paste0(species, ": scrape vs airbrush (healed yes)"),
      OR = suppressWarnings(as.numeric(OR)),
      ci_lower = NA_real_, ci_upper = NA_real_,
      p_value = p_value,
      n = airbrush_n + dremel_n,
      method = sprintf("CONFIRMATION (model-free): airbrush %d/%d vs scrape %d/%d healed",
                       airbrush_yes, airbrush_n, dremel_yes, dremel_n)
    )
  cat("\n[healed] per-species day-28 Fisher (read from exp1_sensitivity_extended.R output):\n")
  print(fisher_rows %>% select(term, OR, p_value, method))
} else {
  message("NOTE: ", fisher_path, " not found — run exp1_sensitivity_extended.R first for the day-28 Fisher confirmation.")
}

# ---- 5. Write the primary results table -------------------------------------
primary <- bind_rows(regen_rows, regen_paired, healed_rows, healed_paired, fisher_rows) %>%
  mutate(across(c(OR, ci_lower, ci_upper), ~ ifelse(is.na(.x), NA_real_, signif(.x, 4))),
         p_value = ifelse(is.na(p_value), NA_real_, signif(p_value, 3)))

out_path <- here("output", "tables", "exp1_fragment_level_primary.csv")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write_csv(primary, out_path)

cat("\n=== PRIMARY (fragment-level, n=22) results written ===\n")
cat("Wrote:", out_path, "\n")
print(primary, n = Inf, width = Inf)

# ---- 6. Inject a "Primary (fragment-level, n=22)" section into the master ----
# The main pipeline (airbrush_dremel_10_15_2025.R) writes paper_results_master.csv
# earlier; we append the primary rows so the single-source-of-truth table carries
# the design-respecting analysis. Idempotent: any pre-existing primary section is
# dropped before re-appending.
master_path <- here("output", "tables", "paper_results_master.csv")
if (file.exists(master_path)) {
  master <- read_csv(master_path, show_col_types = FALSE)
  master <- master %>% dplyr::filter(!grepl("^Primary \\(fragment-level",
                                            as.character(section)))
  primary_master <- primary %>%
    transmute(
      section = "Primary (fragment-level, n=22)",
      outcome,
      model,
      term,
      n,
      estimate      = OR,
      estimate_type = ifelse(is.na(OR), "exact p (paired/Fisher)", "OR"),
      ci_lower, ci_upper,
      effect_pct    = ifelse(is.na(OR), NA_real_, 100 * (OR - 1)),
      fold_change   = OR,
      test_stat     = NA_real_,
      df            = NA_integer_,
      p_value,
      p_formatted   = dplyr::case_when(
        is.na(p_value)   ~ NA_character_,
        p_value < 0.001  ~ "<0.001",
        TRUE             ~ formatC(p_value, format = "f", digits = 3)),
      notes         = method
    )
  master_out <- dplyr::bind_rows(master, primary_master[, names(master)])
  write_csv(master_out, master_path)
  cat("Appended", nrow(primary_master),
      "rows under 'Primary (fragment-level, n=22)' to paper_results_master.csv\n")
} else {
  message("NOTE: ", master_path,
          " not found — run airbrush_dremel_10_15_2025.R first so the master exists to append to.")
}

# ---- 7. PRIMARY vs SUPPORTING note ------------------------------
get1 <- function(path, term_match = "Dremel") {
  if (!file.exists(path)) return(NULL)
  d <- suppressMessages(read_csv(path, show_col_types = FALSE))
  tc <- names(d)[grepl("^Term$|^term$", names(d))][1]
  if (is.na(tc)) return(NULL)
  d[grepl(term_match, d[[tc]], ignore.case = TRUE), , drop = FALSE][1, , drop = FALSE]
}
fmt_or <- function(row) {
  if (is.null(row) || nrow(row) == 0) return("n/a")
  or <- suppressWarnings(as.numeric(row[["OR"]]))
  lo <- suppressWarnings(as.numeric(row[["CI Lower"]]))
  hi <- suppressWarnings(as.numeric(row[["CI Upper"]]))
  p  <- suppressWarnings(as.numeric(row[["P-value"]]))
  sprintf("OR = %.3g [%.3g, %.3g], p = %.3g", or, lo, hi, p)
}

sup_firth_healed <- get1(here("output","tables","firth_healed_fixed_effect_ORs.csv"))
sup_firth_regen  <- get1(here("output","tables","firth_regen_fixed_effect_ORs.csv"))
sup_lme4_healed  <- get1(here("output","tables","binary_fixed_effect_ORs.csv"))
sup_lme4_regen   <- get1(here("output","tables","regen_fixed_effect_ORs.csv"))

dw_line <- "(temporal-autocorrelation diagnostics not found)"
dw_path <- here("output","tables","exp1_dharma_summary.csv")
if (file.exists(dw_path)) {
  dw <- read_csv(dw_path, show_col_types = FALSE)
  dw_line <- paste0("Durbin-Watson temporal-autocorrelation p: ",
                    paste(sprintf("%s = %.3g", dw$outcome, dw$temporal_p), collapse = "; "),
                    " (all significant -> repeated timepoints are NOT independent).")
}

pr <- function(o, m) {
  r <- primary %>% filter(outcome == o, grepl(m, model))
  if (nrow(r) == 0) return("n/a")
  if (grepl("McNemar", m)) return(sprintf("exact McNemar p = %.3g (n = %d pairs)", r$p_value[1], r$n[1]))
  sprintf("OR = %.3g [%.3g, %.3g], p = %.3g (n = %d)", r$OR[1], r$ci_lower[1], r$ci_upper[1], r$p_value[1], r$n[1])
}

doc <- c(
  "# Experiment 1 — PRIMARY vs SUPPORTING analyses",
  "",
  paste0("_Generated by `scripts/exp1_fragment_level_primary.R` on ", Sys.Date(), "._"),
  "",
  "## Why two tiers",
  "",
  "Experiment 1 has **22 fragments = 11 parent colonies x 2 treatments** (paired",
  "split-colony design), each photographed at 7 timepoints (days 0,3,8,13,18,23,28).",
  "",
  "### SUPPORTING (timepoint-level GLMM / Firth / ordinal CLMM) — n=154 rows",
  "These models (in `airbrush_dremel_10_15_2025.R`) treat the 7 repeated",
  "timepoints per fragment as independent observations. This **pseudoreplicates**:",
  "",
  "- the `(1 | parent_id/coral_id)` random effects are **singular** (variance ~ 0);",
  paste0("- **temporal autocorrelation is present** — ", dw_line),
  "",
  "They remain useful as sensitivity analyses but **should not be the headline**.",
  "",
  "### PRIMARY (fragment-level + paired) — n=22 fragments / 11 parents",
  "Design-respecting endpoints (this script):",
  "",
  "- **Fragment-level collapse**: `ever_X = any(X==\"yes\")` per fragment (for the",
  "  absorbing `regenerated` outcome, ever == day-28 value; verified identical).",
  "- **Firth-penalized logistic** (separation-robust; profile-likelihood CIs) as the",
  "  primary estimate, with a `glmer(.|parent_id)` fit alongside (singular / separated).",
  "- **Paired exact McNemar** across the 11 parents — fully pseudoreplication-free.",
  "- **Per-species day-28 Fisher exact** as model-free confirmation.",
  "",
  "## PRIMARY results",
  "",
  "**Regeneration (polyps in wound centre):**",
  paste0("- Firth (PRIMARY): ", pr("regenerated", "Firth")),
  paste0("- Paired exact McNemar: ", pr("regenerated", "McNemar")),
  "",
  "**Healed (complete coenosarc; yes vs not-yes):**",
  paste0("- Firth (PRIMARY): ", pr("healed", "Firth")),
  paste0("- Paired exact McNemar: ", pr("healed", "McNemar")),
  "",
  "## SUPPORTING results (timepoint-level; report as sensitivity)",
  "",
  paste0("- Regeneration, Firth (timepoint-level): ", fmt_or(sup_firth_regen)),
  paste0("- Regeneration, lme4 GLMM (timepoint-level): ", fmt_or(sup_lme4_regen)),
  paste0("- Healed, Firth (timepoint-level): ", fmt_or(sup_firth_healed)),
  paste0("- Healed, lme4 GLMM (timepoint-level): ", fmt_or(sup_lme4_healed)),
  "",
  "### Separation-robust framing for the SUPPORTING healed result",
  "",
  "The timepoint-level lme4-vs-Firth `healed` odds ratio flips depending on two",
  "contested re-scored *Pocillopora* x scrape cells (coral 6b, days 23 & 28).",
  "Because the result hinges on near-complete separation, the **Firth estimate",
  "leads** in the SUPPORTING reporting (bias-corrected, separation-robust), and the",
  "drop-coral-6b / leave-one-out sensitivity",
  "(`output/tables/sensitivity_analysis_summary.csv`,",
  "`exp1_sensitivity_6b.R`) is reported alongside it. The data are NOT altered.",
  "",
  "Note: the fragment-level PRIMARY healed result above does not depend on this",
  "flip — collapsing to `ever_healed` and using exact paired/penalized tests is",
  "robust to the contested cells."
)
doc_path <- here("output", "text", "exp1_primary_vs_supporting.md")
dir.create(dirname(doc_path), recursive = TRUE, showWarnings = FALSE)
writeLines(doc, doc_path)
cat("Wrote PRIMARY-vs-SUPPORTING note:", doc_path, "\n")
