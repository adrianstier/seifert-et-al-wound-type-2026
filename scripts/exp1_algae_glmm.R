# =============================================================================
# Experiment 1 — Algae (debris) outcome GLMM
# -----------------------------------------------------------------------------
# Added 2026-05-19. Per coauthor decision, Experiment 1 is analyzed on
# two response variables: `regenerated` (handled in airbrush_dremel_10_15_2025.R
# section 5c) and ALGAE colonization (`debris`). This script fits the algae
# model with EXACTLY the parallel structure used for `regenerated` (sections
# 5c + 5e of airbrush_dremel_10_15_2025.R): additive vs interaction GLMM (lme4),
# LRT, drop1 main-effect LRTs, emmeans cell probabilities and the
# treatment-within-species contrast, plus a Firth (brglm2) fixed-effects OR
# table. Outputs are written as algae_*/firth_algae_* in output/tables/,
# mirroring the regen_*/firth_regen_* filenames. Nothing else is modified.
#
#   Rscript scripts/exp1_algae_glmm.R
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(readr); library(dplyr); library(stringr)
  library(lme4); library(broom.mixed); library(emmeans); library(brglm2)
})

allowed_treatment <- c("airbrush", "dremel")              # ref = airbrush
allowed_species   <- c("acropora", "pocillopora", "porites")  # ref = acropora
allowed_yn        <- c("no", "yes")

df <- read_csv(here("data", "airbrush_dremel.csv"), show_col_types = FALSE) %>%
  mutate(
    # Defensive trim/lower (matches main script): guards the known Pocillopora
    # trailing-whitespace issue that otherwise drops rows to NA on factor().
    treatment = factor(str_trim(tolower(treatment)), levels = allowed_treatment),
    species   = factor(str_trim(tolower(species)),   levels = allowed_species),
    coral_id  = factor(coral_id),
    debris    = factor(str_trim(tolower(debris)),    levels = allowed_yn),
    parent_id = factor(sub("[[:alpha:]]+$", "", as.character(coral_id))),
    debris01  = if_else(debris == "yes", 1L,
                        if_else(debris == "no", 0L, NA_integer_))
  ) %>%
  filter(!is.na(debris01))

cat("\n--- Algae (debris) GLMM ---\n")
cat("N observations (non-NA): ", nrow(df), "\n", sep = "")
cat("Yes / No counts: ", sum(df$debris01 == 1L), " / ",
    sum(df$debris01 == 0L), "\n", sep = "")
cat("By treatment:\n"); print(with(df, table(treatment, debris01)))

mod_add <- glmer(
  debris01 ~ treatment + species + (1 | parent_id/coral_id),
  data = df, family = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
  na.action = na.exclude
)
mod_int <- glmer(
  debris01 ~ treatment * species + (1 | parent_id/coral_id),
  data = df, family = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
  na.action = na.exclude
)

# LRT interaction vs additive
lrt <- anova(mod_add, mod_int, test = "Chisq")
write_csv(
  tibble(Term = "Treatment × Species Interaction",
         Chi_square = lrt$Chisq[2], DF = lrt$Df[2],
         P_value = lrt$`Pr(>Chisq)`[2]),
  here("output", "tables", "algae_lrt_interaction_vs_additive.csv")
)

# drop1 main effects (additive)
d1 <- as_tibble(drop1(mod_add, test = "Chisq"), rownames = "Term") %>%
  filter(Term %in% c("treatment", "species")) %>%
  transmute(Term, DF = npar, `Chi-square` = LRT, `P-value` = `Pr(Chi)`)
write_csv(d1, here("output", "tables", "algae_glmm_drop1_additive.csv"))

# emmeans: marginal Pr(algae) by species × treatment
emm_cells <- emmeans(mod_add, ~ treatment * species, type = "response")
write_csv(as.data.frame(emm_cells),
          here("output", "tables", "algae_emmeans_cell_probabilities.csv"))

# treatment-within-species contrast (Dremel vs Airbrush), probability scale
emm_t <- emmeans(mod_add, ~ treatment | species, type = "response")
write_csv(as.data.frame(contrast(emm_t, "revpairwise")),
          here("output", "tables", "algae_emmeans_treatment_within_species.csv"))

# Firth-penalized fixed-effects logistic (parallel to firth_regen)
firth_or_table <- function(mod, label) {
  co <- summary(mod)$coefficients
  rows <- intersect(c("treatmentdremel", "speciespocillopora", "speciesporites"),
                     rownames(co))
  # PROFILE-LIKELIHOOD CIs (penalized-likelihood profiling via
  # confint() on the brglm2 fit), replacing the old Wald exp(est +/- z*SE).
  # Estimates and p-values UNCHANGED; only CI bounds differ.
  ci <- suppressMessages(suppressWarnings(stats::confint(mod)))
  tibble(
    Term = dplyr::case_match(rows,
                             "treatmentdremel"    ~ "Dremel vs Airbrush",
                             "speciespocillopora" ~ "Pocillopora vs Acropora",
                             "speciesporites"     ~ "Porites vs Acropora"),
    OR         = exp(co[rows, "Estimate"]),
    `CI Lower` = exp(ci[rows, 1]),
    `CI Upper` = exp(ci[rows, 2]),
    `P-value`  = co[rows, "Pr(>|z|)"],
    outcome    = label
  )
}
mod_firth <- glm(debris01 ~ treatment + species, data = df,
                 family = binomial, method = "brglmFit")
firth_algae <- firth_or_table(mod_firth, "algae")
write_csv(firth_algae, here("output", "tables", "firth_algae_fixed_effect_ORs.csv"))

cat("\n--- Firth (brglm2) — algae ---\n"); print(firth_algae)
cat("Convergence:", mod_firth$converged,
    "| singular(add):", lme4::isSingular(mod_add, tol = 1e-5), "\n")

# -----------------------------------------------------------------------------
# Per-genus Firth ORs (Dremel vs Airbrush within each species).
# Added 2026-05-26 to address reviewer concern that the additive Firth
# OR (which collapses across genera) hides per-genus heterogeneity for algae.
# The treatment × species interaction LRT was significant (χ² = 11.0, df = 2,
# P = 0.004), so we report per-genus Firth ORs as a sensitivity confirming
# directionality is preserved across all three genera even if magnitudes differ.
# Output: firth_algae_per_genus_ORs.csv  (columns: Species, Term, OR, CI Lower,
# CI Upper, P-value, N)
# -----------------------------------------------------------------------------
firth_per_genus <- function(species_name) {
  df_sub <- df %>% filter(species == species_name)
  mod <- glm(debris01 ~ treatment, data = df_sub,
             family = binomial, method = "brglmFit")
  co <- summary(mod)$coefficients
  if ("treatmentdremel" %in% rownames(co)) {
    est <- co["treatmentdremel", "Estimate"]
    pv  <- co["treatmentdremel", "Pr(>|z|)"]
    # profile-likelihood CI (penalized) instead of Wald.
    ci  <- suppressMessages(suppressWarnings(stats::confint(mod)))["treatmentdremel", ]
    tibble(
      Species    = species_name,
      Term       = "Dremel vs Airbrush",
      OR         = exp(est),
      `CI Lower` = exp(ci[1]),
      `CI Upper` = exp(ci[2]),
      `P-value`  = pv,
      N          = nrow(df_sub)
    )
  } else NULL
}
firth_algae_per_genus <- bind_rows(
  firth_per_genus("acropora"),
  firth_per_genus("pocillopora"),
  firth_per_genus("porites")
)
write_csv(firth_algae_per_genus,
          here("output", "tables", "firth_algae_per_genus_ORs.csv"))
cat("\n--- Firth (brglm2) — algae, per-genus (Dremel vs Airbrush) ---\n")
print(firth_algae_per_genus)

cat("\nWrote: algae_lrt_interaction_vs_additive.csv, algae_glmm_drop1_additive.csv,\n",
    "       algae_emmeans_cell_probabilities.csv, algae_emmeans_treatment_within_species.csv,\n",
    "       firth_algae_fixed_effect_ORs.csv, firth_algae_per_genus_ORs.csv\n", sep = "")
