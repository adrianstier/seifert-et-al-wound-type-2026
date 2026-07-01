# =============================================================================
# exp1_multiplicity.R — multiple-testing sensitivity for Exp 1 main effects
# -----------------------------------------------------------------------------
# Exp 1 hypotheses are a-priori and directional (scrape > airbrush; genera
# differ), so unadjusted p-values are the primary inference. This produces a
# defensive BH-FDR / Bonferroni table over the primary main-effect LRTs so the
# manuscript can state the treatment effect survives correction.
# =============================================================================
suppressMessages({library(tidyverse); library(here)})
rd <- function(p) read_csv(here("output","tables",p), show_col_types = FALSE)

grab <- function(file, outcome) {
  d <- rd(file)
  tibble(outcome = outcome,
         term = d$Term,
         chi2 = d$`Chi-square`,
         p_raw = d$`P-value`)
}
fam <- bind_rows(
  grab("binary_glmm_drop1_additive.csv","healed"),
  grab("regen_glmm_drop1_additive.csv","regenerated"),
  grab("algae_glmm_drop1_additive.csv","algae/debris")
) %>%
  mutate(p_BH = p.adjust(p_raw, "BH"),
         p_Bonferroni = p.adjust(p_raw, "bonferroni"),
         survives_Bonferroni = p_Bonferroni < 0.05,
         survives_BH = p_BH < 0.05)
write_csv(fam, here("output","tables","multiplicity_summary_all.csv"))
cat("=== Exp 1 main-effect multiple-testing sensitivity (family of", nrow(fam), "tests) ===\n")
fam %>% mutate(across(c(p_raw,p_BH,p_Bonferroni), ~signif(.x,3))) %>% print(n=Inf, width=Inf)
cat("\nWrote: multiplicity_summary_all.csv\n")
