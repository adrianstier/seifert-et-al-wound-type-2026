# =============================================================================
# exp1_sensitivity_extended.R — robustness battery for the Exp 1 healing effect
# -----------------------------------------------------------------------------
# Complements exp1_sensitivity_6b.R. Adds:
#   (a) coral-7b exclusion (incomplete/nonuniform day-0 wound);
#   (b) missing-data scenarios for the 3 healed NAs (complete-case / best / worst);
#   (c) day-28 Fisher exact tests per species (model-free; the day-28 GLMM is
#       degenerate under separation, so exact tests are the correct tool);
#   (d) a master sensitivity summary stacking 6b + 7b + missing-data scenarios.
# =============================================================================
suppressMessages({library(tidyverse); library(here); library(lme4); library(brglm2)})

prep <- function(df) df %>%
  mutate(species   = factor(str_trim(tolower(species)), levels=c("acropora","pocillopora","porites")),
         treatment = factor(str_trim(tolower(treatment)), levels=c("airbrush","dremel")),
         parent_id = factor(sub("[[:alpha:]]+$","", as.character(coral_id))),
         coral_id  = factor(coral_id))
raw <- prep(read_csv(here("data","airbrush_dremel.csv"), show_col_types=FALSE))

# healed treatment effect for a given healed01 vector + data
fit_healed <- function(d) {
  d <- d %>% filter(!is.na(healed01))
  ctrl <- glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
  m   <- glmer(healed01 ~ treatment + species + (1|parent_id/coral_id), d, family=binomial, control=ctrl)
  mnt <- glmer(healed01 ~ species + (1|parent_id/coral_id), d, family=binomial, control=ctrl)
  mns <- glmer(healed01 ~ treatment + (1|parent_id/coral_id), d, family=binomial, control=ctrl)
  fm  <- glm(healed01 ~ treatment + species, d, family=binomial, method="brglmFit", type="AS_mean")
  tibble(lme4_treatment_OR = exp(fixef(m)[["treatmentdremel"]]),
         firth_treatment_OR = exp(coef(fm)[["treatmentdremel"]]),
         treatment_LRT_p = anova(mnt,m)$`Pr(>Chisq)`[2],
         species_LRT_chi2 = anova(mns,m)$Chisq[2], species_LRT_p = anova(mns,m)$`Pr(>Chisq)`[2],
         n = nrow(d))
}

# ---- (a)+(b) scenario set on healed ----
scen <- list(
  `full (complete-case)` = raw %>% mutate(healed01 = if_else(healed=="yes",1L,0L)),
  `drop coral 6b`        = raw %>% filter(coral_id!="6b") %>% mutate(healed01 = if_else(healed=="yes",1L,0L)),
  `drop coral 7b`        = raw %>% filter(coral_id!="7b") %>% mutate(healed01 = if_else(healed=="yes",1L,0L)),
  `missing healed = yes (best)`  = raw %>% mutate(healed01 = if_else(is.na(healed),1L, if_else(healed=="yes",1L,0L))),
  `missing healed = no (worst)`  = raw %>% mutate(healed01 = if_else(is.na(healed),0L, if_else(healed=="yes",1L,0L)))
)
summary_tbl <- imap_dfr(scen, ~ fit_healed(.x) %>% mutate(scenario = .y)) %>% select(scenario, everything())
write_csv(summary_tbl, here("output","tables","sensitivity_analysis_summary.csv"))

# ---- (c) day-28 Fisher exact per species (scrape vs airbrush; healed yes vs not) ----
d28 <- raw %>% filter(day==28, !is.na(healed)) %>%
  mutate(healed_yes = healed=="yes")
fisher_tbl <- d28 %>% group_split(species) %>% map_dfr(function(g){
  tab <- table(factor(g$treatment, levels=c("airbrush","dremel")),
               factor(g$healed_yes, levels=c(FALSE,TRUE)))
  ft <- fisher.test(tab)
  tibble(species = as.character(g$species[1]),
         airbrush_yes = tab["airbrush","TRUE"], airbrush_n = sum(tab["airbrush",]),
         dremel_yes = tab["dremel","TRUE"],     dremel_n = sum(tab["dremel",]),
         OR = unname(ft$estimate), p_value = ft$p.value)
})
write_csv(fisher_tbl, here("output","tables","endpoint_fisher_by_species_day28.csv"))

cat("\n=== Master sensitivity summary (healed treatment effect) ===\n")
summary_tbl %>% mutate(across(where(is.numeric), ~signif(.x,3))) %>% print(n=Inf, width=Inf)
cat("\n=== Day-28 Fisher exact by species (model-free corroboration) ===\n")
print(fisher_tbl, n=Inf)
cat("\nWrote: sensitivity_analysis_summary.csv, endpoint_fisher_by_species_day28.csv\n")
