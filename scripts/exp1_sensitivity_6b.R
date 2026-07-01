# =============================================================================
# exp1_sensitivity_6b.R — leave-one-coral-out sensitivity for coral 6b
# -----------------------------------------------------------------------------
# Coral 6b (Pocillopora/scrape) is the ONLY Pocillopora-scrape fragment with any
# healed=yes, and the source of the 2026-06-15 re-score that resolved complete
# separation. This script quantifies how much the Exp 1 healing conclusions
# depend on it: refit the healed + regeneration models on (a) full data and
# (b) data with coral 6b dropped, and report treatment OR (lme4 + Firth),
# treatment/species LRTs, and the Pocillopora x scrape healed=yes count.
# Mirrors the data prep in scripts/airbrush_dremel_10_15_2025.R.
# =============================================================================
suppressMessages({library(tidyverse); library(here); library(lme4); library(brglm2)})

df0 <- read_csv(here("data","airbrush_dremel.csv"), show_col_types = FALSE) %>%
  mutate(species   = str_trim(tolower(species)),
         treatment = factor(str_trim(tolower(treatment)), levels = c("airbrush","dremel")),
         species   = factor(species, levels = c("acropora","pocillopora","porites")),
         parent_id = factor(sub("[[:alpha:]]+$","", as.character(coral_id))),
         coral_id  = factor(coral_id),
         healed01  = if_else(healed == "yes", 1L, 0L),
         regen01   = if_else(regenerated == "yes", 1L, 0L))

fit_set <- function(d, outcome) {
  d2 <- d %>% filter(!is.na(.data[[outcome]]))
  f_add <- as.formula(paste0(outcome, " ~ treatment + species + (1|parent_id/coral_id)"))
  f_no_sp <- as.formula(paste0(outcome, " ~ treatment + (1|parent_id/coral_id)"))
  f_no_tr <- as.formula(paste0(outcome, " ~ species + (1|parent_id/coral_id)"))
  ctrl <- glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
  m   <- glmer(f_add,  data=d2, family=binomial, control=ctrl)
  mns <- glmer(f_no_sp,data=d2, family=binomial, control=ctrl)
  mnt <- glmer(f_no_tr,data=d2, family=binomial, control=ctrl)
  lme4_or <- exp(fixef(m)[["treatmentdremel"]])
  sp_lrt  <- anova(mns, m)                       # adding species
  tr_lrt  <- anova(mnt, m)                        # adding treatment
  # Firth (penalized fixed-effects logistic; handles separation)
  fm <- glm(as.formula(paste0(outcome," ~ treatment + species")),
            data=d2, family=binomial, method="brglmFit", type="AS_mean")
  firth_or <- exp(coef(fm)[["treatmentdremel"]])
  tibble(outcome=outcome,
         lme4_treatment_OR = lme4_or,
         firth_treatment_OR = firth_or,
         treatment_LRT_chi2 = tr_lrt$Chisq[2], treatment_LRT_p = tr_lrt$`Pr(>Chisq)`[2],
         species_LRT_chi2   = sp_lrt$Chisq[2],  species_LRT_p   = sp_lrt$`Pr(>Chisq)`[2],
         poc_scrape_yes = sum(d2$species=="pocillopora" & d2$treatment=="dremel" & d2[[outcome]]==1, na.rm=TRUE),
         n = nrow(d2))
}

scenarios <- list(full = df0, drop_6b = df0 %>% filter(coral_id != "6b"))
res <- imap_dfr(scenarios, function(d, nm)
  bind_rows(fit_set(d,"healed01"), fit_set(d,"regen01")) %>% mutate(scenario = nm)) %>%
  select(scenario, outcome, everything())

cat("\n================ 6b leave-one-out sensitivity ================\n")
res %>% mutate(across(where(is.numeric), ~round(.x,3))) %>% print(n=Inf, width=Inf)
write_csv(res, here("output","tables","sensitivity_6b_exclusion.csv"))
cat("\nSaved: output/tables/sensitivity_6b_exclusion.csv\n")
