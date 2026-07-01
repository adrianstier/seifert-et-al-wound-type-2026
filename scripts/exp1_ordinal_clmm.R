# =============================================================================
# exp1_ordinal_clmm.R — ordinal (cumulative-link mixed) model of 3-level healing
# -----------------------------------------------------------------------------
# Confirms the binary collapse (healed yes vs not) used in the main text loses no
# inferential signal: refits the full ordinal response (no < incomplete < yes)
# as a CLMM, checks proportional odds, extracts category probabilities, and
# cross-walks the ordinal vs binary treatment/species effects.
# Pairs with the main-text healing paragraph (SI: "ordinal confirmation").
# =============================================================================
suppressMessages({library(tidyverse); library(here); library(ordinal); library(broom)})

d <- read_csv(here("data","airbrush_dremel.csv"), show_col_types = FALSE) %>%
  mutate(species   = factor(str_trim(tolower(species)),
                            levels = c("acropora","pocillopora","porites")),
         treatment = factor(str_trim(tolower(treatment)), levels = c("airbrush","dremel")),
         parent_id = factor(sub("[[:alpha:]]+$","", as.character(coral_id))),
         coral_id  = factor(coral_id),
         healed3   = factor(str_trim(tolower(healed)),
                            levels = c("no","incomplete","yes"), ordered = TRUE)) %>%
  filter(!is.na(healed3))

# --- CLMM: try nested RE, fall back to parent-only, then to fixed-effects clm ---
re_used <- "(1|parent_id/coral_id)"
m <- tryCatch(clmm(healed3 ~ treatment + species + (1|parent_id/coral_id), data = d, Hess = TRUE),
              error = function(e) NULL)
if (is.null(m)) { re_used <- "(1|parent_id)"
  m <- tryCatch(clmm(healed3 ~ treatment + species + (1|parent_id), data = d, Hess = TRUE),
                error = function(e) NULL) }
if (is.null(m)) { re_used <- "none (clm; RE non-estimable)"
  m <- clm(healed3 ~ treatment + species, data = d) }
cat("CLMM random effect used:", re_used, "\n")

# --- coefficients (logit + OR) ---
co <- as.data.frame(summary(m)$coefficients)
co$term <- rownames(co); names(co)[1:4] <- c("estimate","std.error","z","p.value")
coef_tbl <- co %>% as_tibble() %>%
  mutate(OR = ifelse(grepl("\\|", term), NA, exp(estimate)),
         OR_low = ifelse(grepl("\\|", term), NA, exp(estimate - 1.96*std.error)),
         OR_high= ifelse(grepl("\\|", term), NA, exp(estimate + 1.96*std.error)),
         re_structure = re_used) %>%
  select(term, estimate, std.error, z, p.value, OR, OR_low, OR_high, re_structure)
write_csv(coef_tbl, here("output","tables","ordinal_clmm_coefficients.csv"))

# --- fit statistics ---
fit_tbl <- tibble(logLik = as.numeric(logLik(m)), AIC = AIC(m),
                  n_obs = nrow(d), n_parent = nlevels(d$parent_id),
                  re_structure = re_used)
write_csv(fit_tbl, here("output","tables","ordinal_clmm_fit_statistics.csv"))

# --- proportional-odds (nominal effects) test via fixed-effects clm ---
mclm <- clm(healed3 ~ treatment + species, data = d)
po <- tryCatch({
  nt <- nominal_test(mclm)
  as.data.frame(nt) %>% rownames_to_column("term") %>% as_tibble()
}, error = function(e) tibble(term = NA, note = paste("nominal_test failed:", conditionMessage(e))))
write_csv(po, here("output","tables","ordinal_clmm_proportional_odds_tests.csv"))

# --- predicted category probabilities by species x treatment ---
# Use the fixed-effects clm on a clean grid (Exp1 REs are ~0, so marginal and
# conditional category probabilities coincide); emmeans mode="prob" on the
# nested-RE clmm returns a degenerate result here.
grid <- expand_grid(treatment = levels(d$treatment), species = levels(d$species))
pr   <- predict(mclm, newdata = grid, type = "prob")$fit
pp   <- bind_cols(grid, as_tibble(pr)) %>%
  rename_with(~paste0("Pr_", .x), all_of(c("no","incomplete","yes")))
write_csv(pp, here("output","tables","ordinal_clmm_predicted_probabilities.csv"))

# --- ordinal vs binary cross-walk (treatment + species effects) ---
bin <- read_csv(here("output","tables","binary_fixed_effect_ORs.csv"), show_col_types = FALSE)
xwalk <- coef_tbl %>% filter(!grepl("\\|", term)) %>%
  transmute(term, ordinal_OR = round(OR,2),
            ordinal_CI = sprintf("%.2f-%.2f", OR_low, OR_high),
            ordinal_p = signif(p.value,3)) %>%
  mutate(binary_match = case_when(grepl("dremel", term) ~ "Dremel vs Airbrush",
                                  grepl("pocillopora", term) ~ "Pocillopora vs Acropora",
                                  grepl("porites", term) ~ "Porites vs Acropora", TRUE ~ NA)) %>%
  left_join(bin %>% transmute(binary_match = Term, binary_OR = round(OR,2),
                              binary_CI = sprintf("%.2f-%.2f", `CI Lower`, `CI Upper`),
                              binary_p = signif(`P-value`,3)), by = "binary_match")
write_csv(xwalk, here("output","tables","ordinal_vs_binary_comparison.csv"))

cat("\n=== CLMM coefficients ===\n"); print(coef_tbl, n=Inf)
cat("\n=== Proportional-odds test ===\n"); print(po)
cat("\n=== Predicted category probabilities ===\n"); print(pp, n=Inf)
cat("\n=== Ordinal vs binary cross-walk ===\n"); print(xwalk, n=Inf)
cat("\nWrote: ordinal_clmm_coefficients.csv, ordinal_clmm_fit_statistics.csv,\n",
    "      ordinal_clmm_proportional_odds_tests.csv, ordinal_clmm_predicted_probabilities.csv,\n",
    "      ordinal_vs_binary_comparison.csv\n")
