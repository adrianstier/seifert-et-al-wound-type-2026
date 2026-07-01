# =============================================================================
# Wound Healing Analysis: Publication-Quality Script (Figures -> output/figures)
# Date: 2025-10-27 (refactor); 2026-04-07 (regenerated column + day 3 update);
#       2026-04-21 (independent re-score of photos); 2026-06-15 (an independent re-score for NA-backtrack)
# =============================================================================
# Data version: Wound_Type_exp1_2022 Google Sheet (2026-06-15),
#   ingested 2026-06-18. Backtracked missing observations using a t-1/t+1 neighbour-cell rule
#   (regeneration absorbing); NA now healed 3 / debris 4 /
#   regenerated 0. Two coral-6b re-scores (D23 healed incomplete->yes; D28 healed
#   NA->yes) place 2 healed=yes in Pocillopora/dremel, so the prior COMPLETE
#   SEPARATION is RESOLVED: the lme4 healed GLMM is now estimable (treatment
#   OR ~26); Firth is retained as a sensitivity refit (OR ~17.3), not the
#   forced headline. Regeneration still uses Firth (airbrush cells ~all-zero).
# 2026-04-07 changes: (1) `regenerated` column added between day and healed;
#   (2) first post-wounding timepoint corrected from day 4 to day 3.

# ---- 1. Load Libraries & Setup ----
suppressPackageStartupMessages({
  library(tidyverse)      # dplyr, tidyr, ggplot2, purrr, readr
  library(here)           # here()
  library(janitor)        # clean_names()
  library(lme4)           # glmer()
  library(broom)          # tidy()
  library(broom.mixed)    # tidy() for mixed models
  library(ggpubr)         # theme_pubr()
  library(gt)             # gt() and gtsave()
  library(emmeans)        # emmeans/contrast
  library(scales)         # percent_format()
  library(tidyselect)     # all_of()
  library(ordinal)        # clmm()
})

# Output dirs
dir.create(here("output","figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output","tables"),  recursive = TRUE, showWarnings = FALSE)

# Aesthetic constants
palette_heal <- c("no"="#999999","incomplete"="#F0E442","yes"="#009E73")
theme_fig <- theme_classic(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 9, color = "grey30"),
        plot.margin = margin(10, 10, 10, 10, "mm"))

# ---- 2. Read & Clean Data & QC ----
data_path <- here("data","airbrush_dremel.csv")
if (!file.exists(data_path)) stop("Missing file: ", data_path)

df_raw <- read_csv(data_path, show_col_types = FALSE) %>%
  clean_names()

# Optional metadata join for parent_id if available (safe to skip if file absent)
meta_path <- here("data","airbrush_dremel_metadata.csv")
if (file.exists(meta_path)) {
  meta <- readr::read_csv(meta_path, show_col_types = FALSE) %>% clean_names()
  if (all(c("coral_id","parent_id") %in% names(meta))) {
    df_raw <- df_raw %>% left_join(meta %>% select(coral_id, parent_id), by = "coral_id")
  }
}
# If parent_id missing, derive a fallback (not used in modeling if identical to coral_id)
if (!("parent_id" %in% names(df_raw))) {
  df_raw <- df_raw %>% mutate(parent_id = gsub("[[:alpha:]]+$", "", coral_id))
}

# Normalize optional columns if present
if ("use_image" %in% names(df_raw)) {
  df_raw <- df_raw %>%
    mutate(use_image = if_else(str_detect(tolower(use_image), "yes"), "yes", "no"))
}

# Normalize regenerated (yes/no, with NA allowed)
if ("regenerated" %in% names(df_raw)) {
  df_raw <- df_raw %>%
    mutate(regenerated = if_else(is.na(regenerated), NA_character_,
                                 str_trim(tolower(regenerated))))
}

# Basic structure expectations
required_cols <- c("species","treatment","coral_id","day","healed")
missing_cols <- setdiff(required_cols, names(df_raw))
if (length(missing_cols)) stop("Missing required columns: ", paste(missing_cols, collapse=", "))

# Trim whitespace and normalize case for key text cols
df_raw <- df_raw %>%
  mutate(
    across(c(species, treatment, coral_id, healed,
             intersect(c("debris","use_image","regenerated"), names(df_raw))),
           ~ if (is.character(.x)) str_trim(tolower(.x)) else .x)
  )

# Coerce day to integer safely; warn if non-integer values present
if (any(!is.na(df_raw$day) & df_raw$day != round(df_raw$day))) {
  warning("Some 'day' values are not integers; coercing with round().")
}
df_raw <- df_raw %>% mutate(day = as.integer(round(day)))

# QC: allowed value checks (warn only)
allowed_species   <- c("acropora","pocillopora","porites")
allowed_treatment <- c("airbrush","dremel")
allowed_healed    <- c("no","incomplete","yes")
allowed_yn        <- c("yes","no")

warn_oov <- function(vec, allowed, field) {
  bad <- sort(unique(setdiff(na.omit(vec), allowed)))
  if (length(bad)) warning("Out-of-vocabulary values in '", field, "': ",
                           paste(bad, collapse=", "),
                           " | Allowed: ", paste(allowed, collapse=", "))
}
warn_oov(df_raw$species,   allowed_species,   "species")
warn_oov(df_raw$treatment, allowed_treatment, "treatment")
warn_oov(df_raw$healed,    allowed_healed,    "healed")
if ("debris"      %in% names(df_raw)) warn_oov(df_raw$debris,      allowed_yn, "debris")
if ("use_image"   %in% names(df_raw)) warn_oov(df_raw$use_image,   allowed_yn, "use_image")
if ("regenerated" %in% names(df_raw)) warn_oov(df_raw$regenerated, allowed_yn, "regenerated")

# QC: duplicates per coral_id x day x treatment?
dup_idx <- df_raw %>%
  count(species, treatment, coral_id, day, name = "n") %>%
  filter(n > 1)
if (nrow(dup_idx)) {
  warning("Found duplicate rows for species:treatment:coral_id:day; keeping first instance.")
  df_raw <- df_raw %>% distinct(species, treatment, coral_id, day, .keep_all = TRUE)
}

# ---- 2b. Pairing and repeated-measures structure ----
# coral_id uniquely identifies each coral (sibling) and serves as the random-effect grouping variable.
# A separate pair_id is not used because coral_id and pair_id are identical in this dataset.
df_raw <- df_raw %>% mutate(coral_id = factor(coral_id))

# Tidy analysis df with factors
df <- df_raw %>%
  transmute(
    day         = as.integer(day),
    treatment   = factor(treatment, levels = allowed_treatment),
    species     = factor(species,   levels = allowed_species),
    coral_id    = factor(coral_id),
    healed      = factor(healed,    levels = allowed_healed),
    debris      = if ("debris"      %in% names(df_raw)) factor(debris,      levels = allowed_yn) else factor(NA),
    use_image   = if ("use_image"   %in% names(df_raw)) factor(use_image,   levels = allowed_yn) else factor(NA),
    regenerated = if ("regenerated" %in% names(df_raw)) factor(regenerated, levels = allowed_yn) else factor(NA)
  )

# QC: missingness summary
missing_summary <- df %>%
  summarise(across(everything(), ~sum(is.na(.))), .groups = "drop") %>%
  pivot_longer(everything(), names_to = "column", values_to = "n_missing") %>%
  arrange(desc(n_missing))
print(missing_summary)
readr::write_csv(missing_summary, here("output","tables","qc_missingness_summary.csv"))

# Replicate summary (distinct coral IDs per species × treatment)
rep_summary <- df %>%
  group_by(species, treatment) %>%
  summarise(n_replicates = n_distinct(coral_id), .groups = "drop")
rep_tbl <- rep_summary %>%
  gt() %>%
  cols_label(species = "Species", treatment = "Wound Type", n_replicates = "No. of Coral Replicates") %>%
  tab_header(title = "Replicates by Species & Wound Type") %>%
  fmt_number(columns = "n_replicates", decimals = 0) %>%
  tab_source_note(md("Distinct coral IDs per species and treatment"))
print(rep_tbl)
readr::write_csv(rep_summary, here("output","tables","qc_replicate_summary.csv"))

# ---- 3. Binary outcome prep & descriptives ----
df <- df %>%
  mutate(
    healed_bin = factor(if_else(healed == "yes", "yes", "no"), levels = c("no","yes")),
    healed01   = if_else(healed == "yes", 1L, 0L)
  )



# ---- Derive pairing IDs from coral_id (e.g., "11a" -> parent_id "11") ----
df <- df %>%
  mutate(
    coral_id  = factor(coral_id),                         # keep the original child ID (11a, 11b, ...)
    parent_id = factor(sub("[[:alpha:]]+$", "", as.character(coral_id))),  # strip trailing letters
    pair_ab   = factor(sub("^[0-9]+", "", as.character(coral_id)))         # the suffix: a/b/etc.
  )

# Quick sanity checks (optional but helpful)
# 1) Parent has multiple children?
pairing_check <- df %>%
  distinct(parent_id, coral_id) %>%
  count(parent_id, name = "n_children") %>%
  arrange(parent_id)
readr::write_csv(pairing_check, here("output","tables","qc_parent_children_counts.csv"))

# 2) Any coral_id without a numeric parent extracted?
no_parent <- df %>%
  filter(is.na(parent_id) | parent_id == "")
if (nrow(no_parent) > 0) {
  warning("Some coral_id values did not yield a numeric parent_id. See 'no_parent_preview.csv'.")
  readr::write_csv(no_parent %>% head(20), here("output","tables","no_parent_preview.csv"))
}



df_counts <- df %>%
  group_by(day, treatment, species) %>%
  summarise(prop_healed = mean(healed_bin == "yes", na.rm = TRUE), .groups = "drop")
readr::write_csv(df_counts, here("output","tables","binary_proportion_by_day_species_treatment.csv"))


# ---- 4. Binary GLMMs (additive & interaction; REs for pairing + repeats) ----
# Notes:
# - (1 | parent_id)  : accounts for parent-colony pairing/blocking
# - (1 | coral_id)   : accounts for repeated measures on each coral
# - No fixed 'day' term; time is absorbed via the coral_id random intercept
# - Keep OLRE option separately if you diagnose overdispersion

mod_add <- glmer(
  healed01 ~ treatment + species +
    (1 | parent_id/coral_id),
  data = df, family = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
  na.action = na.exclude
)

mod_int <- glmer(
  healed01 ~ treatment * species +
    (1 | parent_id/coral_id),
  data = df, family = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
  na.action = na.exclude
)

# LRT interaction vs additive
lrt_sp <- anova(mod_add, mod_int, test = "Chisq")
tbl_lrt_bin <- tibble(
  Term         = "Treatment × Species Interaction",
  `Chi-square` = lrt_sp$Chisq[2],
  DF           = lrt_sp$Df[2],
  `P-value`    = lrt_sp$`Pr(>Chisq)`[2]
) %>%
  gt() %>%
  fmt_number(columns = all_of(c("Chi-square","P-value")), decimals = 3) %>%
  tab_header(title = "LRT: Interaction vs Additive GLMM (Binary)") %>%
  tab_source_note(md("Compare additive model against treatment×species"))
print(tbl_lrt_bin)
gtsave(tbl_lrt_bin, here("output","tables","binary_lrt_interaction_vs_additive.html"))
readr::write_csv(
  tibble(Term = "Treatment × Species Interaction",
         Chi_square = lrt_sp$Chisq[2],
         DF = lrt_sp$Df[2],
         P_value = lrt_sp$`Pr(>Chisq)`[2]),
  here("output","tables","binary_lrt_interaction_vs_additive.csv")
)

# Single-term deletion on additive model (main effects)
d1_tbl <- as_tibble(drop1(mod_add, test = "Chisq"), rownames = "Term")
drop_df <- d1_tbl %>%
  filter(Term %in% c("treatment","species")) %>%
  transmute(
    Term,
    DF           = npar,
    `Chi-square` = LRT,
    `P-value`    = `Pr(Chi)`
  )
tbl_drop <- drop_df %>%
  gt() %>%
  fmt_number(columns = all_of(c("Chi-square","P-value")), decimals = 3) %>%
  tab_header(title = "Main-Effect LRTs (Additive GLMM, Binary)") %>%
  tab_source_note(md("Single‐term deletion tests for treatment and species"))
print(tbl_drop)
gtsave(tbl_drop, here("output","tables","binary_glmm_drop1_additive.html"))
readr::write_csv(drop_df, here("output","tables","binary_glmm_drop1_additive.csv"))

# Pairwise species contrasts (binary model)
emm_spp <- emmeans(mod_add, ~ species)
spp_pairs <- contrast(emm_spp, method = "tukey", type = "response")
ci_spp  <- confint(spp_pairs)
pv_df   <- as.data.frame(summary(spp_pairs)) %>% dplyr::select(contrast, p.value)
spp_tbl <- as.data.frame(ci_spp) %>%
  left_join(pv_df, by = "contrast") %>%
  transmute(
    Comparison       = contrast,
    OR               = odds.ratio,
    `% Δ Odds`       = (odds.ratio - 1) * 100,
    `CI Lower (95%)` = asymp.LCL,
    `CI Upper (95%)` = asymp.UCL,
    `P-value`        = p.value
  )
tbl_pairs <- spp_tbl %>%
  gt() %>%
  fmt_number(columns = all_of(c("OR", "% Δ Odds", "CI Lower (95%)", "CI Upper (95%)")), decimals = 2) %>%
  fmt_number(columns = all_of("P-value"), decimals = 3, use_seps = FALSE) %>%
  tab_header(title = "Pairwise Species Contrasts (Tukey-adjusted, Binary)") %>%
  tab_source_note(md("Odds ratios (OR) with Wald 95% CIs; % Δ Odds = (OR–1)×100"))
print(tbl_pairs)
readr::write_csv(spp_tbl, here("output","tables","binary_emmeans_species_pairs.csv"))

# Fixed-effect ORs (binary model)
fixed_df <- tidy(mod_add, effects = "fixed", conf.int = TRUE, conf.method = "Wald") %>%
  filter(term %in% c("treatmentdremel","speciespocillopora","speciesporites")) %>%
  transmute(
    Term = case_when(
      term == "treatmentdremel"    ~ "Dremel vs Airbrush",
      term == "speciespocillopora" ~ "Pocillopora vs Acropora",
      term == "speciesporites"     ~ "Porites vs Acropora"
    ),
    OR         = exp(estimate),
    `CI Lower` = exp(conf.low),
    `CI Upper` = exp(conf.high),
    `P-value`  = p.value
  )
tbl_fixed <- fixed_df %>%
  gt() %>%
  fmt_scientific(columns = all_of(c("OR", "CI Lower", "CI Upper")), decimals = 2) %>%
  fmt_number(columns = all_of("P-value"), decimals = 3, use_seps = FALSE) %>%
  tab_header(title = "Fixed‐Effect Odds Ratios (Additive GLMM, Binary)") %>%
  tab_source_note(md("Wald 95% CIs; random intercept for coral_id"))
print(tbl_fixed)
readr::write_csv(fixed_df, here("output","tables","binary_fixed_effect_ORs.csv"))

# Binary trajectories + population-level predictions
df_plot <- df %>% mutate(healed01 = if_else(healed_bin == "yes", 1L, 0L))
pred_df <- expand_grid(
  treatment = levels(df$treatment),
  species   = levels(df$species),
  day       = sort(unique(df$day))
) %>% mutate(
  pred = predict(mod_add, newdata = ., type = "response", re.form = NA)
)
readr::write_csv(pred_df, here("output","tables","binary_predicted_probabilities_by_day_species_treatment.csv"))


# ---- 5. Ordinal outcome (no / incomplete / yes) ----
df3 <- df %>%
  mutate(
    healed3 = factor(as.character(healed),
                     levels = c("no","incomplete","yes"),
                     ordered = TRUE)
  ) %>% filter(!is.na(healed3))


# ---- 5a. Debris proportion & composition overlay ----

# Proportion of colonies with debris by day × species × treatment
obs_debris <- df3 %>%
  filter(!is.na(debris)) %>%
  group_by(day, species, treatment) %>%
  summarise(prop_debris = mean(debris == "yes", na.rm = TRUE), .groups = "drop") %>%
  arrange(day, species, treatment)

# ---- 5b. Regeneration descriptive summary ----
# New column from the 2026-04-07 data refresh: polyps in center of wound indicate
# complete structural regeneration (distinct from pigment-based `healed`).
# Export a simple proportion table; full modelling is future work.
if ("regenerated" %in% names(df3) && any(!is.na(df3$regenerated))) {
  obs_regen <- df3 %>%
    filter(!is.na(regenerated)) %>%
    group_by(day, species, treatment) %>%
    summarise(
      n             = dplyr::n(),
      n_regenerated = sum(regenerated == "yes"),
      prop_regen    = mean(regenerated == "yes"),
      .groups = "drop"
    ) %>%
    arrange(day, species, treatment)
  readr::write_csv(
    obs_regen,
    here("output","tables","regeneration_proportion_by_day_species_treatment.csv")
  )

  # Cross-tab of regenerated vs healed (overall)
  regen_x_healed <- df3 %>%
    filter(!is.na(regenerated)) %>%
    count(regenerated, healed, name = "n") %>%
    tidyr::pivot_wider(names_from = healed, values_from = n, values_fill = 0L)
  readr::write_csv(
    regen_x_healed,
    here("output","tables","regeneration_by_healed_crosstab.csv")
  )

  # Marginal proportions by treatment
  regen_by_trt <- df3 %>%
    filter(!is.na(regenerated)) %>%
    group_by(treatment) %>%
    summarise(
      n              = dplyr::n(),
      n_regenerated  = sum(regenerated == "yes"),
      prop_regen     = mean(regenerated == "yes"),
      .groups = "drop"
    )
  readr::write_csv(
    regen_by_trt,
    here("output","tables","regeneration_marginal_by_treatment.csv")
  )
  message("Exported regeneration descriptive summaries (3 CSVs in output/tables/).")
}

# Composition plots
df3_props <- df3 %>%
  count(day, species, treatment, healed3, name = "n") %>%
  tidyr::complete(
    day, species, treatment,
    healed3 = factor(c("no","incomplete","yes"),
                     levels = c("no","incomplete","yes"),
                     ordered = TRUE),
    fill = list(n = 0)
  ) %>%
  group_by(day, species, treatment) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  distinct(day, species, treatment, healed3, .keep_all = TRUE)


days_fac <- factor(sort(unique(df3_props$day)))  # ensures the same level order as the bars

# Relabel healing states for publication
df3_props <- df3_props %>%
  mutate(healed3 = factor(healed3,
                          levels = c("no", "incomplete", "yes"),
                          labels = c("Not healed", "Incomplete", "Fully healed"),
                          ordered = TRUE))

palette_heal_pub <- c("Not healed" = "#D55E00", "Incomplete" = "#E69F00", "Fully healed" = "#009E73")

# Add a dummy aesthetic for the debris line legend entry
obs_debris <- obs_debris %>%
  mutate(legend_label = "Algae prevalence")

p_stack <- ggplot(df3_props, aes(x = factor(day), y = prop, fill = healed3)) +
  geom_col(position = "fill", width = 0.9, color = "white", linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_hline(yintercept = 1, linewidth = 0.3, color = "grey50") +
  geom_line(data = obs_debris,
            aes(x = as.numeric(factor(day, levels = levels(days_fac))),
                y = prop_debris, group = 1, linetype = legend_label),
            inherit.aes = FALSE, color = "black", linewidth = 0.8) +
  geom_point(data = obs_debris,
             aes(x = as.numeric(factor(day, levels = levels(days_fac))), y = prop_debris),
             inherit.aes = FALSE, shape = 21, size = 2.0, fill = "white", color = "black") +
  facet_grid(species ~ treatment,
             labeller = labeller(
               species = c(acropora = "Acropora", pocillopora = "Pocillopora", porites = "Porites"),
               treatment = c(airbrush = "Airbrush", dremel = "Scrape")
             )) +
  scale_fill_manual(values = palette_heal_pub, name = "Healing state") +
  scale_linetype_manual(values = "solid", name = NULL) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(-0.08, 1.08),
                     breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Day post-wounding", y = "Percent of fragments") +
  theme_fig +
  theme(strip.text.y = element_text(face = "bold.italic"),
        panel.spacing.y = unit(8, "mm"))

#make and write a figure in pdf and png

# (bare print() removed: it spilled a stray Rplots.pdf in non-interactive
#  Rscript; the figure is persisted by the ggsave() calls just below.)
fig_date <- format(Sys.Date(), "%Y_%m_%d")
ggsave(here("output","figures", paste0("ordinal_composition_with_debris_overlay_", fig_date, ".png")),
       p_stack, width = 170, height = 150, units = "mm", dpi = 300)
ggsave(here("output","figures", paste0("ordinal_composition_with_debris_overlay_", fig_date, ".pdf")),
       p_stack, width = 170, height = 150, units = "mm")


# =============================================================================
# 5b-2. Binary REGENERATION composition figure (revised Figure 2)
# -----------------------------------------------------------------------------
# Added 2026-05-15 per coauthor request: collapse the 3-level `healed`
# composition (no/incomplete/yes) to a BINARY `regenerated` (no/yes) outcome.
# Encoding: white fill + black outline = Not regenerated; green fill =
# Regenerated. Keeps the black algae-prevalence overlay (same `obs_debris`
# layer as p_stack). This block is ADDITIVE — it does not modify or replace
# the 3-level `p_stack` figure above.
# =============================================================================
df_regen_props <- df3 %>%
  filter(!is.na(regenerated)) %>%
  count(day, species, treatment, regenerated, name = "n") %>%
  tidyr::complete(
    day, species, treatment,
    regenerated = factor(c("no", "yes"), levels = c("no", "yes")),
    fill = list(n = 0)
  ) %>%
  group_by(day, species, treatment) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  distinct(day, species, treatment, regenerated, .keep_all = TRUE) %>%
  mutate(regen_lab = factor(regenerated,
                            levels = c("no", "yes"),
                            labels = c("Not regenerated", "Regenerated")))

palette_regen_pub <- c("Not regenerated" = "#FFFFFF", "Regenerated" = "#009E73")

p_regen_stack <- ggplot(df_regen_props, aes(x = factor(day), y = prop, fill = regen_lab)) +
  geom_col(position = "fill", width = 0.9, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_hline(yintercept = 1, linewidth = 0.3, color = "grey50") +
  geom_line(data = obs_debris,
            aes(x = as.numeric(factor(day, levels = levels(days_fac))),
                y = prop_debris, group = 1, linetype = legend_label),
            inherit.aes = FALSE, color = "black", linewidth = 0.8) +
  geom_point(data = obs_debris,
             aes(x = as.numeric(factor(day, levels = levels(days_fac))), y = prop_debris),
             inherit.aes = FALSE, shape = 21, size = 2.0, fill = "white", color = "black") +
  facet_grid(species ~ treatment,
             labeller = labeller(
               species = c(acropora = "Acropora", pocillopora = "Pocillopora", porites = "Porites"),
               treatment = c(airbrush = "Airbrush", dremel = "Scrape")
             )) +
  scale_fill_manual(values = palette_regen_pub, name = "Regeneration") +
  scale_linetype_manual(values = "solid", name = NULL) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(-0.08, 1.08),
                     breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Day post-wounding", y = "Percent of fragments") +
  theme_fig +
  theme(strip.text.y = element_text(face = "bold.italic"),
        panel.spacing.y = unit(8, "mm"))

# (bare print() removed — see note above; ggsave() below persists the figure)
ggsave(here("output","figures", paste0("regenerated_composition_with_debris_overlay_", fig_date, ".png")),
       p_regen_stack, width = 170, height = 150, units = "mm", dpi = 300)
ggsave(here("output","figures", paste0("regenerated_composition_with_debris_overlay_", fig_date, ".pdf")),
       p_regen_stack, width = 170, height = 150, units = "mm")


# =============================================================================
# 5b-2-alt. BINARY HEALED composition figure  (non-manuscript alternate)
# -----------------------------------------------------------------------------
# Added 2026-05-21 as a candidate Fig 2B variant (binary healed yes/not-yes +
# algae overlay, parallel to §5b-2's binary regenerated + algae). On the same
# day the swap was REVERTED: the canonical manuscript Fig 2B is the binary
# `regenerated` figure (§5b-2 above), because the manuscript text references
# Fig 2B in the context of regeneration. This block is kept as a sensitivity
# / alternate figure (do not delete -- it's still produced by `make all`).
#
# Parallel structure to §5b-2 above; same palette (white outlined / green
# #009E73), same 170 x 150 mm size, same facet grid (species x wound type).
# §5a (3-level ordinal healed) remains a separate diagnostic figure.
#
# Filename: healed_composition_with_debris_overlay_<DATE>.{pdf,png}.
# See output/FIGURE_CROSSWALK.md for the canonical Fig 2B target.
# =============================================================================
df_heal_props <- df3 %>%
  mutate(healed01 = factor(if_else(healed == "yes", "yes", "no",
                                   missing = NA_character_),
                           levels = c("no","yes"))) %>%
  filter(!is.na(healed01)) %>%
  count(day, species, treatment, healed01, name = "n") %>%
  tidyr::complete(
    day, species, treatment,
    healed01 = factor(c("no", "yes"), levels = c("no", "yes")),
    fill = list(n = 0)
  ) %>%
  group_by(day, species, treatment) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  distinct(day, species, treatment, healed01, .keep_all = TRUE) %>%
  mutate(heal_lab = factor(healed01,
                           levels = c("no","yes"),
                           labels = c("Not healed", "Healed")))

palette_heal_pub <- c("Not healed" = "#FFFFFF", "Healed" = "#009E73")

p_heal_stack <- ggplot(df_heal_props,
                       aes(x = factor(day), y = prop, fill = heal_lab)) +
  geom_col(position = "fill", width = 0.9, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_hline(yintercept = 1, linewidth = 0.3, color = "grey50") +
  geom_line(data = obs_debris,
            aes(x = as.numeric(factor(day, levels = levels(days_fac))),
                y = prop_debris, group = 1, linetype = legend_label),
            inherit.aes = FALSE, color = "black", linewidth = 0.8) +
  geom_point(data = obs_debris,
             aes(x = as.numeric(factor(day, levels = levels(days_fac))),
                 y = prop_debris),
             inherit.aes = FALSE, shape = 21, size = 2.0,
             fill = "white", color = "black") +
  facet_grid(species ~ treatment,
             labeller = labeller(
               species = c(acropora = "Acropora", pocillopora = "Pocillopora",
                           porites = "Porites"),
               treatment = c(airbrush = "Airbrush", dremel = "Scrape")
             )) +
  scale_fill_manual(values = palette_heal_pub, name = "Healing") +
  scale_linetype_manual(values = "solid", name = NULL) +
  scale_y_continuous(labels = scales::percent_format(),
                     limits = c(-0.08, 1.08),
                     breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Day post-wounding", y = "Percent of fragments") +
  theme_fig +
  theme(strip.text.y = element_text(face = "bold.italic"),
        panel.spacing.y = unit(8, "mm"))

ggsave(here("output","figures", paste0("healed_composition_with_debris_overlay_", fig_date, ".png")),
       p_heal_stack, width = 170, height = 150, units = "mm", dpi = 300)
ggsave(here("output","figures", paste0("healed_composition_with_debris_overlay_", fig_date, ".pdf")),
       p_heal_stack, width = 170, height = 150, units = "mm")


# =============================================================================
# 5b-3. Multi-outcome time-series panel (manuscript Figure 3D rebuild)
# -----------------------------------------------------------------------------
# Added 2026-05-15 per coauthor request: rebuild the Fig 3D-style multi-panel
# time series from the REPRODUCIBLE Exp 1 data in this repo. Definitions agreed
# with coauthors:
#   Healed      = % fragments with complete coenosarc coverage (healed == "yes")
#   Regenerated = % fragments with polyps in wound center      (regenerated == "yes")
#   Algal plug  = % fragments with algae / necrotic tissue     (debris == "yes")
# A "Yellow aggregation" panel was requested but that variable is NOT scored in
# the Exp 1 dataset and no source data exists in-repo -> OMITTED, not fabricated.
# Additive block; does not modify any figure above.
# =============================================================================
outcome_ts <- dplyr::bind_rows(
  df %>% filter(!is.na(healed)) %>%
    group_by(day, species, treatment) %>%
    summarise(n = dplyr::n(), k = sum(healed == "yes"), .groups = "drop") %>%
    mutate(outcome = "Healed (complete coenosarc)"),
  df %>% filter(!is.na(regenerated)) %>%
    group_by(day, species, treatment) %>%
    summarise(n = dplyr::n(), k = sum(regenerated == "yes"), .groups = "drop") %>%
    mutate(outcome = "Regenerated (polyps in center)"),
  df %>% filter(!is.na(debris)) %>%
    group_by(day, species, treatment) %>%
    summarise(n = dplyr::n(), k = sum(debris == "yes"), .groups = "drop") %>%
    mutate(outcome = "Algal plug")
) %>%
  mutate(prop = k / n,
         outcome = factor(outcome,
                          levels = c("Healed (complete coenosarc)",
                                     "Regenerated (polyps in center)",
                                     "Algal plug")))

readr::write_csv(outcome_ts,
                 here("output","tables", "exp1_outcome_timeseries_props.csv"))

p_outcomes <- ggplot(outcome_ts,
                     aes(x = day, y = prop, color = treatment, group = treatment)) +
  geom_line(linewidth = 0.9) +
  geom_point(aes(size = n), alpha = 0.85) +
  facet_grid(outcome ~ species,
             labeller = labeller(
               species = c(acropora = "Acropora", pocillopora = "Pocillopora",
                           porites = "Porites"),
               outcome = ggplot2::label_wrap_gen(18))) +
  scale_color_manual(values = c(airbrush = "#D55E00", dremel = "#009E73"),
                     labels = c(airbrush = "Airbrush", dremel = "Scrape"),
                     name = "Wound type") +
  scale_size_continuous(name = "n fragments", range = c(1.5, 3.5),
                        breaks = c(1, 3, 5)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(-0.02, 1.02), breaks = seq(0, 1, 0.5)) +
  scale_x_continuous(breaks = sort(unique(outcome_ts$day))) +
  labs(x = "Day post-wounding", y = "Percent of fragments") +
  theme_fig +
  theme(strip.text.x = element_text(face = "bold.italic"),
        strip.text.y = element_text(face = "bold"),
        panel.spacing = unit(6, "mm"),
        legend.position = "bottom")

# (bare print() removed — see note above; ggsave() below persists the figure)
ggsave(here("output","figures", paste0("exp1_outcome_timeseries_panels_", fig_date, ".png")),
       p_outcomes, width = 180, height = 175, units = "mm", dpi = 300)
ggsave(here("output","figures", paste0("exp1_outcome_timeseries_panels_", fig_date, ".pdf")),
       p_outcomes, width = 180, height = 175, units = "mm")


# =============================================================================
# 5c. Regenerated outcome GLMMs (parallel structure to Section 4)
# =============================================================================
# `regenerated` codes whether polyps in the wound center indicate complete
# structural regeneration (yes) vs not (no). It is a binary outcome added in the
# 2026-04-07 data refresh. Treatment skew is severe (airbrush 1/72 = 1.4%
# yes; dremel 38/71 = 53.5%), so expect convergence warnings (singular fit /
# inflated SE on the airbrush slope) — analogous to Section 4. Inferences from
# fixed-effect ORs should be interpreted alongside the descriptive proportions
# in `regeneration_marginal_by_treatment.csv`.

if ("regenerated" %in% names(df) && any(!is.na(df$regenerated))) {

  df_regen <- df %>%
    mutate(regen01 = if_else(regenerated == "yes", 1L,
                             if_else(regenerated == "no", 0L, NA_integer_))) %>%
    filter(!is.na(regen01))

  cat("\n--- Regenerated GLMM ---\n")
  cat("N observations (non-NA): ", nrow(df_regen), "\n", sep = "")
  cat("Yes / No counts: ", sum(df_regen$regen01 == 1L), " / ",
      sum(df_regen$regen01 == 0L), "\n", sep = "")

  mod_regen_add <- glmer(
    regen01 ~ treatment + species + (1 | parent_id/coral_id),
    data = df_regen, family = binomial(link = "logit"),
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
    na.action = na.exclude
  )

  mod_regen_int <- glmer(
    regen01 ~ treatment * species + (1 | parent_id/coral_id),
    data = df_regen, family = binomial(link = "logit"),
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
    na.action = na.exclude
  )

  # LRT interaction vs additive
  lrt_regen <- anova(mod_regen_add, mod_regen_int, test = "Chisq")
  readr::write_csv(
    tibble(Term       = "Treatment × Species Interaction",
           Chi_square = lrt_regen$Chisq[2],
           DF         = lrt_regen$Df[2],
           P_value    = lrt_regen$`Pr(>Chisq)`[2]),
    here("output","tables","regen_lrt_interaction_vs_additive.csv")
  )
  tbl_regen_lrt <- tibble(
    Term         = "Treatment × Species Interaction",
    `Chi-square` = lrt_regen$Chisq[2],
    DF           = lrt_regen$Df[2],
    `P-value`    = lrt_regen$`Pr(>Chisq)`[2]
  ) %>%
    gt() %>%
    fmt_number(columns = all_of(c("Chi-square","P-value")), decimals = 3) %>%
    tab_header(title = "LRT: Interaction vs Additive GLMM (Regenerated)")
  gtsave(tbl_regen_lrt, here("output","tables","regen_lrt_interaction_vs_additive.html"))

  # drop1 main effects on additive model
  d1_regen <- as_tibble(drop1(mod_regen_add, test = "Chisq"), rownames = "Term") %>%
    filter(Term %in% c("treatment","species")) %>%
    transmute(Term,
              DF           = npar,
              `Chi-square` = LRT,
              `P-value`    = `Pr(Chi)`)
  readr::write_csv(d1_regen, here("output","tables","regen_glmm_drop1_additive.csv"))
  tbl_regen_drop <- d1_regen %>%
    gt() %>%
    fmt_number(columns = all_of(c("Chi-square","P-value")), decimals = 3) %>%
    tab_header(title = "Main-Effect LRTs (Additive GLMM, Regenerated)")
  gtsave(tbl_regen_drop, here("output","tables","regen_glmm_drop1_additive.html"))

  # Fixed-effect ORs (additive)
  fixed_regen <- broom.mixed::tidy(mod_regen_add, effects = "fixed",
                                   conf.int = TRUE, conf.method = "Wald") %>%
    filter(term %in% c("treatmentdremel","speciespocillopora","speciesporites")) %>%
    transmute(
      Term = case_when(
        term == "treatmentdremel"    ~ "Dremel vs Airbrush",
        term == "speciespocillopora" ~ "Pocillopora vs Acropora",
        term == "speciesporites"     ~ "Porites vs Acropora"
      ),
      OR         = exp(estimate),
      `CI Lower` = exp(conf.low),
      `CI Upper` = exp(conf.high),
      `P-value`  = p.value
    )
  readr::write_csv(fixed_regen, here("output","tables","regen_fixed_effect_ORs.csv"))

  # emmeans: marginal probability of regeneration by species × treatment
  emm_regen_cells <- emmeans(mod_regen_add, ~ treatment * species, type = "response")
  cells_regen_df  <- as.data.frame(emm_regen_cells)
  readr::write_csv(cells_regen_df,
                   here("output","tables","regen_emmeans_cell_probabilities.csv"))

  # Treatment-within-species contrast (Dremel vs Airbrush) on probability scale
  emm_regen_t <- emmeans(mod_regen_add, ~ treatment | species, type = "response")
  t_in_s <- as.data.frame(contrast(emm_regen_t, "revpairwise"))
  readr::write_csv(t_in_s,
                   here("output","tables","regen_emmeans_treatment_within_species.csv"))

  # Species-within-treatment contrast (Tukey-adjusted)
  emm_regen_s <- emmeans(mod_regen_add, ~ species | treatment, type = "response")
  s_in_t <- as.data.frame(contrast(emm_regen_s, "tukey"))
  readr::write_csv(s_in_t,
                   here("output","tables","regen_emmeans_species_within_treatment.csv"))

  # Diagnostics + R²
  od_regen   <- sqrt(sum(residuals(mod_regen_add, type = "pearson")^2, na.rm = TRUE) /
                       df.residual(mod_regen_add))
  sing_regen <- lme4::isSingular(mod_regen_add, tol = 1e-5)
  r2_regen   <- suppressWarnings(MuMIn::r.squaredGLMM(mod_regen_add))
  diag_regen <- tibble(
    model        = "regen_add",
    n_obs        = nrow(df_regen),
    n_yes        = sum(df_regen$regen01 == 1L),
    overdisp     = od_regen,
    singular_fit = sing_regen,
    R2_marginal  = if (is.matrix(r2_regen)) r2_regen[1, 1] else NA_real_,
    R2_conditional = if (is.matrix(r2_regen)) r2_regen[1, 2] else NA_real_,
    AIC          = AIC(mod_regen_add),
    BIC          = BIC(mod_regen_add)
  )
  readr::write_csv(diag_regen, here("output","tables","regen_glmm_diagnostics.csv"))

  cat("Regenerated GLMM tables written (regen_*.csv in output/tables/).\n")
  cat("  Treatment × species LRT: chi2 =", round(lrt_regen$Chisq[2], 3),
      ", p =", signif(lrt_regen$`Pr(>Chisq)`[2], 3), "\n")
  cat("  Singular fit (additive):", sing_regen, "\n")

  # ---- 5d. Regeneration trajectory figure (new findings) ----
  # Observed proportion of fragments with polyp regeneration over time, faceted
  # by species and treatment. Companion to the existing healing composition figure.
  regen_traj <- df_regen %>%
    group_by(day, species, treatment) %>%
    summarise(
      n          = dplyr::n(),
      n_regen    = sum(regen01 == 1L),
      prop_regen = mean(regen01 == 1L),
      .groups    = "drop"
    )

  p_regen <- ggplot(regen_traj,
                    aes(x = day, y = prop_regen, color = treatment, group = treatment)) +
    geom_line(linewidth = 0.9) +
    geom_point(aes(size = n), alpha = 0.85) +
    facet_wrap(~ species, nrow = 1,
               labeller = labeller(species = c(acropora = "Acropora",
                                               pocillopora = "Pocillopora",
                                               porites = "Porites"))) +
    scale_color_manual(values = c(airbrush = "#D55E00", dremel = "#009E73"),
                       labels = c(airbrush = "Airbrush",
                                  dremel   = "Scrape (Dremel)"),
                       name   = "Wound type") +
    scale_size_continuous(name = "n fragments", range = c(1.5, 4),
                          breaks = c(1, 3, 5)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(-0.02, 1.02),
                       breaks = seq(0, 1, 0.25)) +
    scale_x_continuous(breaks = sort(unique(df_regen$day))) +
    labs(x = "Day post-wounding",
         y = "Fragments with polyp regeneration") +
    theme_fig +
    theme(strip.text = element_text(face = "bold.italic"))

  # (bare print() removed — see note above; ggsave() below persists the figure)
  ggsave(here("output","figures",
              paste0("regeneration_trajectory_", fig_date, ".png")),
         p_regen, width = 180, height = 90, units = "mm", dpi = 300)
  ggsave(here("output","figures",
              paste0("regeneration_trajectory_", fig_date, ".pdf")),
         p_regen, width = 180, height = 90, units = "mm")

} else {
  message("Skipping regenerated GLMM: column not present or all NA.")
}


# =============================================================================
# 5e. Firth-penalized sensitivity analysis (brglm2)
# =============================================================================
# Both binomial GLMMs above have singular fits (random-effect variances ≈ 0).
# That collapses the GLMM to a fixed-effects logistic regression in practice,
# but maximum-likelihood standard errors remain inflated under quasi-separation
# (especially the regenerated outcome: only 1/72 airbrush observations are yes).
#
# Firth's bias-corrected logistic regression (Firth 1993) gives finite, well-
# calibrated estimates and CIs even when the MLE is at the boundary. We use
# brglm2::brglmFit which implements the Firth penalty as a glm fitting method.
#
# These rows feed a "Sensitivity (Firth)" block in the master stats table so the
# manuscript can report both the GLMM ORs and the bias-corrected ORs side by side.

suppressPackageStartupMessages({
  library(brglm2)
})

firth_or_table <- function(mod, label) {
  co <- summary(mod)$coefficients
  rows <- intersect(c("treatmentdremel", "speciespocillopora", "speciesporites"),
                    rownames(co))
  # PROFILE-LIKELIHOOD CIs, not Wald exp(est +/- z*SE).
  # brglm2 fits supply a penalized-likelihood profiling method, so confint()
  # returns finite, separation-robust bounds. Estimates and p-values are
  # UNCHANGED (still from summary()); only the CI bounds differ from the old
  # Wald interval.
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

# ---- Healed: Firth-penalized fixed-effects logistic ----
df_healed_firth <- df %>% filter(!is.na(healed01))
mod_healed_firth <- glm(
  healed01 ~ treatment + species,
  data = df_healed_firth, family = binomial,
  method = "brglmFit"
)
firth_healed <- firth_or_table(mod_healed_firth, "healed")
readr::write_csv(firth_healed,
                 here("output","tables","firth_healed_fixed_effect_ORs.csv"))

cat("\n--- Firth (brglm2) — healed ---\n")
print(firth_healed)
cat("Convergence:", mod_healed_firth$converged, "\n")

# ---- Regenerated: Firth-penalized fixed-effects logistic ----
if (exists("df_regen") && nrow(df_regen) > 0) {
  mod_regen_firth <- glm(
    regen01 ~ treatment + species,
    data = df_regen, family = binomial,
    method = "brglmFit"
  )
  firth_regen <- firth_or_table(mod_regen_firth, "regenerated")
  readr::write_csv(firth_regen,
                   here("output","tables","firth_regen_fixed_effect_ORs.csv"))

  cat("\n--- Firth (brglm2) — regenerated ---\n")
  print(firth_regen)
  cat("Convergence:", mod_regen_firth$converged, "\n")
}



# =============================================================================
# 6. Paper-ready summaries: diagnostics, effects, contrasts, master tables
# =============================================================================

suppressPackageStartupMessages({
  library(glue)
  library(MuMIn)     # r.squaredGLMM
  library(DHARMa)    # quick GLMM residual checks (optional)
})

# ---- 6.1 Diagnostics (export) ----
# Overdispersion (Pearson-based)
overdisp_fun <- function(model) {
  rp <- residuals(model, type = "pearson")
  sqrt(sum(rp^2, na.rm = TRUE) / df.residual(model))
}
od_add  <- overdisp_fun(mod_add)
od_int  <- overdisp_fun(mod_int)
sing_add <- lme4::isSingular(mod_add, tol = 1e-5)
sing_int <- lme4::isSingular(mod_int, tol = 1e-5)
r2_add   <- suppressWarnings(MuMIn::r.squaredGLMM(mod_add))
r2_int   <- suppressWarnings(MuMIn::r.squaredGLMM(mod_int))

diag_tbl <- tibble(
  model     = c("mod_add", "mod_int"),
  overdisp  = c(od_add, od_int),
  is_singular = c(sing_add, sing_int),
  R2_marginal = c(r2_add[1, "R2m"], r2_int[1, "R2m"]),
  R2_conditional = c(r2_add[1, "R2c"], r2_int[1, "R2c"]),
  AIC       = c(AIC(mod_add), AIC(mod_int)),
  BIC       = c(BIC(mod_add), BIC(mod_int)),
  logLik    = c(as.numeric(logLik(mod_add)), as.numeric(logLik(mod_int)))
)
readr::write_csv(diag_tbl, here("output", "tables", "diagnostics_glmm_binary.csv"))

# DHARMa residual check for the interaction model. Robust by construction:
# simulate FIRST (no graphics device open yet); only open the PNG once we have
# residuals, and guarantee the device closes even if plot() errors. That is
# what prevents a stray repo-root Rplots.pdf (the old silent try() let the
# default device catch the plot when simulateResiduals failed). A failure is
# now logged LOUDLY instead of silently swallowed, so a missing diagnostic
# is auditable rather than indistinguishable from success.
dharma_png <- here("output", "figures", "dharma_residuals_mod_int.png")
dir.create(dirname(dharma_png), showWarnings = FALSE, recursive = TRUE)
dharma_ok <- tryCatch({
  set.seed(20260515)  # reproducible DHARMa simulation (matches Exp2 DHARMA_SEED)
  simres <- DHARMa::simulateResiduals(mod_int, n = 1000)
  png(dharma_png, width = 1800, height = 1200, res = 180)
  tryCatch(plot(simres), finally = dev.off())
  TRUE
}, error = function(e) {
  warning("DHARMa diagnostic for mod_int FAILED: ", conditionMessage(e),
          " -- ", basename(dharma_png), " was NOT written.", call. = FALSE)
  FALSE
})
if (isTRUE(dharma_ok)) {
  cat("  [OK] DHARMa residual plot:", basename(dharma_png), "\n")
} else {
  cat("  [x] DHARMa residual plot NOT produced (see warning above)\n")
}

# ---- 6.2 Fixed effects (log-odds & OR) with CIs/p-values ----
fixef_add <- broom.mixed::tidy(mod_add, effects = "fixed", conf.int = TRUE, conf.method = "Wald") %>%
  mutate(OR = exp(estimate),
         OR_low = exp(conf.low),
         OR_high = exp(conf.high)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high, OR, OR_low, OR_high)

readr::write_csv(fixef_add, here("output","tables","binary_glmm_additive_fixed_effects_logit_and_OR.csv"))

# ---- 6.3 Variance components + ICCs (latent-scale) ----
# For binomial logit, residual variance is pi^2/3 on the latent scale
latent_var <- (pi^2) / 3
vc <- as.data.frame(VarCorr(mod_add))
# Random terms present: parent_id and coral_id:parent_id
var_parent <- vc %>% filter(grp == "parent_id") %>% pull(vcov)
if (length(var_parent) == 0) var_parent <- 0
var_coral  <- vc %>% filter(grp == "parent_id:coral_id") %>% pull(vcov)
if (length(var_coral) == 0) var_coral <- 0

icc_tbl <- tibble(
  component = c("parent_id", "coral_id_within_parent", "residual_logit"),
  variance  = c(var_parent, var_coral, latent_var),
  ICC       = c(var_parent, var_coral, latent_var) /
    (var_parent + var_coral + latent_var)
)
readr::write_csv(icc_tbl, here("output","tables","binary_glmm_additive_icc_latent.csv"))

# ---- 6.4 Cell means (marginal probabilities) & key contrasts via emmeans ----
# Marginal Pr(healed) by species × treatment (response scale)
emm_cell <- emmeans::emmeans(mod_add, ~ species * treatment, type = "response")
emm_cell_df <- as.data.frame(emm_cell) %>%
  transmute(
    species, treatment,
    prob = prob, SE = SE, lower = asymp.LCL, upper = asymp.UCL
  )
readr::write_csv(emm_cell_df, here("output","tables","binary_emmeans_cell_probabilities.csv"))

# ---- 6.4a Treatment effect within each species (Dremel vs Airbrush) ----
# Goal: BOTH (i) probability difference = Dremel − Airbrush and (ii) odds ratio = Dremel/Airbrush

# 1) Probability difference on the response (probability) scale
# Force the EMMs onto the response scale first, THEN do the contrast.
emm_tr_link  <- emmeans::emmeans(mod_add, ~ treatment | species)           # link scale grid
emm_tr_resp  <- emmeans::regrid(emm_tr_link, transform = "response")       # now on probability scale
tp           <- emmeans::contrast(emm_tr_resp, method = "revpairwise")     # Dremel - Airbrush
tp_df        <- as.data.frame(summary(tp))                                  # has columns like estimate, SE, lower.CL, upper.CL, p.value
tp_ci        <- as.data.frame(confint(tp))                                  # ensure CIs exist regardless of summary() defaults

# Join summary + CI (robust to minor naming diffs)
treat_probdiff_df <- tp_df %>%
  dplyr::select(species, contrast,
                prob_diff = dplyr::any_of(c("estimate","emmean","response","prob","effect")),
                prob_SE   = dplyr::any_of(c("SE","SE.dif","SE.diff")),
                prob_p    = dplyr::any_of(c("p.value","p.value.mixed","p"))) %>%
  dplyr::left_join(
    tp_ci %>% dplyr::select(species, contrast,
                            prob_low  = dplyr::any_of(c("lower.CL","asymp.LCL","LCL")),
                            prob_high = dplyr::any_of(c("upper.CL","asymp.UCL","UCL"))),
    by = c("species","contrast")
  )

# 2) Odds ratio on response scale summary (pairwise on link, summarized as OR)
to            <- emmeans::contrast(emm_tr_link, method = "pairwise", reverse = TRUE)   # Dremel/Airbrush
to_df         <- as.data.frame(summary(to, type = "response"))                          # has odds ratio in 'ratio' or 'odds.ratio'
to_ci         <- as.data.frame(confint(to, type = "response"))

# Make robust to column-name differences (ratio vs odds.ratio; lower.CL vs asymp.LCL)
or_col <- intersect(c("ratio","odds.ratio"), names(to_df))[1]
lo_col <- intersect(c("lower.CL","asymp.LCL","LCL"), names(to_ci))[1]
hi_col <- intersect(c("upper.CL","asymp.UCL","UCL"), names(to_ci))[1]

treat_or_df <- to_df %>%
  dplyr::select(species, contrast,
                OR   = dplyr::all_of(or_col),
                OR_p = dplyr::any_of(c("p.value","p.value.mixed","p"))) %>%
  dplyr::left_join(
    to_ci %>% dplyr::select(species, contrast,
                            OR_low  = dplyr::all_of(lo_col),
                            OR_high = dplyr::all_of(hi_col)),
    by = c("species","contrast")
  )

# 3) Harmonize labels and join both measures
treat_probdiff_df$contrast <- gsub("-", "vs", treat_probdiff_df$contrast, fixed = TRUE)
treat_or_df$contrast       <- gsub("/", "vs", treat_or_df$contrast,       fixed = TRUE)

treat_contr_df <- dplyr::left_join(
  treat_probdiff_df, treat_or_df, by = c("species","contrast")
) %>%
  dplyr::select(
    species, contrast,
    prob_diff, prob_low, prob_high, prob_SE, prob_p,
    OR, OR_low, OR_high, OR_p
  )

# Helpful peek so we don't stay stuck
cat("\n6.4a preview — treatment within species:\n")
print(treat_contr_df)

readr::write_csv(treat_contr_df, here("output","tables","binary_emmeans_treatment_within_species.csv"))

# ---- 6.4.3 Species effect within each treatment (pairwise Tukey) ----
# Self-contained helpers (robust to emmeans version/name differences)
pick_col <- function(df, candidates, label) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) stop(paste0("No column found for ", label,
                                    ". Saw: ", paste(names(df), collapse=", ")))
  hit[1]
}
get_by_col_safe <- function(df, candidates = c("treatment","by1","by")) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) stop("No suitable 'by' column. Saw: ", paste(names(df), collapse=", "))
  hit[1]
}

# A) Probability differences on the response (probability) scale with Tukey adjustment
es_link <- emmeans::emmeans(mod_add, ~ species | treatment)          # link-scale grid
es_resp <- emmeans::regrid(es_link, transform = "response")          # regrid -> prob scale
sp_resp <- emmeans::contrast(es_resp, method = "tukey")              # pairwise species diffs

sp_df <- as.data.frame(summary(sp_resp))
sp_ci <- as.data.frame(confint(sp_resp))

by_col_s <- get_by_col_safe(sp_df, c("treatment","by1","by"))

sp_est <- pick_col(sp_df, c("estimate","emmean","prob","response","effect"), "estimate")
sp_se  <- pick_col(sp_df, c("SE","SE.dif","SE.diff"), "SE")
sp_p   <- pick_col(sp_df, c("p.value","p.value.mixed","p"), "p-value")
sp_l   <- pick_col(sp_ci, c("lower.CL","asymp.LCL","LCL"), "lower CI")
sp_u   <- pick_col(sp_ci, c("upper.CL","asymp.UCL","UCL"), "upper CI")

species_pairs_resp <- sp_df %>%
  dplyr::transmute(
    treatment = .data[[by_col_s]],
    contrast  = contrast,                      # e.g., "pocillopora - acropora"
    prob_diff = .data[[sp_est]],
    prob_SE   = .data[[sp_se]],
    prob_p    = .data[[sp_p]]
  ) %>%
  dplyr::left_join(
    sp_ci %>% dplyr::transmute(
      treatment = .data[[by_col_s]],
      contrast,
      prob_low  = .data[[sp_l]],
      prob_high = .data[[sp_u]]
    ),
    by = c("treatment","contrast")
  )

# B) Odds ratios among species (Tukey): pairwise on link, summarized as ORs on response
sp_or    <- emmeans::contrast(es_link, method = "tukey")
sp_or_df <- as.data.frame(summary(sp_or, type = "response"))
sp_or_ci <- as.data.frame(confint(sp_or,  type = "response"))

by_col_so <- get_by_col_safe(sp_or_df, c("treatment","by1","by"))
or_c      <- pick_col(sp_or_df, c("ratio","odds.ratio"), "odds ratio")
or_p      <- pick_col(sp_or_df, c("p.value","p.value.mixed","p"), "OR p-value")
or_l      <- pick_col(sp_or_ci, c("lower.CL","asymp.LCL","LCL"), "OR lower CI")
or_u      <- pick_col(sp_or_ci, c("upper.CL","asymp.UCL","UCL"), "OR upper CI")

species_pairs_OR <- sp_or_df %>%
  dplyr::transmute(
    treatment = .data[[by_col_so]],
    contrast,
    OR   = .data[[or_c]],
    OR_p = .data[[or_p]]
  ) %>%
  dplyr::left_join(
    sp_or_ci %>% dplyr::transmute(
      treatment = .data[[by_col_so]],
      contrast,
      OR_low  = .data[[or_l]],
      OR_high = .data[[or_u]]
    ),
    by = c("treatment","contrast")
  )

# C) Harmonize labels and join -> species_contr_df
species_pairs_resp$contrast <- gsub("-", "vs", species_pairs_resp$contrast, fixed = TRUE)
species_pairs_OR$contrast   <- gsub("/", "vs", species_pairs_OR$contrast,   fixed = TRUE)

species_contr_df <- dplyr::left_join(
  species_pairs_resp, species_pairs_OR, by = c("treatment","contrast")
)

# Quick preview & save
cat("\n6.4.3 preview — species within treatment:\n")
print(utils::head(species_contr_df, 6))
readr::write_csv(species_contr_df, here("output","tables","binary_emmeans_species_within_treatment.csv"))
# ---- 6.5 Endpoint-only robustness (OPTIONAL; comment out if not needed) ----
# If you want a clearly-time-agnostic estimand, re-fit at the final day only and export same summaries.
# 2026-06-30 (review-restats): a silent try(..., silent=TRUE) previously swallowed
# any failure here and still exited 0, so a broken endpoint block was
# indistinguishable from success. The analysis is unchanged; we now capture the
# error object and re-raise it LOUDLY so a regression is auditable.
endpoint_attempt <- try({
  end_day <- max(df$day, na.rm = TRUE)
  df_end  <- df %>% filter(day == end_day)
  mod_end <- glmer(
    healed01 ~ treatment * species + (1 | parent_id),
    data = df_end, family = binomial, control = glmerControl(optimizer="bobyqa")
  )
  # Export fixed effects (logit + OR)
  fixef_end <- broom.mixed::tidy(mod_end, effects = "fixed", conf.int = TRUE) %>%
    mutate(OR = exp(estimate), OR_low = exp(conf.low), OR_high = exp(conf.high))
  readr::write_csv(fixef_end, here("output","tables","endpoint_glmm_fixed_effects_logit_and_OR.csv"))
  # Export emmeans cell probs + treatment-within-species contrasts
  emm_cell_end <- as.data.frame(emmeans::emmeans(mod_end, ~ species * treatment, type = "response"))
  readr::write_csv(emm_cell_end, here("output","tables","endpoint_emmeans_cell_probabilities.csv"))
  emm_tr_end <- emmeans::emmeans(mod_end, ~ treatment | species, type = "response")
  tr_end_df  <- as.data.frame(summary(contrast(emm_tr_end, method = "revpairwise")))
  readr::write_csv(tr_end_df, here("output","tables","endpoint_emmeans_treatment_within_species.csv"))
}, silent = TRUE)
if (inherits(endpoint_attempt, "try-error")) {
  warning("Endpoint-only (day-28) robustness block FAILED: ",
          conditionMessage(attr(endpoint_attempt, "condition")), call. = FALSE)
  stop("Endpoint-only robustness block failed — see warning above (was silently swallowed before).")
}

# ---- 6.6 Master manuscript-stats table (single source of truth) ----
# Goal: ONE tidy table where every reportable statistic in the manuscript Results
# section can be looked up. Schema:
#   section       — grouping label (Sample size / Observed counts / Model comparison /
#                   Main effects / Fixed-effect ORs / Cell probabilities / Contrasts /
#                   Diagnostics)
#   outcome       — healed / regenerated / debris / design
#   model         — "Binomial GLMM (additive)" / "Binomial GLMM (interaction)" / "Descriptive"
#   term          — parameter / contrast label
#   n             — sample size relevant to that row
#   estimate      — numeric primary estimate (point value reported in text)
#   estimate_type — "OR" / "Pr(healed)" / "Pr(regenerated)" / "Δ probability" / "logit β" /
#                   "χ²" / "count" / "proportion" / "AIC" / "R² (marginal)" etc.
#   ci_lower / ci_upper — 95% Wald CI (NA when not applicable)
#   effect_pct    — effect size as percent: 100 × Δ probability for prob diffs;
#                   100 × (OR − 1) for ORs; signed
#   fold_change   — multiplicative fold (= OR for OR rows; NA otherwise)
#   test_stat     — χ² or z value (NA when not applicable)
#   df            — degrees of freedom (NA when not applicable)
#   p_value       — numeric p
#   p_formatted   — pretty string ("<0.001", "=0.042")
#   notes         — context (additive model, Tukey-adjusted, posterior marginal, etc.)

# Helper: format p-values for manuscript text
fmt_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("<0.001")
  if (p < 0.01)  return(formatC(p, format = "f", digits = 3))
  formatC(p, format = "f", digits = 3)
}
fmt_p_v <- function(p) vapply(p, fmt_p, character(1))

# Helper: build a single row of the master table
mk_row <- function(section, outcome, model, term, n = NA_integer_,
                   estimate = NA_real_, estimate_type = NA_character_,
                   ci_lower = NA_real_, ci_upper = NA_real_,
                   effect_pct = NA_real_, fold_change = NA_real_,
                   test_stat = NA_real_, df = NA_integer_,
                   p_value = NA_real_, notes = NA_character_) {
  tibble(
    section, outcome, model, term, n,
    estimate, estimate_type, ci_lower, ci_upper,
    effect_pct, fold_change, test_stat, df,
    p_value, p_formatted = fmt_p_v(p_value), notes
  )
}

# ---- A. Design / sample-size rows ----
n_rows <- bind_rows(
  mk_row("Sample size", "design", "Descriptive", "Parent colonies (total)",
         n = n_distinct(df$parent_id), estimate = n_distinct(df$parent_id),
         estimate_type = "count",
         notes = "Each parent split into airbrush + scrape fragments (paired design)"),
  mk_row("Sample size", "design", "Descriptive", "Coral fragments (total)",
         n = n_distinct(df$coral_id), estimate = n_distinct(df$coral_id),
         estimate_type = "count"),
  mk_row("Sample size", "design", "Descriptive", "Acropora parent colonies",
         n = n_distinct(df$parent_id[df$species == "acropora"]),
         estimate = n_distinct(df$parent_id[df$species == "acropora"]),
         estimate_type = "count"),
  mk_row("Sample size", "design", "Descriptive", "Pocillopora parent colonies",
         n = n_distinct(df$parent_id[df$species == "pocillopora"]),
         estimate = n_distinct(df$parent_id[df$species == "pocillopora"]),
         estimate_type = "count"),
  mk_row("Sample size", "design", "Descriptive", "Porites parent colonies",
         n = n_distinct(df$parent_id[df$species == "porites"]),
         estimate = n_distinct(df$parent_id[df$species == "porites"]),
         estimate_type = "count"),
  mk_row("Sample size", "design", "Descriptive", "Photographic timepoints per fragment",
         n = length(unique(df$day)), estimate = length(unique(df$day)),
         estimate_type = "count",
         notes = paste("days =", paste(sort(unique(df$day)), collapse = ", "))),
  mk_row("Sample size", "healed", "Descriptive", "Observations with non-NA healed",
         n = sum(!is.na(df$healed01)), estimate = sum(!is.na(df$healed01)),
         estimate_type = "count",
         notes = paste0(sum(is.na(df$healed01)), " NA from photo errors")),
  mk_row("Sample size", "regenerated", "Descriptive", "Observations with non-NA regenerated",
         n = if (exists("df_regen")) nrow(df_regen) else NA_integer_,
         estimate = if (exists("df_regen")) nrow(df_regen) else NA_real_,
         estimate_type = "count",
         notes = "6 NA from images where polyp structure could not be assessed")
)

# ---- B. Observed-count rows (fill manuscript placeholders like "XX/11") ----
# Algae prevalence: did each fragment EVER show debris/algae across the timecourse?
algae_ever <- df %>%
  group_by(treatment, coral_id) %>%
  summarise(ever_algae = any(debris == "yes", na.rm = TRUE), .groups = "drop") %>%
  group_by(treatment) %>%
  summarise(n_fragments = dplyr::n(),
            n_with_algae = sum(ever_algae),
            prop = n_with_algae / n_fragments, .groups = "drop")

algae_rows <- algae_ever %>%
  rowwise() %>%
  do(mk_row(
    section = "Observed counts",
    outcome = "debris",
    model   = "Descriptive",
    term    = paste0(stringr::str_to_title(.$treatment),
                     " fragments with algae/debris at any timepoint"),
    n       = .$n_fragments,
    estimate = .$n_with_algae,
    estimate_type = "count of fragments",
    notes   = sprintf("%d / %d = %.0f%%", .$n_with_algae, .$n_fragments, 100 * .$prop)
  )) %>%
  bind_rows()

# Fully healed at day 28
heal28 <- df %>%
  filter(day == max(day)) %>%
  group_by(treatment) %>%
  summarise(n_fragments = sum(!is.na(healed01)),
            n_healed = sum(healed01 == 1L, na.rm = TRUE),
            prop = n_healed / n_fragments, .groups = "drop")

heal28_rows <- heal28 %>%
  rowwise() %>%
  do(mk_row(
    section = "Observed counts",
    outcome = "healed",
    model   = "Descriptive",
    term    = paste0(stringr::str_to_title(.$treatment), " fragments fully healed at day 28"),
    n       = .$n_fragments,
    estimate = .$n_healed,
    estimate_type = "count of fragments",
    effect_pct = 100 * .$prop,
    notes   = sprintf("%d / %d = %.0f%%", .$n_healed, .$n_fragments, 100 * .$prop)
  )) %>%
  bind_rows()

# Regenerated ever (per fragment)
if (exists("df_regen")) {
  regen_ever <- df_regen %>%
    group_by(treatment, coral_id) %>%
    summarise(ever_regen = any(regen01 == 1L), .groups = "drop") %>%
    group_by(treatment) %>%
    summarise(n_fragments = dplyr::n(),
              n_regen = sum(ever_regen),
              prop = n_regen / n_fragments, .groups = "drop")

  regen_ever_rows <- regen_ever %>%
    rowwise() %>%
    do(mk_row(
      section = "Observed counts",
      outcome = "regenerated",
      model   = "Descriptive",
      term    = paste0(stringr::str_to_title(.$treatment),
                       " fragments with polyp regeneration at any timepoint"),
      n       = .$n_fragments,
      estimate = .$n_regen,
      estimate_type = "count of fragments",
      effect_pct = 100 * .$prop,
      notes   = sprintf("%d / %d = %.0f%%", .$n_regen, .$n_fragments, 100 * .$prop)
    )) %>%
    bind_rows()
} else {
  regen_ever_rows <- tibble()
}

obs_rows <- bind_rows(algae_rows, heal28_rows, regen_ever_rows)

# ---- C. HEALED model rows ----

# C1. LRT interaction vs additive
healed_lrt_row <- mk_row(
  section = "Model comparison",
  outcome = "healed",
  model   = "Binomial GLMM",
  term    = "Treatment × Species interaction",
  test_stat = lrt_sp$Chisq[2],
  df      = lrt_sp$Df[2],
  p_value = lrt_sp$`Pr(>Chisq)`[2],
  estimate = lrt_sp$Chisq[2],
  estimate_type = "χ²",
  notes   = "anova(mod_add, mod_int); additive preferred if NS"
)

# C2. Main effects (drop1 on additive)
healed_drop_rows <- drop_df %>%
  rowwise() %>%
  do(mk_row(
    section = "Main effects",
    outcome = "healed",
    model   = "Binomial GLMM (additive)",
    term    = .$Term,
    test_stat = .$`Chi-square`,
    df      = .$DF,
    p_value = .$`P-value`,
    estimate = .$`Chi-square`,
    estimate_type = "χ²",
    notes   = "drop1(mod_add)"
  )) %>%
  bind_rows()

# C3. Fixed-effect ORs (additive)  — these are the canonical OR lines in the manuscript
healed_or_rows <- fixed_df %>%
  rowwise() %>%
  do(mk_row(
    section = "Fixed-effect ORs",
    outcome = "healed",
    model   = "Binomial GLMM (additive)",
    term    = .$Term,
    estimate = .$OR,
    estimate_type = "OR",
    ci_lower = .$`CI Lower`,
    ci_upper = .$`CI Upper`,
    effect_pct = 100 * (.$OR - 1),
    fold_change = .$OR,
    p_value = .$`P-value`,
    notes   = "SUPPORTING (timepoint-level, pseudoreplicated; singular RE, temporal autocorrelation). Primary = fragment-level (exp1_fragment_level_primary.csv). For healed, lead with the Firth row (separation-robust)."
  )) %>%
  bind_rows()

# C4. Cell probabilities (Pr healed; emmeans on additive)
healed_cell_rows <- emm_cell_df %>%
  rowwise() %>%
  do(mk_row(
    section = "Cell probabilities",
    outcome = "healed",
    model   = "Binomial GLMM (additive)",
    term    = paste0(.$species, " × ", .$treatment),
    estimate = .$prob,
    estimate_type = "Pr(healed)",
    ci_lower = .$lower,
    ci_upper = .$upper,
    effect_pct = 100 * .$prob,
    notes   = "marginal probability, emmeans response scale"
  )) %>%
  bind_rows()

# C5. Treatment within species (probability difference)
healed_trt_rows <- treat_contr_df %>%
  rowwise() %>%
  do(mk_row(
    section = "Contrast: treatment within species",
    outcome = "healed",
    model   = "Binomial GLMM (additive)",
    term    = paste0(.$species, ": ", .$contrast),
    estimate = .$prob_diff,
    estimate_type = "Δ probability",
    ci_lower = .$prob_low,
    ci_upper = .$prob_high,
    effect_pct = 100 * .$prob_diff,
    p_value = .$prob_p,
    notes   = "Dremel − Airbrush; emmeans response scale"
  )) %>%
  bind_rows()

# C6. Species within treatment (Tukey, probability difference)
healed_spp_rows <- species_contr_df %>%
  rowwise() %>%
  do(mk_row(
    section = "Contrast: species within treatment",
    outcome = "healed",
    model   = "Binomial GLMM (additive)",
    term    = paste0(.$treatment, ": ", .$contrast),
    estimate = .$prob_diff,
    estimate_type = "Δ probability",
    ci_lower = .$prob_low,
    ci_upper = .$prob_high,
    effect_pct = 100 * .$prob_diff,
    p_value = .$prob_p,
    notes   = "Tukey-adjusted; emmeans response scale"
  )) %>%
  bind_rows()

# C7. Diagnostics
healed_diag_rows <- bind_rows(
  mk_row("Diagnostics", "healed", "Binomial GLMM (additive)", "AIC",
         estimate = AIC(mod_add), estimate_type = "AIC", notes = "lower is better"),
  mk_row("Diagnostics", "healed", "Binomial GLMM (interaction)", "AIC",
         estimate = AIC(mod_int), estimate_type = "AIC"),
  mk_row("Diagnostics", "healed", "Binomial GLMM (additive)", "Marginal R²",
         estimate = r2_add[1, "R2m"], estimate_type = "R² (marginal)",
         notes = "fixed effects only"),
  mk_row("Diagnostics", "healed", "Binomial GLMM (additive)", "Conditional R²",
         estimate = r2_add[1, "R2c"], estimate_type = "R² (conditional)",
         notes = "fixed + random"),
  mk_row("Diagnostics", "healed", "Binomial GLMM (additive)", "Singular fit",
         estimate = as.numeric(sing_add), estimate_type = "logical",
         notes = ifelse(sing_add, "TRUE — RE variances ≈ 0", "FALSE")),
  mk_row("Diagnostics", "healed", "Binomial GLMM (additive)", "Overdispersion",
         estimate = od_add, estimate_type = "ratio",
         notes = "Pearson; ~1 = OK")
)

# ---- D. REGENERATED model rows ----
if (exists("mod_regen_add")) {

  regen_lrt_row <- mk_row(
    section = "Model comparison", outcome = "regenerated",
    model   = "Binomial GLMM",
    term    = "Treatment × Species interaction",
    test_stat = lrt_regen$Chisq[2], df = lrt_regen$Df[2],
    p_value = lrt_regen$`Pr(>Chisq)`[2],
    estimate = lrt_regen$Chisq[2], estimate_type = "χ²",
    notes   = "anova(mod_regen_add, mod_regen_int)"
  )

  regen_drop_rows <- d1_regen %>%
    rowwise() %>%
    do(mk_row(
      section = "Main effects", outcome = "regenerated",
      model   = "Binomial GLMM (additive)",
      term    = .$Term,
      test_stat = .$`Chi-square`, df = .$DF, p_value = .$`P-value`,
      estimate = .$`Chi-square`, estimate_type = "χ²",
      notes   = "drop1(mod_regen_add)"
    )) %>%
    bind_rows()

  regen_or_rows <- fixed_regen %>%
    rowwise() %>%
    do(mk_row(
      section = "Fixed-effect ORs", outcome = "regenerated",
      model   = "Binomial GLMM (additive)",
      term    = .$Term,
      estimate = .$OR, estimate_type = "OR",
      ci_lower = .$`CI Lower`, ci_upper = .$`CI Upper`,
      effect_pct = 100 * (.$OR - 1),
      fold_change = .$OR,
      p_value = .$`P-value`,
      notes = "SUPPORTING (timepoint-level, pseudoreplicated); CI width reflects near-separation in airbrush group. Primary = fragment-level (exp1_fragment_level_primary.csv)."
    )) %>%
    bind_rows()

  regen_cell_rows <- cells_regen_df %>%
    rowwise() %>%
    do(mk_row(
      section = "Cell probabilities", outcome = "regenerated",
      model   = "Binomial GLMM (additive)",
      term    = paste0(.$species, " × ", .$treatment),
      estimate = .$prob, estimate_type = "Pr(regenerated)",
      ci_lower = .$asymp.LCL, ci_upper = .$asymp.UCL,
      effect_pct = 100 * .$prob,
      notes = "marginal probability, emmeans response scale"
    )) %>%
    bind_rows()

  regen_diag_rows <- bind_rows(
    mk_row("Diagnostics", "regenerated", "Binomial GLMM (additive)", "AIC",
           estimate = diag_regen$AIC, estimate_type = "AIC"),
    mk_row("Diagnostics", "regenerated", "Binomial GLMM (additive)", "Marginal R²",
           estimate = diag_regen$R2_marginal, estimate_type = "R² (marginal)"),
    mk_row("Diagnostics", "regenerated", "Binomial GLMM (additive)", "Conditional R²",
           estimate = diag_regen$R2_conditional, estimate_type = "R² (conditional)"),
    mk_row("Diagnostics", "regenerated", "Binomial GLMM (additive)", "Singular fit",
           estimate = as.numeric(diag_regen$singular_fit), estimate_type = "logical",
           notes = ifelse(diag_regen$singular_fit, "TRUE", "FALSE")),
    mk_row("Diagnostics", "regenerated", "Binomial GLMM (additive)", "Overdispersion",
           estimate = diag_regen$overdisp, estimate_type = "ratio")
  )

  regen_block <- bind_rows(regen_lrt_row, regen_drop_rows, regen_or_rows,
                           regen_cell_rows, regen_diag_rows)
} else {
  regen_block <- tibble()
}

# ---- D2. Firth-penalized sensitivity rows (both outcomes) ----
firth_to_rows <- function(firth_df, outcome_label) {
  firth_df %>%
    rowwise() %>%
    do(mk_row(
      section = "Sensitivity (Firth)",
      outcome = outcome_label,
      model   = "Firth-penalized logistic (brglm2)",
      term    = .$Term,
      estimate = .$OR,
      estimate_type = "OR",
      ci_lower = .$`CI Lower`,
      ci_upper = .$`CI Upper`,
      effect_pct = 100 * (.$OR - 1),
      fold_change = .$OR,
      p_value = .$`P-value`,
      notes = "SUPPORTING (timepoint-level), separation-robust; profile-likelihood CI; fixed-effects only (RE variance ≈ 0 in main GLMM). For healed this Firth row LEADS the supporting reporting. Primary = fragment-level (exp1_fragment_level_primary.csv)."
    )) %>%
    bind_rows()
}

firth_rows <- bind_rows(
  if (exists("firth_healed")) firth_to_rows(firth_healed, "healed") else tibble(),
  if (exists("firth_regen"))  firth_to_rows(firth_regen,  "regenerated") else tibble()
)

# ---- D3. Manuscript-sourced rows for Experiments 2–4 ----
# Stats from `data/manuscript_external_stats.csv` — these come from the
# manuscript draft (2026-02-16) and represent experiments NOT in this repo:
#   Exp 2: longitudinal Porites (glmmTMB with cubic time splines)
#   Exp 3: single-polyp regeneration timing (qualitative)
#   Exp 4: histology (qualitative)
# The underlying datasets live elsewhere; provenance is preserved in the
# `notes` field so the master table is honest about which numbers we computed.
ext_path <- here("data","manuscript_external_stats.csv")
if (file.exists(ext_path)) {
  ext_raw <- readr::read_csv(ext_path, show_col_types = FALSE)
  # 2026-06-30 (review-restats): DROP the Experiment 2 external rows. They were
  # transcribed from a stale 2026-02-16 manuscript draft and CONTRADICT the
  # in-repo Exp 2 re-analysis (exp2_*; penalized binomial GLMM). Keep only the
  # genuinely external, qualitative Exp 3 / Exp 4 rows (no in-repo dataset).
  n_before <- nrow(ext_raw)
  ext_raw <- ext_raw %>% dplyr::filter(!grepl("^2$|exp.?2|experiment.?2",
                                              tolower(as.character(experiment))))
  cat("External stats: dropped", n_before - nrow(ext_raw),
      "stale Exp 2 row(s); kept", nrow(ext_raw), "Exp 3/4 qualitative row(s).\n")
  ext_rows <- ext_raw %>%
    rowwise() %>%
    do(mk_row(
      section = paste0("External: ", .$experiment),
      outcome = .$outcome,
      model   = .$model,
      term    = paste0("[", .$species, "] ", .$term),
      n       = NA_integer_,
      estimate = suppressWarnings(as.numeric(.$estimate)),
      estimate_type = .$estimate_type,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      effect_pct = NA_real_,
      fold_change = NA_real_,
      test_stat = suppressWarnings(as.numeric(.$test_stat)),
      df       = suppressWarnings(as.integer(.$df)),
      p_value  = suppressWarnings(as.numeric(.$p_value)),
      notes    = .$notes
    )) %>%
    bind_rows()
  cat("Loaded", nrow(ext_rows), "manuscript-external rows for Exp 2-4.\n")
} else {
  ext_rows <- tibble()
  message("No external manuscript-stats CSV found at: ", ext_path)
}

# ---- E. Assemble master table ----
results_master <- bind_rows(
  n_rows,
  obs_rows,
  healed_lrt_row,
  healed_drop_rows,
  healed_or_rows,
  healed_cell_rows,
  healed_trt_rows,
  healed_spp_rows,
  healed_diag_rows,
  regen_block,
  firth_rows,
  ext_rows
) %>%
  mutate(
    estimate    = ifelse(is.na(estimate),    NA_real_, signif(estimate,    4)),
    ci_lower    = ifelse(is.na(ci_lower),    NA_real_, signif(ci_lower,    4)),
    ci_upper    = ifelse(is.na(ci_upper),    NA_real_, signif(ci_upper,    4)),
    effect_pct  = ifelse(is.na(effect_pct),  NA_real_, signif(effect_pct,  4)),
    fold_change = ifelse(is.na(fold_change), NA_real_, signif(fold_change, 4)),
    test_stat   = ifelse(is.na(test_stat),   NA_real_, signif(test_stat,   4)),
    p_value     = ifelse(is.na(p_value),     NA_real_, signif(p_value,     3))
  )

readr::write_csv(results_master, here("output","tables","paper_results_master.csv"))

# Pretty HTML for SI / collaborator share
results_master %>%
  gt(groupname_col = "section") %>%
  cols_hide(columns = c("p_value")) %>%   # show p_formatted instead
  cols_label(
    outcome       = "Outcome",
    model         = "Model",
    term          = "Parameter / Contrast",
    n             = "n",
    estimate      = "Estimate",
    estimate_type = "Type",
    ci_lower      = "CI Lower",
    ci_upper      = "CI Upper",
    effect_pct    = "Effect (%)",
    fold_change   = "Fold change",
    test_stat     = html("Test stat (χ² / z)"),
    df            = "df",
    p_formatted   = "p",
    notes         = "Notes"
  ) %>%
  fmt_number(columns = c("estimate","ci_lower","ci_upper","effect_pct",
                         "fold_change","test_stat"),
             decimals = 3, drop_trailing_zeros = TRUE) %>%
  tab_header(
    title = md("**Master Statistics Table — Coral Wound Healing (Experiment 1)**"),
    subtitle = md("All test statistics, effect sizes, and key values reported in the manuscript")
  ) %>%
  tab_source_note(md(paste(
    "Effect (%) for probabilities = 100 × Δ probability.",
    "Effect (%) for odds ratios = 100 × (OR − 1).",
    "Fold change = OR for OR rows."
  ))) %>%
  gtsave(here("output","tables","paper_results_master.html"))

cat("\n=== MASTER STATS TABLE ===\n")
cat("Wrote:", here("output","tables","paper_results_master.csv"), "\n")
cat("Wrote:", here("output","tables","paper_results_master.html"), "\n")
cat("Total rows:", nrow(results_master), "\n")
cat("Sections:", paste(unique(results_master$section), collapse = " | "), "\n")
cat("Outcomes:", paste(unique(results_master$outcome), collapse = " | "), "\n\n")



# ---- 7. Session info (for reproducibility) ----
sink(here("output","tables","sessionInfo.txt"))
print(sessionInfo())
sink()


# =============================================================================
# 8. Narrative summary of key findings (auto-generated prose)
# =============================================================================

suppressPackageStartupMessages({ library(glue); library(dplyr); library(stringr) })

dir.create(here("output","text"), recursive = TRUE, showWarnings = FALSE)

# --- Helpers -----------------------------------------------------------------
fmt_p <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 1e-4, "< 1e-4",
                ifelse(p < 0.001, "< 0.001",
                       sprintf("= %.3f", p))))
}
fmt_ci  <- function(lo, hi, digits = 2) glue("[{round(lo, digits)}, {round(hi, digits)}]")
fmt_pct <- function(x, digits = 1) ifelse(is.na(x), "NA", paste0(round(100*x, digits), "%"))
fmt_num <- function(x, d=3) ifelse(is.na(x), "NA", format(round(x, d), nsmall = d))
nz <- function(x) ifelse(is.na(x), 0, x)

safe_slice1 <- function(df) {
  if (nrow(df) == 0) return(df)
  dplyr::slice_head(df, n = 1)
}

# --- Pull key objects (defensively) ------------------------------------------
diag_add <- diag_tbl %>% dplyr::filter(model == "mod_add") %>% safe_slice1()
diag_int <- diag_tbl %>% dplyr::filter(model == "mod_int") %>% safe_slice1()

# LRT interaction
lrt_stat <- nz(lrt_sp$Chisq[2]); lrt_df <- nz(lrt_sp$Df[2]); lrt_p <- lrt_sp$`Pr(>Chisq)`[2]

# Cell means: best & worst
best_cell <- emm_cell_df %>%
  arrange(desc(prob)) %>%
  mutate(line = glue("{species} × {treatment}: {fmt_pct(prob)} (95% CI {fmt_ci(lower, upper)})")) %>%
  safe_slice1()

worst_cell <- emm_cell_df %>%
  arrange(prob) %>%
  mutate(line = glue("{species} × {treatment}: {fmt_pct(prob)} (95% CI {fmt_ci(lower, upper)})")) %>%
  safe_slice1()

# Treatment within species (Dremel vs Airbrush)
treat_sig <- treat_contr_df %>%
  mutate(sig = !is.na(prob_p) & prob_p < 0.05,
         dir = case_when(prob_diff > 0 ~ "higher", prob_diff < 0 ~ "lower", TRUE ~ "no difference")) %>%
  arrange(prob_p)

treat_lines_sig <- if (nrow(treat_sig)) {
  treat_sig %>%
    transmute(line = glue(
      "- Within {species}, Dremel had {dir} Pr(healed) vs Airbrush by {fmt_pct(abs(prob_diff))} ",
      "(95% CI {fmt_pct(abs(prob_low))}, {fmt_pct(abs(prob_high))}); OR = {round(OR,2)} ",
      "(95% CI {round(OR_low,2)}, {round(OR_high,2)}), p {fmt_p(prob_p)}"
    )) %>%
    pull(line)
} else character(0)

# Species within treatment (Tukey)
spp_sig <- species_contr_df %>%
  mutate(sig = !is.na(prob_p) & prob_p < 0.05,
         dir = case_when(prob_diff > 0 ~ "higher", prob_diff < 0 ~ "lower", TRUE ~ "no difference")) %>%
  arrange(prob_p)

spp_lines_sig <- if (nrow(spp_sig)) {
  spp_sig %>%
    transmute(line = glue(
      "- Under {treatment}, {contrast} in Pr(healed) was {fmt_pct(abs(prob_diff))} ",
      "(95% CI {fmt_pct(abs(prob_low))}, {fmt_pct(abs(prob_high))}); OR = {round(OR,2)} ",
      "(95% CI {round(OR_low,2)}, {round(OR_high,2)}), p {fmt_p(prob_p)}"
    )) %>%
    pull(line)
} else character(0)

# Main effects (drop1) — sourced from drop_df (the section-6.6 master rewrite
# replaced the old `main_effects` intermediate tibble)
main_lines <- if (exists("drop_df") && nrow(drop_df)) {
  drop_df %>%
    transmute(line = glue("- {Term}: χ²({DF}) = {fmt_num(`Chi-square`,3)}, p {fmt_p(`P-value`)}")) %>%
    pull(line)
} else character(0)

# Diagnostics text
diag_text <- if (nrow(diag_add) && nrow(diag_int)) {
  glue(
    "- Additive model (mod_add): overdispersion = {fmt_num(diag_add$overdisp,3)}, ",
    "singular = {diag_add$is_singular}, R²_m = {fmt_num(diag_add$R2_marginal,3)}, R²_c = {fmt_num(diag_add$R2_conditional,3)}, ",
    "AIC = {fmt_num(diag_add$AIC,2)}, BIC = {fmt_num(diag_add$BIC,2)}.\n",
    "- Interaction model (mod_int): overdispersion = {fmt_num(diag_int$overdisp,3)}, ",
    "singular = {diag_int$is_singular}, R²_m = {fmt_num(diag_int$R2_marginal,3)}, R²_c = {fmt_num(diag_int$R2_conditional,3)}, ",
    "AIC = {fmt_num(diag_int$AIC,2)}, BIC = {fmt_num(diag_int$BIC,2)}."
  )
} else "- Diagnostics table not available."

# ICCs text
icc_parent <- icc_tbl %>% dplyr::filter(component == "parent_id") %>% pull(ICC) %>% nz()
icc_coral  <- icc_tbl %>% dplyr::filter(component == "coral_id_within_parent") %>% pull(ICC) %>% nz()
icc_text <- glue("- Latent-scale ICCs: parent_id = {fmt_num(icc_parent,3)}, coral_id-within-parent = {fmt_num(icc_coral,3)}.")

# --- Build narrative ----------------------------------------------------------
header <- "# Summary of Key Findings\n"
model_cmp <- glue(
  "## Model comparison\n",
  "- Treatment × species interaction: χ²({lrt_df}) = {fmt_num(lrt_stat,3)}, p {fmt_p(lrt_p)}.\n",
  "  (From `anova(mod_add, mod_int)`.)\n"
)

main_eff <- if (length(main_lines)) {
  paste0("## Main effects (drop1 on mod_add)\n", paste(main_lines, collapse = "\n"), "\n")
} else ""

cells <- if (nrow(best_cell) && nrow(worst_cell)) {
  glue(
    "## Marginal probabilities (emmeans; response scale)\n",
    "- Highest Pr(healed): {best_cell$line}\n",
    "- Lowest  Pr(healed): {worst_cell$line}\n"
  )
} else "## Marginal probabilities\n- Not available.\n"

treat_sec <- if (length(treat_lines_sig)) {
  paste0("## Treatment contrast within species (Dremel vs Airbrush)\n",
         paste(treat_lines_sig, collapse = "\n"), "\n")
} else {
  "## Treatment contrast within species (Dremel vs Airbrush)\n- No species showed a statistically significant difference at α = 0.05.\n"
}

spp_sec <- if (length(spp_lines_sig)) {
  paste0("## Species contrasts within treatment (Tukey-adjusted)\n",
         paste(spp_lines_sig, collapse = "\n"), "\n")
} else {
  "## Species contrasts within treatment (Tukey-adjusted)\n- No significant pairwise species differences at α = 0.05 within any treatment.\n"
}

diag_sec <- glue("## Model diagnostics\n{diag_text}\n{icc_text}\n")

headline <- if (nrow(best_cell)) {
  glue("**Headline:** The highest healing probability was observed for {best_cell$species} treated with {best_cell$treatment} ",
       "({fmt_pct(best_cell$prob)}; 95% CI {fmt_ci(best_cell$lower, best_cell$upper)}).")
} else "**Headline:** Summary not available (missing emmeans cell estimates)."

report_md <- paste(
  header,
  headline, "\n",
  model_cmp,
  main_eff,
  cells,
  treat_sec,
  spp_sec,
  diag_sec,
  sep = "\n"
)

# --- Write files (Markdown + plain text) -------------------------------------
md_path  <- here("output","text","paper_results_summary.md")
txt_path <- here("output","text","paper_results_summary.txt")
writeLines(report_md, md_path)
writeLines(report_md, txt_path)  # plain text = same content

# Console preview
cat("\n===== SUMMARY (first ~40 lines) =====\n")
cat(paste0(paste(readLines(md_path, n = 40), collapse = "\n"), "\n"))
cat("\nSaved:\n  - ", md_path, "\n  - ", txt_path, "\n", sep = "")
