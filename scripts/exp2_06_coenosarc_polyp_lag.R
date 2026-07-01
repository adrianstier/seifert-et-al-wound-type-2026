# =============================================================================
# Experiment 2 — Coenosarc-to-polyp lag analysis
# -----------------------------------------------------------------------------
# Added 2026-05-26 to address reviewer concern that the manuscript's
# "polyp regeneration approximately one week after coenosarc closure" claim
# was unfalsifiable (presented as a quantitative result but supported only by
# qualitative imaging timing).
#
# This script computes, for each Porites fragment that achieved both coenosarc
# closure AND polyp regeneration during the experiment, the per-fragment lag
# between the first day the wound was scored as having complete coenosarc
# coverage and the first day a polyp was scored in the wound center.
#
# Outputs:
#   output/exp2_tables/coenosarc_to_polyp_lag_per_fragment.csv
#   output/exp2_tables/coenosarc_to_polyp_lag_summary.csv  (n, median, IQR, range)
#
# Sources: data/Porites microscope characterization - complete.xlsx (sheet "data")
# Outcomes used: coenosarc_coverage (yes/no) and polyp_in_center_of_wound (yes/no)
#
#   Rscript scripts/exp2_06_coenosarc_polyp_lag.R
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(readxl); library(dplyr); library(readr); library(tidyr)
})

xlsx_path <- here("data", "Porites microscope characterization - complete.xlsx")
df <- read_excel(xlsx_path, sheet = "data") %>%
  mutate(
    id    = as.character(id),
    day   = as.integer(day),
    coenosarc_yn = tolower(trimws(as.character(coenosarc_coverage))),
    polyp_yn     = tolower(trimws(as.character(polyp_in_center_of_wound)))
  ) %>%
  filter(!is.na(id), !is.na(day))

cat("Fragments in data:", n_distinct(df$id), "\n")
cat("Day range:", paste(range(df$day, na.rm = TRUE), collapse = " – "), "\n")

# First day each fragment achieved coenosarc coverage = yes
first_coen <- df %>%
  filter(coenosarc_yn == "yes") %>%
  group_by(id, treatment) %>%
  summarise(first_coenosarc_day = min(day, na.rm = TRUE), .groups = "drop")

# First day each fragment had a polyp in center of wound
first_polyp <- df %>%
  filter(polyp_yn == "yes") %>%
  group_by(id, treatment) %>%
  summarise(first_polyp_day = min(day, na.rm = TRUE), .groups = "drop")

# Per-fragment lag: only fragments achieving both
lag_per_fragment <- inner_join(first_coen, first_polyp,
                                by = c("id", "treatment")) %>%
  mutate(lag_days = first_polyp_day - first_coenosarc_day) %>%
  arrange(treatment, id)

write_csv(lag_per_fragment,
          here("output", "exp2_tables", "coenosarc_to_polyp_lag_per_fragment.csv"))

cat("\n--- Per-fragment lag (coenosarc → polyp) ---\n")
print(lag_per_fragment)

# Summary: overall and by treatment
summary_tbl <- lag_per_fragment %>%
  group_by(treatment) %>%
  summarise(
    n              = n(),
    median_lag     = median(lag_days, na.rm = TRUE),
    iqr_low        = quantile(lag_days, 0.25, na.rm = TRUE),
    iqr_high       = quantile(lag_days, 0.75, na.rm = TRUE),
    min_lag        = min(lag_days, na.rm = TRUE),
    max_lag        = max(lag_days, na.rm = TRUE),
    mean_lag       = mean(lag_days, na.rm = TRUE),
    sd_lag         = sd(lag_days, na.rm = TRUE),
    .groups = "drop"
  )

overall <- lag_per_fragment %>%
  summarise(
    treatment      = "ALL",
    n              = n(),
    median_lag     = median(lag_days, na.rm = TRUE),
    iqr_low        = quantile(lag_days, 0.25, na.rm = TRUE),
    iqr_high       = quantile(lag_days, 0.75, na.rm = TRUE),
    min_lag        = min(lag_days, na.rm = TRUE),
    max_lag        = max(lag_days, na.rm = TRUE),
    mean_lag       = mean(lag_days, na.rm = TRUE),
    sd_lag         = sd(lag_days, na.rm = TRUE)
  )

summary_tbl <- bind_rows(overall, summary_tbl)
write_csv(summary_tbl,
          here("output", "exp2_tables", "coenosarc_to_polyp_lag_summary.csv"))

cat("\n--- Lag summary (median + IQR) ---\n")
print(summary_tbl)

cat("\nWrote: output/exp2_tables/coenosarc_to_polyp_lag_per_fragment.csv\n")
cat("Wrote: output/exp2_tables/coenosarc_to_polyp_lag_summary.csv\n")
