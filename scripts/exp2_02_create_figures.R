# =============================================================================
# MAIN FIGURE SCRIPT: CREATE PUBLICATION-READY FIGURES (Experiment 2)
#
# PURPOSE:
#   Multi-outcome time-series small multiples for Porites spp. (7 outcomes)
#   -> main text.
#
#   Faceting (not 7 hand-assembled patchwork panels) is used deliberately:
#   per the project figure standard, >4/>6 groups should be faceted, which
#   gives equal panel widths, natively shared axes (labelled once), and a
#   single collected legend. Theme/palette/recoding come from 00_setup.R
#   (single source of truth shared with exp2_01/04/05).
#
# OUTPUT (PDF via cairo_pdf + PNG, sizes in mm, 300 dpi):
#   - output/exp2_figures_main/figure2_porites_all_outcomes.{pdf,png}
#   - companion legend text: *_legend.md next to the figure
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(here)
  library(scales)
})

# theme_pub(), TREATMENT_COLORS (Okabe-Ito), recode_outcome()
source(here("scripts", "00_setup.R"))

OUT_DIR_MAIN <- here("output", "exp2_figures_main")
dir.create(OUT_DIR_MAIN, showWarnings = FALSE, recursive = TRUE)

LINE_SIZE  <- 0.9
POINT_SIZE <- 2.0
ERROR_ALPHA <- 0.15
BASE_SIZE  <- 10  # double-column small multiples (project standard: 10-11)

# -------------------- Load data ---------------------------------------------
cat("\n", strrep("=", 80), "\n LOADING DATA\n", strrep("=", 80), "\n\n", sep = "")

load_species <- function(file, sp) {
  read_excel(here("data", file), sheet = "data") %>%
    clean_names() %>%
    mutate(
      species = sp,
      treatment = tolower(trimws(as.character(treatment))),
      treatment_label = factor(case_when(
        treatment == "airbrush" ~ "Airbrush",
        treatment %in% c("drem", "dremel") ~ "Scrape",
        treatment == "air_drem" ~ "Airbrush-Scrape"
      ), levels = c("Airbrush", "Scrape", "Airbrush-Scrape")),
      day = as.numeric(day),
      id = factor(id)
    ) %>%
    filter(!is.na(treatment_label), !is.na(day))
}

dat_porites  <- load_species("Porites microscope characterization - complete.xlsx",
                             "Porites spp.")

cat("[OK] Loaded", nrow(dat_porites), "Porites observations\n\n")

# -------------------- Long summary across outcomes --------------------------
# n = 1 cells get se = NA (no ribbon drawn) rather than se = 0. A zero-width
# band would falsely imply perfect certainty from a single binary observation.
build_outcome_summary <- function(data, outcomes) {
  labels <- vapply(outcomes, `[[`, character(1), "label")
  purrr::imap_dfr(outcomes, function(info, nm) {
    data %>%
      mutate(value = recode_outcome(.data[[info$col]])) %>%
      filter(!is.na(value)) %>%
      group_by(day, treatment_label) %>%
      summarise(n = dplyr::n(),
                mean = mean(value),
                se = if (dplyr::n() < 2) NA_real_ else sd(value) / sqrt(dplyr::n()),
                .groups = "drop") %>%
      mutate(outcome = info$label)
  }) %>%
    mutate(
      outcome = factor(outcome, levels = labels),
      lower = pmax(mean - se, 0),
      upper = pmin(mean + se, 1)
    )
}

make_species_fig <- function(df, ncol) {
  ggplot(df, aes(day, mean * 100,
                 color = treatment_label, fill = treatment_label)) +
    geom_ribbon(aes(ymin = lower * 100, ymax = upper * 100),
                alpha = ERROR_ALPHA, color = NA, na.rm = TRUE) +
    geom_line(linewidth = LINE_SIZE, alpha = 0.9, na.rm = TRUE) +
    geom_point(size = POINT_SIZE, shape = 21, fill = "white",
               stroke = 0.9, na.rm = TRUE) +
    scale_color_manual(values = setNames(unname(TREATMENT_COLORS), c("Airbrush", "Scrape", "Airbrush-Scrape"))) +
    scale_fill_manual(values = setNames(unname(TREATMENT_COLORS), c("Airbrush", "Scrape", "Airbrush-Scrape"))) +
    scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100),
                       expand = expansion(mult = c(0, 0.02))) +
    facet_wrap(~ outcome, ncol = ncol,
               labeller = label_wrap_gen(18)) +
    labs(x = "Day", y = "Corals with feature (%)") +
    theme_pub(base_size = BASE_SIZE)
}

# -------------------- Figure 2: Porites (main text) -------------------------
cat(strrep("=", 80), "\n CREATING MAIN TEXT FIGURE: PORITES\n",
    strrep("=", 80), "\n\n", sep = "")

PORITES_OUTCOMES <- list(
  list(col = "coenosarc_coverage",       label = "Coenosarc coverage"),
  list(col = "polyp_in_center_of_wound", label = "Regenerated"),
  list(col = "algal_plug",               label = "Algal plug"),
  list(col = "yellow_aggregations",      label = "Yellow aggregations"),
  list(col = "pink",                     label = "Pink aggregations"),
  list(col = "rfp",                      label = "RFP")
)

por_df  <- build_outcome_summary(dat_porites, PORITES_OUTCOMES)
fig_por <- make_species_fig(por_df, ncol = 3)

ggsave(file.path(OUT_DIR_MAIN, "figure2_porites_all_outcomes.pdf"),
       plot = fig_por, width = 180, height = 170, units = "mm",
       dpi = 300, device = cairo_pdf, bg = "white")
ggsave(file.path(OUT_DIR_MAIN, "figure2_porites_all_outcomes.png"),
       plot = fig_por, width = 180, height = 170, units = "mm",
       dpi = 300, bg = "white")

writeLines(c(
  "**Manuscript correspondence:** this is the data panel for **manuscript",
  "Figure 3B** (NOT manuscript Figure 2). The selected **6** Porites",
  "outcomes (coenosarc coverage, regenerated, algal plug, yellow",
  "aggregations, pink aggregations, RFP) match the manuscript Fig 3B set",
  "(using coenosarc coverage rather than composite 'healed'; 2026-05-23",
  "selection).",
  "",
  "**Repo figure: Wound-healing outcomes in *Porites* spp. over time.**",
  "",
  "Each panel shows the percentage of corals exhibiting a given wound",
  "feature (mean across corals) by treatment (Airbrush, Scrape,",
  "Airbrush-Scrape) across post-wounding days. Shaded bands are ±1 SE",
  "across corals; single-observation day/treatment cells are shown without",
  "a band. Statistics: separation-robust binomial GLMMs with a",
  "treatment x natural-spline(day) interaction and a coral random intercept",
  "(see Methods / PUBLICATION_STATISTICS_TABLE.csv)."
), file.path(OUT_DIR_MAIN, "figure2_porites_all_outcomes_legend.md"))

cat("[OK] Saved figure2_porites_all_outcomes.{pdf,png} + legend.md\n\n")

cat(strrep("=", 80), "\n MANUSCRIPT FIGURES COMPLETE\n",
    strrep("=", 80), "\n\n", sep = "")
