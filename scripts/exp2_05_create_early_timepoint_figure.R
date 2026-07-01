################################################################################
# EXPERIMENT 2: EARLY TIMEPOINT FIGURE (DAYS 0-5)
#
# Focused supplement figure: early wound dynamics in Porites for three
# outcomes (Algal plug, Coenosarc coverage, Yellow aggregations), days 0-5.
#
# Faceted (shared axes labelled once, equal panel widths, one collected
# legend) using the shared 00_setup.R theme/palette/recoding — consistent
# with exp2_02 and the project figure standard.
#
# Output (mm, 300 dpi, cairo_pdf + png + companion legend md):
#   - output/exp2_figures_supplement/figureS2_porites_early_dynamics.{pdf,png}
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(here)
  library(scales)
})

# theme_pub(), TREATMENT_COLORS, recode_outcome()
source(here("scripts", "00_setup.R"))

OUT_DIR <- here("output", "exp2_figures_supplement")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

LINE_SIZE   <- 0.9
POINT_SIZE  <- 2.0
ERROR_ALPHA <- 0.15
BASE_SIZE   <- 10

cat("\n", strrep("=", 80),
    "\n CREATING EARLY TIMEPOINT FIGURE (DAYS 0-5)\n",
    strrep("=", 80), "\n\n", sep = "")

dat_porites <- read_excel(
    here("data", "Porites microscope characterization - complete.xlsx"),
    sheet = "data") %>%
  clean_names() %>%
  mutate(
    species = "Porites spp.",
    treatment = tolower(trimws(as.character(treatment))),
    treatment_label = factor(case_when(
      treatment == "airbrush" ~ "Airbrush",
      treatment %in% c("drem", "dremel") ~ "Scrape",
      treatment == "air_drem" ~ "Airbrush + Scrape"
    ), levels = c("Airbrush", "Scrape", "Airbrush + Scrape")),
    day = as.numeric(day),
    id = factor(id)
  ) %>%
  filter(!is.na(treatment_label), !is.na(day), day >= 0, day <= 5)

cat("[OK] Loaded", nrow(dat_porites), "Porites observations (days 0-5)\n")
cat("  Timepoints:", paste(sort(unique(dat_porites$day)), collapse = ", "),
    "\n\n")

OUTCOMES <- list(
  list(col = "algal_plug",          label = "Algal plug"),
  list(col = "coenosarc_coverage",  label = "Coenosarc coverage"),
  list(col = "yellow_aggregations", label = "Yellow aggregations")
)
labels <- vapply(OUTCOMES, `[[`, character(1), "label")

# n = 1 cells -> se = NA (no ribbon) rather than a misleading zero-width band.
early_df <- purrr::imap_dfr(OUTCOMES, function(info, nm) {
  dat_porites %>%
    mutate(value = recode_outcome(.data[[info$col]])) %>%
    filter(!is.na(value)) %>%
    group_by(day, treatment_label) %>%
    summarise(n = dplyr::n(),
              mean = mean(value),
              se = if (dplyr::n() < 2) NA_real_ else sd(value) / sqrt(dplyr::n()),
              .groups = "drop") %>%
    mutate(outcome = info$label)
}) %>%
  mutate(outcome = factor(outcome, levels = labels),
         lower = pmax(mean - se, 0),
         upper = pmin(mean + se, 1))

fig_early <- ggplot(early_df, aes(day, mean * 100,
                                  color = treatment_label,
                                  fill = treatment_label)) +
  geom_ribbon(aes(ymin = lower * 100, ymax = upper * 100),
              alpha = ERROR_ALPHA, color = NA, na.rm = TRUE) +
  geom_line(linewidth = LINE_SIZE, alpha = 0.9, na.rm = TRUE) +
  geom_point(size = POINT_SIZE, shape = 21, fill = "white",
             stroke = 0.9, na.rm = TRUE) +
  scale_color_manual(values = setNames(unname(TREATMENT_COLORS), c("Airbrush", "Scrape", "Airbrush + Scrape"))) +
  scale_fill_manual(values = setNames(unname(TREATMENT_COLORS), c("Airbrush", "Scrape", "Airbrush + Scrape"))) +
  scale_x_continuous(breaks = 0:5, limits = c(0, 5),
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
                     expand = expansion(mult = c(0, 0.02))) +
  facet_wrap(~ outcome, ncol = 3, labeller = label_wrap_gen(18)) +
  labs(x = "Day", y = "Corals with feature (%)") +
  theme_pub(base_size = BASE_SIZE)

ggsave(file.path(OUT_DIR, "figureS2_porites_early_dynamics.pdf"),
       plot = fig_early, width = 180, height = 80, units = "mm",
       dpi = 300, device = cairo_pdf, bg = "white")
ggsave(file.path(OUT_DIR, "figureS2_porites_early_dynamics.png"),
       plot = fig_early, width = 180, height = 80, units = "mm",
       dpi = 300, bg = "white")

writeLines(c(
  "**Manuscript correspondence:** early-timepoint data supplement to",
  "**manuscript Figure 3** (and the Fig 4 early-cellular timeline). It is",
  "**NOT confirmed to be manuscript Figure S2** — do not infer the manuscript",
  "supplement number from the filename.",
  "",
  "**Repo figure: Early wound dynamics in *Porites* spp. (days 0-5).**",
  "",
  "Percentage of corals exhibiting each early-response feature (mean across",
  "corals) by treatment over the first five post-wounding days. Shaded bands",
  "are ±1 SE across corals; single-observation cells are shown without a",
  "band. Complements the full *Porites* outcome time series."
), file.path(OUT_DIR, "figureS2_porites_early_dynamics_legend.md"))

cat("[OK] Saved figureS2_porites_early_dynamics.{pdf,png} + legend.md\n\n")
cat(strrep("=", 80), "\n EARLY TIMEPOINT FIGURE COMPLETE\n",
    strrep("=", 80), "\n\n", sep = "")
