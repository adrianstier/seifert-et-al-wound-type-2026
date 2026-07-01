# =============================================================================
# exp1_dharma_diagnostics.R — residual, temporal, and separation diagnostics
# for the reported Experiment 1 additive GLMMs (healed, regenerated, algae).
# -----------------------------------------------------------------------------
# Adds what the main pipeline lacked: DHARMa simulated-residual tests for the
# ADDITIVE (reported) models — not the interaction model — for all three binary
# outcomes; a temporal-autocorrelation test on day-aggregated residuals (the
# repeated-measures check reviewers asked for); and a formal complete-separation
# test (detectseparation) that justifies the Firth refits.
# Outputs:
#   output/tables/exp1_dharma_summary.csv
#   output/tables/exp1_separation_detection.csv
#   output/figures/exp1_dharma_residuals.{png,pdf}
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse); library(here); library(lme4); library(DHARMa); library(detectseparation)
})
SEED <- 20260515  # matches the Exp 2 DHARMA_SEED for reproducibility

d <- read_csv(here("data","airbrush_dremel.csv"), show_col_types = FALSE) %>%
  mutate(species   = factor(str_trim(tolower(species)), levels = c("acropora","pocillopora","porites")),
         treatment = factor(str_trim(tolower(treatment)), levels = c("airbrush","dremel")),
         parent_id = factor(sub("[[:alpha:]]+$","", as.character(coral_id))),
         coral_id  = factor(coral_id),
         healed01  = if_else(healed == "yes", 1L, 0L),
         regen01   = if_else(regenerated == "yes", 1L, 0L),
         debris01  = if_else(debris == "yes", 1L, 0L))

outcomes <- list(
  Healing            = "healed01",
  Regeneration       = "regen01",
  `Algal colonization` = "debris01")

fit_add <- function(col) {
  dd <- d %>% filter(!is.na(.data[[col]]))
  m  <- glmer(as.formula(paste0(col, " ~ treatment + species + (1|parent_id/coral_id)")),
              data = dd, family = binomial,
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  list(model = m, data = dd)
}
fits <- map(outcomes, fit_add)

# ---- DHARMa residual + temporal-autocorrelation tests --------------------
dharma_rows <- imap_dfr(fits, function(f, nm) {
  set.seed(SEED)
  sr <- simulateResiduals(f$model, n = 1000)
  uni <- testUniformity(sr, plot = FALSE)
  dsp <- testDispersion(sr, plot = FALSE)
  out <- testOutliers(sr, plot = FALSE)
  # temporal autocorrelation on residuals aggregated within day (one value/day)
  srd <- recalculateResiduals(sr, group = f$data$day)
  tac <- tryCatch(testTemporalAutocorrelation(srd, time = sort(unique(f$data$day)), plot = FALSE),
                  error = function(e) list(p.value = NA_real_, statistic = NA_real_))
  tibble(outcome = nm, n = nrow(f$data),
         uniformity_KS_p = unname(uni$p.value),
         dispersion_p    = unname(dsp$p.value),
         outlier_p       = unname(out$p.value),
         temporal_DW     = unname(tac$statistic[1]),
         temporal_p      = unname(tac$p.value))
})
write_csv(dharma_rows, here("output","tables","exp1_dharma_summary.csv"))

# ---- formal complete-separation test (fixed-effects logistic) ------------
sep_rows <- imap_dfr(outcomes, function(col, nm) {
  dd <- d %>% filter(!is.na(.data[[col]]))
  ds <- glm(as.formula(paste0(col, " ~ treatment + species")), data = dd,
            family = binomial, method = "detect_separation")
  b   <- coef(ds)              # 0 = finite MLE; Inf/-Inf = separated coefficient
  sep <- names(b)[is.infinite(b)]
  tibble(outcome = nm,
         complete_separation = any(is.infinite(b)),
         separated_terms = if (length(sep)) paste(sep, collapse = "; ") else "none")
})
write_csv(sep_rows, here("output","tables","exp1_separation_detection.csv"))

# ---- composite residual figure: clean ggplot panels (no DHARMa base-graphics
#      title/label clutter) — QQ uniformity + residual-vs-predicted per outcome ---
suppressPackageStartupMessages(library(patchwork))
theme_pub <- function(base_size = 9) ggplot2::theme_bw(base_size = base_size) +
  ggplot2::theme(panel.grid.minor = element_blank(), strip.background = element_blank(),
                 plot.title = element_text(size = rel(0.95), face = "bold"),
                 plot.margin = margin(4, 6, 4, 6))
panel_pair <- function(f, nm) {
  set.seed(SEED); sr <- simulateResiduals(f$model, n = 1000)
  r <- sr$scaledResiduals; n <- length(r)
  qq <- tibble(expected = qunif(ppoints(n)), observed = sort(r))
  rp <- tibble(pred = sr$fittedPredictedResponse, res = r)
  p1 <- ggplot(qq, aes(expected, observed)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = 2) +
    geom_point(size = 0.5, alpha = 0.5, colour = "#0072B2") +
    labs(x = "Expected (uniform)", y = "Observed residual", title = paste0(nm, ": QQ uniformity")) +
    theme_pub()
  p2 <- ggplot(rp, aes(pred, res)) +
    geom_hline(yintercept = c(0.25, 0.5, 0.75), colour = "grey85", linetype = 2) +
    geom_point(size = 0.5, alpha = 0.4, colour = "#0072B2") +
    geom_smooth(method = "loess", se = FALSE, colour = "#D55E00", linewidth = 0.6, na.rm = TRUE) +
    labs(x = "Predicted probability", y = "Scaled residual", title = paste0(nm, ": residual vs. predicted")) +
    theme_pub()
  p1 | p2
}
dharma_fig <- wrap_plots(imap(fits, panel_pair), ncol = 1) +
  plot_annotation(title = "Experiment 1 additive GLMMs: DHARMa simulated residuals (seeded, n = 1000)",
                  theme = theme(plot.title = element_text(size = 11, face = "bold")))
ggsave(here("output","figures","exp1_dharma_residuals.pdf"), dharma_fig,
       width = 180, height = 215, units = "mm", dpi = 300, device = cairo_pdf, bg = "white")
ggsave(here("output","figures","exp1_dharma_residuals.png"), dharma_fig,
       width = 180, height = 215, units = "mm", dpi = 300, bg = "white")

cat("\n=== DHARMa summary (Exp 1 additive models) ===\n"); print(as.data.frame(dharma_rows), digits = 3)
cat("\n=== Separation detection ===\n"); print(as.data.frame(sep_rows))
cat("\nWrote: exp1_dharma_summary.csv, exp1_separation_detection.csv, exp1_dharma_residuals.{png,pdf}\n")
