# =============================================================================
# SHARED SETUP — Experiment 2 (sourced by exp2_01..exp2_05)
#
# Single source of truth for:
#   - outcome recoding
#   - the separation-robust GLMM fitting ladder (so the model used for
#     INFERENCE in exp2_01 is byte-identical to the model DIAGNOSED in exp2_04)
#   - the reproducibility seed for stochastic diagnostics
#   - the publication figure theme + palette (figure scripts)
#
# Rationale: the code review found exp2_01 and exp2_04 fitted *different*
# models (exp2_04 always fitted the naive standard GLMM, even for outcomes
# where exp2_01 escalated to a penalized fit), and that recode_outcome /
# theme / palette were copy-pasted and already drifting across scripts.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(glmmTMB)
  library(splines)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# Reproducibility: fixed seed for DHARMa Monte-Carlo residual simulation.
# Set immediately before every simulateResiduals() call.
# -----------------------------------------------------------------------------
DHARMA_SEED <- 20260515L

# -----------------------------------------------------------------------------
# Outcome recoding: yes/no -> 1/0, everything else -> NA (never 0.5).
# -----------------------------------------------------------------------------
recode_outcome <- function(x) {
  x_lower <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    x_lower %in% c("yes", "y", "true", "t", "1") ~ 1,
    x_lower %in% c("no", "n", "false", "f", "0") ~ 0,
    TRUE ~ NA_real_
  )
}

# -----------------------------------------------------------------------------
# Separation-robust GLMM ladder
#
# Weakly-informative prior for complete separation (Gelman et al. 2008,
# Stat Med 27:2865-2873): Normal(0,3) on log-odds slopes, wider Normal(0,10)
# on the intercept. This is a penalized-likelihood / MAP estimate fitted with
# a PROPER binomial likelihood on the real 0/1 response — it is NOT the former
# (invalid) trick of squeezing the response toward 0.5 and refitting binomial.
# -----------------------------------------------------------------------------
SEPARATION_PRIORS <- data.frame(
  prior = c("normal(0, 3)", "normal(0, 10)"),
  class = c("fixef", "fixef"),
  coef  = c("", "(Intercept)")
)
METHOD_STANDARD  <- "Standard GLMM"
METHOD_PENALIZED <- "Penalized GLMM (weakly-informative Normal prior)"
METHOD_OLRE      <- "OLRE"

# glmmTMB-native health check. The previous code tested
# summary(m)$optinfo$conv$lme4$code — an lme4 structure glmmTMB never
# populates, so the convergence guard was silently always TRUE. Here:
# finite & non-extreme SEs, optimizer convergence == 0, positive-definite
# Hessian (sdr$pdHess).
model_is_healthy <- function(m) {
  if (is.null(m)) return(FALSE)
  ct <- tryCatch(summary(m)$coefficients$cond, error = function(e) NULL)
  if (is.null(ct)) return(FALSE)
  se <- ct[, "Std. Error"]
  all(is.finite(se)) && all(se <= 100) &&
    isTRUE(m$fit$convergence == 0) && isTRUE(m$sdr$pdHess)
}

# Per-outcome modelling frame from a species data frame.
prepare_outcome_data <- function(data, outcome_col) {
  data %>%
    dplyr::mutate(
      value  = recode_outcome(.data[[outcome_col]]),
      day_c  = day - mean(day, na.rm = TRUE),
      obs_id = factor(dplyr::row_number())
    ) %>%
    dplyr::filter(!is.na(value))
}

# Separation-robust fitting ladder. Returns list(method, model, success).
# Approach 1: standard binomial GLMM.
# Approach 2: same model, penalized with the weakly-informative prior.
# Approach 3: observation-level random effect (last resort).
fit_best_glmm <- function(dat_outcome, verbose = TRUE) {
  say <- function(...) if (verbose) cat(...)
  f <- value ~ treatment_label * ns(day_c, 3) + (1 | id)

  std <- tryCatch({
    m <- glmmTMB(f, data = dat_outcome, family = binomial())
    if (model_is_healthy(m)) {
      say("    [OK] Standard GLMM converged\n")
      list(method = METHOD_STANDARD, model = m, success = TRUE)
    } else {
      say("    [x] Standard GLMM unstable (separation / non-convergence)\n")
      list(method = METHOD_STANDARD, model = NULL, success = FALSE)
    }
  }, error = function(e) {
    say("    [x] Standard GLMM failed:", conditionMessage(e), "\n")
    list(method = METHOD_STANDARD, model = NULL, success = FALSE)
  })
  if (isTRUE(std$success)) return(std)

  say("  Trying penalized GLMM (weakly-informative Normal prior)...\n")
  pen <- tryCatch({
    m <- glmmTMB(f, data = dat_outcome, family = binomial(),
                 priors = SEPARATION_PRIORS)
    if (model_is_healthy(m)) {
      say("    [OK] Penalized GLMM converged\n")
      list(method = METHOD_PENALIZED, model = m, success = TRUE)
    } else {
      say("    [x] Penalized GLMM unstable\n")
      list(method = METHOD_PENALIZED, model = NULL, success = FALSE)
    }
  }, error = function(e) {
    say("    [x] Penalized GLMM failed:", conditionMessage(e), "\n")
    list(method = METHOD_PENALIZED, model = NULL, success = FALSE)
  })
  if (isTRUE(pen$success)) return(pen)

  say("  Trying observation-level random effects...\n")
  tryCatch({
    m <- glmmTMB(
      value ~ treatment_label * ns(day_c, 3) + (1 | id) + (1 | obs_id),
      data = dat_outcome, family = binomial()
    )
    if (model_is_healthy(m)) {
      say("    [OK] OLRE model converged\n")
      list(method = METHOD_OLRE, model = m, success = TRUE)
    } else {
      say("    [x] OLRE model unstable\n")
      list(method = METHOD_OLRE, model = NULL, success = FALSE)
    }
  }, error = function(e) {
    say("    [x] OLRE model failed:", conditionMessage(e), "\n")
    list(method = METHOD_OLRE, model = NULL, success = FALSE)
  })
}

# -----------------------------------------------------------------------------
# JOINT treatment x time interaction test (PRIMARY interaction result).
#
# The full interaction block is treatment_label : ns(day_c, 3) = the set of all
# treatment-contrast x spline-basis coefficients. With 3 treatment levels and a
# cubic natural spline this block is 2 contrasts x 3 bases = 6 coefficients
# (6 df) -- NOT a single coefficient. The previous "min of the individual Wald
# p-values" summary was not a valid single test and is parameterization-
# dependent (centred vs raw day changes which basis-coefficient is smallest).
#
# Method per outcome:
#   - Standard binomial GLMM (no separation): a clean likelihood-ratio test of
#     full vs reduced (treatment + ns(day_c,3), no interaction) via anova().
#   - Penalized GLMM (weakly-informative prior) or OLRE: the penalized objective
#     is not a plain likelihood, so an LRT is not clean. Use a joint Wald test
#     on the interaction-block coefficients:  W = b' [V]^-1 b  ~  chi^2(df),
#     where b and V are the fixed-effect estimates and covariance for the block.
#
# Returns list(joint_p, joint_method, joint_stat, joint_df).
# -----------------------------------------------------------------------------
joint_interaction_test <- function(dat_outcome, best_model) {
  if (is.null(best_model) || !isTRUE(best_model$success)) {
    return(list(joint_p = NA_real_, joint_method = "model failed",
                joint_stat = NA_real_, joint_df = NA_integer_))
  }
  m_full <- best_model$model
  method <- best_model$method

  beta <- glmmTMB::fixef(m_full)$cond
  inter_idx <- grep("treatment_label.*ns\\(day_c", names(beta))
  k <- length(inter_idx)
  if (k == 0) {
    return(list(joint_p = NA_real_, joint_method = "no interaction terms",
                joint_stat = NA_real_, joint_df = 0L))
  }

  is_standard <- identical(method, METHOD_STANDARD)

  # --- Clean LRT for the standard (un-penalized) GLMM only -------------------
  if (is_standard) {
    f_red <- value ~ treatment_label + ns(day_c, 3) + (1 | id)
    lrt <- tryCatch({
      m_red <- glmmTMB(f_red, data = dat_outcome, family = binomial())
      an <- stats::anova(m_red, m_full)
      list(p   = an[["Pr(>Chisq)"]][2],
           chi = an[["Chisq"]][2],
           df  = an[["Chi Df"]][2])
    }, error = function(e) NULL)
    if (!is.null(lrt) && is.finite(lrt$p)) {
      return(list(joint_p = lrt$p,
                  joint_method = sprintf("LRT, %d df", lrt$df),
                  joint_stat = lrt$chi, joint_df = lrt$df))
    }
    # else fall through to joint Wald
  }

  # --- Joint Wald on the interaction block (penalized / OLRE / LRT failure) ---
  V    <- as.matrix(stats::vcov(m_full)$cond)
  b_I  <- beta[inter_idx]
  V_I  <- V[inter_idx, inter_idx, drop = FALSE]
  W    <- tryCatch(as.numeric(t(b_I) %*% solve(V_I) %*% b_I),
                   error = function(e) NA_real_)
  p    <- if (is.finite(W)) stats::pchisq(W, df = k, lower.tail = FALSE) else NA_real_
  lbl  <- if (is_standard) sprintf("Joint Wald, %d df (LRT unavailable)", k)
          else if (identical(method, METHOD_PENALIZED))
            sprintf("Joint Wald, %d df (penalized)", k)
          else sprintf("Joint Wald, %d df", k)
  list(joint_p = p, joint_method = lbl, joint_stat = W, joint_df = k)
}

# -----------------------------------------------------------------------------
# Publication figure theme + palette (used by exp2_02 / exp2_05).
# Okabe-Ito qualitative palette; text sizes flow from base_size via rel().
# -----------------------------------------------------------------------------
TREATMENT_COLORS <- c(
  "Airbrush"        = "#0072B2",
  "Dremel"          = "#D55E00",
  "Airbrush-Dremel" = "#009E73"
)

theme_pub <- function(base_size = 10) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(size = ggplot2::rel(0.9),
                                                face = "bold"),
      axis.text        = ggplot2::element_text(size = ggplot2::rel(0.85)),
      legend.position  = "bottom",
      legend.title     = ggplot2::element_blank(),
      plot.tag         = ggplot2::element_text(size = ggplot2::rel(1.1),
                                                face = "bold"),
      plot.margin      = ggplot2::margin(8, 8, 8, 8, "mm")
    )
}
