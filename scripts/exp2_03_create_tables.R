# =============================================================================
# CREATE PUBLICATION-READY STATISTICAL TABLES (Experiment 2)
#
# Reads the tidy fixed-effect coefficient CSVs written by exp2_01
# (output/exp2_models/*_tidy.csv) -- the authoritative machine-readable model
# output. This REPLACES the former approach of regex-scraping printed
# summary() text, which mis-parsed glmmTMB coefficient layouts (wrong line
# offset), did not handle the OLRE method, and silently reported parse
# failures as "No significant interactions". The hardcoded KEY FINDINGS block
# (which could contradict the regenerated CSV) is also removed; findings are
# now derived from the table.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# METHOD_* label constants, shared with exp2_01/02/04/05
source(here("scripts", "00_setup.R"))

MODELS_DIR <- here("output", "exp2_models")
TABLES_DIR <- here("output", "exp2_tables")
dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)

tidy_files <- list.files(MODELS_DIR, pattern = "_tidy\\.csv$",
                         full.names = TRUE)

if (length(tidy_files) == 0) {
  stop("No *_tidy.csv files in ", MODELS_DIR,
       " -- run scripts/exp2_01_fit_glmm_models.R first.")
}
cat("Found", length(tidy_files), "tidy coefficient files\n\n")

INTERACTION_RX <- "treatment_label.*ns\\(day_c"

summarise_outcome <- function(file_path) {
  d <- readr::read_csv(file_path, show_col_types = FALSE)

  needed <- c("Species", "Outcome", "Method", "term", "estimate",
              "p.value", "conf.low", "conf.high")
  if (!all(needed %in% names(d))) {
    return(tibble(
      Species = NA_character_, Outcome = basename(file_path),
      Method = NA_character_, N_Interactions = 0L, N_Significant = 0L,
      Min_P = NA_real_, Status = "PARSE FAILURE: missing columns",
      Significant_Interactions = "PARSE FAILURE"
    ))
  }

  inter <- d %>% filter(grepl(INTERACTION_RX, term))

  if (nrow(inter) == 0) {
    return(tibble(
      Species = d$Species[1], Outcome = d$Outcome[1], Method = d$Method[1],
      N_Interactions = 0L, N_Significant = 0L, Min_P = NA_real_,
      Status = "No interaction terms", Significant_Interactions = "None"
    ))
  }

  p <- inter$p.value
  if (all(is.na(p))) {
    # Explicit failure -- never silently fall through to a null result.
    return(tibble(
      Species = d$Species[1], Outcome = d$Outcome[1], Method = d$Method[1],
      N_Interactions = nrow(inter), N_Significant = 0L, Min_P = NA_real_,
      Status = "PARSE FAILURE: all interaction p-values NA",
      Significant_Interactions = "PARSE FAILURE"
    ))
  }

  sig <- inter$term[!is.na(p) & p < 0.05]
  tibble(
    Species = d$Species[1],
    Outcome = d$Outcome[1],
    Method  = d$Method[1],
    N_Interactions = nrow(inter),
    N_Significant  = sum(p < 0.05, na.rm = TRUE),
    Min_P = min(p, na.rm = TRUE),
    Status = "OK",
    Significant_Interactions = if (length(sig) > 0)
      paste(sig, collapse = "; ") else "None"
  )
}

all_stats <- map_dfr(tidy_files, summarise_outcome)

# ---- PRIMARY interaction result: the JOINT treatment x time test -----------
# exp2_01 writes the joint test (LRT for standard GLMMs, joint Wald on the full
# interaction block for penalized/OLRE fits) to ALL_OUTCOMES_SUMMARY.csv. Join
# it in by Species + Outcome. The per-coefficient min Wald p (Min_P) is kept
# only as a SECONDARY, parameterization-dependent descriptor.
summary_csv <- file.path(MODELS_DIR, "ALL_OUTCOMES_SUMMARY.csv")
if (!file.exists(summary_csv)) {
  stop("Missing ", summary_csv,
       " -- run scripts/exp2_01_fit_glmm_models.R first (it writes the joint test).")
}
joint_tbl <- readr::read_csv(summary_csv, show_col_types = FALSE) %>%
  select(Species, Outcome, Joint_Interaction_P, Joint_Method, Joint_DF)

summary_table <- all_stats %>%
  left_join(joint_tbl, by = c("Species", "Outcome")) %>%
  mutate(
    Statistical_Result = case_when(
      grepl("PARSE FAILURE", Status) ~ Status,
      is.finite(Joint_Interaction_P) & Joint_Interaction_P < 0.05 ~
        paste0("Significant treatment x time interaction (",
               Joint_Method, ", p = ", signif(Joint_Interaction_P, 3), ")"),
      is.finite(Joint_Interaction_P) & Joint_Interaction_P < 0.10 ~
        paste0("Marginal treatment x time interaction (",
               Joint_Method, ", p = ", signif(Joint_Interaction_P, 3), ")"),
      is.finite(Joint_Interaction_P) ~
        paste0("No significant treatment x time interaction (",
               Joint_Method, ", p = ", signif(Joint_Interaction_P, 3), ")"),
      TRUE ~ "Joint test unavailable"
    ),
    Method_Note = case_when(
      Method == METHOD_PENALIZED ~
        "Penalized GLMM (weakly-informative Normal prior; Gelman 2008; separation)",
      Method == METHOD_STANDARD ~ "Standard binomial GLMM",
      Method == METHOD_OLRE     ~ "GLMM + observation-level random effect",
      TRUE ~ Method
    )
  ) %>%
  rename(Secondary_Min_Wald_P = Min_P) %>%
  select(Species, Outcome, Method_Note,
         Joint_Interaction_P, Joint_Method, Statistical_Result,
         Secondary_Min_Wald_P, Status,
         N_Significant, Significant_Interactions) %>%
  arrange(Species, Joint_Interaction_P)

# Canonical filename is UPPERCASE to match README / MANIFEST / docs and to
# avoid a case-divergent duplicate on case-sensitive filesystems (Linux/CI).
out_csv <- file.path(TABLES_DIR, "PUBLICATION_STATISTICS_TABLE.csv")
write_csv(summary_table, out_csv)

# Drop any stale lowercase duplicate left by older runs.
old_lc <- file.path(TABLES_DIR, "publication_statistics_table.csv")
if (file.exists(old_lc) &&
    !identical(normalizePath(old_lc, mustWork = FALSE),
               normalizePath(out_csv, mustWork = FALSE))) {
  file.remove(old_lc)
}

# ---- Multiple-testing correction (regenerated, NOT a static orphan) --------
# Family = the SIX a-priori Porites outcomes, each contributing ONE joint
# treatment x time p-value (not the ~36 individual basis coefficients the old
# min-p procedure implicitly screened). Sensitivity only: the a-priori directed
# hypotheses use the unadjusted JOINT p as the primary inference (author
# decision, see RESULTS_MASTER). BH-FDR & Bonferroni are computed WITHIN each
# species' outcome family from the regenerated joint p, so this file can never
# drift from PUBLICATION_STATISTICS_TABLE.csv.
mt <- summary_table %>%
  filter(is.finite(Joint_Interaction_P)) %>%
  group_by(Species) %>%
  mutate(
    p_raw        = Joint_Interaction_P,
    p_BH_FDR     = p.adjust(Joint_Interaction_P, method = "BH"),
    p_Bonferroni = p.adjust(Joint_Interaction_P, method = "bonferroni"),
    sig_raw        = p_raw < 0.05,
    sig_BH         = p_BH_FDR < 0.05,
    sig_Bonferroni = p_Bonferroni < 0.05,
    family_size  = dplyr::n()
  ) %>%
  ungroup() %>%
  select(Species, Outcome, p_raw, p_BH_FDR, p_Bonferroni,
         sig_raw, sig_BH, sig_Bonferroni, family_size) %>%
  arrange(Species, p_raw)
write_csv(mt, file.path(TABLES_DIR, "MULTIPLE_TESTING_CORRECTION.csv"))

# ---- A-priori multiplicity statement (methods note, regenerated) -----------
# Written to output/text/ so the manuscript Methods can quote it verbatim and it
# can never drift from the regenerated correction table.
TEXT_DIR <- here("output", "text")
dir.create(TEXT_DIR, showWarnings = FALSE, recursive = TRUE)
fam_n <- mt %>% count(Species, name = "n_outcomes")
apriori_lines <- c(
  "Experiment 2 -- multiplicity framing (a-priori directed hypotheses)",
  paste0("Generated by scripts/exp2_03_create_tables.R on ", Sys.Date()),
  "",
  "The six Porites outcomes (RFP, polyps in centre of wound, coenosarc",
  "coverage, yellow aggregations, pink aggregations, algal plug) were specified",
  "a priori as directed hypotheses about how skeletal vs tissue-only wounding",
  "alters the regeneration trajectory. The PRIMARY treatment x time result for",
  "each outcome is a single JOINT test of the full treatment_label : ns(day_c,3)",
  "interaction block (likelihood-ratio test for the standard binomial GLMM;",
  "joint Wald test on the interaction-block coefficients for the penalized /",
  "separation fits). Under this a-priori directed-hypothesis design the",
  "unadjusted joint p is primary.",
  "",
  paste0("The true testing family is therefore SIX joint tests (one per outcome),",
         " not the ~36 individual spline-basis coefficients that the former"),
  "min-of-Wald-p summary implicitly screened. The multiple-comparison correction",
  "OF RECORD is Benjamini-Hochberg (controls the false-discovery rate across the",
  "six-outcome family, the appropriate target for an a-priori directed-hypothesis",
  "design); Bonferroni (family-wise error) is reported only as a more conservative",
  "sensitivity. Both are tabulated over the six joint p-values in",
  "MULTIPLE_TESTING_CORRECTION.csv (see sig_BH and sig_Bonferroni).",
  ""
)
for (i in seq_len(nrow(fam_n))) {
  apriori_lines <- c(apriori_lines,
    sprintf("Family size for %s: %d outcomes.", fam_n$Species[i], fam_n$n_outcomes[i]))
}
writeLines(apriori_lines, file.path(TEXT_DIR, "exp2_multiplicity_apriori_note.txt"))
cat("[OK] A-priori multiplicity note saved to: output/text/exp2_multiplicity_apriori_note.txt\n\n")

# ---- console report --------------------------------------------------------
cat(strrep("=", 80), "\n PUBLICATION-READY STATISTICAL SUMMARY\n",
    strrep("=", 80), "\n\n", sep = "")
print(summary_table, n = Inf, width = Inf)

fails <- summary_table %>% filter(grepl("PARSE FAILURE", Status))
if (nrow(fails) > 0) {
  cat("\n*** ", nrow(fails),
      " outcome(s) FAILED to parse -- investigate before using this table:\n",
      sep = "")
  print(fails %>% select(Species, Outcome, Status))
}

# ---- KEY FINDINGS: derived from the regenerated table, never hardcoded -----
# Based on the PRIMARY joint treatment x time interaction p (not min-p).
cat("\nKEY FINDINGS (joint treatment x time interaction, derived from models):\n")
sig_tbl <- summary_table %>%
  filter(is.finite(Joint_Interaction_P) & Joint_Interaction_P < 0.05) %>%
  arrange(Joint_Interaction_P)
if (nrow(sig_tbl) == 0) {
  cat("  - No outcomes show a significant joint treatment x time interaction.\n")
} else {
  for (i in seq_len(nrow(sig_tbl))) {
    cat(sprintf(
      "  - %s %s: significant joint treatment x time interaction (%s, p = %.4g)\n",
      sig_tbl$Species[i], sig_tbl$Outcome[i],
      sig_tbl$Joint_Method[i], sig_tbl$Joint_Interaction_P[i]))
  }
}
ns_tbl <- summary_table %>%
  filter(is.finite(Joint_Interaction_P) & Joint_Interaction_P >= 0.05) %>%
  arrange(Joint_Interaction_P)
for (i in seq_len(nrow(ns_tbl))) {
  cat(sprintf("  - %s %s: no significant joint interaction (%s, p = %.4g)\n",
              ns_tbl$Species[i], ns_tbl$Outcome[i],
              ns_tbl$Joint_Method[i], ns_tbl$Joint_Interaction_P[i]))
}

cat("\n[OK] Publication table saved to:",
    file.path("output", "exp2_tables", "PUBLICATION_STATISTICS_TABLE.csv"),
    "\n\n")
