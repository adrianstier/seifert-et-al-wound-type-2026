# =============================================================================
# verify_pipeline.R — reproducibility / integrity self-check
#
# Asserts the invariants that must hold for the analysis to be "working":
# data integrity, Exp 1 / Exp 2 numeric invariants, model-selection coherence,
# doc/crosswalk-traces-to-CSV, and repo hygiene. Exits non-zero on any failure.
#
# Run:  Rscript scripts/verify_pipeline.R          (or: make verify)
# Assumes the pipeline outputs exist (run `make all` first for a clean check).
# =============================================================================

suppressPackageStartupMessages({ library(readr); library(here) })

.PASS <- 0L; .FAIL <- 0L
ok <- function(label, cond, detail = "") {
  cond <- isTRUE(tryCatch(cond, error = function(e) { detail <<- conditionMessage(e); FALSE }))
  if (cond) { .PASS <<- .PASS + 1L; cat(sprintf("  PASS  %s\n", label)) }
  else      { .FAIL <<- .FAIL + 1L; cat(sprintf("  FAIL  %s%s\n", label,
                                       if (nzchar(detail)) paste0("  [", detail, "]") else "")) }
}
near <- function(x, target, rtol = 0.02) is.finite(x) && abs(x - target) <= rtol * abs(target)
section <- function(s) cat(sprintf("\n== %s ==\n", s))
rd <- function(p) suppressWarnings(read_csv(here(p), show_col_types = FALSE))

# ---- A. scripts parse -------------------------------------------------------
section("A. scripts parse")
for (f in list.files(here("scripts"), pattern = "\\.R$", full.names = TRUE)) {
  ok(paste("parses:", basename(f)),
     { invisible(parse(f)); TRUE })
}

# ---- B. data integrity ------------------------------------------------------
section("B. data integrity (2026-06-15 NA-backtrack re-score)")
d <- rd("data/airbrush_dremel.csv")
ok("airbrush_dremel.csv is 154 x 9", nrow(d) == 154 && ncol(d) == 9)
ok("NA healed=3 debris=4 regenerated=0 (t-1/t+1 neighbour-cell backtrack)",
   sum(is.na(d$healed)) == 3 && sum(is.na(d$debris)) == 4 && sum(is.na(d$regenerated)) == 0)
ok("regenerated 36 yes / 118 no",
   sum(d$regenerated == "yes", na.rm = TRUE) == 36 &&
   sum(d$regenerated == "no",  na.rm = TRUE) == 118)
cellv <- function(df, sp, tr, id, dy, col) {
  r <- df[tolower(trimws(df$species))==sp & tolower(df$treatment)==tr &
          tolower(as.character(df$coral_id))==id & as.integer(df$day)==dy, ]
  if (nrow(r)==1) as.character(r[[col]]) else NA_character_
}
ok("backtrack re-score acropora/dremel/10b/d18 regenerated == no", cellv(d,"acropora","dremel","10b",18,"regenerated") == "no")
ok("backtrack re-score pocillopora/dremel/6b/d23 healed == yes",   cellv(d,"pocillopora","dremel","6b",23,"healed") == "yes")
ok("backtrack NA-fill pocillopora/dremel/6b/d28 healed == yes",    cellv(d,"pocillopora","dremel","6b",28,"healed") == "yes")
old <- rd("data/archive_exp1_2022/airbrush_dremel_pre-backtrack-2026-06-18.csv")
ok("archived pre-backtrack data present & differs (6b/d23 healed == incomplete)", cellv(old,"pocillopora","dremel","6b",23,"healed") == "incomplete")
ok("separation resolved: Pocillopora x dremel healed=yes == 2",
   sum(tolower(trimws(d$species))=="pocillopora" & d$treatment=="dremel" & d$healed=="yes", na.rm=TRUE) == 2)

# ---- C. Exp 1 numeric invariants -------------------------------------------
section("C. Exp 1 invariants")
h <- rd("output/tables/firth_healed_fixed_effect_ORs.csv")
r <- rd("output/tables/firth_regen_fixed_effect_ORs.csv")
hor <- function(t, col="OR") h[[col]][grepl(t, h$Term)][1]
ror <- function(t, col="OR") r[[col]][grepl(t, r$Term)][1]
hb  <- rd("output/tables/binary_fixed_effect_ORs.csv")            # lme4 healed (now estimable)
hbor <- function(t, col="OR") hb[[col]][grepl(t, hb$Term)][1]
ok("Healed lme4 Scrape/Airbrush OR ~ 26.0 (finite; separation resolved)", near(hbor("Dremel|Scrape|treatment"), 25.96, 0.05))
ok("Healed Firth Scrape/Airbrush OR ~ 17.3 (sensitivity refit)",          near(hor("Dremel|Scrape|treatment"), 17.30, 0.03))
ok("Healed Firth Pocillopora vs Acropora ~ 0.236",                        near(hor("Pocillopora"), 0.2355, 0.05))
ok("Healed Firth Porites vs Acropora ~ 0.498",                            near(hor("Porites"), 0.498, 0.03))
ok("Regen Firth Scrape/Airbrush OR ~ 52.2",                               near(ror("Dremel|Scrape|treatment"), 52.21, 0.03))
ok("Regen Firth Porites vs Acropora ~ 1.07",                              near(ror("Porites"), 1.07, 0.05))
dr <- rd("output/tables/regen_glmm_drop1_additive.csv")
ok("Regen drop1 treatment chi2 ~ 31.79", near(dr$`Chi-square`[grepl("treatment", dr$Term)][1], 31.79, 0.03))
ok("Regen drop1 species chi2 ~ 14.13",   near(dr$`Chi-square`[grepl("species", dr$Term)][1],   14.13, 0.03))
hdr <- rd("output/tables/binary_glmm_drop1_additive.csv")
ok("Healed drop1 species now NS (chi2 ~ 4.80, p > 0.05)",
   near(hdr$`Chi-square`[grepl("species", hdr$Term)][1], 4.80, 0.05) &&
   suppressWarnings(as.numeric(hdr$`P-value`[grepl("species", hdr$Term)][1])) > 0.05)
lr <- rd("output/tables/regen_lrt_interaction_vs_additive.csv")
ok("Regen trt x species interaction NS -> additive (p > 0.05)",
   suppressWarnings(as.numeric(lr$P_value[1])) > 0.05)

# ---- D. Exp 2 invariants ----------------------------------------------------
section("D. Exp 2 invariants")
# Exp 2 is Porites-only -> 6 Porites outcomes (composite "Healed" dropped
# 2026-05-25; redundant with the three constituent outcomes plotted in Fig 3
# and was inflating the multiple-testing denominator).
p <- rd("output/exp2_tables/PUBLICATION_STATISTICS_TABLE.csv")
ok("PUBLICATION_STATISTICS_TABLE has 6 rows (Porites-only, m=6)",
   nrow(p) == 6)
por <- p[grepl("Porites", p$Species), ]
ok("all Exp 2 publication rows are Porites",
   all(grepl("Porites", p$Species)))
pink <- por[por$Outcome == "Pink", ]
ok("Porites Pink uses Standard GLMM (no separation)", grepl("Standard", pink$Method_Note))
ok("Porites non-Pink outcomes use Penalized prior",
   all(grepl("Penalized", por$Method_Note[por$Outcome != "Pink"])))
ok("Porites RFP joint interaction p < 1e-4 (~1.5e-6)", por$Joint_Interaction_P[por$Outcome == "RFP"] < 1e-4)
ok("Porites 4 of 6 outcomes significant (joint treatment x time p<0.05)",
   sum(por$Joint_Interaction_P < 0.05, na.rm = TRUE) == 4 && nrow(por) == 6)
mt <- rd("output/exp2_tables/MULTIPLE_TESTING_CORRECTION.csv")
mtp <- mt[grepl("Porites", mt$Species), ]
survBH <- mtp$Outcome[mtp$p_BH_FDR < 0.05]
ok("BH-FDR (correction of record): RFP/Pink/Polyps/Yellow survive; Coenosarc does not",
   all(c("RFP") %in% survBH) && any(grepl("Pink", survBH)) &&
   any(grepl("Polyps", survBH)) && any(grepl("Yellow", survBH)) &&
   !any(grepl("Coenosarc", survBH)))
survBonf <- mtp$Outcome[mtp$p_Bonferroni < 0.05]
ok("Bonferroni (conservative sensitivity): RFP & Pink survive",
   all(c("RFP") %in% survBonf) && any(grepl("Pink", survBonf)))

# ---- D2. Exp 1 design-respecting PRIMARY (fragment-level + paired, n=22) -----
section("D2. Exp 1 fragment-level primary (pseudoreplication-free)")
fp <- rd("output/tables/exp1_fragment_level_primary.csv")
ok("fragment-level primary carries regenerated + healed results",
   any(fp$outcome == "regenerated") && any(fp$outcome == "healed"))
mcn <- fp[grepl("McNemar", fp$model), ]
ok("paired McNemar significant for regenerated (p<0.05)",
   any(mcn$outcome == "regenerated" & as.numeric(mcn$p_value) < 0.05))
ok("paired McNemar significant for healed (p<0.05)",
   any(mcn$outcome == "healed" & as.numeric(mcn$p_value) < 0.05))
frth <- fp[grepl("Firth", fp$model), ]
ok("fragment-level Firth regen OR > 1 (scrape favoured over airbrush)",
   any(frth$outcome == "regenerated" & as.numeric(frth$OR) > 1))

# ---- E. model-selection coherence (exp2_04 diagnoses inference model) -------
section("E. exp2_04 diagnoses the exp2_01 inference model")
allm <- rd("output/exp2_models/ALL_OUTCOMES_SUMMARY.csv")
dh   <- rd("output/exp2_diagnostics/dharma_summary.csv")
key  <- function(df) paste(df$Species, df$Outcome)
m1 <- setNames(allm$Method, key(allm)); m4 <- setNames(dh$Method, key(dh))
common <- intersect(names(m1), names(m4))
ok("DHARMa Method == exp2_01 Method for every outcome",
   length(common) > 0 && all(m1[common] == m4[common]))

# ---- F. canonical CSV values ------------------------------------------------
section("F. canonical CSV values")
ok("regen OR matches canonical 52.2 in CSV", near(ror("Dremel|Scrape|treatment"), 52.21, 0.03))
ok("healed Firth OR matches canonical 17.3 in CSV", near(hor("Dremel|Scrape|treatment"), 17.30, 0.03))
ok("RFP joint-interaction P matches canonical ~1.5e-6 in CSV", por$Joint_Interaction_P[por$Outcome=="RFP"] < 1e-4)

# ---- G. repo hygiene --------------------------------------------------------
section("G. repo hygiene")
ok("no stray Rplots.pdf at repo root", !file.exists(here("Rplots.pdf")))
ok("UPPERCASE canonical Exp2 tables exist",
   file.exists(here("output/exp2_tables/PUBLICATION_STATISTICS_TABLE.csv")) &&
   file.exists(here("output/exp2_models/ALL_OUTCOMES_SUMMARY.csv")))
cff <- readLines(here("CITATION.cff"))
ok("CITATION.cff: attributed (>=4 authors with ORCIDs, no XXXX placeholder)",
   sum(grepl("family-names:", cff)) >= 4 &&
   sum(grepl("orcid:", cff)) >= 4 &&
   !any(grepl("XXXX", cff)))
ok("LICENSE names the copyright holder",
   { lic <- readLines(here("LICENSE"))
     any(grepl("Copyright \\(c\\) [0-9]{4} .+[A-Za-z]", lic)) &&
     !any(grepl("information removed", lic)) })
# (manuscript/ intentionally absent — code-only public analysis repo)

# ---- H. model diagnostics (DHARMa, separation, spline-df) -------------------
section("H. model diagnostics & supplement")
dh <- rd("output/tables/exp1_dharma_summary.csv")
ok("Exp1 DHARMa: 3 outcomes, residuals clean (uniformity/dispersion/outlier p > 0.05)",
   nrow(dh) == 3 && all(dh$uniformity_KS_p > 0.05) && all(dh$dispersion_p > 0.05) && all(dh$outlier_p > 0.05))
ok("Exp1 DHARMa: day-aggregated temporal autocorrelation present (p < 0.05; documented)",
   all(dh$temporal_p < 0.05))
sp <- rd("output/tables/exp1_separation_detection.csv")
ok("No complete separation in any Exp1 model after re-score (Firth = sensitivity)",
   all(sp$complete_separation == FALSE))
ss <- rd("output/exp2_tables/spline_df_sensitivity_all_outcomes.csv")
ok("Exp2 spline-df sensitivity covers 6 outcomes x 3 df (18 rows)", nrow(ss) == 18)
ok("Exp2 RFP time interaction robust across spline df (p < 0.05 at df 3/4/5)",
   all(ss$min_interaction_p[ss$outcome == "RFP"] < 0.05))
# (Word supplement builder intentionally excluded — analysis code only)

# ---- summary ----------------------------------------------------------------
cat(sprintf("\n================  %d passed, %d failed  ================\n", .PASS, .FAIL))
if (.FAIL > 0L) quit(status = 1L, save = "no")
cat("ALL CHECKS PASSED\n")
