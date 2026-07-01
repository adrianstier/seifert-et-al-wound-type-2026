#!/usr/bin/env Rscript
# Install the packages required to reproduce this repo.
# This is a lightweight stand-in for a pinned renv.lock (see README "Known
# open items" — a true renv snapshot is the recommended next step).
# Exact environment of the last successful Exp 1 run: output/text/sessionInfo.txt

# Package list reconciled against actual library()/require()/pkg:: calls in
# all scripts/*.R (excluding scripts/archive/) on 2026-05-15.
# Note: dplyr, ggplot2, tidyr, stringr, readr, purrr are attached via
# tidyverse; broom and tidyselect are listed explicitly because they are
# loaded/used directly. stats, utils, and splines are base R (splines is
# retained in the vector for the requireNamespace check but skipped on install).
pkgs <- c(
  # core
  "tidyverse", "here", "janitor", "readxl", "glue", "scales",
  "broom", "tidyselect",
  # Experiment 1 (authoritative)
  "lme4", "brglm2", "ordinal", "emmeans", "broom.mixed", "MuMIn", "DHARMa",
  "gt", "ggpubr",
  # Experiment 2
  "glmmTMB", "patchwork", "splines"           # splines ships with base R
)
need <- setdiff(pkgs, rownames(installed.packages()))
if (length(need)) {
  message("Installing: ", paste(need, collapse = ", "))
  install.packages(need, repos = "https://cloud.r-project.org")
} else {
  message("All required packages already installed.")
}
invisible(lapply(setdiff(pkgs, "splines"), requireNamespace, quietly = TRUE))
message("Dependency check complete. R ", getRversion(),
        " (last validated run: R 4.5.2).")
