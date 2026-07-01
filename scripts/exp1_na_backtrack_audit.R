# =============================================================================
# exp1_na_backtrack_audit.R — auditable record of the 2026-06-15 re-score
# -----------------------------------------------------------------------------
# Diffs the pre-backtrack raw data against the current data and classifies every
# changed analytical cell (regenerated/healed/debris) against a t-1/t+1
# neighbour-cell rule: for each cell it records the t-1 and t+1 observed states, the
# value the rule would infer, the re-scored value, and a classification
# (NA-fill matches rule / NA-fill beyond rule / re-score / left NA).
# Backs Supplementary Tables S1.5 (backtrack audit) and S2.6 (separation log).
# =============================================================================
suppressMessages({library(tidyverse); library(here)})

norm <- function(p) read_csv(p, show_col_types = FALSE) |>
  mutate(across(c(species,treatment,coral_id,regenerated,healed,debris),
                ~na_if(str_squish(as.character(.x)),"NA")))
base <- norm(here("data","archive_exp1_2022","airbrush_dremel_pre-backtrack-2026-06-18.csv"))
cur  <- norm(here("data","airbrush_dremel.csv"))
key  <- c("species","treatment","coral_id","day")

# bracketing prediction from the PRE-backtrack data (regeneration absorbing)
predict_fragment <- function(v, day, absorbing = NULL) {
  o <- order(day); v <- v[o]; res <- rep(NA_character_, length(v))
  for (i in which(is.na(v))) {
    p <- if (i>1  && any(!is.na(v[seq_len(i-1)])))      v[max(which(!is.na(v[seq_len(i-1)])))] else NA
    n <- if (i<length(v) && any(!is.na(v[(i+1):length(v)]))) v[i+which(!is.na(v[(i+1):length(v)]))[1]] else NA
    if (!is.null(absorbing) && !is.na(p) && p==absorbing) { res[i] <- absorbing; next }
    if (!is.na(p) && !is.na(n) && p==n)                   { res[i] <- p; next }
  }
  res[order(o)]
}
# neighbour states for annotation
neighbours <- function(v, day) {
  o <- order(day); v <- v[o]; pv <- nv <- rep(NA_character_, length(v))
  for (i in seq_along(v)) {
    pv[i] <- if (i>1  && any(!is.na(v[seq_len(i-1)])))      v[max(which(!is.na(v[seq_len(i-1)])))] else NA
    nv[i] <- if (i<length(v) && any(!is.na(v[(i+1):length(v)]))) v[i+which(!is.na(v[(i+1):length(v)]))[1]] else NA
  }
  list(prev=pv[order(o)], nxt=nv[order(o)])
}

enrich <- function(df) df |> arrange(coral_id,day) |> group_by(coral_id) |> mutate(
  pred_regenerated = predict_fragment(regenerated,day,absorbing="yes"),
  pred_healed = predict_fragment(healed,day), pred_debris = predict_fragment(debris,day),
  prev_regenerated = neighbours(regenerated,day)$prev, nxt_regenerated = neighbours(regenerated,day)$nxt,
  prev_healed = neighbours(healed,day)$prev, nxt_healed = neighbours(healed,day)$nxt,
  prev_debris = neighbours(debris,day)$prev, nxt_debris = neighbours(debris,day)$nxt) |> ungroup()
be <- enrich(base)

audit <- map_dfr(c("regenerated","healed","debris"), function(col) {
  b <- be |> select(all_of(key), base=all_of(col),
                    prev=all_of(paste0("prev_",col)), nxt=all_of(paste0("nxt_",col)),
                    pred=all_of(paste0("pred_",col)))
  m <- cur |> select(all_of(key), now=all_of(col))
  full_join(b,m,by=key) |> mutate(outcome=col)
}) |>
  filter((is.na(base)!=is.na(now)) | (!is.na(base)&!is.na(now)&base!=now)) |>
  mutate(classification = case_when(
    !is.na(base) & !is.na(now)            ~ "RE-SCORE (already-scored cell)",
    is.na(base) & !is.na(now) & !is.na(pred) & pred==now ~ "NA-fill: matches bracketing rule",
    is.na(base) & !is.na(now) & !is.na(pred) & pred!=now ~ "NA-fill: DISAGREES with rule",
    is.na(base) & !is.na(now) & is.na(pred)              ~ "NA-fill: beyond rule (unresolved)",
    is.na(base) & is.na(now)                              ~ "left NA",
    TRUE ~ "other"),
    separation_relevant = (outcome=="healed" & species=="pocillopora" & treatment=="dremel" & now=="yes")) |>
  select(species,treatment,coral_id,day,outcome,
         t_minus1=prev, t_plus1=nxt, prior_value=base, rule_predicts=pred,
         rescored_value=now, classification, separation_relevant) |>
  arrange(classification,coral_id,day)

cat("\n================ NA-backtrack + separation audit ================\n")
print(audit, n=Inf, width=Inf)
cat("\n--- tally ---\n"); audit |> count(outcome, classification) |> print(n=Inf)
cat("\n--- separation-relevant cells (Pocillopora x scrape -> healed=yes) ---\n")
audit |> filter(separation_relevant) |> select(coral_id,day,prior_value,rescored_value,classification) |> print()
write_csv(audit, here("output","tables","na_backtrack_audit_2026-06-18.csv"))
cat("\nSaved: output/tables/na_backtrack_audit_2026-06-18.csv\n")
