# =============================================================================
# 03_lock_manuscript_values.R -- derive and lock the display values of the
# manuscript from the outputs of 01 and 02 (machine-readable value lock).
#
# Produces (all under output/):
#   latitude_support_grouped.csv        Section 3.6 (latitude gradient, all 47 vs supported)
#   subset_declines_grouped.csv         pre-to-post declines by prefecture subset
#   severity_declines_grouped.csv       42.9% / 43.3% (Section 3.7)
#   blup_extremes_grouped.csv, blup_change_extremes_grouped.csv   Sections 3.5 / 3.8
#   pooled_national_declines_grouped.csv
#   MANUSCRIPT_VALUE_LOCK.csv           every display value, with its source file
#   lock_mismatches.csv                 only if the rerun differs from results/
#
# If results/MANUSCRIPT_VALUE_LOCK.csv (the archived lock of the manuscript
# run) is present, the rerun is compared with it key by key.
# Prerequisites: 01_core_reanalysis.R and 02_attributable_reanalysis.R.
# =============================================================================

source("R/00_setup.R")
start_log("03_lock_manuscript_values")

obj <- readRDS(file.path(OUT_DIR, "grouped_stage1_main_severity.rds"))
MAIN <- obj$MAIN; SEV <- obj$SEV; NON <- obj$NON

pm <- read.csv(file.path(DATA_DIR, "pref_metadata.csv"))
support <- read_out("p95_support.csv")
subsets <- read_out("subset_sensitivity_grouped.csv")
sev <- read_out("severity_pooled_grouped.csv")
af <- read_out("attributable_summary_grouped.csv")
death <- read_out("pooled_national_death_grouped.csv")
blup_df <- read_out("prefecture_blup_p95_grouped.csv")
chg <- read_out("prefecture_change_grouped.csv")

cat("=== LATITUDE SUPPORT SENSITIVITY ===\n")
lat_support <- list()
for (p in periods) for (set_name in c("all_47", "observed_p95_only")) {
  ids <- if (set_name == "all_47") prefectures else
    support$pref_code[support$period == p & support$supported]
  ids <- intersect(ids, converged_ids(MAIN, p))
  yy <- do.call(rbind, MAIN$coef[[p]][ids]); ss <- MAIN$vcov[[p]][ids]
  x <- pm$latitude[match(ids, pm$pref_code)]
  f0 <- mixmeta(yy ~ 1, S = ss, method = "reml")
  f1 <- mixmeta(yy ~ x, S = ss, method = "reml")
  cc <- coef(f1); VV <- vcov(f1); ix <- grep("[.]x$", names(cc))
  ww <- wald_block(cc, VV, ix)
  t0 <- sum(diag(as.matrix(f0$Psi))); t1 <- sum(diag(as.matrix(f1$Psi)))
  lat_support[[length(lat_support) + 1]] <- data.frame(
    period = p, set = set_name, n = length(ids), chi2 = ww["chi2"],
    df = ww["df"], p = ww["p"], pct_explained = max(0, (t0 - t1) / t0) * 100)
}
lat_support <- do.call(rbind, lat_support)
print(lat_support); write_csv(lat_support, "latitude_support_grouped.csv")

cat("\n=== DECLINES BY SET ===\n")
set_declines <- do.call(rbind, lapply(unique(subsets$set), function(s) {
  z <- subsets[subsets$set == s, ]
  pre <- z$rr[z$period == "pre_pandemic"]
  post <- z$rr[z$period == "post_pandemic"]
  data.frame(set = s, decline_pct = 100 * (1 - post / pre))
}))
print(set_declines); write_csv(set_declines, "subset_declines_grouped.csv")

cat("\n=== SEVERITY DECLINES ===\n")
sev_declines <- do.call(rbind, lapply(unique(sev$outcome), function(s) {
  z <- sev[sev$outcome == s, ]
  pre <- z$rr[z$period == "pre_pandemic"]
  post <- z$rr[z$period == "post_pandemic"]
  data.frame(outcome = s, decline_pct = 100 * (1 - post / pre))
}))
print(sev_declines); write_csv(sev_declines, "severity_declines_grouped.csv")

cat("\n=== BLUP EXTREMES ===\n")
extremes <- list()
for (p in periods) {
  z <- blup_df[blup_df$period == p, ]
  extremes[[length(extremes) + 1]] <- data.frame(
    period = p, min_rr = min(z$rr), min_pref = z$pref_name[which.min(z$rr)],
    max_rr = max(z$rr), max_pref = z$pref_name[which.max(z$rr)])
}
extremes <- do.call(rbind, extremes)
change_extremes <- data.frame(
  min_change = min(chg$pct_change), min_pref = chg$pref_name[which.min(chg$pct_change)],
  max_change = max(chg$pct_change), max_pref = chg$pref_name[which.max(chg$pct_change)],
  n_declined = sum(chg$pct_change < 0))
print(extremes); print(change_extremes)
write_csv(extremes, "blup_extremes_grouped.csv")
write_csv(change_extremes, "blup_change_extremes_grouped.csv")

cat("\n=== DEATH DECLINE ===\n")
death_declines <- do.call(rbind, lapply(unique(death$outcome), function(o) {
  z <- death[death$outcome == o, ]
  data.frame(outcome = o,
             decline_pct = 100 * (1 - z$rr[z$period == "post_pandemic"] /
                                         z$rr[z$period == "pre_pandemic"]))
}))
print(death_declines); write_csv(death_declines, "pooled_national_declines_grouped.csv")

cat("\n=== ATTRIBUTABLE DISPLAY VALUES ===\n")
af_main <- af[af$outcome == "all_severity" & af$reference == "P50", ]
af_sev <- af[af$outcome == "severe" & af$reference == "P50", ]
print(af_main); print(af_sev)

# Long-form, machine-readable value lock. Each row names the source file.
lock <- list()
add <- function(key, value, display, source) {
  lock[[length(lock) + 1]] <<- data.frame(key = key, value = as.numeric(value),
                                          display = display, source = source)
}

t2 <- read_out("table2_grouped_exact.csv")
for (j in seq_len(nrow(t2))) {
  z <- t2[j, ]
  stem <- paste("table2", tolower(z$percentile), z$period, sep = ".")
  add(paste0(stem, ".rr"), z$rr, sprintf("%.2f", z$rr), "table2_grouped_exact.csv")
  add(paste0(stem, ".lo"), z$lo, sprintf("%.2f", z$lo), "table2_grouped_exact.csv")
  add(paste0(stem, ".hi"), z$hi, sprintf("%.2f", z$hi), "table2_grouped_exact.csv")
}

pc <- read_out("period_contrasts_independent.csv")
for (j in seq_len(nrow(pc))) {
  z <- pc[j, ]; stem <- gsub(" ", "_", z$contrast)
  add(paste0("contrast.", stem, ".ratio"), z$ratio, sprintf("%.3f", z$ratio), "period_contrasts_independent.csv")
  add(paste0("contrast.", stem, ".lo"), z$lo, sprintf("%.3f", z$lo), "period_contrasts_independent.csv")
  add(paste0("contrast.", stem, ".hi"), z$hi, sprintf("%.3f", z$hi), "period_contrasts_independent.csv")
  add(paste0("contrast.", stem, ".decline"), z$decline_pct, sprintf("%.1f%%", z$decline_pct), "period_contrasts_independent.csv")
}

for (j in seq_len(nrow(sev))) {
  z <- sev[j, ]; stem <- paste("severity", z$outcome, z$period, sep = ".")
  add(paste0(stem, ".rr"), z$rr, sprintf("%.2f", z$rr), "severity_pooled_grouped.csv")
  add(paste0(stem, ".lo"), z$lo, sprintf("%.2f", z$lo), "severity_pooled_grouped.csv")
  add(paste0(stem, ".hi"), z$hi, sprintf("%.2f", z$hi), "severity_pooled_grouped.csv")
}

for (j in seq_len(nrow(af_main))) {
  z <- af_main[j, ]; stem <- paste("af", z$period, sep = ".")
  add(paste0(stem, ".an"), z$an, format(round(z$an), big.mark = ","), "attributable_summary_grouped.csv")
  add(paste0(stem, ".af"), z$af, sprintf("%.1f%%", z$af), "attributable_summary_grouped.csv")
  add(paste0(stem, ".annual_an"), z$an_annual, format(round(z$an_annual), big.mark = ","), "attributable_summary_grouped.csv")
}

lock <- do.call(rbind, lock); rownames(lock) <- NULL
write_csv(lock, "MANUSCRIPT_VALUE_LOCK.csv")

stopifnot(nrow(lock) > 50, change_extremes$n_declined == 47,
          all(sev_declines$decline_pct > 0))

# ---- Comparison with the archived lock of the manuscript run -----------------
cat("\n=== LOCK CHECK AGAINST results/MANUSCRIPT_VALUE_LOCK.csv ===\n")
status <- "PASS"
archived <- file.path(RESULTS_DIR, "MANUSCRIPT_VALUE_LOCK.csv")
if (file.exists(archived)) {
  ref <- read.csv(archived, stringsAsFactors = FALSE)
  cmp <- merge(ref[, c("key", "value", "display")], lock[, c("key", "value", "display")],
               by = "key", suffixes = c("_archived", "_rerun"), all = TRUE)
  cmp$display_match <- !is.na(cmp$display_archived) & !is.na(cmp$display_rerun) &
    trimws(cmp$display_archived) == trimws(cmp$display_rerun)
  cmp$value_match <- !is.na(cmp$value_archived) & !is.na(cmp$value_rerun) &
    abs(cmp$value_archived - cmp$value_rerun) <= 1e-6 * pmax(1, abs(cmp$value_archived))
  mism <- cmp[!(cmp$display_match & cmp$value_match), ]
  if (nrow(mism) == 0) {
    cat(sprintf("LOCK_CHECK|PASS|%d keys identical to the archived manuscript lock\n", nrow(cmp)))
  } else {
    print(mism[, c("key", "display_archived", "display_rerun", "value_archived", "value_rerun")])
    write_csv(mism, "lock_mismatches.csv")
    cat(sprintf("LOCK_CHECK|FAIL|%d of %d keys differ (see output/lock_mismatches.csv)\n",
                nrow(mism), nrow(cmp)))
    status <- "FAIL"
    warning("Rerun differs from results/MANUSCRIPT_VALUE_LOCK.csv; see output/lock_mismatches.csv")
  }
} else {
  cat("LOCK_CHECK|SKIPPED|results/MANUSCRIPT_VALUE_LOCK.csv not present\n")
}

end_log("03_lock_manuscript_values", status)
