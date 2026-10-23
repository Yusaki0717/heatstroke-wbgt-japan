# =============================================================================
# 05_supplementary_descriptives.R -- descriptive tables and supplementary items
# that are derived directly from the analysis data and the fitted objects of 01.
#
# NOTE ON PROVENANCE. Unlike 01-04, which are the archived scripts of the
# manuscript run (2026-08-25) with paths adapted for this repository, this
# script was written for the public release to make the remaining descriptive
# items of the manuscript reproducible from the same data objects. Its values
# can be checked against the manuscript directly.
#
# Produces (all under output/):
#   table1_descriptives.csv              Table 1 (incl. convergence row, strata counts)
#   tableS3_days_above_30.csv            Table S3
#   alert_threshold_proxy.csv            33 C level in the observed distribution
#   alert_threshold_by_year.csv          prefecture-days at/above 33 C by year
#   alert_coverage_prefecture_days.csv   policy-coverage prefecture-days (Response, Comment 7)
#   attributable_share_ge33.csv          share of heat-attributable EAD on days >= 33 C
#   figures/FigureS1_grouped.png/.pdf    Figure S1 (RR vs latitude, fitted meta-regression)
# Prerequisites: 01_core_reanalysis.R (RDS + CSV outputs).
# =============================================================================

source("R/00_setup.R")
load_attrdl()
suppressMessages(library(ggplot2))
start_log("05_supplementary_descriptives")

obj <- readRDS(file.path(OUT_DIR, "grouped_stage1_main_severity.rds"))
MAIN <- obj$MAIN

# ---- Table 1: descriptive statistics by period ------------------------------
cat("=== TABLE 1 ===\n")
describe_block <- function(d) {
  d %>% summarise(
    prefecture_days          = n(),
    prefecture_days_obs_wbgt = sum(!is.na(wbgt_max)),
    year_month_strata        = n_distinct(paste(pref_code, year, month)),
    total_ead                = sum(total),
    daily_mean_ead           = mean(total),
    daily_median_ead         = median(total),
    daily_p95_ead            = as.numeric(quantile(total, 0.95)),
    daily_max_ead            = max(total),
    pct_elderly              = 100 * sum(age_elderly) / sum(total),
    pct_severe_plus_death    = 100 * sum(severe_total) / sum(total),
    pct_death                = 100 * sum(death_total) / sum(total),
    wbgt_mean                = mean(wbgt_max, na.rm = TRUE),
    wbgt_sd                  = sd(wbgt_max, na.rm = TRUE),
    wbgt_median              = median(wbgt_max, na.rm = TRUE),
    wbgt_p25                 = as.numeric(quantile(wbgt_max, 0.25, na.rm = TRUE)),
    wbgt_p75                 = as.numeric(quantile(wbgt_max, 0.75, na.rm = TRUE)),
    wbgt_p90                 = as.numeric(quantile(wbgt_max, 0.90, na.rm = TRUE)),
    wbgt_p95                 = as.numeric(quantile(wbgt_max, 0.95, na.rm = TRUE)),
    wbgt_min                 = min(wbgt_max, na.rm = TRUE),
    wbgt_max_obs             = max(wbgt_max, na.rm = TRUE),
    .groups = "drop")
}
t1_periods <- dat %>% group_by(period) %>% describe_block() %>%
  mutate(period = as.character(period))
t1_total <- dat %>% mutate(period = "total") %>% group_by(period) %>% describe_block()
t1 <- bind_rows(t1_periods, t1_total)
conv <- vapply(periods, function(p) length(converged_ids(MAIN, p)), integer(1))
t1$dlnm_converged_prefectures <- c(conv[t1_periods$period], NA_integer_)
print(as.data.frame(t1))
write_csv(t1, "table1_descriptives.csv")
cat(sprintf("TABLE1|prefecture_days=%d|observed_wbgt=%d|missing_wbgt_days=%d|ead_on_missing=%d|strata=%s\n",
            nrow(dat), sum(!is.na(dat$wbgt_max)), sum(is.na(dat$wbgt_max)),
            sum(dat$total[is.na(dat$wbgt_max)]),
            paste(t1_periods$year_month_strata, collapse = "/")))
stopifnot(nrow(t1) == 4L, all(conv == 47L))

# ---- Table S3: days above 30 C for the three most truncated prefectures ------
cat("\n=== TABLE S3 ===\n")
s3 <- dat %>% filter(pref_code %in% c(1, 2, 35)) %>%
  group_by(pref_code, pref_name, period) %>%
  summarise(days_above_30 = sum(wbgt_max > 30, na.rm = TRUE),
            summer_days   = sum(!is.na(wbgt_max)),
            pct_summer_days = 100 * days_above_30 / summer_days,
            max_wbgt = max(wbgt_max, na.rm = TRUE),
            days_above_upper_knot = sum(wbgt_max > knots_var[2], na.rm = TRUE),
            .groups = "drop")
print(as.data.frame(s3)); write_csv(s3, "tableS3_days_above_30.csv")

# ---- 33 C alert threshold as a proxy (Section 2.3.6, Response Comment 7) ----
cat("\n=== 33 C THRESHOLD PROXY ===\n")
THR <- 33
obs <- dat %>% filter(!is.na(wbgt_max))
proxy_period <- obs %>% group_by(period) %>%
  summarise(prefecture_days = n(), days_ge_33 = sum(wbgt_max >= THR),
            pct_ge_33 = 100 * days_ge_33 / prefecture_days, .groups = "drop") %>%
  mutate(period = as.character(period))
proxy_all <- data.frame(period = "total", prefecture_days = nrow(obs),
                        days_ge_33 = sum(obs$wbgt_max >= THR),
                        pct_ge_33 = 100 * mean(obs$wbgt_max >= THR))
proxy <- bind_rows(proxy_period, proxy_all)
proxy$percentile_of_33C <- 100 * mean(wbgt_all < THR)
print(proxy); write_csv(proxy, "alert_threshold_proxy.csv")
by_year <- obs %>% group_by(year) %>%
  summarise(prefecture_days = n(), days_ge_33 = sum(wbgt_max >= THR),
            prefectures_ge_33 = n_distinct(pref_code[wbgt_max >= THR]), .groups = "drop")
print(by_year); write_csv(by_year, "alert_threshold_by_year.csv")

# Heatstroke Alert coverage: Kanto-Koshin pilot from 1 July 2020 (Ibaraki,
# Tochigi, Gunma, Saitama, Chiba, Tokyo, Kanagawa, Yamanashi, Nagano),
# nationwide from April 2021.
pilot <- c(8, 9, 10, 11, 12, 13, 14, 19, 20)
coverage <- dat %>%
  mutate(coverage = case_when(
    year <= 2019 ~ "none (2015-2019)",
    year == 2020 & pref_code %in% pilot & month >= 7 ~ "Kanto-Koshin pilot (Jul-Sep 2020)",
    year == 2020 ~ "none (2020 remainder)",
    year == 2021 ~ "nationwide (2021)",
    TRUE ~ "nationwide (2022-2023)")) %>%
  group_by(period, coverage) %>%
  summarise(calendar_prefecture_days = n(),
            observed_wbgt_prefecture_days = sum(!is.na(wbgt_max)), .groups = "drop")
print(as.data.frame(coverage)); write_csv(coverage, "alert_coverage_prefecture_days.csv")

# Share of heat-attributable dispatches (median reference) falling on days with
# WBGTmax >= 33 C, from the day-level attributable numbers of the first-stage fits.
cat("\n=== ATTRIBUTABLE SHARE ON DAYS >= 33 C ===\n")
grouped_lag_matrix <- function(x, group, lag_days) {
  ans <- matrix(NA_real_, nrow = length(x), ncol = lag_days + 1)
  for (g in unique(group)) {
    ii <- which(group == g); z <- x[ii]
    for (ell in 0:lag_days) {
      ans[ii, ell + 1] <- if (ell == 0) z else c(rep(NA_real_, ell), head(z, -ell))
    }
  }
  ans
}
share_rows <- list()
for (p in periods) for (i in prefectures) {
  d <- dat %>% filter(pref_code == i, period == p, !is.na(wbgt_max)) %>% arrange(date)
  if (nrow(d) < 60 || sum(d$total) < 10) next
  cb <- crossbasis(d$wbgt_max, lag = 5,
                   argvar = list(fun = "ns", knots = knots_var, Boundary.knots = WBGT_RANGE),
                   arglag = list(fun = "ns", knots = knots_lag),
                   group = factor(d$year))
  mo <- gnm(total ~ cb + factor(dow) + holiday, eliminate = d$stratum,
            family = quasipoisson(), data = d)
  lagmat <- grouped_lag_matrix(d$wbgt_max, d$year, 5)
  an_day <- attrdl(lagmat, cb, d$total, mo, type = "an", tot = FALSE,
                   cen = wbgt_cen, range = c(wbgt_cen, max(d$wbgt_max)))
  ok <- !is.na(an_day)
  # Day-level AN is undefined on the first lag days of each summer (NA). As in
  # attrdl(tot = TRUE), rescale each prefecture's sum to its full dispatch
  # denominator so that the totals equal those of 02_attributable_reanalysis.R.
  scale <- sum(d$total) / sum(d$total[ok])
  share_rows[[length(share_rows) + 1]] <- data.frame(
    period = p, pref_code = i,
    an_total = sum(an_day[ok]) * scale,
    an_ge_33 = sum(an_day[ok & d$wbgt_max >= THR]) * scale,
    days = sum(ok), days_ge_33 = sum(ok & d$wbgt_max >= THR))
}
share <- do.call(rbind, share_rows)
share_summary <- bind_rows(
  share %>% group_by(period) %>%
    summarise(across(c(an_total, an_ge_33, days, days_ge_33), sum), .groups = "drop") %>%
    mutate(period = as.character(period)),
  share %>% summarise(across(c(an_total, an_ge_33, days, days_ge_33), sum)) %>%
    mutate(period = "total")) %>%
  mutate(pct_days_ge_33 = 100 * days_ge_33 / days,
         pct_attributable_on_ge_33 = 100 * an_ge_33 / an_total)
# The national totals must equal the attributable numbers of script 02.
af02 <- read_out("attributable_summary_grouped.csv")
af02 <- af02[af02$outcome == "all_severity" & af02$reference == "P50", ]
stopifnot(isTRUE(all.equal(share_summary$an_total[match(af02$period, share_summary$period)],
                           af02$an, tolerance = 1e-8)))
print(as.data.frame(share_summary))
write_csv(share_summary, "attributable_share_ge33.csv")
stopifnot(nrow(share) == 141L)

# ---- Figure S1: prefecture RR at P95 against latitude ------------------------
cat("\n=== FIGURE S1 ===\n")
pm <- read.csv(file.path(DATA_DIR, "pref_metadata.csv"))
blup_df <- read_out("prefecture_blup_p95_grouped.csv")
mr <- read_out("metaregression_grouped.csv")
bd95 <- basis_difference(Q_REPORT[["P95"]])
lat_grid <- seq(min(pm$latitude), max(pm$latitude), length.out = 100)
period_titles <- c(pre_pandemic = "Pre-pandemic (2015-2019)",
                   pandemic = "Pandemic (2020-2021)",
                   post_pandemic = "Post-pandemic (2022-2023)")
fit_rows <- list(); ann_rows <- list()
for (p in periods) {
  ids <- converged_ids(MAIN, p)
  yy <- do.call(rbind, MAIN$coef[[p]][ids]); ss <- MAIN$vcov[[p]][ids]
  x <- pm$latitude[match(ids, pm$pref_code)]
  f1 <- mixmeta(yy ~ x, S = ss, method = "reml")
  pred <- as.matrix(predict(f1, newdata = data.frame(x = lat_grid)))
  fit_rows[[p]] <- data.frame(period = p, latitude = lat_grid,
                              rr = exp(as.numeric(pred %*% bd95)))
  z <- mr[mr$period == p & mr$variable == "latitude", ]
  p_txt <- if (z$p < 0.001) "P < 0.001" else sprintf("P = %.3f", z$p)
  ann_rows[[p]] <- data.frame(period = p,
    label = sprintf("Wald %s\n%.1f%% of between-prefecture\nvariance explained", p_txt, z$pct_explained))
}
fit_df <- bind_rows(fit_rows); ann_df <- bind_rows(ann_rows)
pts <- blup_df %>% left_join(pm[, c("pref_code", "latitude")], by = "pref_code")
for (d_ in list("pts", "fit_df", "ann_df")) {
  assign(d_, within(get(d_), period <- factor(period, levels = periods,
                                              labels = period_titles[periods])))
}
figS1 <- ggplot(pts, aes(latitude, rr)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, colour = "grey60", linewidth = 0.4) +
  geom_point(size = 1.8) +
  geom_line(data = fit_df, aes(latitude, rr), colour = "#B2182B", linewidth = 1) +
  geom_text(data = ann_df, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.05, vjust = 1.2, size = 3.2, inherit.aes = FALSE) +
  scale_y_log10() +
  facet_wrap(~ period, nrow = 1) +
  labs(x = expression("Latitude of prefectural capital ("*degree*"N)"),
       y = expression("Cumulative RR at WBGT P95 (31.9"*degree*"C) vs 27.5"*degree*"C, log scale")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
ggsave(file.path(FIG_OUT, "FigureS1_grouped.png"), figS1, width = 14, height = 5, dpi = 300)
ggsave(file.path(FIG_OUT, "FigureS1_grouped.pdf"), figS1, width = 14, height = 5)
cat("FIGS1|PASS\n")

end_log("05_supplementary_descriptives")
