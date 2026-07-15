# =============================================================================
# Temporal changes in the association between wet-bulb globe temperature and
# heat-related emergency ambulance dispatches across Japan, 2015-2023
#
# Two-stage distributed lag non-linear model (DLNM) + multivariate meta-analysis
# comparing the WBGTmax-EAD association across three periods (pre-pandemic,
# pandemic, post-pandemic).
#
# Author:  Pengyu Huang
# Paper:   Huang P, Zhang J. Environmental Health Perspectives (under review),
#          manuscript hp-2026-00307n.
# License: MIT (see LICENSE)
#
# ------------------------------ DATA (not included) --------------------------
# This script requires two publicly available datasets that are NOT distributed
# in this repository. Download them from the official sources and place them in
# ./data/ (see README.md for the exact file list and column specifications):
#
#   1. FDMA heat-related emergency ambulance dispatch data (9 annual .xlsx files)
#      https://www.fdma.go.jp/disaster/heatstroke/post3.html
#
#   2. Ministry of the Environment daily WBGT data, saved as
#      data/wbgt_daily_2015_2023.csv  (columns: date, pref_code, wbgt_max)
#      https://www.wbgt.env.go.jp/record_data.php
#
# ------------------------------ HOW TO RUN -----------------------------------
# 1. Place all data files in ./data/  (see README.md)
# 2. Set the working directory to the repository root, or open the .Rproj
# 3. source("R/heatstroke_dlnm_analysis.R")
# Outputs (tables + figures) are written to ./output/
# =============================================================================

################################################################################
# FDMA Heatstroke × WBGT DLNM Analysis
# Multi-period comparison: pre-pandemic / pandemic / post-pandemic
# 
# Structure:
#   Part 0: Setup & dependencies
#   Part 1: Data loading & cleaning (FDMA)
#   Part 2: WBGT data loading & merging
#   Part 3: Descriptive statistics
#   Part 4: First-stage prefecture-level DLNM
#   Part 5: Second-stage meta-analysis
#   Part 6: Period comparison & visualization
#   Part 7: Sensitivity analyses (notes; see also R/sensitivity_knots_S5.R)
#   Part 8: Meta-regression (notes)
################################################################################

# ==============================================================================
# PART 0: SETUP
# ==============================================================================

# Install if needed (run once)
# install.packages(c("readxl", "dplyr", "tidyr", "lubridate", "dlnm", "gnm",
#                     "mixmeta", "ggplot2", "patchwork", "sf", "viridis"))

library(readxl)
library(dplyr)
library(tidyr)
library(lubridate)
library(dlnm)
library(gnm)
library(mixmeta)
library(ggplot2)
library(patchwork)
library(splines)

# Create output directories if they do not exist
dir.create("output", showWarnings = FALSE)
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

# Set working directory to where your data files are
# setwd("path/to/your/data")

# ==============================================================================
# PART 1: LOAD & CLEAN FDMA DATA
# ==============================================================================

# --- 1.1 File mapping ---
file_map <- tibble(
  file = c(
    "data/heatstroke003_data_h27.xlsx",  # 2015
    "data/heatstroke003_data_h28.xlsx",  # 2016
    "data/heatstroke003_data_h29.xlsx",  # 2017
    "data/heatstroke003_data_h30.xlsx",  # 2018
    "data/heatstroke003_data_r1.xlsx",   # 2019
    "data/heatstroke003_data_r2.xlsx",   # 2020
    "data/heatstroke003_data_r3.xlsx",   # 2021
    "data/heatstroke003_data_r4.xlsx",   # 2022
    "data/heatstroke003_data_r5.xlsx"    # 2023
  ),
  year = 2015:2023
)

# --- 1.2 Read all files ---
read_fdma_file <- function(filepath, yr) {
  sheets <- excel_sheets(filepath)
  # Only keep June-September sheets for comparability (2020 lacks May)
  # Sheet names are like "2017_06", "2022_07", etc.
  
  dfs <- lapply(sheets, function(sh) {
    df <- read_excel(filepath, sheet = sh)
    # Standardize column names to English
    colnames(df) <- c(
      "date", "pref_code", "total",
      "age_neonate", "age_infant", "age_child", "age_adult", "age_elderly", "age_unknown",
      "sev_death", "sev_severe", "sev_moderate", "sev_mild", "sev_other",
      "loc_home", "loc_work1", "loc_work2", "loc_school",
      "loc_indoor_public", "loc_outdoor_public", "loc_road", "loc_other"
    )
    df$date <- as.Date(df$date)
    df
  })
  
  bind_rows(dfs)
}

cat("Loading FDMA data...\n")
fdma_all <- bind_rows(
  lapply(1:nrow(file_map), function(i) {
    cat("  Reading:", file_map$file[i], "\n")
    read_fdma_file(file_map$file[i], file_map$year[i])
  })
)

# --- 1.3 Filter to June-September only (统一分析窗口) ---
fdma <- fdma_all %>%
  filter(month(date) >= 6 & month(date) <= 9) %>%
  mutate(
    year  = year(date),
    month = month(date),
    doy   = yday(date),
    dow   = wday(date, label = FALSE),  # 1=Sun, 7=Sat
    # Period classification
    period = case_when(
      year %in% 2015:2019 ~ "pre_pandemic",
      year %in% 2020:2021 ~ "pandemic",
      year %in% 2022:2023 ~ "post_pandemic"
    ),
    period = factor(period, levels = c("pre_pandemic", "pandemic", "post_pandemic"))
  )

# --- 1.4 Prefecture name mapping ---
pref_names <- c(
  "Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima",
  "Ibaraki", "Tochigi", "Gunma", "Saitama", "Chiba", "Tokyo", "Kanagawa",
  "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano",
  "Gifu", "Shizuoka", "Aichi", "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo",
  "Nara", "Wakayama", "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi",
  "Tokushima", "Kagawa", "Ehime", "Kochi", "Fukuoka", "Saga", "Nagasaki",
  "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa"
)
pref_lookup <- tibble(pref_code = 1:47, pref_name = pref_names)

fdma <- fdma %>% left_join(pref_lookup, by = "pref_code")

# --- 1.5 Data quality checks ---
cat("\n=== DATA QUALITY CHECKS ===\n")
cat("Total rows:", nrow(fdma), "\n")
cat("Date range:", as.character(min(fdma$date)), "to", as.character(max(fdma$date)), "\n")
cat("Years:", paste(sort(unique(fdma$year)), collapse = ", "), "\n")
cat("Prefectures:", n_distinct(fdma$pref_code), "\n")
cat("Missing total:", sum(is.na(fdma$total)), "\n")
cat("Negative total:", sum(fdma$total < 0, na.rm = TRUE), "\n")

# Cross-check: age subtotals should equal total
fdma <- fdma %>%
  mutate(
    age_sum = age_neonate + age_infant + age_child + age_adult + age_elderly + age_unknown,
    age_check = abs(total - age_sum)
  )
cat("Age sum mismatches (>0):", sum(fdma$age_check > 0, na.rm = TRUE), "\n")

# Summary by period
cat("\nSummary by period:\n")
fdma %>%
  group_by(period) %>%
  summarise(
    years = paste(range(year), collapse = "-"),
    n_days = n_distinct(date),
    total_cases = sum(total),
    daily_mean = round(mean(total), 2),
    .groups = "drop"
  ) %>%
  print()


# ==============================================================================
# PART 2: WBGT DATA
# ==============================================================================

# -------------------------------------------------------------------------
# NOTE: WBGT data must be downloaded separately from:
#   https://www.wbgt.env.go.jp/record_data.php
# 
# You need daily WBGTmax for each prefecture capital, 2015-2023, June-Sept.
# Save as CSV with columns: date, pref_code, wbgt_max
# 
# If WBGT is not yet available, you can use JMA Tmax as a substitute
# for initial testing. The code below shows both options.
# -------------------------------------------------------------------------

#--- Option A: Load WBGT data (primary) ---
#Uncomment when you have the WBGT CSV ready:

wbgt <- read.csv("data/wbgt_daily_2015_2023.csv")
wbgt$date <- as.Date(wbgt$date)

dat <- fdma %>%
  left_join(wbgt, by = c("date", "pref_code"))

# --- Option B: Placeholder with simulated WBGT for testing code ---
# REMOVE THIS BLOCK once you have real WBGT data
#cat("\n*** WBGT data not yet loaded. Using simulated data for code testing. ***\n")
#cat("*** Replace with real WBGT data before running actual analysis. ***\n\n")

#set.seed(42)
#wbgt_sim <- fdma %>%
  #distinct(date, pref_code) %>%
 # mutate(
    # Simulate WBGT roughly correlated with latitude and season
    #lat_proxy = 45 - pref_code * 0.4,  # rough N-S gradient
   # season_effect = sin((yday(date) - 152) / 61 * pi) * 5,  # peak in late July
   # wbgt_max = 25 + season_effect - lat_proxy * 0.1 + rnorm(n(), 0, 2)
  #) %>%
  #select(date, pref_code, wbgt_max)

#dat <- fdma %>%
  #left_join(wbgt_sim, by = c("date", "pref_code"))

cat("Missing WBGT after merge:", sum(is.na(dat$wbgt_max)), "\n")

# --- Japanese public holidays (June-September) ---
# Marine Day, Mountain Day, and other holidays during summer months
holidays_jp <- as.Date(c(
  # 2015
  "2015-07-20", "2015-09-21", "2015-09-22", "2015-09-23",
  # 2016
  "2016-07-18", "2016-08-11", "2016-09-19", "2016-09-22",
  # 2017
  "2017-07-17", "2017-08-11", "2017-09-18", "2017-09-23",
  # 2018
  "2018-07-16", "2018-08-11", "2018-09-17", "2018-09-23", "2018-09-24",
  # 2019
  "2019-07-15", "2019-08-11", "2019-08-12", "2019-09-16", "2019-09-23",
  # 2020
  "2020-07-23", "2020-07-24", "2020-08-10", "2020-09-21", "2020-09-22",
  # 2021
  "2021-07-22", "2021-07-23", "2021-08-08", "2021-08-09", "2021-09-20", "2021-09-23",
  # 2022
  "2022-07-18", "2022-08-11", "2022-09-19", "2022-09-23",
  # 2023
  "2023-07-17", "2023-08-11", "2023-09-18", "2023-09-23"
))

dat <- dat %>%
  mutate(
    holiday = as.integer(date %in% holidays_jp),
    # Year-month stratum for conditional model
    stratum = paste(year, month, sep = "-"),
    stratum = factor(stratum)
  )


# ==============================================================================
# PART 3: DESCRIPTIVE STATISTICS
# ==============================================================================

cat("\n=== DESCRIPTIVE STATISTICS ===\n")

# --- Table 1: Summary by period ---
table1 <- dat %>%
  group_by(period) %>%
  summarise(
    n_pref_days  = n(),
    total_cases  = sum(total),
    daily_mean   = round(mean(total), 2),
    daily_median = median(total),
    daily_p95    = quantile(total, 0.95),
    daily_max    = max(total),
    wbgt_mean    = round(mean(wbgt_max, na.rm = TRUE), 1),
    wbgt_sd      = round(sd(wbgt_max, na.rm = TRUE), 1),
    wbgt_p50     = round(quantile(wbgt_max, 0.50, na.rm = TRUE), 1),
    wbgt_p95     = round(quantile(wbgt_max, 0.95, na.rm = TRUE), 1),
    pct_elderly  = round(sum(age_elderly) / sum(total) * 100, 1),
    pct_severe   = round(sum(sev_death + sev_severe) / sum(total) * 100, 1),
    .groups = "drop"
  )
print(table1)

# --- Table 1 by prefecture (for supplementary) ---
table1_pref <- dat %>%
  group_by(pref_code, pref_name, period) %>%
  summarise(
    total_cases = sum(total),
    daily_mean  = round(mean(total), 2),
    wbgt_mean   = round(mean(wbgt_max, na.rm = TRUE), 1),
    .groups = "drop"
  )

# Save Table 1
write.csv(table1, "output/table1_summary.csv", row.names = FALSE)
write.csv(table1_pref, "output/table1_by_prefecture.csv", row.names = FALSE)


# ==============================================================================
# PART 4: FIRST-STAGE PREFECTURE-LEVEL DLNM
# ==============================================================================

cat("\n=== FIRST-STAGE DLNM ===\n")

# --- 4.1 Define WBGT percentiles from FULL study period (for comparability) ---
wbgt_all <- dat$wbgt_max[!is.na(dat$wbgt_max)]
knots_var <- quantile(wbgt_all, c(0.50, 0.90))
knots_lag <- logknots(5, 2)  # 2 internal knots on log scale for lag 0-5

cat("Exposure knots (50th, 90th percentile of full period):",
    round(knots_var, 1), "\n")
cat("Lag knots:", round(knots_lag, 2), "\n")

# Reference value: WBGT at minimum morbidity (will estimate from data)
# For now use median as centering value
wbgt_cen <- median(wbgt_all)

# Store percentiles for later comparison
wbgt_p95  <- quantile(wbgt_all, 0.95)
wbgt_p75  <- quantile(wbgt_all, 0.75)
wbgt_p50  <- quantile(wbgt_all, 0.50)
wbgt_range <- range(wbgt_all)

cat("WBGT centering value:", round(wbgt_cen, 1), "\n")
cat("WBGT 95th percentile:", round(wbgt_p95, 1), "\n")
cat("WBGT range:", round(wbgt_range, 1), "\n")

# --- 4.2 First-stage: run DLNM for each prefecture × period ---
prefectures <- sort(unique(dat$pref_code))
periods <- levels(dat$period)

# Storage for coefficients and (co)variances
coef_list <- list()
vcov_list <- list()
conv_list <- list()  # convergence tracker

# Degrees of freedom for cross-basis
# Exposure: ns with 2 internal knots → df = 2 + 1 = 3
# Lag: ns with 2 internal knots on log scale → df = 2 + 1 = 3
# Total cb parameters: 3 × 3 = 9

for (p in periods) {
  cat("\n--- Period:", p, "---\n")
  
  coef_list[[p]] <- list()
  vcov_list[[p]] <- list()
  conv_list[[p]] <- list()
  
  for (i in prefectures) {
    pref_dat <- dat %>%
      filter(pref_code == i, period == p) %>%
      arrange(date) %>%
      # Need to drop NAs in wbgt for crossbasis
      filter(!is.na(wbgt_max))
    
    pref_name_i <- pref_lookup$pref_name[i]
    
    # Skip if insufficient data
    if (nrow(pref_dat) < 60 || sum(pref_dat$total) < 10) {
      cat("  Skipping", pref_name_i, "(n=", nrow(pref_dat),
          ", total=", sum(pref_dat$total), ")\n")
      coef_list[[p]][[i]] <- NULL
      vcov_list[[p]][[i]] <- NULL
      conv_list[[p]][[i]] <- FALSE
      next
    }
    
    tryCatch({
      # Cross-basis with FIXED knots (from full period)
      cb <- crossbasis(
        pref_dat$wbgt_max,
        lag = 5,
        argvar = list(fun = "ns", knots = knots_var),
        arglag = list(fun = "ns", knots = knots_lag)
      )
      
      # Conditional quasi-Poisson, stratified by year-month
      model <- gnm(
        total ~ cb + factor(dow) + holiday,
        eliminate = stratum,
        family = quasipoisson(),
        data = pref_dat
      )
      
      # Extract reduced coefficients (for meta-analysis)
      # Reduce to overall cumulative association
      red <- crossreduce(cb, model, cen = wbgt_cen)
      
      coef_list[[p]][[i]] <- coef(red)
      vcov_list[[p]][[i]] <- vcov(red)
      conv_list[[p]][[i]] <- TRUE
      
    }, error = function(e) {
      cat("  ERROR for", pref_name_i, ":", conditionMessage(e), "\n")
      coef_list[[p]][[i]] <<- NULL
      vcov_list[[p]][[i]] <<- NULL
      conv_list[[p]][[i]] <<- FALSE
    })
  }
  
  n_ok <- sum(sapply(conv_list[[p]], isTRUE))
  cat("  Converged:", n_ok, "/", length(prefectures), "\n")
}


# ==============================================================================
# PART 5: SECOND-STAGE META-ANALYSIS
# ==============================================================================

cat("\n=== SECOND-STAGE META-ANALYSIS ===\n")

# Storage for meta-analysis results
meta_results <- list()
blup_results <- list()

for (p in periods) {
  cat("\n--- Period:", p, "---\n")
  
  # Collect converged prefectures
  ok_idx <- which(sapply(conv_list[[p]], isTRUE))
  
  if (length(ok_idx) < 5) {
    cat("  Too few converged prefectures. Skipping meta-analysis.\n")
    next
  }
  
  # Stack coefficients and covariances
  ymat <- do.call(rbind, coef_list[[p]][ok_idx])
  Slist <- vcov_list[[p]][ok_idx]
  
  cat("  Pooling", nrow(ymat), "prefectures\n")
  
  # Random-effects meta-analysis
  meta_fit <- mixmeta(ymat ~ 1, S = Slist, method = "reml")
  
  cat("  Meta-analysis log-likelihood:", round(logLik(meta_fit), 2), "\n")
  
  meta_results[[p]] <- meta_fit
  
  # BLUP (Best Linear Unbiased Predictions) for each prefecture
  blup_results[[p]] <- blup(meta_fit, vcov = TRUE)
}


# ==============================================================================
# PART 6: VISUALIZATION & PERIOD COMPARISON
# ==============================================================================

cat("\n=== GENERATING FIGURES ===\n")

dir.create("output", showWarnings = FALSE)
dir.create("output/figures", showWarnings = FALSE)

# --- 6.1 Pooled exposure-response curves by period ---

# Create prediction basis (same structure as first-stage)
wbgt_pred <- seq(wbgt_range[1], wbgt_range[2], length.out = 100)
basis_pred <- ns(wbgt_pred, knots = knots_var,
                 Boundary.knots = wbgt_range)

# For each period, predict pooled curve
plot_data <- list()

for (p in periods) {
  if (is.null(meta_results[[p]])) next
  
  # Get pooled coefficients
  pooled_coef <- coef(meta_results[[p]])
  pooled_vcov <- vcov(meta_results[[p]])
  
  # Predict log-RR
  log_rr <- basis_pred %*% pooled_coef
  
  # Confidence intervals
  se <- sqrt(diag(basis_pred %*% pooled_vcov %*% t(basis_pred)))
  
  # Center at reference (wbgt_cen)
  ref_idx <- which.min(abs(wbgt_pred - wbgt_cen))
  log_rr_cen <- log_rr - log_rr[ref_idx]
  
  plot_data[[p]] <- tibble(
    wbgt    = wbgt_pred,
    rr      = exp(log_rr_cen),
    rr_low  = exp(log_rr_cen - 1.96 * se),
    rr_high = exp(log_rr_cen + 1.96 * se),
    period  = p
  )
}

plot_df <- bind_rows(plot_data) %>%
  mutate(period = factor(period, levels = periods,
                         labels = c("Pre-pandemic\n(2015-2019)",
                                    "Pandemic\n(2020-2021)",
                                    "Post-pandemic\n(2022-2023)")))

# Main figure: pooled exposure-response
fig1 <- ggplot(plot_df, aes(x = wbgt, y = rr)) +
  geom_ribbon(aes(ymin = rr_low, ymax = rr_high, fill = period), alpha = 0.2) +
  geom_line(aes(color = period), linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = wbgt_p95, linetype = "dotted", color = "grey40") +
  annotate("text", x = wbgt_p95 + 0.3, y = max(plot_df$rr) * 0.95,
           label = "95th\npercentile", size = 3, hjust = 0) +
  scale_color_manual(values = c("#2166AC", "#B2182B", "#4DAF4A")) +
  scale_fill_manual(values = c("#2166AC", "#B2182B", "#4DAF4A")) +
  labs(
    x = expression("Daily maximum WBGT ("*degree*"C)"),
    y = "Cumulative Relative Risk",
    title = "Pooled exposure-response: WBGT and heatstroke ambulance dispatches",
    subtitle = paste0("Reference: WBGT = ", round(wbgt_cen, 1), "°C (median). ",
                      "Random-effects meta-analysis across 47 prefectures"),
    color = "Period", fill = "Period"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("output/figures/fig1_pooled_exposure_response.pdf", fig1,
       width = 10, height = 7)
ggsave("output/figures/fig1_pooled_exposure_response.png", fig1,
       width = 10, height = 7, dpi = 300)

cat("  Fig 1 saved.\n")


# --- 6.2 Period comparison: RR at key percentiles ---

rr_comparison <- list()

for (p in periods) {
  if (is.null(meta_results[[p]])) next
  
  pooled_coef <- coef(meta_results[[p]])
  pooled_vcov <- vcov(meta_results[[p]])
  
  # 用和first-stage相同的basis重建预测
  basis_var <- ns(wbgt_pred, knots = knots_var,
                  Boundary.knots = wbgt_range)
  
  for (pctl_name in c("p75", "p90", "p95")) {
    pctl_val <- switch(pctl_name,
                       p75 = as.numeric(wbgt_p75),
                       p90 = as.numeric(quantile(wbgt_all, 0.90)),
                       p95 = as.numeric(wbgt_p95)
    )
    
    # 在预测网格上找最近的index
    idx_val <- which.min(abs(wbgt_pred - pctl_val))
    idx_ref <- which.min(abs(wbgt_pred - wbgt_cen))
    
    b_diff <- as.numeric(basis_var[idx_val, ] - basis_var[idx_ref, ])
    
    log_rr <- sum(b_diff * pooled_coef)
    se <- sqrt(as.numeric(t(b_diff) %*% pooled_vcov %*% b_diff))
    
    rr_comparison <- c(rr_comparison, list(data.frame(
      period = p,
      percentile = pctl_name,
      wbgt_value = round(pctl_val, 1),
      rr = round(exp(log_rr), 3),
      rr_low = round(exp(log_rr - 1.96 * se), 3),
      rr_high = round(exp(log_rr + 1.96 * se), 3)
    )))
  }
}

rr_table <- do.call(rbind, rr_comparison)
cat("\nCumulative RR at key WBGT percentiles (vs median):\n")
print(rr_table)
write.csv(rr_table, "output/rr_comparison_by_period.csv", row.names = FALSE)

# --- 6.3 Lag-response at 95th percentile of WBGT ---

lag_plot_data <- list()

for (p in periods) {
  if (is.null(meta_results[[p]])) next
  
  # Re-run one representative prefecture to get lag structure
  # Or better: use full cross-basis prediction from meta
  # Here we take a simpler approach: re-fit first-stage for a large prefecture
  # (e.g., Tokyo, pref_code=13) and extract lag-response
  
  pref_dat <- dat %>%
    filter(pref_code == 13, period == p, !is.na(wbgt_max)) %>%
    arrange(date)
  
  if (nrow(pref_dat) < 60) next
  
  cb <- crossbasis(
    pref_dat$wbgt_max,
    lag = 5,
    argvar = list(fun = "ns", knots = knots_var),
    arglag = list(fun = "ns", knots = knots_lag)
  )
  
  model <- gnm(
    total ~ cb + factor(dow) + holiday,
    eliminate = stratum,
    family = quasipoisson(),
    data = pref_dat
  )
  
  # Predict lag-response at 95th percentile
  cp <- crosspred(cb, model, at = wbgt_p95, cen = wbgt_cen)
  
  lag_plot_data[[p]] <- tibble(
    lag    = 0:5,
    rr     = cp$matRRfit[1, ],
    rr_low = cp$matRRlow[1, ],
    rr_high = cp$matRRhigh[1, ],
    period = p
  )
}

lag_df <- bind_rows(lag_plot_data) %>%
  mutate(period = factor(period, levels = periods,
                         labels = c("Pre-pandemic\n(2015-2019)",
                                    "Pandemic\n(2020-2021)",
                                    "Post-pandemic\n(2022-2023)")))

fig2 <- ggplot(lag_df, aes(x = lag, y = rr, color = period, fill = period)) +
  geom_ribbon(aes(ymin = rr_low, ymax = rr_high), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  scale_x_continuous(breaks = 0:5) +
  scale_color_manual(values = c("#2166AC", "#B2182B", "#4DAF4A")) +
  scale_fill_manual(values = c("#2166AC", "#B2182B", "#4DAF4A")) +
  labs(
    x = "Lag (days)",
    y = "Relative Risk",
    title = paste0("Lag-response at WBGT 95th percentile (",
                   round(wbgt_p95, 1), "°C) — Tokyo"),
    subtitle = paste0("Reference: WBGT = ", round(wbgt_cen, 1), "°C"),
    color = "Period", fill = "Period"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave("output/figures/fig2_lag_response_tokyo.pdf", fig2,
       width = 8, height = 6)
ggsave("output/figures/fig2_lag_response_tokyo.png", fig2,
       width = 8, height = 6, dpi = 300)

cat("  Fig 2 saved.\n")


# --- 6.4 Descriptive time series plot ---

daily_national <- dat %>%
  group_by(date, period, year) %>%
  summarise(
    total = sum(total),
    wbgt_mean = mean(wbgt_max, na.rm = TRUE),
    .groups = "drop"
  )

fig3a <- ggplot(daily_national, aes(x = date, y = total, color = period)) +
  geom_line(alpha = 0.7, linewidth = 0.3) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  scale_color_manual(values = c("#2166AC", "#B2182B", "#4DAF4A")) +
  labs(x = NULL, y = "Daily national dispatches",
       title = "Heatstroke ambulance dispatches, Japan 2015-2023") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

fig3b <- ggplot(daily_national, aes(x = date, y = wbgt_mean, color = period)) +
  geom_line(alpha = 0.7, linewidth = 0.3) +
  scale_color_manual(values = c("#2166AC", "#B2182B", "#4DAF4A")) +
  labs(x = NULL, y = "Mean daily WBGTmax (°C)",
       color = "Period") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

fig3 <- fig3a / fig3b + plot_layout(heights = c(2, 1))

ggsave("output/figures/fig3_timeseries.pdf", fig3,
       width = 12, height = 7)
ggsave("output/figures/fig3_timeseries.png", fig3,
       width = 12, height = 7, dpi = 300)

cat("  Fig 3 saved.\n")


# ==============================================================================
# PART 7: SENSITIVITY ANALYSES
# ==============================================================================

cat("\n=== SENSITIVITY ANALYSIS NOTES ===\n")
cat("
The following sensitivity analyses should be run:
1. Alternative lag: re-run Part 4-5 with lag=7 and lag=10
   - Change: lag = 7 in crossbasis(), knots_lag <- logknots(7, 2)
   
2. Alternative exposure: replace wbgt_max with Tmax from JMA
   - Requires separate JMA data download
   
3. Exclude 2020 from pandemic period (most COVID-affected year)
   - Change: filter(year != 2020) in pandemic subset
   
4. Two-period comparison (pre vs post only, excluding pandemic)
   - Simpler comparison: 2015-2019 vs 2022-2023
   
5. Subgroup by severity: re-run using sev_severe + sev_death as outcome
   
6. Subgroup by age: re-run using age_elderly as outcome

To run any sensitivity analysis, modify the relevant parameters
and re-execute Parts 4-6.
\n")


# ==============================================================================
# PART 8: META-REGRESSION (for Environment International upgrade)
# ==============================================================================

cat("=== META-REGRESSION (OPTIONAL) ===\n")
cat("
To examine spatial heterogeneity drivers, add prefecture-level
meta-regression variables to the second-stage model:

  Prefecture-level predictors (obtain from e-Stat):
  - latitude of capital city
  - urbanization_rate (DID人口比率)
  - pct_elderly (65歳以上人口比率)
  - population_density
  
  Example code:
  meta_fit <- mixmeta(ymat ~ latitude + pct_elderly + urban_rate,
                      S = Slist, method = 'reml')
  
  This answers 'WHY do effects differ across prefectures?'
  — critical for Environment International vs Environmental Health.
\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Outputs saved to ./output/\n")
cat("Next steps:\n")
cat("  1. Download and integrate real WBGT data\n")
cat("  2. Run the full pipeline\n")
cat("  3. Check convergence and diagnostic plots\n")
cat("  4. Run sensitivity analyses\n")
source("table1_precise.R")