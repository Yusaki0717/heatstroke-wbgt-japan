################################################################################
# Sensitivity analyses reported in manuscript Table 4
#   - Lag 0-7               (max lag extended from 5 to 7 days)
#   - Excluding 2020        (pandemic period restricted to 2021)
#   - Elderly (>=65) only   (outcome: age_elderly)
# (Severe-case-only row is implemented in revision_sensitivity_analyses.R, C2.)
#
# Prerequisite: run heatstroke_dlnm_analysis.R first
# (needs dat, periods, prefectures, knots_var, wbgt_all, wbgt_cen in workspace).
# RR extraction follows the SAME 100-point grid path as Part 6.2 of the main
# script, so results are directly comparable with Table 2 / Table 4.
################################################################################

library(dlnm)
library(gnm)
library(mixmeta)
library(dplyr)
library(splines)

# --- Shared prediction basis: identical to main script Part 6.2 ----------------
wbgt_range_s <- range(wbgt_all)
wbgt_pred    <- seq(wbgt_range_s[1], wbgt_range_s[2], length.out = 100)
basis_var    <- ns(wbgt_pred, knots = knots_var, Boundary.knots = wbgt_range_s)
wbgt_p95_val <- as.numeric(quantile(wbgt_all, 0.95))
idx_val <- which.min(abs(wbgt_pred - wbgt_p95_val))
idx_ref <- which.min(abs(wbgt_pred - wbgt_cen))
bdiff   <- as.numeric(basis_var[idx_val, ] - basis_var[idx_ref, ])

# --- Generic two-stage runner -------------------------------------------------
# outcome:    column name of the outcome in dat ("total" or "age_elderly")
# lag:        maximum lag for crossbasis (5 = main, 7 = sensitivity)
# drop_years: calendar years to exclude before fitting (NULL = none)
run_variant <- function(outcome = "total", lag = 5, drop_years = NULL, label = "") {
  kl <- logknots(lag, 2)
  rows <- list()
  for (p in periods) {
    coefs <- list(); vcovs <- list(); n_conv <- 0L
    for (i in prefectures) {
      pref_dat <- dat %>%
        filter(pref_code == i, period == p, !is.na(wbgt_max))
      if (!is.null(drop_years)) {
        pref_dat <- pref_dat %>% filter(!(year %in% drop_years))
      }
      pref_dat <- pref_dat %>% arrange(date) %>% droplevels()
      if (nrow(pref_dat) < 60 || sum(pref_dat[[outcome]], na.rm = TRUE) < 10) next
      fit <- tryCatch({
        cb <- crossbasis(
          pref_dat$wbgt_max, lag = lag,
          argvar = list(fun = "ns", knots = knots_var),
          arglag = list(fun = "ns", knots = kl)
        )
        fml <- as.formula(paste(outcome, "~ cb + factor(dow) + holiday"))
        model <- gnm(fml, eliminate = stratum,
                     family = quasipoisson(), data = pref_dat)
        red <- crossreduce(cb, model, cen = wbgt_cen)
        list(coef = coef(red), vcov = vcov(red))
      }, error = function(e) NULL)
      if (!is.null(fit)) {
        coefs[[length(coefs) + 1L]] <- fit$coef
        vcovs[[length(vcovs) + 1L]] <- fit$vcov
        n_conv <- n_conv + 1L
      }
    }
    cat("  ", label, "|", p, ": converged", n_conv, "/", length(prefectures), "
")
    if (n_conv < 5) next
    ymat <- do.call(rbind, coefs)
    mfit <- mixmeta(ymat ~ 1, S = vcovs, method = "reml")
    est <- sum(bdiff * coef(mfit))
    se  <- sqrt(as.numeric(t(bdiff) %*% vcov(mfit) %*% bdiff))
    rows[[p]] <- data.frame(
      analysis = label, period = p, n_pref = n_conv,
      rr      = round(exp(est), 2),
      rr_low  = round(exp(est - 1.96 * se), 2),
      rr_high = round(exp(est + 1.96 * se), 2)
    )
  }
  do.call(rbind, rows)
}

# --- Run the three Table 4 variants -------------------------------------------
cat("\n=== TABLE 4 SENSITIVITY ANALYSES ===\n")
res <- rbind(
  run_variant(lag = 7,                  label = "Lag 0-7"),
  run_variant(drop_years = 2020,        label = "Excluding 2020"),
  run_variant(outcome = "age_elderly",  label = "Elderly only (>=65)")
)
res$period <- factor(res$period, levels = periods)
res <- res %>% arrange(analysis, period)

cat("\nPooled RR at P95 by variant and period:\n")
print(as.data.frame(res), row.names = FALSE)
write.csv(res, "output/revision/table4_sensitivity.csv", row.names = FALSE)

# --- Self-check against manuscript Table 4 ------------------------------------
# These are the values reported in the revised manuscript (code-archive audit,
# single-pipeline rerun). This script must reproduce them exactly.
cat("\nExpected (Table 4, revised manuscript):\n")
cat("  Lag 0-7:        5.18 (4.90-5.47) / 3.82 (3.54-4.13) / 3.48 (3.21-3.77)\n")
cat("  Excluding 2020: 5.06 (4.79-5.35) / 3.93 (3.61-4.28) / 3.24 (3.03-3.46)\n")
cat("  Elderly only:   4.92 (4.65-5.21) / 3.70 (3.43-3.99) / 3.32 (3.08-3.59)\n")
cat("Note: pre-pandemic and post-pandemic periods contain no 2020 data, so\n")
cat("the 'Excluding 2020' row reproduces the main-analysis values there;\n")
cat("the pandemic-period estimate uses 2021 data only (43/47 converged).\n")

# --- Fig. 5: forest plot of Table 4 -------------------------------------------
# Assembles the published figure from pipeline outputs only:
#   main row   <- output/rr_comparison_by_period.csv   (main script, Part 6.2)
#   3 variants <- res (this script, just written to CSV)
#   severe row <- output/revision/severe_trajectory_full.csv (revision_sensitivity_analyses.R, C2)
library(ggplot2)

main_csv   <- "output/rr_comparison_by_period.csv"
severe_csv <- "output/revision/severe_trajectory_full.csv"
if (file.exists(main_csv) && file.exists(severe_csv)) {
  tab2 <- read.csv(main_csv)
  main_row <- tab2[tolower(tab2$percentile) == "p95",
                   c("period", "rr", "rr_low", "rr_high")]
  if (nrow(main_row) == 0) {
    stop("No P95 rows found in ", main_csv,
         " (percentile values present: ",
         paste(unique(tab2$percentile), collapse = ", "), ")")
  }
  main_row$analysis <- "Main analysis"

  sev <- read.csv(severe_csv)
  sev_row <- sev[, c("period", "rr", "rr_low", "rr_high")]
  sev_row$analysis <- "Severe cases"

  fig5_df <- rbind(main_row, res[, c("analysis", "period", "rr", "rr_low", "rr_high")], sev_row)
  fig5_df$analysis <- factor(fig5_df$analysis,
    levels = c("Main analysis", "Lag 0-7", "Excluding 2020",
               "Elderly only (>=65)", "Severe cases"),
    labels = c("Main analysis", "Lag 0-7", "Exclude 2020",
               "Elderly only", "Severe cases"))
  fig5_df$period <- factor(fig5_df$period,
    levels = c("pre_pandemic", "pandemic", "post_pandemic"),
    labels = c("Pre-pandemic", "Pandemic", "Post-pandemic"))

  period_cols <- c("Pre-pandemic" = "#2166AC", "Pandemic" = "#B2182B",
                   "Post-pandemic" = "#4DAF4A")
  p5 <- ggplot(fig5_df, aes(x = analysis, y = rr, colour = period,
                            group = period)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(ymin = rr_low, ymax = rr_high),
                  position = position_dodge(width = 0.6), width = 0.25,
                  linewidth = 0.9) +
    geom_point(position = position_dodge(width = 0.6), size = 3) +
    scale_colour_manual(values = period_cols, name = "Period") +
    scale_x_discrete(limits = rev) +
    coord_flip() +
    labs(y = "Cumulative RR at WBGT P95 (31.9\u00b0C)", x = NULL) +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom",
          panel.grid.major.x = element_line(colour = "grey90"),
          panel.grid.minor = element_blank())
  ggsave("output/revision/fig5_sensitivity.png", p5,
         width = 2400, height = 1500, units = "px", dpi = 200)
  cat("\nFig. 5 written to output/revision/fig5_sensitivity.png\n")
} else {
  cat("\n(Fig. 5 skipped: run heatstroke_dlnm_analysis.R and",
      "revision_sensitivity_analyses.R first to generate the input CSVs.)\n")
}
