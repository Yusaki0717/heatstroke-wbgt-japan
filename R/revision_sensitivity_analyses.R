################################################################################
# Response-to-reviewers analyses (all correspond to published content)
#
# Manuscript hp-2026-00307n. Each block below is labelled with the specific
# reviewer comment it answers and the manuscript location it supports:
#   B1  -> Table S4  (AF sensitivity to reference/centring choice)
#   B2  -> main text Section 4.2 (16.3% vs 20.3% burden comparison)
#   B3  -> Fig. 7a annotation values (verification only; no separate table)
#   C2  -> Table 4 "Severe cases only" row (8.14 / 4.97 / 5.04) and the
#          severe-case attributable numbers in Table 5
#   C3  -> Table S3 (convergence diagnostic: high-WBGT-day counts for
#          Hokkaido, Aomori, Yamaguchi), cited in main text Section 3.2
#
# PREREQUISITE: run, in order,
#   1. R/heatstroke_dlnm_analysis.R   (produces dat, conv_list, knots_var,
#                                       knots_lag, wbgt_cen, wbgt_all, periods,
#                                       prefectures, pref_lookup)
#   2. R/spatial_maps_and_af.R        (produces af_df, af_severe_df -- B2, B3,
#                                       and C2b depend on these)
#   3. this script
#
################################################################################


library(dplyr)
library(tidyr)
library(dlnm)
library(gnm)
library(splines)

dir.create("output/revision", showWarnings = FALSE, recursive = TRUE)

# Period year counts (for annualization)
period_years <- c(pre_pandemic = 5, pandemic = 2, post_pandemic = 2)

# ==============================================================================
# B1. AF SENSITIVITY TO REFERENCE (CENTERING) CHOICE
#     Reviewer 2: "justify the choice of the median as the reference and show
#      how sensitive the AF is to plausible alternative references"
#
#     We recompute national AF by period under several reference values.
#     Because the WBGT-EAD relationship is monotonic, the AF magnitude is
#     largely set by the reference; this quantifies that dependence.
# ==============================================================================
cat("\n=== B1. AF SENSITIVITY TO REFERENCE CHOICE ===\n")

ref_candidates <- c(
  p10 = as.numeric(quantile(wbgt_all, 0.10)),
  p25 = as.numeric(quantile(wbgt_all, 0.25)),
  p50 = as.numeric(quantile(wbgt_all, 0.50)),   # original (27.5C)
  p75 = as.numeric(quantile(wbgt_all, 0.75))
)
cat("Reference candidates (WBGT C):\n"); print(round(ref_candidates, 1))

af_ref_list <- list()

for (ref_name in names(ref_candidates)) {
  ref_val <- as.numeric(ref_candidates[ref_name])

  for (p in periods) {
    for (i in prefectures) {
      pref_dat <- dat %>%
        filter(pref_code == i, period == p, !is.na(wbgt_max)) %>%
        arrange(date)

      if (nrow(pref_dat) < 60 || sum(pref_dat$total) < 10) next

      tryCatch({
        cb <- crossbasis(
          pref_dat$wbgt_max, lag = 5,
          argvar = list(fun = "ns", knots = knots_var),
          arglag = list(fun = "ns", knots = knots_lag)
        )
        model <- gnm(
          total ~ cb + factor(dow) + holiday,
          eliminate = stratum, family = quasipoisson(), data = pref_dat
        )
        # attributable number for WBGT above the reference value only
        an <- attrdl(
          pref_dat$wbgt_max, cb, pref_dat$total, model,
          type = "an", cen = ref_val,
          range = c(ref_val, max(pref_dat$wbgt_max))
        )
        af_ref_list <- c(af_ref_list, list(data.frame(
          ref_name = ref_name,
          ref_val  = round(ref_val, 1),
          pref_code = i,
          period = p,
          total_ead = sum(pref_dat$total),
          an_heat = as.numeric(an),
          stringsAsFactors = FALSE
        )))
      }, error = function(e) {})
    }
  }
  cat("  Reference", ref_name, "(", round(ref_val,1), "C) done\n")
}

af_ref_df <- do.call(rbind, af_ref_list)

# National AF by reference x period (summed over prefectures that converged)
af_ref_summary <- af_ref_df %>%
  group_by(ref_name, ref_val, period) %>%
  summarise(
    n_pref = n_distinct(pref_code),
    total_ead = sum(total_ead),
    total_an = round(sum(an_heat)),
    national_af = round(sum(an_heat) / sum(total_ead) * 100, 1),
    .groups = "drop"
  ) %>%
  mutate(period = factor(period, levels = periods)) %>%
  arrange(ref_val, period)

cat("\nNational AF by reference and period:\n")
print(as.data.frame(af_ref_summary))
write.csv(af_ref_summary, "output/revision/af_reference_sensitivity.csv", row.names = FALSE)

# Wide format for easy manuscript table
af_ref_wide <- af_ref_summary %>%
  select(ref_name, ref_val, period, national_af) %>%
  pivot_wider(names_from = period, values_from = national_af)
cat("\nAF (%) wide format (rows = reference, cols = period):\n")
print(as.data.frame(af_ref_wide))
write.csv(af_ref_wide, "output/revision/af_reference_sensitivity_wide.csv", row.names = FALSE)


# ==============================================================================
# B2. BURDEN RESTRICTED TO COMMONLY-CONVERGED PREFECTURES
#     Reviewer 2: "recompute the burden comparison restricted to the prefectures
#      that converged in all three periods" (pre=44, pandemic=46, post=47;
#      excluded ones are northern high-RR prefectures, so the three totals are
#      not built from the same geographic base).
# ==============================================================================
cat("\n=== B2. BURDEN RESTRICTED TO COMMON CONVERGED PREFECTURES ===\n")

converged_by_period <- lapply(periods, function(p) {
  which(sapply(conv_list[[p]], isTRUE))
})
names(converged_by_period) <- periods

for (p in periods) {
  cat("  ", p, ": ", length(converged_by_period[[p]]), " converged\n", sep = "")
}

common_pref <- Reduce(intersect, converged_by_period)
cat("Prefectures converged in ALL three periods: ", length(common_pref), "\n")
cat("Codes: ", paste(common_pref, collapse = ", "), "\n")

excluded_from_common <- setdiff(sort(unique(dat$pref_code)), common_pref)
cat("Excluded (not converged in >=1 period): ",
    paste(pref_lookup$pref_name[excluded_from_common], collapse = ", "), "\n")

# --- Recompute all-severity burden using ONLY common prefectures ---
# Uses af_df produced by incremental_analysis.R (columns: pref_code, period,
# total_ead, an_heat, ...). We restrict to common_pref for ALL periods.

if (!exists("af_df")) stop("af_df not found. Run incremental_analysis.R first.")

burden_common <- af_df %>%
  filter(pref_code %in% common_pref) %>%
  group_by(period) %>%
  summarise(
    n_pref = n_distinct(pref_code),
    total_ead = sum(total_ead),
    total_an = round(sum(an_heat)),
    af = round(sum(an_heat) / sum(total_ead) * 100, 1),
    .groups = "drop"
  ) %>%
  mutate(
    n_years = period_years[as.character(period)],
    annual_an = round(total_an / n_years)
  ) %>%
  mutate(period = factor(period, levels = periods)) %>%
  arrange(period)

cat("\n--- Burden using COMMON prefectures only (n =", length(common_pref), ") ---\n")
print(as.data.frame(burden_common))

# --- For comparison: original burden using all converged prefectures per period ---
burden_original <- af_df %>%
  group_by(period) %>%
  summarise(
    n_pref = n_distinct(pref_code),
    total_ead = sum(total_ead),
    total_an = round(sum(an_heat)),
    af = round(sum(an_heat) / sum(total_ead) * 100, 1),
    .groups = "drop"
  ) %>%
  mutate(
    n_years = period_years[as.character(period)],
    annual_an = round(total_an / n_years)
  ) %>%
  mutate(period = factor(period, levels = periods)) %>%
  arrange(period)

cat("\n--- Burden using ALL converged prefectures (original) ---\n")
print(as.data.frame(burden_original))

# --- Percent change pre -> post for both, to show robustness of the ~20% claim ---
pct_change <- function(tbl) {
  pre  <- tbl$annual_an[tbl$period == "pre_pandemic"]
  post <- tbl$annual_an[tbl$period == "post_pandemic"]
  round((post - pre) / pre * 100, 1)
}
cat("\nAnnualized AN change pre->post:\n")
cat("  Common prefectures:  ", pct_change(burden_common), "%\n")
cat("  All converged (orig):", pct_change(burden_original), "%\n")

burden_compare <- bind_rows(
  burden_common   %>% mutate(set = "common_prefectures"),
  burden_original %>% mutate(set = "all_converged")
)
write.csv(burden_compare, "output/revision/burden_common_prefectures.csv", row.names = FALSE)


# ==============================================================================
# B3. FIG 7a ANNOTATION CHECK
#     Reviewer 2: "annualized attributable counts annotated in Fig. 7a do not
#      appear to match Table 5... pre-pandemic seems to be 5-yr total / 2 rather
#      than / 5". Verify the correct annualized numbers here so the figure can
#      be corrected. (Main text values 42,672 / 35,974 / 51,322 are correct.)
# ==============================================================================
cat("\n=== B3. FIG 7a ANNOTATION CHECK ===\n")

fig7a_check <- af_df %>%
  group_by(period) %>%
  summarise(
    total_ead = sum(total_ead),
    total_an = round(sum(an_heat)),
    af = round(sum(an_heat) / sum(total_ead) * 100, 1),
    .groups = "drop"
  ) %>%
  mutate(
    n_years = period_years[as.character(period)],
    annual_an_CORRECT = round(total_an / n_years),
    # the erroneous version R2 suspected (dividing everything by 2):
    annual_an_IF_DIV2 = round(total_an / 2)
  ) %>%
  mutate(period = factor(period, levels = periods)) %>%
  arrange(period)

cat("\nCorrect vs suspected-erroneous annualized AN:\n")
print(as.data.frame(fig7a_check))
cat("\n--> Use the 'annual_an_CORRECT' column for Fig 7a annotations.\n")
cat("    Expected: pre ~42,672 | pandemic ~35,974 | post ~51,322\n")
write.csv(fig7a_check, "output/revision/fig7a_annotation_check.csv", row.names = FALSE)


# ==============================================================================
# C2. FULL 3-PERIOD SEVERE-CASE TRAJECTORY + ATTRIBUTABLE NUMBERS
#     Reviewer 2: severe RR falls 8.14 -> 4.97 (pre->pandemic, ~39%, CIs do not
#      overlap) then plateaus (4.97 -> 5.04). The manuscript compared only the
#      last two periods. Report the full trajectory and give severe-case AN
#      symmetrically with all-severity (Table 5 currently omits severe AN detail).
#
#     RR values come from the second-stage meta on severe cases. If you already
#     have severe-case meta objects, reuse them. Below we (re)fit severe-case
#     first stage + meta to extract the pooled RR at P95 for all three periods
#     with CIs, and also summarize severe-case AN from af_severe_df.
# ==============================================================================
cat("\n=== C2. FULL SEVERE-CASE TRAJECTORY ===\n")

# Ensure severe outcome exists
if (!"severe_total" %in% names(dat)) {
  dat$severe_total <- dat$sev_death + dat$sev_severe
}

# --- C2a: pooled severe-case RR at P95 for each period (with CI) ---
wbgt_p95_val <- as.numeric(quantile(wbgt_all, 0.95))

severe_coef <- list(); severe_vcov <- list(); severe_conv <- list()

for (p in periods) {
  severe_coef[[p]] <- list(); severe_vcov[[p]] <- list(); severe_conv[[p]] <- list()
  for (i in prefectures) {
    pref_dat <- dat %>%
      filter(pref_code == i, period == p, !is.na(wbgt_max)) %>%
      arrange(date)
    if (nrow(pref_dat) < 60 || sum(pref_dat$severe_total) < 5) {
      severe_conv[[p]][[i]] <- FALSE; next
    }
    tryCatch({
      cb <- crossbasis(
        pref_dat$wbgt_max, lag = 5,
        argvar = list(fun = "ns", knots = knots_var),
        arglag = list(fun = "ns", knots = knots_lag)
      )
      model <- gnm(
        severe_total ~ cb + factor(dow) + holiday,
        eliminate = stratum, family = quasipoisson(), data = pref_dat
      )
      red <- crossreduce(cb, model, cen = wbgt_cen)
      severe_coef[[p]][[i]] <- coef(red)
      severe_vcov[[p]][[i]] <- vcov(red)
      severe_conv[[p]][[i]] <- TRUE
    }, error = function(e) {
      severe_conv[[p]][[i]] <<- FALSE
    })
  }
  n_ok <- sum(sapply(severe_conv[[p]], isTRUE))
  cat("  ", p, ": severe-case converged ", n_ok, "/", length(prefectures), "\n", sep = "")
}

# Prediction basis on the same natural-spline structure used for reduction
wbgt_range_loc <- range(wbgt_all)
severe_rr_rows <- list()

for (p in periods) {
  ok_idx <- which(sapply(severe_conv[[p]], isTRUE))
  if (length(ok_idx) < 5) {
    cat("  ", p, ": too few severe-converged prefectures for meta\n", sep = ""); next
  }
  ymat <- do.call(rbind, severe_coef[[p]][ok_idx])
  Slist <- severe_vcov[[p]][ok_idx]
  mfit <- mixmeta(ymat ~ 1, S = Slist, method = "reml")

  pooled_coef <- coef(mfit)
  pooled_vcov <- vcov(mfit)

  # basis at P95 and at reference, same ns() basis as crossreduce output
  basis_at <- ns(c(wbgt_cen, wbgt_p95_val), knots = knots_var,
                 Boundary.knots = wbgt_range_loc)
  bdiff <- as.numeric(basis_at[2, ] - basis_at[1, ])
  logrr <- sum(bdiff * pooled_coef)
  se    <- sqrt(as.numeric(t(bdiff) %*% pooled_vcov %*% bdiff))

  severe_rr_rows[[p]] <- data.frame(
    period = p,
    n_pref = length(ok_idx),
    rr = round(exp(logrr), 2),
    rr_low = round(exp(logrr - 1.96 * se), 2),
    rr_high = round(exp(logrr + 1.96 * se), 2)
  )
}
severe_rr_tbl <- do.call(rbind, severe_rr_rows)
severe_rr_tbl$period <- factor(severe_rr_tbl$period, levels = periods)
severe_rr_tbl <- severe_rr_tbl %>% arrange(period)

cat("\nPooled severe-case RR at P95 (should reproduce ~8.14 / 4.97 / 5.04):\n")
print(as.data.frame(severe_rr_tbl))

# --- C2b: severe-case AN by period (from af_severe_df in incremental script) ---
if (exists("af_severe_df")) {
  severe_an_tbl <- af_severe_df %>%
    group_by(period) %>%
    summarise(
      n_pref = n_distinct(pref_code),
      total_severe = sum(total_severe),
      severe_an = round(sum(an_heat)),
      severe_af = round(sum(an_heat) / sum(total_severe) * 100, 1),
      .groups = "drop"
    ) %>%
    mutate(
      n_years = period_years[as.character(period)],
      severe_an_annual = round(severe_an / n_years)
    ) %>%
    mutate(period = factor(period, levels = periods)) %>%
    arrange(period)
  cat("\nSevere-case AN by period (for symmetric Table 5):\n")
  print(as.data.frame(severe_an_tbl))
} else {
  severe_an_tbl <- NULL
  cat("\n(af_severe_df not found; run incremental_analysis.R severe block for AN.)\n")
}

# --- C2c: pre->pandemic percent change with CI-overlap flag ---
pre  <- severe_rr_tbl[severe_rr_tbl$period == "pre_pandemic", ]
pan  <- severe_rr_tbl[severe_rr_tbl$period == "pandemic", ]
post <- severe_rr_tbl[severe_rr_tbl$period == "post_pandemic", ]
cat("\nSevere-case trajectory summary:\n")
cat(sprintf("  pre->pandemic: %.2f -> %.2f  (%.1f%% change); CIs overlap: %s\n",
            pre$rr, pan$rr, (pan$rr - pre$rr)/pre$rr*100,
            ifelse(pre$rr_low <= pan$rr_high & pan$rr_low <= pre$rr_high, "YES", "NO")))
cat(sprintf("  pandemic->post: %.2f -> %.2f  (%.1f%% change); CIs overlap: %s\n",
            pan$rr, post$rr, (post$rr - pan$rr)/pan$rr*100,
            ifelse(pan$rr_low <= post$rr_high & post$rr_low <= pan$rr_high, "YES", "NO")))

# Combine and save
severe_full <- severe_rr_tbl
if (!is.null(severe_an_tbl)) {
  severe_full <- severe_rr_tbl %>%
    left_join(severe_an_tbl %>% select(period, total_severe, severe_an,
                                       severe_af, severe_an_annual),
              by = "period")
}
write.csv(severe_full, "output/revision/severe_trajectory_full.csv", row.names = FALSE)


# ==============================================================================
# C3. CONVERGENCE DIAGNOSTIC (Reviewer 1, point 1)
#     "improved convergence likely reflects the secular increase in WBGT —
#      testing this hypothesis would strengthen the study."
#     Show, for the prefectures that failed in the pre-pandemic period
#     (Hokkaido, Aomori, Yamaguchi), how the count of high-WBGT days rose
#     across periods.
# ==============================================================================
cat("\n=== C3. CONVERGENCE DIAGNOSTIC ===\n")

noncov_pre <- which(!sapply(conv_list[["pre_pandemic"]], isTRUE))
cat("Non-converged in pre-pandemic:",
    paste(pref_lookup$pref_name[noncov_pre], collapse = ", "), "\n")

# thresholds of interest for "high" WBGT exposure
thr30 <- 30
thr_p90 <- as.numeric(quantile(wbgt_all, 0.90))

conv_diag <- dat %>%
  filter(pref_code %in% noncov_pre) %>%
  group_by(pref_code, pref_name, period) %>%
  summarise(
    n_days = sum(!is.na(wbgt_max)),
    days_over_30 = sum(wbgt_max > thr30, na.rm = TRUE),
    pct_over_30 = round(mean(wbgt_max > thr30, na.rm = TRUE) * 100, 1),
    days_over_p90 = sum(wbgt_max > thr_p90, na.rm = TRUE),
    max_wbgt = round(max(wbgt_max, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  mutate(period = factor(period, levels = periods)) %>%
  arrange(pref_code, period)

cat("\nHigh-WBGT exposure by period for pre-pandemic non-converged prefectures:\n")
cat("(threshold 30C and pooled P90 =", round(thr_p90,1), "C)\n")
print(as.data.frame(conv_diag))
write.csv(conv_diag, "output/revision/convergence_diagnostic.csv", row.names = FALSE)


# ==============================================================================
cat("\n=== ALL REVISION ANALYSES COMPLETE ===\n")
cat("Outputs in output/revision/:\n")
cat("  B1  af_reference_sensitivity.csv / _wide.csv\n")
cat("  B2  burden_common_prefectures.csv\n")
cat("  B3  fig7a_annotation_check.csv\n")
cat("  C2  severe_trajectory_full.csv\n")
cat("  C3  convergence_diagnostic.csv\n")
