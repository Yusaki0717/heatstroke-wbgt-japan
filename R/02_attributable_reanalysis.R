# =============================================================================
# 02_attributable_reanalysis.R -- attributable number / fraction (Section 2.3.5)
#
# Uses the attrdl() function of Gasparrini & Leone (2014) with lagged exposure
# histories built separately within each prefecture-summer, so no lag window
# crosses an off-season boundary. Monte Carlo 95% CIs (1,000 draws) for the
# national all-severity totals at the median reference.
#
# Produces (all under output/):
#   attributable_prefecture_all_grouped.csv      per prefecture-period, 4 references
#   attributable_prefecture_severe_grouped.csv   per prefecture-period, median reference
#   attributable_summary_grouped.csv             Table 5, Table S4, Section 3.8
#   attributable_pre_post_change.csv             18.1% / 26.1% comparison (Section 3.8)
# Prerequisite: none beyond 00_setup.R (first-stage fits are redone here).
# =============================================================================

source("R/00_setup.R")
load_attrdl()
start_log("02_attributable_reanalysis")

NSIM <- 1000
set.seed(20260825)
n_summers <- c(pre_pandemic = 5, pandemic = 2, post_pandemic = 2)
REFS <- as.numeric(quantile(wbgt_all, c(0.10, 0.25, 0.50, 0.75)))
names(REFS) <- c("P10", "P25", "P50", "P75")

grouped_lag_matrix <- function(x, group, lag_days) {
  ans <- matrix(NA_real_, nrow = length(x), ncol = lag_days + 1)
  for (g in unique(group)) {
    ii <- which(group == g)
    z <- x[ii]
    for (ell in 0:lag_days) {
      ans[ii, ell + 1] <- if (ell == 0) z else c(rep(NA_real_, ell), head(z, -ell))
    }
  }
  ans
}

fit_attr_outcome <- function(outcome, min_events, simulate = FALSE,
                             refs_use = names(REFS)) {
  detail <- list(); sim_by_period <- setNames(vector("list", length(periods)), periods)
  for (p in periods) {
    sim_by_period[[p]] <- if (simulate) rep(0, NSIM) else NULL
    for (i in prefectures) {
      d <- dat %>% filter(pref_code == i, period == p, !is.na(wbgt_max)) %>% arrange(date)
      if (nrow(d) < 60 || sum(d[[outcome]]) < min_events) next
      y <- d[[outcome]]
      ans <- tryCatch({
        cb <- crossbasis(
          d$wbgt_max, lag = 5,
          argvar = list(fun = "ns", knots = knots_var, Boundary.knots = WBGT_RANGE),
          arglag = list(fun = "ns", knots = knots_lag),
          group = factor(d$year)
        )
        fml <- reformulate(c("cb", "factor(dow)", "holiday"), response = outcome)
        mo <- gnm(fml, eliminate = d$stratum, family = quasipoisson(), data = d)
        lagmat <- grouped_lag_matrix(d$wbgt_max, d$year, 5)
        an_ref <- setNames(rep(NA_real_, length(REFS)), names(REFS))
        an_ref[refs_use] <- vapply(REFS[refs_use], function(ref) {
          attrdl(lagmat, cb, y, mo, type = "an", cen = ref,
                 range = c(ref, max(d$wbgt_max)))
        }, numeric(1))
        sims <- if (simulate) {
          attrdl(lagmat, cb, y, mo, type = "an", cen = REFS[["P50"]],
                 range = c(REFS[["P50"]], max(d$wbgt_max)),
                 sim = TRUE, nsim = NSIM)
        } else NULL
        list(an = an_ref, sims = sims, total = sum(y))
      }, error = function(e) {
        message(sprintf("FAILED AF %s / %s / pref %d: %s", outcome, p, i, e$message))
        NULL
      })
      if (is.null(ans)) next
      detail[[length(detail) + 1]] <- data.frame(
        outcome = outcome, period = p, pref_code = i, total_ead = ans$total,
        an_p10 = ans$an["P10"], an_p25 = ans$an["P25"],
        an_p50 = ans$an["P50"], an_p75 = ans$an["P75"])
      if (simulate) sim_by_period[[p]] <- sim_by_period[[p]] + ans$sims
    }
    cat(sprintf("AF_FIT|%s|%s|prefectures=%d\n", outcome, p,
                sum(vapply(detail, function(x) x$period[1] == p, logical(1)))))
  }
  list(detail = do.call(rbind, detail), sims = sim_by_period)
}

ALL <- fit_attr_outcome("total", 10, simulate = TRUE)
# The manuscript reports severe-case AF only at the median reference.  Do not
# extrapolate the sparse severe model to unrequested low reference values.
SEV <- fit_attr_outcome("severe_total", 5, simulate = FALSE, refs_use = "P50")

write_csv(ALL$detail, "attributable_prefecture_all_grouped.csv")
write_csv(SEV$detail, "attributable_prefecture_severe_grouped.csv")

aggregate_attr <- function(obj, outcome) {
  out <- list()
  for (p in periods) {
    z <- obj$detail[obj$detail$period == p, ]
    total <- sum(z$total_ead)
    for (rn in names(REFS)) {
      if (all(is.na(z[[paste0("an_", tolower(rn))]]))) next
      an <- sum(z[[paste0("an_", tolower(rn))]])
      row <- data.frame(outcome = outcome, period = p, reference = rn,
                        reference_wbgt = REFS[[rn]], total_ead = total,
                        an = an, af = 100 * an / total,
                        summers = n_summers[[p]], an_annual = an / n_summers[[p]],
                        ead_annual = total / n_summers[[p]])
      if (rn == "P50" && !is.null(obj$sims[[p]])) {
        ci <- quantile(obj$sims[[p]], c(0.025, 0.975))
        row$an_lo <- ci[1]; row$an_hi <- ci[2]
        row$af_lo <- 100 * ci[1] / total; row$af_hi <- 100 * ci[2] / total
      } else {
        row$an_lo <- NA_real_; row$an_hi <- NA_real_
        row$af_lo <- NA_real_; row$af_hi <- NA_real_
      }
      out[[length(out) + 1]] <- row
    }
  }
  do.call(rbind, out)
}

res <- rbind(aggregate_attr(ALL, "all_severity"),
             aggregate_attr(SEV, "severe"))
rownames(res) <- NULL
write_csv(res, "attributable_summary_grouped.csv")
print(res)

m <- res[res$outcome == "all_severity" & res$reference == "P50", ]
pre <- m[m$period == "pre_pandemic", ]
post <- m[m$period == "post_pandemic", ]
summary_change <- data.frame(
  metric = c("annualized_total_ead", "annualized_attributable_number", "AF_percentage_points"),
  value = c(100 * (post$ead_annual / pre$ead_annual - 1),
            100 * (post$an_annual / pre$an_annual - 1),
            post$af - pre$af))
write_csv(summary_change, "attributable_pre_post_change.csv")
print(summary_change)

stopifnot(nrow(ALL$detail) == 141L, nrow(SEV$detail) == 141L,
          all(is.finite(res$an)), all(is.finite(res$af)))
end_log("02_attributable_reanalysis")
