################################################################################
# PART 7A: SENSITIVITY ANALYSIS — PERIOD-SPECIFIC vs POOLED KNOTS
# Response to Reviewer 2, Comment 2.2 ("moving basis" concern)
# SELF-CONTAINED / ROBUST VERSION
# --------------------------------------------------------------------------
# HOW TO RUN
#   Run the main script (Parts 0-6) first so that `dat` and the packages
#   (dplyr, dlnm, gnm, mixmeta, splines) are in memory, THEN:
#       source("sensitivity_knots_S5.R")
#
# WHY THIS VERSION
#   The previous version relied on the objects `periods`, `dat$period` and
#   `meta_results`. If any of those has been re-labelled in the session (e.g.
#   `periods` reset to plotting labels), the period subsets come back empty and
#   every fit fails. This version depends ONLY on `dat` and defines the three
#   periods directly from `dat$year`, and it recomputes BOTH the pooled-knot
#   (main) and period-specific-knot fits from scratch with the identical
#   pipeline (conditional quasi-Poisson via gnm, year x month strata, lag 0-5,
#   ns argvar/arglag, REML meta, same fixed reference 27.5 C and read-point).
#
# SELF-CHECK
#   The RR_pooled_knots column should reproduce your Table 2 values 5.06 / 3.71
#   / 3.24. If it does, RR_period_knots is the sensitivity result.
#
# OUTPUT
#   output/tableS5_knot_sensitivity.csv        <- primary (fill into Table S5)
#   output/tableS5B_fully_period_specific.csv  <- optional secondary check
################################################################################

stopifnot(exists("dat"))
dir.create("output", showWarnings = FALSE)

cat("\n=== SENSITIVITY: period-specific vs pooled knots (Reviewer 2.2) ===\n")

## ---- 0. Recompute all quantities from `dat` (robust to session state) --------
wbgt_all_  <- dat$wbgt_max[!is.na(dat$wbgt_max)]
knots_pool <- quantile(wbgt_all_, c(0.50, 0.90))   # pooled whole-study knots (main)
wbgt_cen_  <- as.numeric(median(wbgt_all_))         # fixed centring (27.5)
wbgt_rng_  <- range(wbgt_all_)
wbgt_p95_  <- as.numeric(quantile(wbgt_all_, 0.95)) # read-point (31.9)
prefs_     <- sort(unique(dat$pref_code))
klag_      <- logknots(5, 2)

# three periods defined directly from year (label -> years), in order
period_defs <- list(
  "Pre-pandemic (2015-2019)"  = 2015:2019,
  "Pandemic (2020-2021)"      = 2020:2021,
  "Post-pandemic (2022-2023)" = 2022:2023
)

cat("Pooled knots (main):", round(as.numeric(knots_pool), 1),
    "| centring:", round(wbgt_cen_, 1), "| read-point (P95):", round(wbgt_p95_, 1), "\n")

## ---- 1. Helper: fit first + second stage for a set of years & knots ---------
fit_stage12 <- function(yrs, knots_used, cen_used) {
  cl <- list(); vl <- list()
  for (i in prefs_) {
    pd <- dat %>%
      filter(pref_code == i, year %in% yrs, !is.na(wbgt_max)) %>%
      arrange(date)
    if (nrow(pd) < 60 || sum(pd$total) < 10) next
    res <- tryCatch(
      suppressWarnings({
        cb <- crossbasis(pd$wbgt_max, lag = 5,
                         argvar = list(fun = "ns", knots = knots_used),
                         arglag = list(fun = "ns", knots = klag_))
        m  <- gnm(total ~ cb + factor(dow) + holiday,
                  eliminate = stratum, family = quasipoisson(), data = pd)
        rd <- crossreduce(cb, m, cen = cen_used)
        list(coef = coef(rd), vcov = vcov(rd))
      }),
      error = function(e) NULL
    )
    if (!is.null(res)) { cl[[as.character(i)]] <- res$coef; vl[[as.character(i)]] <- res$vcov }
  }
  fit <- if (length(cl) >= 5) mixmeta(do.call(rbind, cl) ~ 1, S = vl, method = "reml") else NULL
  list(fit = fit, n_ok = length(cl))
}

## ---- 2. Helper: pooled cumulative RR at the 95th percentile (31.9 C) ---------
extract_rr_p95 <- function(meta_fit, knots_used) {
  if (is.null(meta_fit)) return(c(rr = NA, low = NA, high = NA))
  wp   <- seq(wbgt_rng_[1], wbgt_rng_[2], length.out = 100)
  bas  <- ns(wp, knots = knots_used, Boundary.knots = wbgt_rng_)
  i_x  <- which.min(abs(wp - wbgt_p95_))
  i_r  <- which.min(abs(wp - wbgt_cen_))
  bdif <- as.numeric(bas[i_x, ] - bas[i_r, ])
  cf   <- coef(meta_fit); vc <- vcov(meta_fit)
  logrr <- sum(bdif * cf)
  se    <- sqrt(as.numeric(t(bdif) %*% vc %*% bdif))
  c(rr = exp(logrr), low = exp(logrr - 1.96 * se), high = exp(logrr + 1.96 * se))
}
fmt <- function(v) if (is.na(v["rr"])) "NA" else
  sprintf("%.2f (%.2f-%.2f)", v["rr"], v["low"], v["high"])

## ---- 3. Primary sensitivity: pooled vs period-specific knots -----------------
rows <- lapply(names(period_defs), function(lbl) {
  yrs   <- period_defs[[lbl]]
  wper  <- dat$wbgt_max[dat$year %in% yrs & !is.na(dat$wbgt_max)]
  kper  <- quantile(wper, c(0.50, 0.90))              # period-specific knots

  main  <- fit_stage12(yrs, knots_pool, wbgt_cen_)     # pooled knots (= main)
  sens  <- fit_stage12(yrs, kper,       wbgt_cen_)     # period-specific knots

  cat(sprintf("  [%s] converged: pooled-knot %d/47, period-knot %d/47 | period knots = %s / %s\n",
              lbl, main$n_ok, sens$n_ok,
              round(kper[1], 1), round(kper[2], 1)))

  main_rr <- extract_rr_p95(main$fit, knots_pool)
  sens_rr <- extract_rr_p95(sens$fit, kper)
  pdiff   <- if (is.na(sens_rr["rr"]) || is.na(main_rr["rr"])) NA_character_ else
                 sprintf("%+.1f%%", (sens_rr["rr"] / main_rr["rr"] - 1) * 100)
  data.frame(
    period          = lbl,
    knots_pooled_C  = paste0(round(knots_pool[1], 1), " / ", round(knots_pool[2], 1)),
    RR_pooled_knots = fmt(main_rr),
    knots_period_C  = paste0(round(kper[1], 1), " / ", round(kper[2], 1)),
    RR_period_knots = fmt(sens_rr),
    pct_diff_RR     = pdiff,
    stringsAsFactors = FALSE
  )
})
tableS5 <- do.call(rbind, rows)

cat("\n----- Table S5: cumulative RR at the 95th percentile (31.9 C),",
    "pooled vs period-specific knots -----\n")
print(tableS5, row.names = FALSE)
write.csv(tableS5, "output/tableS5_knot_sensitivity.csv", row.names = FALSE)
cat("\nSaved -> output/tableS5_knot_sensitivity.csv\n")
cat("(Self-check: RR_pooled_knots should read 5.06 / 3.71 / 3.24.)\n")

## ---- 4. (OPTIONAL) Fully period-specific basis (knots + centring + read-pt) --
rowsB <- lapply(names(period_defs), function(lbl) {
  yrs   <- period_defs[[lbl]]
  wper  <- dat$wbgt_max[dat$year %in% yrs & !is.na(dat$wbgt_max)]
  cen_p <- as.numeric(median(wper)); p95_p <- as.numeric(quantile(wper, 0.95))
  rng_p <- range(wper); kper <- quantile(wper, c(0.50, 0.90))

  f <- fit_stage12(yrs, kper, cen_p)     # period-specific knots AND centring
  if (is.null(f$fit)) {
    return(data.frame(period = lbl, ref_median_C = round(cen_p, 1),
                      read_p95_C = round(p95_p, 1), RR_fully_period = NA,
                      stringsAsFactors = FALSE))
  }
  wp   <- seq(rng_p[1], rng_p[2], length.out = 100)
  bas  <- ns(wp, knots = kper, Boundary.knots = rng_p)
  i_x  <- which.min(abs(wp - p95_p)); i_r <- which.min(abs(wp - cen_p))
  bdif <- as.numeric(bas[i_x, ] - bas[i_r, ])
  logrr <- sum(bdif * coef(f$fit))
  se    <- sqrt(as.numeric(t(bdif) %*% vcov(f$fit) %*% bdif))
  data.frame(
    period          = lbl,
    ref_median_C    = round(cen_p, 1),
    read_p95_C      = round(p95_p, 1),
    RR_fully_period = sprintf("%.2f (%.2f-%.2f)",
                              exp(logrr), exp(logrr - 1.96 * se), exp(logrr + 1.96 * se)),
    stringsAsFactors = FALSE
  )
})
tableS5B <- do.call(rbind, rowsB)

cat("\n----- (Optional) Fully period-specific basis:",
    "RR at each period's own 95th pct vs own median -----\n")
print(tableS5B, row.names = FALSE)
write.csv(tableS5B, "output/tableS5B_fully_period_specific.csv", row.names = FALSE)
cat("\nSaved -> output/tableS5B_fully_period_specific.csv\n")

cat("\n=== Knot-sensitivity analysis complete ===\n")
