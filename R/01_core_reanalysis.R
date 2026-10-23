# =============================================================================
# 01_core_reanalysis.R -- coefficient-level reanalysis (corrected pipeline)
#
# Corrections relative to the July-2026 public release (manuscript Section
# 2.3.1 and Response to Reviewers): boundary knots fixed to the pooled WBGT
# range, lag histories reset at every summer, exact WBGT read-points.
#
# Produces (all under output/):
#   table2_grouped_exact.csv            Table 2
#   period_contrasts_independent.csv    Table 2 footnote / Section 3.2 contrasts
#   heterogeneity_grouped.csv           I2, Cochran Q, trace of Psi (Table 2 footnote)
#   period_moderator_stacked.csv        stacked period-moderator model (Section 2.3.2)
#   severity_pooled_grouped.csv         Table 4 severe / non-severe rows
#   severity_interaction_fixed.csv      period x severity interaction (Section 3.7)
#   severity_period_tests_fixed.csv     severe-only period tests (Section 3.7)
#   metaregression_grouped.csv          Table 3 (+ BH-adjusted P, variance explained)
#   prefecture_blup_p95_grouped.csv     Table S1
#   prefecture_change_grouped.csv       Table S2
#   p95_support.csv, subset_sensitivity_grouped.csv   Section 3.6 / Table S6 sensitivity
#   table4_grouped_exact.csv            Table 4
#   period_specific_knots_grouped.csv   Table S5
#   within_season_grouped_exact.csv     within-season depletion (Section 4.1)
#   pooled_national_death_grouped.csv   pooled national death model (Section 2.3.6)
#   tokyo_cumulative_p95_grouped.csv    Fig. 2 caption values
#   grouped_stage1_*.rds, sessionInfo.txt
# =============================================================================

source("R/00_setup.R")
start_log("01_core_reanalysis")

cat(sprintf("WBGT|median=%.9f|P75=%.9f|P90=%.9f|P95=%.9f|range=%.9f,%.9f\n",
            wbgt_cen, Q_REPORT["P75"], Q_REPORT["P90"], Q_REPORT["P95"],
            WBGT_RANGE[1], WBGT_RANGE[2]))

# Daily within-summer covariates are defined before any period/prefecture subset.
# "High-WBGT day" for the depletion analysis = WBGTmax above the pooled P95 (31.9 C).
dat <- dat %>% arrange(pref_code, year, date) %>% group_by(pref_code, year) %>%
  mutate(excess = pmax(0, wbgt_max - wbgt_cen),
         cumheat = dplyr::lag(cumsum(ifelse(is.na(excess), 0, excess)), default = 0),
         hot = as.integer(wbgt_max > Q_REPORT[["P95"]]),
         cumhot = dplyr::lag(cumsum(ifelse(is.na(hot), 0, hot)), default = 0)) %>%
  ungroup()

cat("\n=== MAIN AND SEVERITY FIRST STAGE ===\n")
MAIN <- fit_stage1_grouped("total", 10, label = "main_lag0_5")
SEV  <- fit_stage1_grouped("severe_total", 5, label = "severe_lag0_5")
NON  <- fit_stage1_grouped("nonsevere_total", 5, label = "nonsevere_lag0_5")
saveRDS(list(MAIN = MAIN, SEV = SEV, NON = NON),
        file.path(OUT_DIR, "grouped_stage1_main_severity.rds"))

# With the common basis the first stage is estimable for all 47 prefectures in
# every period (Table 1, convergence row). Stop if that is not reproduced.
stopifnot(all(vapply(periods, function(p) length(converged_ids(MAIN, p)), integer(1)) == 47L))
stopifnot(all(vapply(periods, function(p) length(converged_ids(SEV, p)), integer(1)) == 47L))
stopifnot(all(vapply(periods, function(p) length(converged_ids(NON, p)), integer(1)) == 47L))

cat("\n=== TABLE 2 ===\n")
t2 <- list(); pools_main <- list()
for (qname in names(Q_REPORT)) for (p in periods) {
  po <- pool_grouped(MAIN, p, at = Q_REPORT[[qname]])
  if (qname == "P95") pools_main[[p]] <- po
  t2[[length(t2) + 1]] <- data.frame(
    percentile = qname, wbgt = Q_REPORT[[qname]], period = p, n_pref = po$n,
    rr = po$rr["rr"], lo = po$rr["lo"], hi = po$rr["hi"])
  cat(sprintf("TABLE2|%s|%s|rr=%.9f|lo=%.9f|hi=%.9f\n",
              qname, p, po$rr["rr"], po$rr["lo"], po$rr["hi"]))
}
t2 <- do.call(rbind, t2); rownames(t2) <- NULL
write_csv(t2, "table2_grouped_exact.csv")

cat("\n=== INDEPENDENT PERIOD CONTRASTS ===\n")
# The three periods comprise disjoint years, so the contrast variance is the
# sum of the two period variances (Section 2.3.2).
pc <- list()
pair_defs <- list(
  c("pandemic", "pre_pandemic"),
  c("post_pandemic", "pre_pandemic"),
  c("post_pandemic", "pandemic")
)
for (z in pair_defs) {
  a <- pools_main[[z[1]]]$rr; b <- pools_main[[z[2]]]$rr
  d <- a["log_rr"] - b["log_rr"]
  se <- sqrt(a["se"]^2 + b["se"]^2)
  pv <- 2 * pnorm(-abs(d / se))
  row <- data.frame(contrast = paste(z[1], "vs", z[2]), ratio = exp(d),
                    lo = exp(d - 1.96 * se), hi = exp(d + 1.96 * se),
                    decline_pct = 100 * (1 - exp(d)), z = d / se, p = pv)
  pc[[length(pc) + 1]] <- row
  cat(sprintf("PAIR|%s|ratio=%.9f|lo=%.9f|hi=%.9f|decline=%.6f|p=%.12g\n",
              row$contrast, row$ratio, row$lo, row$hi, row$decline_pct, row$p))
}
pc <- do.call(rbind, pc); write_csv(pc, "period_contrasts_independent.csv")

main_rr <- setNames(vapply(periods, function(p) pools_main[[p]]$rr["rr"], numeric(1)), periods)
excess_decline <- 100 * (1 - (main_rr["post_pandemic"] - 1) /
                                  (main_rr["pre_pandemic"] - 1))
cat(sprintf("EXCESS_RISK_DECLINE|%.9f\n", excess_decline))

cat("\n=== HETEROGENEITY ===\n")
het <- do.call(rbind, lapply(periods, function(p) {
  po <- pools_main[[p]]
  data.frame(period = p, Q = po$Q, df = po$df, p = po$pQ,
             I2 = po$I2, tau2_trace = po$tau2)
}))
print(het); write_csv(het, "heterogeneity_grouped.csv")

cat("\n=== STACKED PERIOD MODERATOR ===\n")
st <- stack_one_outcome(MAIN)
fit_st <- mixmeta(st$Y ~ period, S = st$S, data = st$M,
                  random = ~ 1 | pref_code, method = "reml")
cf <- coef(fit_st); V <- vcov(fit_st); nm <- names(cf)
ip <- grep("[.]period", nm)
joint <- wald_block(cf, V, ip)
bd95 <- basis_difference(Q_REPORT[["P95"]])
ipan <- grep("[.]periodpandemic$", nm)
ipos <- grep("[.]periodpost_pandemic$", nm)
sm <- rbind(
  data.frame(contrast = "pandemic vs pre", t(contrast_block(cf, V, ipan, bd = bd95))),
  data.frame(contrast = "post vs pre", t(contrast_block(cf, V, ipos, bd = bd95))),
  data.frame(contrast = "post vs pandemic", t(contrast_block(cf, V, ipos, ipan, bd95)))
)
sm$joint_chi2 <- joint["chi2"]; sm$joint_df <- joint["df"]; sm$joint_p <- joint["p"]
print(sm); write_csv(sm, "period_moderator_stacked.csv")

cat("\n=== SEVERITY POOLED RESULTS ===\n")
sev_rows <- list()
for (lab in c("severe", "nonsevere")) {
  ob <- if (lab == "severe") SEV else NON
  for (p in periods) {
    po <- pool_grouped(ob, p)
    sev_rows[[length(sev_rows) + 1]] <- data.frame(
      outcome = lab, period = p, n_pref = po$n,
      rr = po$rr["rr"], lo = po$rr["lo"], hi = po$rr["hi"])
  }
}
sev_rows <- do.call(rbind, sev_rows); print(sev_rows)
write_csv(sev_rows, "severity_pooled_grouped.csv")

cat("\n=== FIXED-EFFECTS PERIOD BY SEVERITY INTERACTION ===\n")
Y <- list(); S <- list(); M <- list()
for (sv in c("nonsevere", "severe")) {
  ob <- if (sv == "severe") SEV else NON
  for (p in periods) for (i in converged_ids(ob, p)) {
    Y[[length(Y) + 1]] <- ob$coef[[p]][[i]]
    S[[length(S) + 1]] <- ob$vcov[[p]][[i]]
    M[[length(M) + 1]] <- data.frame(pref_code = i, period = p, sev = sv)
  }
}
Y <- do.call(rbind, Y); M <- do.call(rbind, M)
M$period <- factor(M$period, levels = periods)
M$sev <- factor(M$sev, levels = c("nonsevere", "severe"))
fit_int <- mixmeta(Y ~ period * sev, S = S, data = M, method = "fixed")
ci <- coef(fit_int); Vi <- vcov(fit_int); ni <- names(ci)
ii <- grep("period.*:sev|sev.*:period", ni)
int_w <- wald_block(ci, Vi, ii)
int_out <- data.frame(estimator = "fixed", chi2 = int_w["chi2"],
                      df = int_w["df"], p = int_w["p"])
print(int_out); write_csv(int_out, "severity_interaction_fixed.csv")
stopifnot(int_out$p > 0.05)

cat("\n=== FIXED-EFFECTS SEVERE-ONLY PERIOD TESTS ===\n")
ys <- Y[M$sev == "severe", , drop = FALSE]
Ss <- S[M$sev == "severe"]
ms <- M[M$sev == "severe", ]
fit_so <- mixmeta(ys ~ period, S = Ss, data = ms, method = "fixed")
cs <- coef(fit_so); Vs <- vcov(fit_so); ns <- names(cs)
isp <- grep("[.]period", ns)
span <- grep("[.]periodpandemic$", ns)
spos <- grep("[.]periodpost_pandemic$", ns)
sev_tests <- rbind(
  data.frame(test = "overall",
             chi2 = unname(wald_block(cs, Vs, isp)["chi2"]),
             df = unname(wald_block(cs, Vs, isp)["df"]),
             p = unname(wald_block(cs, Vs, isp)["p"])),
  data.frame(test = "pandemic vs pre",
             chi2 = unname(contrast_block(cs, Vs, span, bd = bd95)["chi2"]),
             df = unname(contrast_block(cs, Vs, span, bd = bd95)["df"]),
             p = unname(contrast_block(cs, Vs, span, bd = bd95)["p"])),
  data.frame(test = "post vs pre",
             chi2 = unname(contrast_block(cs, Vs, spos, bd = bd95)["chi2"]),
             df = unname(contrast_block(cs, Vs, spos, bd = bd95)["df"]),
             p = unname(contrast_block(cs, Vs, spos, bd = bd95)["p"])),
  data.frame(test = "post vs pandemic",
             chi2 = unname(contrast_block(cs, Vs, spos, span, bd95)["chi2"]),
             df = unname(contrast_block(cs, Vs, spos, span, bd95)["df"]),
             p = unname(contrast_block(cs, Vs, spos, span, bd95)["p"]))
)
print(sev_tests); write_csv(sev_tests, "severity_period_tests_fixed.csv")

cat("\n=== META-REGRESSION ===\n")
pm <- read.csv(file.path(DATA_DIR, "pref_metadata.csv"))
mr <- list()
for (p in periods) {
  ids <- converged_ids(MAIN, p)
  yy <- do.call(rbind, MAIN$coef[[p]][ids]); ss <- MAIN$vcov[[p]][ids]
  md <- pm[match(ids, pm$pref_code), ]
  for (vn in c("latitude", "elderly_prop_2020", "log_pop_density")) {
    x <- md[[vn]]
    f0 <- mixmeta(yy ~ 1, S = ss, method = "reml")
    f1 <- mixmeta(yy ~ x, S = ss, method = "reml")
    cc <- coef(f1); VV <- vcov(f1); ix <- grep("[.]x$", names(cc))
    ww <- wald_block(cc, VV, ix)
    t0 <- sum(diag(as.matrix(f0$Psi))); t1 <- sum(diag(as.matrix(f1$Psi)))
    mr[[length(mr) + 1]] <- data.frame(
      period = p, variable = vn, chi2 = ww["chi2"], df = ww["df"], p = ww["p"],
      pct_explained = max(0, (t0 - t1) / t0) * 100)
  }
}
mr <- do.call(rbind, mr); mr$p_bh <- p.adjust(mr$p, method = "BH")
print(mr); write_csv(mr, "metaregression_grouped.csv")

cat("\n=== BLUP PREFECTURE RESULTS ===\n")
pr <- list()
for (p in periods) {
  po <- pools_main[[p]]; bl <- blup(po$fit, vcov = TRUE)
  for (k in seq_along(po$ids)) {
    e <- sum(bd95 * bl[[k]]$blup)
    se <- sqrt(as.numeric(t(bd95) %*% bl[[k]]$vcov %*% bd95))
    pr[[length(pr) + 1]] <- data.frame(
      pref_code = po$ids[k], pref_name = PREF_NAMES[po$ids[k]], period = p,
      rr = exp(e), lo = exp(e - 1.96 * se), hi = exp(e + 1.96 * se))
  }
}
pr <- do.call(rbind, pr); write_csv(pr, "prefecture_blup_p95_grouped.csv")
chg <- merge(pr[pr$period == "pre_pandemic", c("pref_code","pref_name","rr")],
             pr[pr$period == "post_pandemic", c("pref_code","rr")],
             by = "pref_code", suffixes = c("_pre", "_post"))
chg$pct_change <- 100 * (chg$rr_post / chg$rr_pre - 1)
write_csv(chg, "prefecture_change_grouped.csv")
cat(sprintf("BLUP|pre_range=%.9f,%.9f|post_range=%.9f,%.9f|change=%.9f,%.9f|declined=%d\n",
            min(pr$rr[pr$period == "pre_pandemic"]), max(pr$rr[pr$period == "pre_pandemic"]),
            min(pr$rr[pr$period == "post_pandemic"]), max(pr$rr[pr$period == "post_pandemic"]),
            min(chg$pct_change), max(chg$pct_change), sum(chg$pct_change < 0)))

cat("\n=== SUPPORT AND STATION-EXCLUSION SENSITIVITY ===\n")
# Yamaguchi (35): anomalous pre-2020 exposure series (Section 2.2.2, Table S6).
# Fukushima (7) and Kochi (39): observation-site changes (Section 2.2.2).
support <- do.call(rbind, lapply(periods, function(p) data.frame(
  period = p, pref_code = prefectures,
  max_wbgt = vapply(prefectures, function(i)
    max(dat$wbgt_max[dat$period == p & dat$pref_code == i], na.rm = TRUE), numeric(1)))))
support$supported <- support$max_wbgt >= Q_REPORT[["P95"]]
write_csv(support, "p95_support.csv")
sets <- list(
  all_47 = prefectures,
  exclude_yamaguchi = setdiff(prefectures, 35),
  exclude_station_changes = setdiff(prefectures, c(7, 39)),
  exclude_yamaguchi_and_station_changes = setdiff(prefectures, c(7, 35, 39))
)
sr <- list()
for (sn in names(sets)) for (p in periods) {
  po <- pool_grouped(MAIN, p, sets[[sn]])
  sr[[length(sr) + 1]] <- data.frame(set = sn, period = p, n = po$n,
                                     rr = po$rr["rr"], lo = po$rr["lo"], hi = po$rr["hi"])
}
for (p in periods) {
  ids <- support$pref_code[support$period == p & support$supported]
  po <- pool_grouped(MAIN, p, ids)
  sr[[length(sr) + 1]] <- data.frame(set = "observed_p95_only", period = p, n = po$n,
                                     rr = po$rr["rr"], lo = po$rr["lo"], hi = po$rr["hi"])
}
sr <- do.call(rbind, sr); write_csv(sr, "subset_sensitivity_grouped.csv")

cat("\n=== TABLE 4 AND PERIOD-SPECIFIC KNOTS ===\n")
LAG7 <- fit_stage1_grouped("total", 10, lag_days = 7, label = "lag0_7")
EX20 <- fit_stage1_grouped("total", 10, drop_years = 2020, label = "exclude_2020")
ELD  <- fit_stage1_grouped("age_elderly", 10, label = "elderly")
kper <- setNames(lapply(periods, function(p)
  as.numeric(quantile(dat$wbgt_max[dat$period == p], c(0.50, 0.90), na.rm = TRUE))), periods)
PK <- fit_stage1_grouped("total", 10, knots_by_period = kper, label = "period_specific_knots")

t4 <- list()
objects <- list("Main analysis (lag 0-5)" = MAIN, "Lag 0-7" = LAG7,
                "Excluding 2020" = EX20, "Elderly only (>=65)" = ELD,
                "Severe cases" = SEV, "Non-severe cases" = NON)
for (an in names(objects)) for (p in periods) {
  # The exclude-2020 pre/post fits equal MAIN by construction; both were refit
  # here to provide an independent archived check.
  po <- pool_grouped(objects[[an]], p)
  t4[[length(t4) + 1]] <- data.frame(analysis = an, period = p, n_pref = po$n,
                                     rr = po$rr["rr"], lo = po$rr["lo"], hi = po$rr["hi"])
}
t4 <- do.call(rbind, t4); write_csv(t4, "table4_grouped_exact.csv")

pk <- do.call(rbind, lapply(periods, function(p) {
  a <- pool_grouped(MAIN, p); b <- pool_grouped(PK, p)
  data.frame(period = p, pooled_knot_50 = knots_var[1], pooled_knot_90 = knots_var[2],
             pooled_rr = a$rr["rr"], period_knot_50 = kper[[p]][1],
             period_knot_90 = kper[[p]][2], period_knot_rr = b$rr["rr"],
             pct_difference = 100 * (b$rr["rr"] / a$rr["rr"] - 1))
}))
write_csv(pk, "period_specific_knots_grouped.csv")
saveRDS(list(LAG7 = LAG7, EX20 = EX20, ELD = ELD, PK = PK),
        file.path(OUT_DIR, "grouped_stage1_sensitivity.rds"))

cat("\n=== WITHIN-SEASON DEPLETION ===\n")
CH <- fit_stage1_grouped("total", 10, extra_terms = "cumheat", label = "plus_cumheat")
CD <- fit_stage1_grouped("total", 10, extra_terms = "cumhot", label = "plus_cumhot")
EA <- fit_stage1_grouped("total", 10, season_half = "early", label = "early_summer")
LA <- fit_stage1_grouped("total", 10, season_half = "late", label = "late_summer")

stack_contrasts <- function(obj, label) {
  zz <- stack_one_outcome(obj)
  ff <- mixmeta(zz$Y ~ period, S = zz$S, data = zz$M,
                random = ~ 1 | pref_code, method = "reml")
  cc <- coef(ff); VV <- vcov(ff); nn <- names(cc)
  a <- grep("[.]periodpandemic$", nn); b <- grep("[.]periodpost_pandemic$", nn)
  rbind(data.frame(model = label, contrast = "pandemic vs pre",
                   t(contrast_block(cc, VV, a, bd = bd95))),
        data.frame(model = label, contrast = "post vs pre",
                   t(contrast_block(cc, VV, b, bd = bd95))))
}
depl <- rbind(stack_contrasts(MAIN, "unadjusted"),
              stack_contrasts(CH, "+ cumulative heat"),
              stack_contrasts(CD, "+ cumulative hot-day count"),
              stack_contrasts(EA, "early summer"),
              stack_contrasts(LA, "late summer"))
print(depl); write_csv(depl, "within_season_grouped_exact.csv")

cat("\n=== POOLED NATIONAL DEATH CHECK ===\n")
# Deaths alone are too sparse for the two-stage design; a single pooled national
# conditional quasi-Poisson model (prefecture x year x month strata, lag groups
# prefecture x year) is reported for internal comparison only (Section 2.3.6).
dpool <- dat %>% filter(!is.na(wbgt_max)) %>% arrange(pref_code, date) %>%
  mutate(strat_pref = factor(paste(pref_code, year, month, sep = "-")))
death_rows <- list()
for (oc in c("death_total", "severe_total", "total")) for (p in periods) {
  dd <- dpool %>% filter(period == p)
  cb <- crossbasis(dd$wbgt_max, lag = 5,
                   argvar = list(fun = "ns", knots = knots_var, Boundary.knots = WBGT_RANGE),
                   arglag = list(fun = "ns", knots = knots_lag),
                   group = interaction(dd$pref_code, dd$year, drop = TRUE))
  fml <- reformulate(c("cb", "factor(dow)", "holiday"), response = oc)
  mo <- gnm(fml, eliminate = dd$strat_pref, family = quasipoisson(), data = dd)
  rd <- crossreduce(cb, mo, cen = wbgt_cen)
  bd <- basis_difference(Q_REPORT[["P95"]])
  e <- sum(bd * coef(rd)); se <- sqrt(as.numeric(t(bd) %*% vcov(rd) %*% bd))
  death_rows[[length(death_rows) + 1]] <- data.frame(
    outcome = oc, period = p, events = sum(dd[[oc]]), rr = exp(e),
    lo = exp(e - 1.96 * se), hi = exp(e + 1.96 * se))
}
death_rows <- do.call(rbind, death_rows); print(death_rows)
write_csv(death_rows, "pooled_national_death_grouped.csv")

cat("\n=== TOKYO CUMULATIVE RR FOR FIGURE 2 CAPTION ===\n")
tokyo <- do.call(rbind, lapply(periods, function(p) {
  cf13 <- MAIN$coef[[p]][[13]]; vc13 <- MAIN$vcov[[p]][[13]]
  bd <- basis_difference(Q_REPORT[["P95"]])
  e <- sum(bd * cf13); se <- sqrt(as.numeric(t(bd) %*% vc13 %*% bd))
  data.frame(period = p, rr = exp(e), lo = exp(e - 1.96 * se), hi = exp(e + 1.96 * se))
}))
print(tokyo); write_csv(tokyo, "tokyo_cumulative_p95_grouped.csv")

stopifnot(main_rr["pre_pandemic"] > main_rr["pandemic"],
          main_rr["pandemic"] > main_rr["post_pandemic"],
          abs(Q_REPORT[["P95"]] - 31.9) < 0.05)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
end_log("01_core_reanalysis")
