# =============================================================================
# C2d. FORMAL PERIOD x EXPOSURE INTERACTION TEST (severe cases)
#
# Motivation: the "plateau" claim (pandemic RR 4.97 vs post-pandemic 5.04)
# currently rests on visual CI overlap. Overlapping CIs are NOT a formal test
# of no difference. This module tests, in a single multivariate meta-regression,
# whether the reduced severe-case exposure-response coefficients differ by
# period — i.e. a Wald test on the period terms — using the SAME machinery as
# the Table 3 meta-regression (mixmeta + Wald).
#
# PREREQUISITE: run the main script and then revision_sensitivity_analyses.R
# first, so that these objects exist in the workspace:
#                 severe_coef, severe_vcov, severe_conv, periods, prefectures,
#                 knots_var, knots_lag, wbgt_cen, wbgt_all
# No new fitting of first-stage models is done here; we reuse the reduced
# coefficients C2 already produced.
# =============================================================================

suppressMessages({ library(mixmeta) })

cat("\n=== C2d. Period x exposure interaction test (severe cases) ===\n")

# ---- 1. Assemble ALL converged severe-case reduced coefficients across periods
# Each prefecture-period contributes one row: a length-3 reduced coefficient
# vector (ns with 2 internal knots -> df = 3) plus its 3x3 (co)variance.
ylist <- list()      # coefficient rows
Slist_all <- list()  # (co)variance matrices
period_vec <- character(0)

for (p in periods) {
  ok_idx <- which(sapply(severe_conv[[p]], isTRUE))
  for (i in ok_idx) {
    ylist[[length(ylist) + 1]]      <- severe_coef[[p]][[i]]
    Slist_all[[length(Slist_all) + 1]] <- severe_vcov[[p]][[i]]
    period_vec <- c(period_vec, p)
  }
}

Y <- do.call(rbind, ylist)                       # (sum n_pref) x 3
period_f <- factor(period_vec, levels = periods) # pre / pandemic / post

cat("Total prefecture-period units entering the pooled model:", nrow(Y), "\n")
cat("By period: ",
    paste(names(table(period_f)), table(period_f), sep = "=", collapse = ", "),
    "\n\n")

# ---- 2. Two multivariate meta-analysis models -------------------------------
# Reduced model: common exposure-response across periods (intercept only)
m_reduced <- mixmeta(Y ~ 1,        S = Slist_all, method = "reml")
# Full model: exposure-response allowed to differ by period
m_full    <- mixmeta(Y ~ period_f, S = Slist_all, method = "fixed")
# NOTE on method: LR tests require ML/REML with identical random structure.
# Here each unit carries its own within-unit S (two-stage design), so we test
# the period terms with a multivariate Wald test on the fixed-effects model,
# which is the standard two-stage approach (Gasparrini & Armstrong 2013).
m_reduced_f <- mixmeta(Y ~ 1,        S = Slist_all, method = "fixed")

# ---- 3. Wald test on the period coefficients --------------------------------
# The "period" effect spans 2 dummy terms x 3 basis dimensions = 6 coefficients.
# Test H0: all period-related coefficients = 0.
cf   <- coef(m_full)
vcov_full <- vcov(m_full)
# names of coefficients associated with period_f
per_terms <- grep("period_f", names(cf), value = TRUE)
cat("Period-related coefficients tested (", length(per_terms), "):\n", sep = "")
print(per_terms)

b_per <- cf[per_terms]
V_per <- vcov_full[per_terms, per_terms, drop = FALSE]
wald_stat <- as.numeric(t(b_per) %*% solve(V_per) %*% b_per)
df_wald   <- length(per_terms)
p_wald    <- pchisq(wald_stat, df = df_wald, lower.tail = FALSE)

cat(sprintf("\nOVERALL period x exposure Wald test:  chi2 = %.2f, df = %d, P = %.4g\n",
            wald_stat, df_wald, p_wald))
cat("  H0: severe-case exposure-response is identical across all three periods.\n")
if (p_wald < 0.05) {
  cat("  -> Reject H0: the severe-case association differs across periods overall.\n")
} else {
  cat("  -> Do not reject H0.\n")
}

# ---- 4. Targeted contrast: pandemic vs post-pandemic (the 'plateau') --------
# Re-level so pandemic is the reference; then the post-pandemic terms isolate
# the pandemic -> post difference that the plateau claim is about.
period_f2 <- relevel(period_f, ref = "pandemic")
m_full2   <- mixmeta(Y ~ period_f2, S = Slist_all, method = "fixed")
cf2 <- coef(m_full2); V2 <- vcov(m_full2)
post_terms <- grep("post_pandemic", names(cf2), value = TRUE)

b_post <- cf2[post_terms]
V_post <- V2[post_terms, post_terms, drop = FALSE]
wald_post <- as.numeric(t(b_post) %*% solve(V_post) %*% b_post)
df_post   <- length(post_terms)
p_post    <- pchisq(wald_post, df = df_post, lower.tail = FALSE)

cat(sprintf("\nPANDEMIC vs POST-PANDEMIC contrast (the 'plateau'):\n"))
cat(sprintf("  Wald chi2 = %.2f, df = %d, P = %.4g\n", wald_post, df_post, p_post))
if (p_post < 0.05) {
  cat("  -> Severe-case association DID change from pandemic to post-pandemic\n")
  cat("     (plateau interpretation would need qualifying).\n")
} else {
  cat("  -> No significant pandemic->post change: consistent with a plateau\n")
  cat("     (now supported by a formal test, not only overlapping CIs).\n")
}

# ---- 5. Targeted contrast: pre vs pandemic (the initial drop) ---------------
pre_terms <- grep("pre_pandemic", names(cf2), value = TRUE)  # ref = pandemic
b_pre <- cf2[pre_terms]
V_pre <- V2[pre_terms, pre_terms, drop = FALSE]
wald_pre <- as.numeric(t(b_pre) %*% solve(V_pre) %*% b_pre)
df_pre   <- length(pre_terms)
p_pre    <- pchisq(wald_pre, df = df_pre, lower.tail = FALSE)
cat(sprintf("\nPRE vs PANDEMIC contrast (the initial drop):\n"))
cat(sprintf("  Wald chi2 = %.2f, df = %d, P = %.4g\n", wald_pre, df_pre, p_pre))
if (p_pre < 0.05) {
  cat("  -> Severe-case association dropped significantly pre -> pandemic\n")
  cat("     (consistent with the 8.14 -> 4.97 decline).\n")
}

# ---- 6. Save a compact summary ----------------------------------------------
interaction_summary <- data.frame(
  contrast = c("overall (3 periods)", "pre vs pandemic", "pandemic vs post"),
  chi2     = round(c(wald_stat, wald_pre, wald_post), 2),
  df       = c(df_wald, df_pre, df_post),
  p_value  = signif(c(p_wald, p_pre, p_post), 3)
)
cat("\n--- Interaction test summary ---\n")
print(interaction_summary, row.names = FALSE)
dir.create("output/revision", showWarnings = FALSE, recursive = TRUE)
write.csv(interaction_summary, "output/revision/severe_interaction_test.csv",
          row.names = FALSE)
cat("\nSaved -> output/revision/severe_interaction_test.csv\n")

# --- Self-check against manuscript Section 3.7 --------------------------------
cat("\nExpected (Section 3.7, revised manuscript):\n")
cat("  overall chi2(6) = 47.22 (p < 0.001);",
    "pre->pandemic chi2(3) = 26.97 (p < 0.001);",
    "pandemic->post chi2(3) = 7.03 (p = 0.071)\n")
