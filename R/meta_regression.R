################################################################################
# Meta-regression (manuscript Table 3)
#
# Univariate random-effects meta-regressions of the pooled prefecture-level
# reduced coefficients on three prefecture covariates (latitude of the capital
# city, % population aged >=65, log population density), with multivariate
# Wald tests on the covariate coefficients.
#
# Prerequisite: run heatstroke_dlnm_analysis.R first
# (needs periods, coef_list, vcov_list, conv_list, meta_results in workspace).
#
# Note: mixmeta returns multivariate coefficients; p-values cannot be
# extracted as in lm. The multivariate Wald statistic is computed as
# chi2 = beta' V^-1 beta over the outcome-dimension blocks.
################################################################################

library(mixmeta)

cat("\n=== META-REGRESSION ===\n\n")

# Prefecture-level covariates (capital-city latitude; % aged >=65; population density)
pref_meta <- data.frame(
  pref_code = 1:47,
  latitude = c(
    43.06, 40.82, 39.70, 38.27, 39.72, 38.24, 37.75,
    36.34, 36.57, 36.39, 35.86, 35.61, 35.69, 35.45,
    37.90, 36.70, 36.59, 36.07, 35.66, 36.65,
    35.39, 34.98, 35.18, 34.73, 35.00, 35.01, 34.69,
    34.69, 34.69, 34.23, 35.50, 35.47, 34.66, 34.40,
    34.19, 34.07, 34.34, 33.84, 33.56, 33.59, 33.25,
    32.75, 32.79, 33.24, 31.91, 31.56, 26.21
  ),
  pct_elderly = c(
    32.2, 34.0, 33.5, 28.3, 37.5, 34.6, 32.6,
    29.8, 28.9, 30.2, 26.8, 28.0, 23.1, 25.6,
    33.0, 33.3, 30.4, 31.5, 30.5, 33.0,
    30.6, 30.3, 25.3, 30.6, 26.5, 29.2, 27.6,
    29.9, 32.0, 33.1, 33.0, 35.0, 29.5, 29.0,
    34.8, 34.0, 32.6, 33.6, 35.5, 28.0, 31.5,
    32.5, 31.5, 34.0, 33.5, 32.5, 22.2
  ),
  pop_density = c(
    66, 130, 80, 315, 83, 117, 133,
    470, 302, 307, 1934, 1218, 6410, 3823,
    178, 245, 271, 186, 183, 152,
    189, 474, 1460, 310, 341, 560, 4631,
    650, 360, 196, 159, 101, 267, 332,
    225, 175, 514, 236, 100, 1024, 333,
    325, 236, 181, 139, 178, 637
  )
)
pref_meta$log_density <- log(pref_meta$pop_density)

# Wald tests on meta-regression coefficients
# (mixmeta returns multivariate coefficients; extract the slope blocks)

mr_results <- list()

for (p in periods) {
  if (is.null(meta_results[[p]])) next
  
  ok_idx <- which(sapply(conv_list[[p]], isTRUE))
  ymat <- do.call(rbind, coef_list[[p]][ok_idx])
  Slist <- vcov_list[[p]][ok_idx]
  
  meta_sub <- pref_meta[pref_meta$pref_code %in% ok_idx, ]
  meta_sub <- meta_sub[match(ok_idx, meta_sub$pref_code), ]
  
  if (nrow(meta_sub) != nrow(ymat)) {
    cat("Dimension mismatch for", p, "- skipping\n")
    next
  }
  
  cat("--- Period:", p, "(n =", nrow(ymat), "prefectures) ---\n")
  
  # Univariate meta-regression, one covariate at a time
  vars <- list(
    latitude = meta_sub$latitude,
    pct_elderly = meta_sub$pct_elderly,
    log_density = meta_sub$log_density
  )
  
  for (vname in names(vars)) {
    tryCatch({
      x <- vars[[vname]]
      mr <- mixmeta(ymat ~ x, S = Slist, method = "reml")
      
      # Extract coefficients: for a multivariate outcome, coef() may
      # return a matrix (rows = predictors, columns = outcome dimensions)
      cf <- coef(mr)
      vc <- vcov(mr)
      
      # cf may be a vector or a matrix, depending on outcome dimension
      if (is.matrix(cf)) {
        # Second row (effect of x), all outcome columns
        beta_x <- cf[2, ]
        # Corresponding variances live in the diagonal block of vcov;
        # vcov is (n_pred * n_out) x (n_pred * n_out)
        n_out <- ncol(ymat)
        # the x coefficients occupy indices n_out+1 .. 2*n_out
        idx_x <- (n_out + 1):(2 * n_out)
        vcov_x <- vc[idx_x, idx_x]
      } else {
        # Single-outcome case
        n_pred <- length(cf) / ncol(ymat)
        beta_x <- cf[(ncol(ymat) + 1):(2 * ncol(ymat))]
        idx_x <- (ncol(ymat) + 1):(2 * ncol(ymat))
        vcov_x <- vc[idx_x, idx_x]
      }
      
      # Multivariate Wald test: chi2 = beta' %*% V^-1 %*% beta
      wald_stat <- as.numeric(t(beta_x) %*% solve(vcov_x) %*% beta_x)
      wald_df <- length(beta_x)
      wald_p <- 1 - pchisq(wald_stat, df = wald_df)
      
      cat(sprintf("  %s: Wald chi2(%d) = %.2f, p = %.4f\n",
                  vname, wald_df, wald_stat, wald_p))
      
      mr_results <- c(mr_results, list(data.frame(
        period = p, variable = vname,
        wald_chi2 = round(wald_stat, 2),
        df = wald_df,
        p_value = round(wald_p, 4)
      )))
      
    }, error = function(e) {
      cat(sprintf("  %s: ERROR - %s\n", vname, e$message))
    })
  }

  cat("\n")
}

# Save results
if (length(mr_results) > 0) {
  mr_table <- do.call(rbind, mr_results)
  cat("Summary of univariate meta-regression:\n")
  print(mr_table)
  write.csv(mr_table, "output/metaregression_results.csv", row.names = FALSE)
}

# --- Self-check against manuscript Table 3 ------------------------------------
# These are the Wald statistics reported in the revised manuscript (code-archive
# audit, single-pipeline rerun). This script must reproduce them.
cat("\nExpected (Table 3, revised manuscript):\n")
cat("  Latitude:          66.75 (<0.001) / 92.71 (<0.001) / 30.80 (<0.001)\n")
cat("  % Population >=65:  2.97 (0.396)   /  2.01 (0.571)   / 10.05 (0.018)\n")
cat("  Log pop. density:  10.50 (0.015)   /  7.67 (0.053)   /  4.88 (0.181)\n")

cat("\n=== META-REGRESSION COMPLETE ===\n")
