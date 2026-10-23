# =============================================================================
# 00_setup.R -- shared setup for the corrected (summer-grouped) pipeline
#
# Paper:  Huang P, Zhang J. Temporal changes in the association between
#         wet-bulb globe temperature and heat-related emergency ambulance
#         dispatches across Japan: a multi-prefecture time-series analysis
#         spanning the COVID-19 pandemic and national Heatstroke Alert
#         implementation (2015-2023). Manuscript hp-2026-00307n.
#
# This file is sourced first by every analysis script (01-05). It
#   1. loads and cleans the FDMA dispatch data and the daily WBGTmax series
#      (data are NOT distributed with this repository; see README.md),
#   2. defines the common exposure basis on the pooled whole-study WBGT
#      distribution: internal knots at the 50th/90th percentiles, boundary
#      knots at the pooled minimum/maximum, centring at the pooled median,
#      reporting percentiles P75/P90/P95 (manuscript Section 2.3.1),
#   3. defines the first-stage / second-stage helper functions.
#
# Provenance. For the public release this file merges the two private helper
# files of the archived manuscript run of 2026-08-25 (`_common.R`, which
# prepared the data exactly as Parts 1-2 of the July-2026 public script did,
# and `grouped_helpers.R`, whose functions are reproduced below unchanged).
#
# Always run from the repository root (the folder that contains run_all.R).
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(dlnm)
  library(gnm)
  library(mixmeta)
  library(splines)
})

# ---- 0. Paths (all relative to the repository root) --------------------------
if (!file.exists(file.path("R", "00_setup.R"))) {
  stop("Run from the repository root: setwd() to the folder containing run_all.R.")
}
DATA_DIR    <- "data"
OUT_DIR     <- "output"                      # everything a rerun produces
LOG_DIR     <- file.path(OUT_DIR, "logs")
FIG_OUT     <- file.path(OUT_DIR, "figures")
RESULTS_DIR <- "results"                     # archived outputs of the manuscript run
for (d in c(OUT_DIR, LOG_DIR, FIG_OUT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

write_csv <- function(x, name) {
  write.csv(x, file.path(OUT_DIR, name), row.names = FALSE, fileEncoding = "UTF-8")
}
read_out <- function(name) read.csv(file.path(OUT_DIR, name), stringsAsFactors = FALSE)

# Each script writes output/logs/<name>.log and must end with "STATUS=PASS".
start_log <- function(name) {
  while (sink.number() > 0) sink()
  sink(file.path(LOG_DIR, paste0(name, ".log")), split = TRUE)
  cat("START|", name, "|", format(Sys.time(), tz = "Asia/Tokyo"), "\n", sep = "")
  invisible(name)
}
end_log <- function(name, status = "PASS") {
  cat("END|", name, "|", format(Sys.time(), tz = "Asia/Tokyo"),
      "|STATUS=", status, "\n", sep = "")
  while (sink.number() > 0) sink()
  invisible(status)
}

# ---- 0b. attrdl() by Gasparrini & Leone (2014) ------------------------------
# attrdl() is not part of the dlnm package and is not redistributed here (no
# licence is attached to the upstream file). It is fetched once into R/attrdl.R
# (git-ignored) from the author's repository; place the file there manually if
# the machine has no internet access.
ATTRDL_URL <- paste0("https://raw.githubusercontent.com/gasparrini/",
                     "2014_gasparrini_BMCmrm_Rcodedata/master/attrdl.R")
load_attrdl <- function(path = file.path("R", "attrdl.R")) {
  if (!file.exists(path)) {
    message("R/attrdl.R not found; downloading from ", ATTRDL_URL)
    ok <- tryCatch({
      utils::download.file(ATTRDL_URL, path, mode = "wb", quiet = TRUE); TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (!ok || !file.exists(path) || file.size(path) < 1000) {
      if (file.exists(path)) unlink(path)
      stop("Could not download attrdl.R. Download it manually from\n  ", ATTRDL_URL,
           "\nand save it as R/attrdl.R, then rerun.")
    }
  }
  source(path, local = globalenv())
  if (!exists("attrdl", mode = "function")) stop("R/attrdl.R did not define attrdl().")
  invisible(TRUE)
}

# ---- 1. Data ----------------------------------------------------------------
# The nine FDMA workbooks and the daily WBGT CSV are loaded once per session;
# set options(hswbgt.reload = TRUE) to force a reload.
if (!exists(".HSWBGT_DATA_LOADED", envir = globalenv()) ||
    isTRUE(getOption("hswbgt.reload"))) {

  find_data <- function(f) {
    p <- file.path(DATA_DIR, f)
    ifelse(file.exists(p), p, f)   # ./data/ first, then the working directory
  }
  file_map <- data.frame(
    file = find_data(c(
      "heatstroke003_data_h27.xlsx",  # 2015
      "heatstroke003_data_h28.xlsx",  # 2016
      "heatstroke003_data_h29.xlsx",  # 2017
      "heatstroke003_data_h30.xlsx",  # 2018
      "heatstroke003_data_r1.xlsx",   # 2019
      "heatstroke003_data_r2.xlsx",   # 2020
      "heatstroke003_data_r3.xlsx",   # 2021
      "heatstroke003_data_r4.xlsx",   # 2022
      "heatstroke003_data_r5.xlsx"    # 2023
    )),
    year = 2015:2023, stringsAsFactors = FALSE
  )
  if (any(!file.exists(file_map$file))) {
    stop("FDMA files not found in ./data/ (see README.md):\n  ",
         paste(file_map$file[!file.exists(file_map$file)], collapse = "\n  "))
  }

  # Column order of the FDMA monthly sheets (one sheet per month, e.g. "2017_06")
  FDMA_COLS <- c(
    "date", "pref_code", "total",
    "age_neonate", "age_infant", "age_child", "age_adult", "age_elderly", "age_unknown",
    "sev_death", "sev_severe", "sev_moderate", "sev_mild", "sev_other",
    "loc_home", "loc_work1", "loc_work2", "loc_school",
    "loc_indoor_public", "loc_outdoor_public", "loc_road", "loc_other"
  )
  read_fdma_file <- function(filepath) {
    sheets <- excel_sheets(filepath)
    dfs <- lapply(sheets, function(sh) {
      # .name_repair = "minimal": the published sheets repeat a header ("その他");
      # names are overwritten by position below, so no repair message is needed.
      df <- read_excel(filepath, sheet = sh, .name_repair = "minimal")
      if (ncol(df) != length(FDMA_COLS)) {
        stop(sprintf("%s / sheet %s: expected %d columns, found %d",
                     basename(filepath), sh, length(FDMA_COLS), ncol(df)))
      }
      colnames(df) <- FDMA_COLS
      df$date <- as.Date(df$date)
      df
    })
    bind_rows(dfs)
  }

  cat("Loading FDMA data...\n")
  fdma_all <- bind_rows(lapply(seq_len(nrow(file_map)), function(i) {
    cat("  Reading:", file_map$file[i], "\n")
    read_fdma_file(file_map$file[i])
  }))

  # June-September only (common analysis window; 2020 lacks May) + period labels
  fdma <- fdma_all %>%
    filter(month(date) >= 6 & month(date) <= 9) %>%
    mutate(
      year  = year(date),
      month = month(date),
      doy   = yday(date),
      dow   = wday(date, label = FALSE),       # 1 = Sunday (reference level), 7 = Saturday
      period = case_when(
        year %in% 2015:2019 ~ "pre_pandemic",
        year %in% 2020:2021 ~ "pandemic",
        year %in% 2022:2023 ~ "post_pandemic"
      ),
      period = factor(period, levels = c("pre_pandemic", "pandemic", "post_pandemic"))
    )

  PREF_NAMES <- c(
    "Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima",
    "Ibaraki", "Tochigi", "Gunma", "Saitama", "Chiba", "Tokyo", "Kanagawa",
    "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano",
    "Gifu", "Shizuoka", "Aichi", "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo",
    "Nara", "Wakayama", "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi",
    "Tokushima", "Kagawa", "Ehime", "Kochi", "Fukuoka", "Saga", "Nagasaki",
    "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa"
  )
  pref_lookup <- data.frame(pref_code = 1:47, pref_name = PREF_NAMES,
                            stringsAsFactors = FALSE)
  fdma <- fdma %>% left_join(pref_lookup, by = "pref_code")

  cat("\n=== DATA QUALITY CHECKS (FDMA) ===\n")
  cat("Rows:", nrow(fdma), "| dates:", as.character(min(fdma$date)), "to",
      as.character(max(fdma$date)), "| prefectures:", n_distinct(fdma$pref_code), "\n")
  cat("Missing total:", sum(is.na(fdma$total)),
      "| negative total:", sum(fdma$total < 0, na.rm = TRUE), "\n")
  age_mismatch <- with(fdma, abs(total - (age_neonate + age_infant + age_child +
                                            age_adult + age_elderly + age_unknown)))
  cat("Age-subtotal mismatches (>0):", sum(age_mismatch > 0, na.rm = TRUE), "\n")
  stopifnot(n_distinct(fdma$pref_code) == 47, all(sort(unique(fdma$pref_code)) == 1:47),
            all(sort(unique(fdma$year)) == 2015:2023), all(fdma$month %in% 6:9))

  # Daily WBGTmax per prefectural capital (date, pref_code, wbgt_max)
  wbgt_csv <- find_data("wbgt_daily_2015_2023.csv")
  if (!file.exists(wbgt_csv)) {
    stop("wbgt_daily_2015_2023.csv not found in ./data/ (see README.md).")
  }
  wbgt <- read.csv(wbgt_csv)
  stopifnot(all(c("date", "pref_code", "wbgt_max") %in% names(wbgt)))
  wbgt$date <- as.Date(wbgt$date)
  if (anyDuplicated(wbgt[, c("date", "pref_code")])) {
    stop("wbgt_daily_2015_2023.csv has duplicated date x pref_code rows.")
  }

  dat <- fdma %>% left_join(wbgt[, c("date", "pref_code", "wbgt_max")],
                            by = c("date", "pref_code"))
  cat("Prefecture-days:", nrow(dat), "| missing WBGT after merge:",
      sum(is.na(dat$wbgt_max)), "| EAD on days without WBGT:",
      sum(dat$total[is.na(dat$wbgt_max)]), "\n")

  # Japanese public holidays within June-September (binary covariate)
  holidays_jp <- as.Date(c(
    "2015-07-20", "2015-09-21", "2015-09-22", "2015-09-23",
    "2016-07-18", "2016-08-11", "2016-09-19", "2016-09-22",
    "2017-07-17", "2017-08-11", "2017-09-18", "2017-09-23",
    "2018-07-16", "2018-08-11", "2018-09-17", "2018-09-23", "2018-09-24",
    "2019-07-15", "2019-08-11", "2019-08-12", "2019-09-16", "2019-09-23",
    "2020-07-23", "2020-07-24", "2020-08-10", "2020-09-21", "2020-09-22",
    "2021-07-22", "2021-07-23", "2021-08-08", "2021-08-09", "2021-09-20", "2021-09-23",
    "2022-07-18", "2022-08-11", "2022-09-19", "2022-09-23",
    "2023-07-17", "2023-08-11", "2023-09-18", "2023-09-23"
  ))

  dat <- dat %>%
    mutate(
      holiday = as.integer(date %in% holidays_jp),
      stratum = factor(paste(year, month, sep = "-")),   # year x month strata (conditional model)
      severe_total    = sev_death + sev_severe,           # "severe cases": death or severe
      nonsevere_total = sev_moderate + sev_mild,          # mutually exclusive complement
      death_total     = sev_death
    )

  cat("\nSummary by period:\n")
  dat %>% group_by(period) %>%
    summarise(years = paste(range(year), collapse = "-"), prefecture_days = n(),
              total_ead = sum(total), daily_mean = round(mean(total), 2),
              .groups = "drop") %>% as.data.frame() %>% print()

  rm(fdma_all, wbgt, age_mismatch)
  assign(".HSWBGT_DATA_LOADED", TRUE, envir = globalenv())
}

# ---- 2. Common exposure basis (Section 2.3.1) --------------------------------
wbgt_all    <- dat$wbgt_max[!is.na(dat$wbgt_max)]
wbgt_cen    <- median(wbgt_all)                     # centring value (pooled median, 27.5 C)
WBGT_RANGE  <- range(wbgt_all)                      # boundary knots, fixed for every prefecture-period
knots_var   <- quantile(wbgt_all, c(0.50, 0.90))    # internal exposure knots (27.5 / 31.4 C)
knots_lag   <- logknots(5, 2)                       # lag knots, lag 0-5
Q_REPORT    <- as.numeric(quantile(wbgt_all, c(0.75, 0.90, 0.95)))
names(Q_REPORT) <- c("P75", "P90", "P95")           # reporting percentiles (30.1 / 31.4 / 31.9 C)
periods     <- levels(dat$period)
prefectures <- sort(unique(dat$pref_code))
stopifnot(length(prefectures) == 47)

cat(sprintf("SETUP|median=%.3f|knots=%.3f/%.3f|P95=%.3f|range=%.1f,%.1f|prefecture_days=%d\n",
            wbgt_cen, knots_var[1], knots_var[2], Q_REPORT[["P95"]],
            WBGT_RANGE[1], WBGT_RANGE[2], nrow(dat)))

# ---- 3. Helpers (reproduced from the archived grouped_helpers.R) -------------
basis_difference <- function(at, var_knots = knots_var) {
  as.numeric(
    ns(at, knots = var_knots, Boundary.knots = WBGT_RANGE) -
      ns(wbgt_cen, knots = var_knots, Boundary.knots = WBGT_RANGE)
  )
}

rr_from_fit <- function(fit, at, var_knots = knots_var) {
  bd <- basis_difference(at, var_knots)
  est <- sum(bd * coef(fit))
  se <- sqrt(as.numeric(t(bd) %*% vcov(fit) %*% bd))
  c(log_rr = est, se = se, rr = exp(est),
    lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se))
}

# First stage: conditional quasi-Poisson DLNM per prefecture x period.
#   * exposure basis anchored to the pooled range (Boundary.knots = WBGT_RANGE),
#   * lag history reset at every summer (group = calendar year),
#   * year x month strata eliminated by gnm, DOW factor (Sunday reference), holiday.
fit_stage1_grouped <- function(outcome, min_events = 10, lag_days = 5,
                               extra_terms = character(), drop_years = integer(),
                               season_half = "all", knots_by_period = NULL,
                               label = outcome) {
  CF <- setNames(vector("list", length(periods)), periods)
  VC <- setNames(vector("list", length(periods)), periods)
  OK <- setNames(vector("list", length(periods)), periods)
  KL <- setNames(vector("list", length(periods)), periods)

  for (p in periods) {
    CF[[p]] <- vector("list", max(prefectures))
    VC[[p]] <- vector("list", max(prefectures))
    OK[[p]] <- as.list(rep(FALSE, max(prefectures)))
    kv <- if (is.null(knots_by_period)) knots_var else knots_by_period[[p]]
    KL[[p]] <- kv

    for (i in prefectures) {
      d <- dat %>%
        filter(pref_code == i, period == p, !(year %in% drop_years))
      if (season_half == "early") d <- d %>% filter(month %in% c(6, 7))
      if (season_half == "late")  d <- d %>% filter(month %in% c(8, 9))
      # Preserve the historical missing-WBGT convention while resetting each
      # summer's lag history explicitly.
      d <- d %>% filter(!is.na(wbgt_max)) %>% arrange(date)
      if (nrow(d) < 60 || sum(d[[outcome]], na.rm = TRUE) < min_events) next

      ans <- tryCatch({
        cb <- crossbasis(
          d$wbgt_max,
          lag = lag_days,
          argvar = list(fun = "ns", knots = kv, Boundary.knots = WBGT_RANGE),
          arglag = list(fun = "ns", knots = logknots(lag_days, 2)),
          group = factor(d$year)
        )
        rhs <- c("cb", "factor(dow)", "holiday", extra_terms)
        fml <- reformulate(rhs, response = outcome)
        mo <- gnm(fml, eliminate = d$stratum, family = quasipoisson(), data = d)
        rd <- crossreduce(cb, mo, cen = wbgt_cen)
        list(coef = coef(rd), vcov = vcov(rd))
      }, error = function(e) {
        message(sprintf("FAILED %s / %s / pref %d: %s", label, p, i, e$message))
        NULL
      })

      if (!is.null(ans)) {
        CF[[p]][[i]] <- ans$coef
        VC[[p]][[i]] <- ans$vcov
        OK[[p]][[i]] <- TRUE
      }
    }
    cat(sprintf("FIT|%s|%s|converged=%d\n", label, p,
                sum(vapply(OK[[p]], isTRUE, logical(1)))))
  }
  list(coef = CF, vcov = VC, conv = OK, knots = KL,
       outcome = outcome, lag = lag_days, label = label)
}

converged_ids <- function(obj, p) {
  which(vapply(obj$conv[[p]], isTRUE, logical(1)))
}

# Second stage: multivariate random-effects meta-analysis (REML) of the reduced
# coefficients; RR read analytically at `at` versus the centring value.
pool_grouped <- function(obj, p, idx = prefectures, at = Q_REPORT[["P95"]],
                         method = "reml") {
  idx <- intersect(idx, converged_ids(obj, p))
  stopifnot(length(idx) >= 5)
  Y <- do.call(rbind, obj$coef[[p]][idx])
  S <- obj$vcov[[p]][idx]
  fit <- mixmeta(Y ~ 1, S = S, method = method)
  rr <- rr_from_fit(fit, at, obj$knots[[p]])
  qs <- summary(fit)$qstat
  Q <- as.numeric(qs$Q)[1]
  df <- as.numeric(qs$df)[1]
  pQ <- as.numeric(qs$pvalue)[1]
  list(n = length(idx), ids = idx, fit = fit, rr = rr,
       Q = Q, df = df, pQ = pQ,
       I2 = max(0, (Q - df) / Q) * 100,
       tau2 = sum(diag(as.matrix(fit$Psi))) )
}

stack_one_outcome <- function(obj) {
  Y <- list(); S <- list(); M <- list()
  for (p in periods) for (i in converged_ids(obj, p)) {
    Y[[length(Y) + 1]] <- obj$coef[[p]][[i]]
    S[[length(S) + 1]] <- obj$vcov[[p]][[i]]
    M[[length(M) + 1]] <- data.frame(pref_code = i, period = p)
  }
  M <- do.call(rbind, M)
  M$period <- factor(M$period, levels = periods)
  list(Y = do.call(rbind, Y), S = S, M = M)
}

wald_block <- function(cf, V, idx) {
  b <- cf[idx]
  Vb <- V[idx, idx, drop = FALSE]
  W <- as.numeric(t(b) %*% solve(Vb) %*% b)
  c(chi2 = W, df = length(idx), p = pchisq(W, length(idx), lower.tail = FALSE))
}

contrast_block <- function(cf, V, idx_a, idx_b = NULL, bd) {
  if (is.null(idx_b)) {
    d <- cf[idx_a]
    Vd <- V[idx_a, idx_a, drop = FALSE]
  } else {
    d <- cf[idx_a] - cf[idx_b]
    Vd <- V[idx_a, idx_a, drop = FALSE] + V[idx_b, idx_b, drop = FALSE] -
      V[idx_a, idx_b, drop = FALSE] - V[idx_b, idx_a, drop = FALSE]
  }
  e <- sum(bd * d)
  se <- sqrt(as.numeric(t(bd) %*% Vd %*% bd))
  W <- as.numeric(t(d) %*% solve(Vd) %*% d)
  c(log_ratio = e, se = se, ratio = exp(e),
    lo = exp(e - 1.96 * se), hi = exp(e + 1.96 * se),
    decline_pct = 100 * (1 - exp(e)), chi2 = W, df = length(d),
    p = pchisq(W, length(d), lower.tail = FALSE))
}
