# Archived outputs of the manuscript run

Machine-readable outputs of the archived manuscript run (25 August 2026;
Windows, R 4.5.3, dlnm 2.4.10, mixmeta 1.2.2, gnm 1.1-5) from which every
value in the revision-2 manuscript was taken. Because the input data cannot
be redistributed, these files are what a reader can check the manuscript
against without rerunning the pipeline, and what `R/03_lock_manuscript_values.R`
compares a rerun with.

Contents:

- `MANUSCRIPT_VALUE_LOCK.csv` — long-form value lock: `key, value, display, source`
  (66 keys: Table 2, period contrasts, severity strata, attributable fractions).
- Result tables with the same names as a rerun writes under `output/`
  (`table2_grouped_exact.csv`, `table4_grouped_exact.csv`, `metaregression_grouped.csv`,
  `prefecture_blup_p95_grouped.csv`, `prefecture_change_grouped.csv`,
  `attributable_summary_grouped.csv`, `severity_*.csv`, `period_contrasts_independent.csv`, …).
- `alert_burden_summary_grouped.csv`, `alert_burden_prefecture_grouped.csv` — the
  33 °C alert-day burden analysis of the archived run (its script is not part of
  this release; `R/05_supplementary_descriptives.R` reproduces the summary as
  `attributable_share_ge33.csv`).
- `grouped_stage1_main_severity.rds`, `grouped_stage1_sensitivity.rds` — the
  first-stage reduced coefficients and covariance matrices (all 47 prefectures,
  three periods; main, severe, non-severe, lag 0–7, excluding 2020, elderly,
  period-specific knots), from which every second-stage result can be recomputed
  without the raw data.
- `logs/` — run logs of scripts 01–04 and of the alert-burden script, each ending
  in `STATUS=PASS`.
- `sessionInfo.txt` — R and package versions of the archived run.

A rerun of this repository on the study data reproduces these files; see
"Verification" in the top-level README for the one documented difference
(Monte Carlo interval bounds, platform-dependent at the 0.1 % level).
