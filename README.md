# Temporal changes in the WBGT–heat ambulance dispatch association across Japan (2015–2023)

Analysis code for a two-stage distributed lag non-linear model (DLNM) with
multivariate meta-analysis, comparing the association between the daily maximum
wet-bulb globe temperature (WBGTmax) and heat-related emergency ambulance
dispatches (EAD) across the 47 Japanese prefectures in three periods
(pre-pandemic 2015–2019, pandemic 2020–2021, post-pandemic 2022–2023).

**Paper:** Huang P, Zhang J. *Temporal changes in the association between
wet-bulb globe temperature and heat-related emergency ambulance dispatches
across Japan: a multi-prefecture time-series analysis spanning the COVID-19
pandemic and national Heatstroke Alert implementation (2015–2023).*
Manuscript hp-2026-00307n, revision 2 (under review).

> **Release status.** This repository holds the corrected ("summer-grouped")
> pipeline described in the revision-2 manuscript and Response to Reviewers.
> Scripts `01`–`04` are the archived scripts of the manuscript run of
> 25 August 2026 with file paths adapted to this layout; `00` and `05` were
> prepared for the public release (see [Provenance](#provenance)). The
> pipeline has **not yet been re-executed from this repository layout**. Before
> the release is tagged and archived with a DOI, the archived machine-readable
> outputs will be placed in `results/` and the lock check in script `03` will
> confirm that a rerun reproduces them.

---

## What changed relative to the July 2026 scripts

The seven scripts published here in July 2026 implemented an earlier
specification. While answering the reviewers (Response to Reviewers,
"The exposure basis, and why all results have changed"; manuscript
Sections 2.3.1 and 3.2) three defects were found and corrected:

| | July 2026 scripts | This release |
|---|---|---|
| Exposure basis | `ns()` boundary knots followed each prefecture's own observed WBGT range, so the reduced coefficients pooled in the second stage were not expressed in a common basis | Boundary knots fixed at the pooled whole-study minimum and maximum for every prefecture and period (`Boundary.knots = WBGT_RANGE`) |
| Lag history | June–September series of successive years were appended without a grouping variable, so lags bridged the off-season gap | `group = year` in every `crossbasis()`: the lag history resets at the start of each summer |
| Reporting points | Nearest point of a 100-point grid (the "31.9 °C" column was read at 31.792 °C) | Evaluated analytically at the exact percentile |
| Estimable prefectures (pre / pandemic / post) | 44 / 46 / 47 | 47 / 47 / 47 |
| Pooled RR at P95 | 5.06 / 3.71 / 3.24 | 10.45 / 8.78 / 5.90 |
| Severe cases, RR at P95 | 8.14 / 4.97 / 5.04 | 21.38 / 16.22 / 12.20 |

The superseded scripts remain in the git history (last at commit `987416d`,
18 July 2026) for the record only. Their hard-coded "self-check" values belong
to the superseded specification, and the `attrdl()` function they called was
never present in the tree, so their attributable-fraction blocks failed
silently (Response to Reviewers, Comment 6).

---

## Data (not included in this repository)

The two datasets are publicly available from Japanese government sources and
are **not redistributed here**. Download them and place them in `data/`.

### 1. Heat-related emergency ambulance dispatch data (FDMA)

Fire and Disaster Management Agency (総務省消防庁), heatstroke ambulance
transport data. Download the nine annual workbooks and save them in `data/`
with these exact names:

| File | Year | Japanese era |
|------|------|--------------|
| `heatstroke003_data_h27.xlsx` | 2015 | H27 |
| `heatstroke003_data_h28.xlsx` | 2016 | H28 |
| `heatstroke003_data_h29.xlsx` | 2017 | H29 |
| `heatstroke003_data_h30.xlsx` | 2018 | H30 |
| `heatstroke003_data_r1.xlsx`  | 2019 | R1  |
| `heatstroke003_data_r2.xlsx`  | 2020 | R2  |
| `heatstroke003_data_r3.xlsx`  | 2021 | R3  |
| `heatstroke003_data_r4.xlsx`  | 2022 | R4  |
| `heatstroke003_data_r5.xlsx`  | 2023 | R5  |

Source: https://www.fdma.go.jp/disaster/heatstroke/post3.html

Each workbook contains one worksheet per month (e.g. `2017_06`). The loader in
`R/00_setup.R` keeps the June–September sheets and renames the columns, in
their published order, to `date, pref_code, total`, six age columns
(`age_neonate … age_unknown`), five severity columns
(`sev_death, sev_severe, sev_moderate, sev_mild, sev_other`) and, from 2017,
eight location columns (the 2015 and 2016 workbooks have 14 columns without
them; the loader accepts both layouts and ignores empty rows and unnamed
columns). The analysis uses `date, pref_code, total, sev_*, age_elderly`. "Severe cases" are `sev_death + sev_severe`; "non-severe cases"
are `sev_moderate + sev_mild` (`sev_other`, 0.53 % of dispatches, belongs to
neither stratum).

### 2. Daily WBGT data (Ministry of the Environment)

Daily maximum WBGT for each prefectural capital, June–September 2015–2023,
saved as a single CSV `data/wbgt_daily_2015_2023.csv` with exactly these
columns:

```
date,pref_code,wbgt_max
2015-06-01,1,23.4
...
```

- `date` — ISO `YYYY-MM-DD`
- `pref_code` — prefecture code 1–47 (matching the FDMA `pref_code`)
- `wbgt_max` — daily maximum WBGT in °C, defined as the maximum of the
  24 hourly "actual values" (実況値) of the MOE Heat Illness Prevention
  Information System at one designated observation site per prefectural
  capital (manuscript Section 2.2.2)

Source: https://www.wbgt.env.go.jp/record_data.php

The conversion of the hourly MOE records into this daily file is not part of
this repository (see [Not included](#not-included-in-this-repository)).
Prefecture-days without a WBGT value (four in the study data) are kept in the
descriptive counts and dropped from the model fits and attributable
calculations, which is the convention used throughout the manuscript.

### 3. Prefecture covariates (included)

`data/pref_metadata.csv` holds the prefecture-level covariates of the
meta-regression (Section 2.3.3): latitude of the prefectural capital,
proportion of residents aged ≥ 65 years (2020), population density (2020) and
its natural logarithm.

---

## Requirements

The archived manuscript run used R 4.5.3 with dlnm 2.4.10, mixmeta 1.2.2,
gnm 1.1-5, dplyr 1.2.0 and tsModel 0.6-2 (manuscript Section 2.3.6);
`output/sessionInfo.txt` records the versions of any rerun.

```r
install.packages(c(
  "readxl", "dplyr", "tidyr", "lubridate", "tibble",   # data handling
  "dlnm", "gnm", "mixmeta", "tsModel", "splines",       # modelling
  "ggplot2", "patchwork", "cowplot", "viridis",          # figures
  "sf", "rnaturalearth", "rnaturalearthdata"            # maps (Figures 6 and 7)
))
# ne_states() in R/04 needs the high-resolution boundaries:
install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev")
```

**`attrdl.R`.** Attributable numbers and fractions are computed with the
`attrdl()` function of Gasparrini & Leone (2014), which is not part of the
dlnm package. The upstream file carries no licence, so it is not redistributed
here: `R/00_setup.R` downloads it once into `R/attrdl.R` (git-ignored) from
https://github.com/gasparrini/2014_gasparrini_BMCmrm_Rcodedata. On a machine
without internet access, download that file manually and save it as
`R/attrdl.R`. The function calls `tsModel::Lag()`, hence the tsModel
dependency.

---

## How to run

Place the data files in `data/`, then from the repository root:

```r
source("run_all.R")        # or:  Rscript run_all.R
```

`run_all.R` runs the five scripts below in order and stops at the first error.
Each script sources `R/00_setup.R` itself, so they can also be run one at a
time (01 → 02 → 03 → 04 → 05). Every script writes `output/logs/<script>.log`
and must end with `STATUS=PASS`; assertions inside the scripts
(`stopifnot`) stop the run if, for example, the first stage is not estimable
for all 47 prefectures in every period or the period ordering of the pooled
estimates is not reproduced.

All outputs go to `output/` (tables as CSV, `output/figures/` for figures,
`output/logs/` for logs, `output/sessionInfo.txt`); `output/` is git-ignored.

### Value lock

`R/03_lock_manuscript_values.R` writes `output/MANUSCRIPT_VALUE_LOCK.csv`,
a long-form table of every display value used in the manuscript
(`key, value, display, source`). If `results/MANUSCRIPT_VALUE_LOCK.csv` — the
archived lock of the manuscript run — is present, the rerun is compared with
it key by key; any difference is written to `output/lock_mismatches.csv` and
the log of script 03 ends in `STATUS=FAIL`. This replaces the hard-coded
self-check values of the July 2026 scripts. As the Response to Reviewers
notes, such checks establish that a rerun reproduces the archived values, not
that the specification is correct.

---

## Scripts and the manuscript elements they produce

| Script | Manuscript elements | Main outputs (`output/`) |
|---|---|---|
| `R/00_setup.R` | Data preparation (Sections 2.2, 2.3.1): June–September window, periods, year × month strata, day of week, holidays, severity strata; common exposure basis (knots 27.5 / 31.4 °C, boundary knots at the pooled range, centring 27.5 °C, read-points P75 / P90 / P95); helper functions | — |
| `R/01_core_reanalysis.R` | Table 2 with heterogeneity (I², Q, trace of Ψ) and pairwise period contrasts (Section 3.2); stacked period-moderator model (Section 2.3.2); Table 3 with Benjamini–Hochberg P values and variance explained; Table 4 (all six rows); Table S1; Table S2; Table S5; support-restricted and station-exclusion sensitivity analyses (Sections 2.2.2, 3.6, Table S6 footnote); severe vs non-severe interaction and severe-only period tests (Section 3.7); within-season depletion (Sections 2.3.6, 4.1); pooled national death model (Section 2.3.6); Tokyo first-stage RR for the Fig. 2 caption | `table2_grouped_exact.csv`, `period_contrasts_independent.csv`, `heterogeneity_grouped.csv`, `period_moderator_stacked.csv`, `metaregression_grouped.csv`, `table4_grouped_exact.csv`, `prefecture_blup_p95_grouped.csv`, `prefecture_change_grouped.csv`, `period_specific_knots_grouped.csv`, `subset_sensitivity_grouped.csv`, `severity_*.csv`, `within_season_grouped_exact.csv`, `pooled_national_death_grouped.csv`, `tokyo_cumulative_p95_grouped.csv`, `grouped_stage1_*.rds` |
| `R/02_attributable_reanalysis.R` | Table 5 with Monte Carlo CIs; Table S4; Section 3.8 (18.1 % vs 26.1 %) | `attributable_summary_grouped.csv`, `attributable_prefecture_*.csv`, `attributable_pre_post_change.csv` |
| `R/03_lock_manuscript_values.R` | Section 3.6 latitude gradient with and without support at 31.9 °C; declines and extremes quoted in Sections 3.5, 3.7 and 3.8; the value lock | `MANUSCRIPT_VALUE_LOCK.csv`, `latitude_support_grouped.csv`, `*_declines_grouped.csv`, `blup_*extremes*.csv` |
| `R/04_generate_grouped_figures.R` | Figures 1–7 and the TOC graphic | `figures/Figure*_grouped.*`, `figures/TOC_grouped.*` |
| `R/05_supplementary_descriptives.R` | Table 1 (incl. convergence row and strata counts); Table S3; Figure S1; the 33 °C alert-threshold proxy and Heatstroke Alert coverage counts (Section 2.3.6; Response to Reviewers, Comment 7); share of heat-attributable dispatches on days ≥ 33 °C | `table1_descriptives.csv`, `tableS3_days_above_30.csv`, `alert_*.csv`, `attributable_share_ge33.csv`, `figures/FigureS1_grouped.*` |

Two definitions that the manuscript states only in words are fixed in code:
the "high-WBGT days" of the within-season depletion analysis are days with
WBGTmax above the pooled 95th percentile (31.9 °C) (`R/01`), and the
Kantō-Kōshin pilot region of the Heatstroke Alert comprises Ibaraki, Tochigi,
Gunma, Saitama, Chiba, Tokyo, Kanagawa, Yamanashi and Nagano (`R/05`).

---

## Not included in this repository

- **Table S6 and the station-code check** (globe-temperature diagnostic for
  the Yamaguchi station; observation-site changes in Fukushima and Kochi).
  Both require the raw hourly MOE files, which include the globe temperature
  and station codes; the daily `wbgt_max` file used here does not. The
  exclusion sensitivity analyses that follow from them *are* included
  (`subset_sensitivity_grouped.csv`).
- **Conversion of the hourly MOE records to daily WBGTmax** (Section 2.2.2).
- **The death-only first-stage eligibility counts** (21 / 13 / 10
  prefectures, Section 2.3.6); the pooled national death model that replaces
  the two-stage design for deaths *is* included.
- **Raw data** (see above).

---

## Provenance

- `R/01`–`R/04` reproduce the archived scripts of the manuscript run
  (`01_core_reanalysis.R`, `02_attributable_reanalysis.R`,
  `03_lock_manuscript_values.R`, `04_generate_grouped_figures.R`, 25 August
  2026) with three kinds of change only: absolute local paths replaced by
  repository-relative ones (`output/`, `data/`), the private helper files
  replaced by `R/00_setup.R`, and the logging calls replaced by the
  `start_log()` / `end_log()` pair. Script 03 additionally gained the
  comparison with `results/MANUSCRIPT_VALUE_LOCK.csv`.
- `R/00_setup.R` merges the two private helper files of that run: the data
  preparation (identical to Parts 1–2 of the July 2026 main script, which is
  the convention the archived helpers refer to) and the modelling helpers,
  reproduced unchanged.
- `R/05_supplementary_descriptives.R` was written for the public release so
  that the descriptive items of the manuscript can be reproduced from the same
  data objects. Its outputs can be checked against the manuscript directly.
- `results/` is reserved for the archived outputs of the manuscript run (see
  `results/README.md`).

---

## Citation

If you use this code, please cite the paper above. A `CITATION.cff` and a
DOI-archived release will be added on acceptance.

## License

MIT — see [LICENSE](LICENSE). `R/attrdl.R`, fetched at run time, is
© Antonio Gasparrini and is not covered by this licence.
