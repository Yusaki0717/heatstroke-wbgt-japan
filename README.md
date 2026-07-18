# Temporal changes in the WBGT–heat ambulance dispatch association across Japan (2015–2023)

Analysis code for a two-stage distributed lag non-linear model (DLNM) with
multivariate meta-analysis, comparing the association between the daily maximum
wet-bulb globe temperature (WBGTmax) and heat-related emergency ambulance
dispatches (EAD) across 47 Japanese prefectures over three periods
(pre-pandemic 2015–2019, pandemic 2020–2021, post-pandemic 2022–2023).

**Paper:** Huang P, Zhang J. *Temporal changes in the association between
wet-bulb globe temperature and heat-related emergency ambulance dispatches
across Japan: a multi-prefecture time-series analysis spanning the COVID-19
pandemic and national Heatstroke Alert implementation (2015–2023).*
*Environmental Health Perspectives* (manuscript hp-2026-00307n; under review).

---

## Data (not included in this repository)

The two datasets used in the analysis are publicly available from Japanese
government sources. They are **not redistributed here**; download them from the
official sources below and place them in the `data/` directory.

### 1. Heat-related emergency ambulance dispatch data (FDMA)

Fire and Disaster Management Agency (総務省消防庁), heatstroke ambulance
transport data. Download the nine annual files and place them in `data/` with
these exact names:

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
`Part 1` keeps the June–September sheets and standardises the columns to:
`date, pref_code, total`, the age columns
(`age_neonate, age_infant, age_child, age_adult, age_elderly, age_unknown`),
the severity columns (`sev_death, sev_severe, sev_moderate, sev_mild, sev_other`),
and the location columns. Only `date, pref_code, total, sev_*, age_elderly` are
used in the published analysis.

### 2. Daily WBGT data (Ministry of the Environment)

Daily maximum WBGT for each prefectural capital, 2015–2023, June–September.
Save as a single CSV at `data/wbgt_daily_2015_2023.csv` with exactly these
columns:

```
date,pref_code,wbgt_max
2015-06-01,1,23.4
...
```

- `date` — ISO `YYYY-MM-DD`
- `pref_code` — prefecture code 1–47 (matching the FDMA `pref_code`)
- `wbgt_max` — daily maximum WBGT in °C

Source: https://www.wbgt.env.go.jp/record_data.php

---

## Requirements

- R (≥ 4.2 recommended)
- R packages:

```r
install.packages(c(
  "readxl", "dplyr", "tidyr", "lubridate",
  "dlnm", "gnm", "mixmeta",
  "ggplot2", "patchwork", "splines",
  "sf", "viridis", "rnaturalearth", "rnaturalearthdata"
))
```

`splines` ships with base R. The core modelling packages are
[`dlnm`](https://cran.r-project.org/package=dlnm) (cross-basis construction and
reduction), [`gnm`](https://cran.r-project.org/package=gnm) (conditional
quasi-Poisson first stage), and
[`mixmeta`](https://cran.r-project.org/package=mixmeta) (multivariate REML
meta-analysis). `sf`, `viridis`, and `rnaturalearth`/`rnaturalearthdata` are
needed only for the choropleth maps (`R/spatial_maps_and_af.R`).

---

## How to run

Run the scripts **in order** from the repository root; each depends on
objects left in the workspace by the previous ones.

```r
source("R/heatstroke_dlnm_analysis.R")        # 1. main pipeline (required first)
source("R/spatial_maps_and_af.R")             # 2. maps + AF/AN (Fig 4, Fig 6, Table 5)
source("R/revision_sensitivity_analyses.R")   # 3. response-to-reviewers analyses
source("R/sensitivity_knots_S5.R")            # 4. knot-placement sensitivity (Table S5)
source("R/meta_regression.R")                 # 5. prefecture-level meta-regression (Table 3)
source("R/sensitivity_table4.R")              # 6. Table 4 sensitivity rows + Fig 5
source("R/severe_interaction_test.R")         # 7. period × exposure interaction test
```

or run everything in one go:

```r
source("run_all.R")   # sources the seven scripts above in order
```

Place all data files in `data/` first (see above). Each script creates its own
`output/` subfolders as needed and writes CSVs and figures there.

> **Note on `R/spatial_maps_and_af.R`:** its opening block ("0. Bridge")
> derives two intermediate files (`output/pref_rr_at_p95.csv`,
> `output/rr_change_post_vs_pre.csv`) from the `blup_results` object produced
> by the main script. This bridge was reconstructed for this repository
> release and has been verified: `pref_rr_at_p95.csv` reproduces the Table S1
> values exactly (e.g. post-pandemic Kagoshima 2.36, Aichi 4.03), and
> `rr_change_post_vs_pre.csv` correctly restricts to the 44 prefectures that
> converged in both the pre- and post-pandemic periods (excluding Hokkaido,
> Aomori, and Yamaguchi), matching the range reported in Table S2
> (−10.2% to −50.6%).

### Built-in self-checks

Every pipeline script verifies its own output against the values reported in
the manuscript:

- `heatstroke_dlnm_analysis.R` compares `rr_table` against **Table 2**
  programmatically (all 9 percentile×period cells) and raises a warning on any
  mismatch.
- `revision_sensitivity_analyses.R` (block C2) prints the expected severe-case
  row of **Table 4** (8.14 / 4.97 / 5.04) alongside its output.
- `meta_regression.R` prints the expected **Table 3** Wald statistics.
- `sensitivity_table4.R` prints the expected **Table 4** sensitivity rows.
- `severe_interaction_test.R` prints the expected interaction-test statistics
  reported in **Section 3.7**.

A successful end-to-end run therefore reproduces every published number; any
drift between code and manuscript is visible immediately in the console output.

---

## Repository contents

```
.
├── README.md
├── LICENSE
├── .gitignore
├── run_all.R                             # one-shot runner for the full pipeline
├── R/
│   ├── heatstroke_dlnm_analysis.R        # 1. full pipeline (Parts 0–6 + notes)
│   ├── spatial_maps_and_af.R             # 2. maps + attributable fraction/number
│   ├── revision_sensitivity_analyses.R   # 3. response-to-reviewers analyses
│   ├── sensitivity_knots_S5.R            # 4. knot-placement sensitivity (Table S5)
│   ├── meta_regression.R                 # 5. prefecture-level meta-regression (Table 3)
│   ├── sensitivity_table4.R              # 6. Table 4 sensitivity rows + Fig 5
│   └── severe_interaction_test.R         # 7. period × exposure interaction test
└── data/
    └── README.md                         # placeholder; put downloaded data here
```

### What the main script does (Parts 0–6)

- **Part 0** — packages and output directories.
- **Part 1** — load and clean the nine FDMA workbooks; restrict to June–September;
  assign each observation to a study period.
- **Part 2** — load the WBGT CSV and merge it with the FDMA data by
  `date` and `pref_code`.
- **Part 3** — descriptive statistics (`table1_*`).
- **Part 4** — first stage: for each prefecture within each period, a conditional
  quasi-Poisson DLNM (`gnm`, eliminating year×month strata) with natural-spline
  cross-bases for exposure and lag (lag 0–5 days). Exposure knots are placed at
  the 50th and 90th percentiles of the **pooled whole-study** WBGT distribution
  and held fixed across periods; the reference is the pooled whole-study median.
- **Part 5** — second stage: multivariate REML meta-analysis (`mixmeta`) of the
  reduced coefficients, with best linear unbiased predictions per prefecture.
- **Part 6** — pooled exposure–response and lag–response curves, cumulative
  relative risks at reference percentiles, and the period comparison
  (`rr_comparison_by_period.csv`, figures), followed by the Table 2 self-check.

The script closes with short **notes** pointing to the dedicated scripts that
implement the sensitivity analyses (lag 0–7, exclusion of 2020, elderly
subgroup: `sensitivity_table4.R`; severe cases: `revision_sensitivity_analyses.R`
block C2) and the meta-regression (`meta_regression.R`).

### What `spatial_maps_and_af.R` does

- **Part 0 (Bridge)** — derives per-prefecture RR-at-P95 tables from
  `blup_results` (see note above).
- **Part 1** — downloads a Japan prefecture shapefile (NaturalEarth) and
  matches it to `pref_code`.
- **Part A** — three choropleth maps: post-pandemic RR (→ Fig 6), all three
  periods side by side, and percent change post- vs pre-pandemic (→ Fig 4).
- **Part B** — attributable fraction/number by prefecture and period, for
  all-severity and severe-case outcomes (→ Table 5), plus an AF map and bar
  chart.

### What `revision_sensitivity_analyses.R` does

Each block answers a specific reviewer comment and supports a specific part
of the manuscript:

| Block | Manuscript location |
|-------|----------------------|
| B1 | Table S4 — AF sensitivity to the reference (centring) value |
| B2 | Section 4.2 — burden comparison restricted to commonly-converged prefectures (16.3% vs 20.3%) |
| B3 | Fig. 7a — verification of the annualized attributable-number annotations |
| C2 | Table 4 "Severe cases only" row (8.14 / 4.97 / 5.04) and the severe-case attributable numbers in Table 5 |
| C3 | Table S3 — convergence diagnostic (high-WBGT-day counts for Hokkaido, Aomori, Yamaguchi), cited in Section 3.2 |

### Sensitivity and meta-regression scripts

- **`R/sensitivity_knots_S5.R`** — knot-placement sensitivity reported as
  **Table S5** (re-estimating each period with knots at that period's own WBGT
  percentiles).
- **`R/meta_regression.R`** — univariate random-effects meta-regressions of the
  reduced coefficients on prefecture-level latitude, % population aged ≥65, and
  log population density, with multivariate Wald tests (**Table 3**). The
  prefecture covariate table is embedded in the script.
- **`R/sensitivity_table4.R`** — the three all-severity sensitivity rows of
  **Table 4** (lag 0–7 days; exclusion of year 2020; elderly ≥65 subgroup),
  using the same two-stage runner as the main analysis, and assembly of
  **Fig. 5** from the pipeline CSV outputs.
- **`R/severe_interaction_test.R`** — formal test of whether the exposure–response
  association differs across periods: reduced coefficients from all
  prefecture–period units are stacked into a single multivariate random-effects
  meta-regression with period as a categorical predictor, and the period
  coefficients are compared with multivariate Wald tests (overall, and
  pre→pandemic / pandemic→post contrasts), reported in Sections 2.3.4 and 3.7.

---

## Citation

If you use this code, please cite the paper above. A `CITATION.cff` will be
added on publication.

## License

MIT — see [LICENSE](LICENSE).
