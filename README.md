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
  "ggplot2", "patchwork", "splines"
))
```

`splines` ships with base R. The core modelling packages are
[`dlnm`](https://cran.r-project.org/package=dlnm) (cross-basis construction and
reduction), [`gnm`](https://cran.r-project.org/package=gnm) (conditional
quasi-Poisson first stage), and
[`mixmeta`](https://cran.r-project.org/package=mixmeta) (multivariate REML
meta-analysis).

---

## How to run

1. Clone the repository and place all data files in `data/` (see above).
2. From the repository root:

```r
source("R/heatstroke_dlnm_analysis.R")
```

The script creates `output/` and `output/figures/` and writes the summary
tables (`table1_summary.csv`, `rr_comparison_by_period.csv`, …) and the main
figures there.

3. To reproduce the knot-placement sensitivity analysis (Table S5), run the
main script first (so that `dat` and the fitted objects are in memory), then:

```r
source("R/sensitivity_knots_S5.R")
```

---

## Repository contents

```
.
├── README.md
├── LICENSE
├── .gitignore
├── R/
│   ├── heatstroke_dlnm_analysis.R   # full pipeline (Parts 0–8)
│   └── sensitivity_knots_S5.R       # knot-placement sensitivity (Table S5)
└── data/
    └── README.md                    # placeholder; put downloaded data here
```

### What the main script does (Parts 0–8)

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
  (`rr_comparison_by_period.csv`, figures).

### Sensitivity and extension code

- **`R/sensitivity_knots_S5.R`** is complete, runnable code for the
  knot-placement sensitivity analysis reported as **Table S5** (re-estimating
  each period with knots at that period's own WBGT percentiles).
- **Parts 7 and 8** of the main script are **documented as notes**, not executable
  blocks. Part 7 records the parameter changes for the other sensitivity analyses
  reported in the paper (alternative maximum lag; exclusion of 2020; severity and
  age subgroups), each of which is produced by re-running Parts 4–6 with the
  stated modification. Part 8 sketches the prefecture-level meta-regression.
  These can be made fully scripted on request.

---

## Citation

If you use this code, please cite the paper above. A `CITATION.cff` will be
added on publication.

## License

MIT — see [LICENSE](LICENSE).
