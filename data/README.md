# Data directory

The two datasets used in the analysis are **not included** in this repository
because they are redistributed from Japanese government sources. Download them
from the official sources and place them here with the exact file names below.
See the top-level `README.md` for the full file list and column specifications.

## Required files (not distributed)

**FDMA heat-related ambulance transport data**, nine annual workbooks
`heatstroke003_data_h27.xlsx` (2015) … `heatstroke003_data_r5.xlsx` (2023),
one worksheet per month (e.g. `2017_06`).
Source: https://www.fdma.go.jp/disaster/heatstroke/post3.html

**Ministry of the Environment WBGT data**, saved as a single CSV
`wbgt_daily_2015_2023.csv` with columns `date, pref_code, wbgt_max`
(daily maximum of the 24 hourly "actual values" at one designated site per
prefectural capital, June–September 2015–2023).
Source: https://www.wbgt.env.go.jp/record_data.php

## Included

**`pref_metadata.csv`** — the prefecture-level covariates used in the
meta-regression (manuscript Section 2.3.3): latitude and longitude of the
prefectural capital, proportion of residents aged 65 years and older (2020),
population density (persons/km², 2020) and its natural logarithm. Columns:
`pref_code, pref_name_en, pref_name_jp, latitude, longitude,
elderly_prop_2020, pop_density_2020, log_pop_density`.
