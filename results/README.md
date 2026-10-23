# Archived outputs of the manuscript run

This folder is reserved for the machine-readable outputs of the archived
manuscript run (25 August 2026) from which every value in the revision-2
manuscript was taken. Because the input data cannot be redistributed, these
files are what a reader can check the manuscript against without rerunning the
pipeline, and what `R/03_lock_manuscript_values.R` compares a rerun with.

Expected contents (same names as produced under `output/` by a rerun):

- `MANUSCRIPT_VALUE_LOCK.csv` — long-form value lock: `key, value, display, source`
- the named result tables, e.g. `table2_grouped_exact.csv`, `table4_grouped_exact.csv`,
  `metaregression_grouped.csv`, `prefecture_blup_p95_grouped.csv`,
  `prefecture_change_grouped.csv`, `attributable_summary_grouped.csv`,
  `severity_interaction_fixed.csv`, `period_contrasts_independent.csv`, …
- `logs/*.log` — run logs ending in `STATUS=PASS`
- `sessionInfo.txt` — R and package versions of the archived run

Status: **not yet populated**. The authors will copy the archived outputs here
before the release is tagged and archived with a DOI.
