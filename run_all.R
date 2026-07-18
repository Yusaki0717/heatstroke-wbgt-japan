# =============================================================================
# run_all.R — full pipeline for hp-2026-00307n
#
# Runs every analysis script in dependency order.
# Prerequisite: the data files must already be in data/ (see README.md).
#
# Usage: source this file from anywhere —
#
#   source("run_all.R")                    # if it is in your working directory
#   source("C:/path/to/repo/run_all.R")    # full path works too
#
# The script locates its own folder and switches the working directory there,
# so relative paths in the analysis scripts (data/, output/, R/) always resolve
# against the repository root, no matter where the R session was started.
#
# Watch for the self-check lines in the console output (see README,
# "Built-in self-checks"); the main script warns if Table 2 is not reproduced.
# =============================================================================

# --- 0. Locate this script's own directory ------------------------------------
.this_file <- tryCatch(
  normalizePath(sys.frame(1)$ofile, mustWork = FALSE),
  error = function(e) NULL
)
if (is.null(.this_file) &&
    requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  .this_file <- rstudioapi::getActiveDocumentContext()$path
}
if (is.null(.this_file) || !nzchar(.this_file)) {
  .script_dir <- getwd()
} else {
  .script_dir <- dirname(normalizePath(.this_file))
}
setwd(.script_dir)
cat("Pipeline root:", getwd(), "\n")

# --- 1. Resolve script locations (R/ subfolder preferred, flat layout ok) -----
script_names <- c(
  "heatstroke_dlnm_analysis.R",       # 1. main pipeline (required first)
  "spatial_maps_and_af.R",            # 2. maps + AF/AN (Fig 4, Fig 6, Table 5)
  "revision_sensitivity_analyses.R",  # 3. response-to-reviewers analyses
  "sensitivity_knots_S5.R",           # 4. knot-placement sensitivity (Table S5)
  "meta_regression.R",                # 5. prefecture-level meta-regression (Table 3)
  "sensitivity_table4.R",             # 6. Table 4 sensitivity rows + Fig 5
  "severe_interaction_test.R"         # 7. period x exposure interaction test
)

resolve_script <- function(f) {
  cand <- c(file.path("R", f), f)
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0) return(NA_character_)
  normalizePath(hit[1])
}
scripts <- vapply(script_names, resolve_script, character(1))

if (any(is.na(scripts))) {
  cat("\nThe following scripts were NOT found under", getwd(), "\n")
  cat(paste0("  - ", script_names[is.na(scripts)], collapse = "\n"), "\n")
  cat("(looked in both ./R/ and ./). Check that run_all.R sits in the\n")
  cat("repository root together with the R/ folder, then try again.\n")
  stop("Missing analysis scripts; aborting before any computation.")
}

# --- 2. Run in order -----------------------------------------------------------
for (s in scripts) {
  cat("\n", strrep("=", 76), "\n", sep = "")
  cat("RUNNING:", s, "\n")
  cat(strrep("=", 76), "\n\n", sep = "")
  source(s, encoding = "UTF-8")
}

cat("\n", strrep("=", 76), "\n", sep = "")
cat("ALL SCRIPTS COMPLETE. Check the self-check lines above before trusting\n")
cat("any downstream use of the outputs.\n")
