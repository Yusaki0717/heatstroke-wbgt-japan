# =============================================================================
# run_all.R -- run the complete corrected pipeline for hp-2026-00307n
#
# Order (each script sources R/00_setup.R itself):
#   R/01_core_reanalysis.R          first stage, meta-analysis, Tables 2-4, S1, S2, S5
#   R/02_attributable_reanalysis.R  attributable numbers / fractions, Table 5, S4
#   R/03_lock_manuscript_values.R   MANUSCRIPT_VALUE_LOCK.csv (+ check against results/)
#   R/04_generate_grouped_figures.R Figures 1-7, TOC graphic
#   R/05_supplementary_descriptives.R Table 1, Table S3, Figure S1, 33 C proxy
#
# Prerequisites: data files in data/ (see README.md). Run from anywhere:
#   source("run_all.R")      or      Rscript run_all.R
# The script switches to its own folder, so relative paths always resolve
# against the repository root. Every script stops on its first error and
# writes output/logs/<script>.log ending in STATUS=PASS.
# =============================================================================

.locate_root <- function() {
  f <- NULL
  m <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(m)) f <- sub("^--file=", "", m[1])
  if (is.null(f)) {
    f <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = FALSE), error = function(e) NULL)
  }
  if ((is.null(f) || !nzchar(f)) &&
      requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    f <- rstudioapi::getActiveDocumentContext()$path
  }
  if (is.null(f) || !nzchar(f)) getwd() else dirname(normalizePath(f))
}
setwd(.locate_root())
cat("Pipeline root:", getwd(), "\n")

scripts <- c(
  "R/01_core_reanalysis.R",
  "R/02_attributable_reanalysis.R",
  "R/03_lock_manuscript_values.R",
  "R/04_generate_grouped_figures.R",
  "R/05_supplementary_descriptives.R"
)
missing <- scripts[!file.exists(scripts)]
if (length(missing)) {
  stop("Missing analysis scripts (is run_all.R in the repository root?):\n  ",
       paste(missing, collapse = "\n  "))
}

t_start <- Sys.time()
for (s in scripts) {
  cat("\n", strrep("=", 76), "\nRUNNING: ", s, "\n", strrep("=", 76), "\n\n", sep = "")
  tryCatch(
    source(s, encoding = "UTF-8"),
    error = function(e) {
      while (sink.number() > 0) sink()
      stop(sprintf("Pipeline stopped in %s:\n  %s", s, conditionMessage(e)), call. = FALSE)
    }
  )
}

# ---- Summary of per-script status lines ------------------------------------
logs <- file.path("output", "logs", paste0(sub("\\.R$", "", basename(scripts)), ".log"))
status <- vapply(logs, function(f) {
  if (!file.exists(f)) return("NO LOG")
  x <- readLines(f, warn = FALSE)
  st <- grep("^END\\|", x, value = TRUE)
  if (length(st)) sub(".*STATUS=", "", st[length(st)]) else "INCOMPLETE"
}, character(1))
cat("\n", strrep("=", 76), "\n", sep = "")
print(data.frame(log = basename(logs), status = unname(status)), row.names = FALSE)
cat(sprintf("Elapsed: %.1f min\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
if (!all(status == "PASS")) {
  stop("One or more scripts did not finish with STATUS=PASS; see output/logs/.")
}
cat("ALL SCRIPTS COMPLETE. Outputs are in output/ (tables, logs) and output/figures/.\n")
