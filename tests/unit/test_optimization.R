# tests/unit/test_optimization.R
library(yaml)
source("scripts/R/logging_utils.R")

setup_strict_logging()
log_info("Starting QC optimization review validation tests...", stage = "test_optimization")

report_path <- "reports/QC_OPTIMIZATION_REPORT.md"
opt_dir <- "reports/qc_optimization"

# Check report presence
if (!file.exists(report_path)) {
  log_error("Missing QC optimization report.", stage = "test_optimization")
  quit(status = 1)
}

# Check TSV files
stats_tsv <- file.path(opt_dir, "sensitivity_summary_stats.tsv")
overlap_tsv <- file.path(opt_dir, "overlap_summary.tsv")

if (!file.exists(stats_tsv)) {
  log_error("Missing sensitivity summary statistics TSV.", stage = "test_optimization")
  quit(status = 1)
}
if (!file.exists(overlap_tsv)) {
  log_error("Missing overlap summary TSV.", stage = "test_optimization")
  quit(status = 1)
}

# Load TSV and verify rows
stats <- read.delim(stats_tsv, sep = "\t")
if (nrow(stats) == 0) {
  log_error("Sensitivity summary statistics TSV is empty.", stage = "test_optimization")
  quit(status = 1)
}

overlaps <- read.delim(overlap_tsv, sep = "\t")
if (nrow(overlaps) == 0) {
  log_error("Overlap summary TSV is empty.", stage = "test_optimization")
  quit(status = 1)
}

# Check sensitivity curves figures
expected_figs <- c("sensitivity_mt.png", "sensitivity_ribo.png", "sensitivity_features.png", "sensitivity_counts.png", "distribution_comparison.png")
for (fig in expected_figs) {
  fpath <- file.path(opt_dir, fig)
  if (!file.exists(fpath)) {
    log_error(sprintf("Missing optimization figure: %s", fig), stage = "test_optimization")
    quit(status = 1)
  }
}

log_info("QC optimization review validation tests PASSED.", stage = "test_optimization")
quit(status = 0)
