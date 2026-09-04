#!/usr/bin/env Rscript
# =============================================================================
# Phase 5 - M40a - build the sensitivity matrices that need their own cNMF run.
#
#   loo_<patient>  malignant cells with that patient REMOVED - tests whether a
#                  program survives when the patient that dominates it is gone
#   highconf       malignancy_confidence == "High" only
#
# Same declared gene universe as the primary run, so any difference is caused
# by the cells that were withheld and not by a different gene selection.
# =============================================================================
source("scripts/phase5/utils/phase5_common.R")
set.seed(42)
OUT <- file.path(P5_PROG, "cnmf_input")
facts <- list(milestone = "M40a", generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))

md <- readRDS(file.path(P5_PROG, "phase5_cell_metadata.rds"))
cnts <- readRDS(file.path(P5_PROG, "phase5_working_counts.rds"))
genes <- readLines(file.path(OUT, "primary", "genes.txt"))
mal <- md[md$malignancy_refined == "Malignant", ]
p5_msg("malignant cells %d | primary gene universe %d", nrow(mal), length(genes))

write_run <- function(cells, tag) {
  d <- file.path(OUT, tag); dir.create(d, showWarnings = FALSE, recursive = TRUE)
  sub <- cnts[genes, cells, drop = FALSE]
  sub <- sub[Matrix::rowSums(sub) > 0, , drop = FALSE]
  Matrix::writeMM(as(sub, "dgCMatrix"), file.path(d, "counts.mtx"))
  writeLines(rownames(sub), file.path(d, "genes.txt"))
  writeLines(colnames(sub), file.path(d, "cells.txt"))
  p5_msg("  [%s] %d genes x %d cells", tag, nrow(sub), ncol(sub))
  list(tag = tag, genes = nrow(sub), cells = ncol(sub))
}

runs <- list()
for (s in SAMPLES) {
  ids <- mal$cell_id[mal$sample_id != s]
  runs[[paste0("loo_", s)]] <- write_run(ids, paste0("loo_", s))
}
hc <- mal$cell_id[mal$malignancy_confidence == "High"]
p5_msg("High-confidence malignant cells: %d (%.1f%% of the compartment)",
       length(hc), 100 * length(hc) / nrow(mal))
print(table(mal$sample_id[mal$malignancy_confidence == "High"]))
runs$highconf <- write_run(hc, "highconf")

facts$runs <- runs
facts$loo_design <- "each run removes ONE patient entirely; a program that vanishes when its dominant patient is withheld is not robust"
facts$highconf_n <- length(hc)
p5_json(facts, file.path(P5_VAL, "m40a_sensitivity_inputs_facts.json"))
p5_sec("M40a complete")
