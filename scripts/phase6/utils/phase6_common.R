# =============================================================================
# Phase 6 shared constants and helpers.
#
# Phase 6 combines the Phase 5 continuous malignant programs with the frozen
# Phase 4 SCEVAN CNA clones and asks whether malignant transcriptional
# heterogeneity is clone-constrained, diverse within clones, or mixed - and
# which TF and pathway activities distinguish the programs.
#
# Two rules dominate every script here:
#   * MPNST_3's SCEVAN clone structure failed the Phase 4 immune sanity gate
#     and is excluded from EVERY clone-based inference.
#   * clone labels are patient-scoped. MPNST_1_clone1 and MPNST_4_clone1 are
#     unrelated names and are never pooled as homologous.
# =============================================================================
suppressPackageStartupMessages({
  library(Matrix); library(dplyr); library(tidyr); library(jsonlite)
})

P5_OBJECT <- "results/phase5/phase5_final_object.rds"
P4_OBJECT <- "results/phase4/phase4_final_object.rds"
P4_MD5    <- "e85ba8486e456917e2483f2773bdbaf3"
P4_SHA256 <- "a658447744e5ee8621fc5a3127492ceb11bce37bac4d6f3f0f7255bdcb647ae4"

P6_CLONE <- "results/phase6/clone_program"
P6_REG   <- "results/phase6/regulatory"
P6_PATH  <- "results/phase6/pathways"
P6_PLAS  <- "results/phase6/plasticity"
P6_CNA   <- "results/phase6/cna_expression"
P6_VAL   <- "results/phase6/validation"
P6_TAB   <- "results/phase6/tables"
P6_TABF  <- "results/phase6/tables/final"
P6_FIG   <- "results/phase6/figures/final"

SAMPLES        <- c("MPNST_1", "MPNST_2", "MPNST_3", "MPNST_4")
CLONE_RELIABLE <- c("MPNST_1", "MPNST_2", "MPNST_4")
CLONE_EXCLUDED <- "MPNST_3"
CLONE_EXCLUSION_REASON <- paste(
  "MPNST_3 failed the Phase 4 SCEVAN immune sanity gate: its primary and",
  "sensitivity runs are inverted (agreement 0.0864 against 0.9663-0.9974),",
  "six of nine canonical immune populations are called ~100% malignant, and",
  "its three inferred clones are T/NK, plasma/pDC/B and plasma-dominated.",
  "Its clone structure carries no tumour-clone inferential weight.")

# declared clone-size thresholds for the diversity analyses
CLONE_SIZE_GRID <- c(20L, 50L, 100L)
CLONE_SIZE_PRIMARY <- 20L

p6_msg <- function(...) cat(sprintf(...), "\n", sep = "")
p6_sec <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n",
                          strrep("=", 78), "\n", sep = "")
p6_json <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  write(toJSON(x, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA), path)
  p6_msg("  [json] %s", path)
}
p6_tsv <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  p6_msg("  [tsv]  %s (%d x %d)", path, nrow(x), ncol(x))
}
p6_md5 <- function(p) unname(tools::md5sum(p))
p6_sha256 <- function(p) sub(" .*$", "", system2("sha256sum", shQuote(p), stdout = TRUE)[1])

# Shannon entropy and its effective-number transform, on a probability vector
shannon <- function(p) { p <- p[p > 0]; -sum(p * log(p)) }
effective_n <- function(p) exp(shannon(p))

P6_CAPTION <- paste(
  "n = 4 patients and sample_id = patient = dataset. MPNST_3's SCEVAN clone structure failed the Phase 4 immune sanity gate",
  "and is EXCLUDED from every clone-based inference; clone labels are patient-scoped and never homologous across patients.",
  "SCEVAN infers broad copy number from expression - it is not DNA sequencing, and CNA-expression association is an internal",
  "consistency analysis, NOT independent validation. Within-clone program diversity is not observed state switching.",
  sep = "\n")
