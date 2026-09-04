# =============================================================================
# Phase 5 shared constants and helpers.
#
# Phase 5 asks whether CONTINUOUS malignant transcriptional programs recur
# across patients even though Phase 4's DISCRETE malignant states did not
# (0 of 8 recurrent; 97.2% of malignant cells in patient-private states).
#
# Nothing in Phase 1-4 is modified. The Phase 4 object is opened read-only and
# every Phase 5 artefact is written under results/phase5/.
# =============================================================================
suppressPackageStartupMessages({
  library(Matrix); library(dplyr); library(tidyr); library(jsonlite)
})

P5_ROOT    <- "/local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST"
P4_OBJECT  <- "results/phase4/phase4_final_object.rds"
P4_MD5     <- "e85ba8486e456917e2483f2773bdbaf3"
P4_SHA256  <- "a658447744e5ee8621fc5a3127492ceb11bce37bac4d6f3f0f7255bdcb647ae4"

P5_PROG <- "results/phase5/programs"
P5_ECM  <- "results/phase5/malignant_ecm"
P5_VAL  <- "results/phase5/validation"
P5_TAB  <- "results/phase5/tables"
P5_TABF <- "results/phase5/tables/final"
P5_FIG  <- "results/phase5/figures/final"
P5_REP  <- "reports/phase5"

SAMPLES  <- c("MPNST_1", "MPNST_2", "MPNST_3", "MPNST_4")
# MPNST_3 failed the Phase 4 SCEVAN immune sanity gate; its clones are immune
# artefacts and carry no tumour-clone inferential weight (Phase 4 handoff A1/A2).
CLONE_RELIABLE <- c("MPNST_1", "MPNST_2", "MPNST_4")

# frozen Phase 4 counts - asserted, never assumed
P4_FROZEN <- list(
  total_cells = 19716L,
  refined = c(Malignant = 6434L, `Non-malignant` = 9078L,
              Ambiguous = 3766L, `Excluded-low-quality` = 438L),
  fibroblast = c(total = 5064L, Malignant = 4036L,
                 `Non-malignant` = 908L, Ambiguous = 120L),
  n_clones = c(MPNST_1 = 7L, MPNST_2 = 4L, MPNST_3 = 3L, MPNST_4 = 8L))

p5_msg <- function(...) cat(sprintf(...), "\n", sep = "")
p5_sec <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n",
                          strrep("=", 78), "\n", sep = "")

p5_json <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  # digits = NA keeps full precision: jsonlite's default (4) silently rounds.
  write(toJSON(x, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA),
        path)
  p5_msg("  [json] %s", path)
}

p5_tsv <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  p5_msg("  [tsv]  %s (%d x %d)", path, nrow(x), ncol(x))
}

# ---- checksum helpers -------------------------------------------------------
p5_md5 <- function(p) unname(tools::md5sum(p))
p5_sha256 <- function(p) {
  out <- system2("sha256sum", shQuote(p), stdout = TRUE)
  sub(" .*$", "", out[1])
}

# ---- canonical cell-cycle genes (Seurat's Tirosh et al. lists) --------------
# Used ONLY for the declared cell-cycle sensitivity analysis. Cell cycle may be
# a real malignant program, so it is never removed from the primary analysis.
p5_cc_genes <- function() unique(c(Seurat::cc.genes.updated.2019$s.genes,
                                   Seurat::cc.genes.updated.2019$g2m.genes))
