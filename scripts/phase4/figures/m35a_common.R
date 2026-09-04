# =============================================================================
# Phase 4 - M35A shared helpers for SCEVAN figure consolidation.
#
# M35A surfaces and consolidates EXISTING Phase 4 SCEVAN evidence. It re-runs
# nothing: no SCEVAN call, no clone assignment, no malignancy rule and no
# threshold is recomputed here. Every number plotted is read back from the
# frozen M29/M30 outputs and asserted against the frozen handoff before use.
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scales)
  library(grid); library(patchwork)
})
source("scripts/phase4/utils/phase4_plot_utils.R")

M35A_MAL     <- "results/phase4/malignancy"
M35A_TAB     <- "results/phase4/tables"
M35A_TABF    <- "results/phase4/tables/final"
M35A_FIG     <- "results/phase4/figures/final"
M35A_SCEVAN  <- "results/phase4/scevan/by_sample"
M35A_SAMPLES <- c("MPNST_1", "MPNST_2", "MPNST_3", "MPNST_4")
M35A_RELIABLE <- c("MPNST_1", "MPNST_2", "MPNST_4")

# ---- frozen values from PHASE4_HANDOFF.md -----------------------------------
# Every one of these is asserted, not assumed. A figure is never written from
# an input that disagrees with the frozen handoff.
M35A_FROZEN <- list(
  total_cells   = 19716L,
  refined       = c(Malignant = 6434L, `Non-malignant` = 9078L,
                    Ambiguous = 3766L, `Excluded-low-quality` = 438L),
  fibroblast    = c(total = 5064L, Malignant = 4036L,
                    `Non-malignant` = 908L, Ambiguous = 120L),
  candidate     = c(Malignant = 836L, Ambiguous = 395L, `Non-malignant` = 0L),
  mpnst_tumor   = c(Malignant = 1405L, Ambiguous = 2015L),
  n_clones      = c(MPNST_1 = 7L, MPNST_2 = 4L, MPNST_3 = 3L, MPNST_4 = 8L),
  agreement     = c(MPNST_1 = 0.99514117, MPNST_2 = 0.96628722,
                    MPNST_3 = 0.08639456, MPNST_4 = 0.99738258)
)

m35a_load_meta <- function() {
  md <- readRDS(file.path(M35A_MAL, "phase4_malignancy_metadata.rds"))
  md$cell_id <- rownames(md)
  md$malignancy_refined <- factor(md$malignancy_refined,
    levels = c("Malignant", "Non-malignant", "Ambiguous", "Excluded-low-quality"))
  md$sample_id <- factor(md$sample_id, levels = M35A_SAMPLES)
  md
}

# ---- validation -------------------------------------------------------------
# Hard stop, not a warning. If the plotting input has drifted from the frozen
# handoff the correct action is to diagnose it, never to quietly plot it.
m35a_validate <- function(md) {
  chk <- list()
  add <- function(name, got, want) {
    ok <- identical(as.integer(got), as.integer(want))
    chk[[name]] <<- list(check = name, expected = as.integer(want),
                         observed = as.integer(got), pass = ok)
    cat(sprintf("  [%s] %-42s expected %-7s observed %s\n",
                if (ok) "OK  " else "FAIL", name,
                paste(want, collapse = "/"), paste(got, collapse = "/")))
    ok
  }
  ok <- TRUE
  ok <- add("total_cells", nrow(md), M35A_FROZEN$total_cells) && ok
  ref <- table(md$malignancy_refined)
  for (k in names(M35A_FROZEN$refined))
    ok <- add(paste0("refined_", k), ref[[k]], M35A_FROZEN$refined[[k]]) && ok
  fb <- md[md$annotation_ccc == "Fibroblast", ]
  ok <- add("fibroblast_total", nrow(fb), M35A_FROZEN$fibroblast[["total"]]) && ok
  for (k in c("Malignant", "Non-malignant", "Ambiguous"))
    ok <- add(paste0("fibroblast_", k), sum(fb$malignancy_refined == k),
              M35A_FROZEN$fibroblast[[k]]) && ok
  cm <- md[md$annotation_ccc == "Candidate-Malignant-Unresolved", ]
  for (k in c("Malignant", "Ambiguous", "Non-malignant"))
    ok <- add(paste0("candidate_", k), sum(cm$malignancy_refined == k),
              M35A_FROZEN$candidate[[k]]) && ok
  mt <- md[md$annotation_ccc == "MPNST-Tumor", ]
  for (k in c("Malignant", "Ambiguous"))
    ok <- add(paste0("mpnst_tumor_", k), sum(mt$malignancy_refined == k),
              M35A_FROZEN$mpnst_tumor[[k]]) && ok
  for (s in M35A_SAMPLES) {
    n <- length(unique(na.omit(md$scevan_clone[md$sample_id == s])))
    ok <- add(paste0("clones_", s), n, M35A_FROZEN$n_clones[[s]]) && ok
  }
  if (!ok) stop("M35A: plotting input contradicts the frozen Phase 4 handoff. ",
                "Figure freeze STOPPED - diagnose before plotting.")
  cat("  all frozen-value checks passed\n")
  invisible(do.call(rbind, lapply(chk, as.data.frame)))
}

# ---- annotation grouping for the composition figures ------------------------
# The three disputed identities stay separate and named; the canonical lineages
# that are not in dispute are collapsed so they cannot dominate the legend.
M35A_ANN_LEVELS <- c("Fibroblast", "Candidate-Malignant-Unresolved", "MPNST-Tumor",
                     "Endothelial", "Pericyte-VSMC", "Plasma-cell",
                     "Other immune", "Other")
M35A_OTHER_IMMUNE <- c("CD4-T", "CD8-T", "NK", "T-cell-other", "B-cell",
                       "Macrophage", "Monocyte", "Dendritic", "Plasmacytoid-DC")

m35a_group_ann <- function(x) {
  x <- as.character(x)
  out <- ifelse(x %in% c("Fibroblast", "Candidate-Malignant-Unresolved",
                         "MPNST-Tumor", "Endothelial", "Pericyte-VSMC",
                         "Plasma-cell"), x,
         ifelse(x %in% M35A_OTHER_IMMUNE, "Other immune", "Other"))
  factor(out, levels = M35A_ANN_LEVELS)
}

M35A_ANN_COL <- c(
  "Fibroblast"                     = "#D95F02",
  "Candidate-Malignant-Unresolved" = "#7570B3",
  "MPNST-Tumor"                    = "#B2182B",
  "Endothelial"                    = "#66A61E",
  "Pericyte-VSMC"                  = "#A6761D",
  "Plasma-cell"                    = "#1F78B4",
  "Other immune"                   = "#9EC9E2",
  "Other"                          = "grey72")

M35A_CAPTION <- paste(
  "SCEVAN infers copy number from gene expression. This is NOT DNA sequencing and provides no DNA-level proof.",
  "Inference is reliable for broad chromosomal, arm-level and large-segment events only; no single-gene CNV call is asserted.",
  "A SCEVAN non-malignant call is not proof of non-malignancy. n = 4 patients; no population-level claim follows.",
  sep = "\n")

# Long subtitles clip at the page edge unless they are wrapped. Wrap on words,
# preserving any newline the caller wrote deliberately.
wrap_sub <- function(x, width = 150) {
  paste(vapply(strsplit(x, "\n", fixed = TRUE)[[1]],
               function(l) paste(strwrap(l, width = width), collapse = "\n"),
               character(1)), collapse = "\n")
}

# hg38 chr1-22 lengths, used only to lay out a shared genome axis so profiles
# from different patients (which retain different numbers of genes) are drawn
# on the same coordinate.
M35A_CHRLEN <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
                 159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
                 114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
                 58617616, 64444167, 46709983, 50818468)
names(M35A_CHRLEN) <- as.character(1:22)
M35A_CHROFF <- c(0, cumsum(as.numeric(M35A_CHRLEN))[-22])
names(M35A_CHROFF) <- names(M35A_CHRLEN)
M35A_CHRMID <- M35A_CHROFF + M35A_CHRLEN / 2
M35A_GTOT   <- sum(as.numeric(M35A_CHRLEN))

m35a_msg <- function(...) cat(sprintf(...), "\n", sep = "")
