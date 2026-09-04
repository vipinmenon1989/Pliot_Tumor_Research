#!/usr/bin/env Rscript
# =============================================================================
# Phase 4 · M35 — the CNA heatmap required by §20 and §52.
#
# A genome-ordered inferred-CNA heatmap per patient, with annotation tracks for
# the SCEVAN call, the refined malignancy call, the SCEVAN clone and the frozen
# Phase 2 annotation. This is the figure that lets a reader see for themselves
# that the promoted "fibroblasts" carry the same CNA architecture as the cells
# Phase 2 already called malignant.
#
# Rendering only - no new analysis (§56).
# =============================================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ComplexHeatmap); library(circlize); library(grid); library(Matrix)
})
set.seed(42)

SC  <- "results/phase4/scevan/by_sample"
MAL <- "results/phase4/malignancy"
FIG <- "results/phase4/figures"
SAMPLES <- c("MPNST_1","MPNST_2","MPNST_3","MPNST_4")
MAX_CELLS <- 2000L      # per patient, for legibility; sampling is seeded
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
log_ <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

calls <- read.delim(file.path(MAL, "PHASE4_MALIGNANCY_CALLS.tsv"), stringsAsFactors = FALSE)
rownames(calls) <- calls$cell_id
rel <- read.delim("results/phase4/tables/SCEVAN_SAMPLE_RELIABILITY.tsv", stringsAsFactors = FALSE)

P4_MAL <- c(Malignant = "#B2182B", `Non-malignant` = "#2166AC",
            Ambiguous = "#F0A202", `Excluded-low-quality` = "grey65")
P4_SCEVAN <- c(malignant = "#B2182B", `non-malignant` = "#2166AC",
               `not-assessed` = "grey75")

panels <- list()
for (s in SAMPLES) {
  f  <- file.path(SC, s, "primary", "output", paste0(s, "_primary_CNAmtx.RData"))
  fa <- file.path(SC, s, "primary", "output", paste0(s, "_primary_count_mtx_annot.RData"))
  if (!file.exists(f) || !file.exists(fa)) { log_("skip ", s); next }
  e <- new.env(); load(f, envir = e)
  nm <- ls(e)[vapply(ls(e), function(n) is.matrix(get(n, e)), logical(1))][1]
  M <- get(nm, e)
  a <- new.env(); load(fa, envir = a); ann <- get("count_mtx_annot", a)
  ann <- ann[!duplicated(ann$gene_name), ]; rownames(ann) <- ann$gene_name
  g <- intersect(rownames(M), rownames(ann))
  M <- M[g, , drop = FALSE]; ann <- ann[g, , drop = FALSE]
  o <- order(as.numeric(as.character(ann$seqnames)), as.numeric(as.character(ann$start)))
  M <- M[o, , drop = FALSE]; chr <- as.numeric(as.character(ann$seqnames))[o]

  cd <- calls[colnames(M), , drop = FALSE]
  # Order cells so structure is visible: refined call, then clone, then Phase 2 label.
  ord <- order(factor(cd$malignancy_refined,
                      levels = c("Malignant","Ambiguous","Non-malignant","Excluded-low-quality")),
               ifelse(is.na(cd$scevan_clone), "zzz", cd$scevan_clone),
               cd$annotation_ccc_phase3)
  M <- M[, ord, drop = FALSE]; cd <- cd[ord, , drop = FALSE]
  if (ncol(M) > MAX_CELLS) {
    keep <- sort(sample(ncol(M), MAX_CELLS))
    M <- M[, keep, drop = FALSE]; cd <- cd[keep, , drop = FALSE]
    sub_note <- sprintf(" (%s of %s cells, seeded random subsample)",
                        format(MAX_CELLS, big.mark=","), format(length(ord), big.mark=","))
  } else sub_note <- sprintf(" (all %s assessed cells)", format(ncol(M), big.mark=","))

  clones <- sort(unique(na.omit(cd$scevan_clone)))
  cpal <- if (length(clones))
    setNames(colorRampPalette(RColorBrewer::brewer.pal(8,"Dark2"))(length(clones)), clones) else character(0)
  anns <- sort(unique(cd$annotation_ccc_phase3))
  apal <- setNames(colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(length(anns)), anns)

  ha <- rowAnnotation(
    `Phase 2 annotation` = cd$annotation_ccc_phase3,
    `SCEVAN call`        = cd$malignancy_scevan,
    `refined call`       = cd$malignancy_refined,
    `SCEVAN clone`       = ifelse(is.na(cd$scevan_clone), "none", cd$scevan_clone),
    col = list(`Phase 2 annotation` = apal,
               `SCEVAN call` = P4_SCEVAN,
               `refined call` = P4_MAL,
               `SCEVAN clone` = c(cpal, none = "grey90")),
    annotation_name_gp = gpar(fontsize = 7),
    simple_anno_size = unit(3.2, "mm"),
    annotation_legend_param = list(labels_gp = gpar(fontsize = 6),
                                   title_gp = gpar(fontsize = 7, fontface = "bold")))

  reliable <- rel$scevan_reliable[rel$sample_id == s]
  ttl <- sprintf("%s — inferred CNA%s%s", s, sub_note,
                 if (length(reliable) && !isTRUE(reliable))
                   "\nSANITY-GATE FAILURE: this sample's malignant/normal partition is unreliable" else "")
  panels[[s]] <- list(M = t(M), chr = chr, ha = ha, title = ttl)
  rm(M, e, a); invisible(gc(FALSE))
  log_("prepared ", s)
}

if (!length(panels)) { log_("no CNA matrices available; figure not produced"); quit(status = 0) }

cf <- colorRamp2(c(-0.25, 0, 0.25), c("#2166AC", "white", "#B2182B"))
for (dev_ in c("pdf","png")) {
  fp <- file.path(FIG, paste0("04_scevan_cna_heatmap.", dev_))
  if (dev_ == "pdf") pdf(fp, width = 15.5, height = 4.6 * length(panels))
  else png(fp, width = 15.5, height = 4.6 * length(panels), units = "in", res = 200, bg = "white")
  pushViewport(viewport(layout = grid.layout(length(panels), 1)))
  for (i in seq_along(panels)) {
    p <- panels[[i]]
    pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1))
    ht <- Heatmap(p$M, name = "inferred\nrelative CNA", col = cf,
      cluster_rows = FALSE, cluster_columns = FALSE,
      show_row_names = FALSE, show_column_names = FALSE,
      column_split = factor(p$chr, levels = sort(unique(p$chr))),
      column_title_gp = gpar(fontsize = 6),
      column_gap = unit(0.4, "mm"), border = TRUE,
      row_title = "cells (ordered by refined call, then clone)",
      row_title_gp = gpar(fontsize = 7),
      left_annotation = p$ha,
      heatmap_legend_param = list(labels_gp = gpar(fontsize = 6),
                                  title_gp = gpar(fontsize = 7, fontface = "bold")),
      use_raster = TRUE, raster_quality = 3)
    draw(ht, column_title = p$title,
         column_title_gp = gpar(fontsize = 10, fontface = "bold"),
         newpage = FALSE, merge_legend = TRUE,
         heatmap_legend_side = "right", annotation_legend_side = "right")
    popViewport()
  }
  popViewport(); dev.off()
  cat(sprintf("  [fig] %s\n", fp))
}
log_("CNA heatmap complete. Columns are chromosomes 1-22 in genomic order;")
log_("chrX/chrY are absent because SCEVAN's annotation does not cover them.")
