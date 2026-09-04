#!/usr/bin/env Rscript
# =============================================================================
# Phase 4 - M35A - FIGURE 1
#   17_scevan_native_cna_MPNST1.pdf   18_..._MPNST2.pdf   19_..._MPNST4.pdf
#   20_scevan_cna_reliable_patients.pdf
#
# The native SCEVAN CNA heatmaps already exist as PNG. They are SURFACED here,
# never redrawn: the CNA matrix, its cell ordering and its subclone track are
# SCEVAN's own output, embedded verbatim. The Phase 2 / clone / malignancy
# annotation the native plot cannot carry is supplied as an ALIGNED COMPANION
# PANEL built from the same frozen clone assignments, so no CNA value is
# touched and no fragile re-plot of SCEVAN internals is attempted.
# =============================================================================
source("scripts/phase4/figures/m35a_common.R")
suppressPackageStartupMessages({ library(png); library(grid) })

md <- m35a_load_meta()
m35a_msg("== M35A figures 17-20: native SCEVAN CNA heatmaps ==")
m35a_validate(md)

# ---- read a native PNG, trim its uniform white margin -----------------------
read_native <- function(path, pad = 6L) {
  stopifnot(file.exists(path))
  x <- png::readPNG(path)
  if (length(dim(x)) == 2L) x <- array(x, c(dim(x), 1L))
  ink <- apply(x[, , 1:min(3L, dim(x)[3]), drop = FALSE], c(1, 2),
               function(v) any(v < 0.985))
  rr <- range(which(rowSums(ink) > 0)); cc <- range(which(colSums(ink) > 0))
  rr <- c(max(1L, rr[1] - pad), min(dim(x)[1], rr[2] + pad))
  cc <- c(max(1L, cc[1] - pad), min(dim(x)[2], cc[2] + pad))
  x[rr[1]:rr[2], cc[1]:cc[2], , drop = FALSE]
}
native_grob <- function(path) {
  img <- read_native(path)
  rasterGrob(img, interpolate = TRUE, width = unit(1, "npc"),
             height = unit(1, "npc"))
}
native_panel <- function(path, title, sub = NULL) {
  wrap_elements(full = native_grob(path)) +
    labs(title = title, subtitle = sub) +
    theme(plot.title = element_text(face = "bold", size = 11,
                                    margin = margin(b = 3)),
          plot.subtitle = element_text(size = 8.5, colour = "grey25",
                                       margin = margin(b = 6)),
          plot.margin = margin(t = 10, r = 4, b = 4, l = 4))
}
nf <- function(sample, what) file.path(M35A_SCEVAN, sample, "primary", "output",
                                       sprintf("%s_primary%s.png", sample, what))

# ---- companion annotation panel, from the same frozen clone assignments -----
cl <- md |> filter(!is.na(scevan_clone), scevan_clone != "") |>
  mutate(ann = m35a_group_ann(annotation_ccc))

companion <- function(sample) {
  d <- cl |> filter(sample_id == sample)
  sz <- d |> count(scevan_clone, name = "clone_size")
  comp <- d |> count(scevan_clone, ann, name = "n") |>
    left_join(sz, by = "scevan_clone") |> mutate(frac = n / clone_size)
  mal <- d |> group_by(scevan_clone) |>
    summarise(fm = mean(malignancy_refined == "Malignant"), .groups = "drop") |>
    left_join(sz, by = "scevan_clone")
  lv <- sz |> arrange(desc(clone_size)) |> pull(scevan_clone)
  lab <- setNames(sub("^MPNST_[0-9]+_clone", "clone ", lv), lv)
  comp$scevan_clone <- factor(comp$scevan_clone, levels = rev(lv))
  mal$scevan_clone  <- factor(mal$scevan_clone,  levels = rev(lv))

  p1 <- ggplot(comp, aes(y = scevan_clone, x = frac, fill = ann)) +
    geom_col(width = 0.78, colour = "white", linewidth = 0.15) +
    scale_fill_manual(values = M35A_ANN_COL, name = "Phase 2 annotation_ccc",
                      drop = FALSE) +
    scale_x_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
    scale_y_discrete(labels = lab) +
    labs(x = "fraction of clone", y = NULL,
         title = sprintf("%s - Phase 2 identity of each CNA-defined clone", sample),
         subtitle = "Same clone assignments as the heatmap above. SCEVAN never saw these labels.") +
    p4_theme(9) + theme(plot.title = element_text(size = 10))
  p2 <- ggplot(mal, aes(y = scevan_clone, x = fm)) +
    geom_col(width = 0.78, fill = "#B2182B", alpha = 0.85) +
    geom_text(aes(x = fm, label = sprintf(" %.0f%%  (n=%s)", 100 * fm,
                  format(clone_size, big.mark = ","))),
              hjust = 0, size = 2.5, colour = "grey20") +
    scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.42),
                       breaks = c(0, 0.5, 1), expand = c(0, 0)) +
    scale_y_discrete(labels = NULL) +
    labs(x = "refined Malignant", y = NULL,
         title = "malignancy_refined", subtitle = "clone size at right") +
    p4_theme(9) + theme(plot.title = element_text(size = 10))
  p1 + p2 + plot_layout(widths = c(2.6, 1))
}

# ---- per-patient figures 17 / 18 / 19 ---------------------------------------
per_patient <- list(MPNST_1 = "17_scevan_native_cna_MPNST1.pdf",
                    MPNST_2 = "18_scevan_native_cna_MPNST2.pdf",
                    MPNST_4 = "19_scevan_native_cna_MPNST4.pdf")
nclone <- tapply(cl$scevan_clone, droplevels(cl$sample_id),
                 function(v) length(unique(v)))

for (s in names(per_patient)) {
  out <- file.path(M35A_FIG, per_patient[[s]])
  pdf(out, width = 13.5, height = 15.5, onefile = TRUE)

  hd <- native_panel(nf(s, "heatmap_subclones"),
    sprintf("%s - native SCEVAN CNA heatmap, tumour cells, subclone-annotated", s),
    wrap_sub(paste0("SCEVAN output, embedded verbatim: rows are cells grouped by CNA-defined subclone (left colour bar), ",
           "columns are genes in genomic order, chr1-22. ",
           sprintf("%d subclones. Red = inferred relative gain, blue = inferred relative loss. ",
                   nclone[[s]]),
           "Chromosome X and Y are not assessable by this method."), 150))
  pg1 <- hd / companion(s) + plot_layout(heights = c(2.5, 1)) +
    plot_annotation(
      title = sprintf("Figure %s. %s - CNA-defined subclone architecture beside its Phase 2 composition",
                      sub("_.*", "", per_patient[[s]]), s),
      caption = M35A_CAPTION,
      theme = theme(plot.title = element_text(face = "bold", size = 13),
                    plot.caption = element_text(size = 7.5, colour = "grey35", hjust = 0)))
  print(pg1)

  pg2 <- (native_panel(nf(s, "heatmap"),
            sprintf("%s - all assessed cells, malignant / normal track", s),
            "SCEVAN's own classification of every assessed cell, primary run (norm_cell = NULL).") /
          native_panel(nf(s, "onlytumorheatmap"),
            sprintf("%s - tumour cells only", s),
            "The malignant compartment on its own, same matrix and ordering.") /
          native_panel(nf(s, "consensus"),
            sprintf("%s - consensus clonal CN profile", s),
            "Segment-level consensus copy number across the malignant compartment.")) +
    plot_layout(heights = c(1.25, 1.1, 0.5)) +
    plot_annotation(
      title = sprintf("%s - supporting native SCEVAN CNA output (primary run)", s),
      caption = M35A_CAPTION,
      theme = theme(plot.title = element_text(face = "bold", size = 13),
                    plot.caption = element_text(size = 7.5, colour = "grey35", hjust = 0)))
  print(pg2)
  dev.off()
  m35a_msg("  [fig] %s (2 pages)", out)
}

# ---- figure 20: the three reliable patients on one page ---------------------
lbl <- c(MPNST_1 = "MPNST_1 - 7 subclones; all 7 mix Candidate-Malignant-Unresolved with MPNST-Tumor, 6 of 7 also carry Fibroblast cells",
         MPNST_2 = "MPNST_2 - 4 subclones; ALL 4 are fibroblast-dominated",
         MPNST_4 = "MPNST_4 - 8 subclones; 7 of 8 are fibroblast-dominated (clone 8 is endothelial)")
p20 <- (native_panel(nf("MPNST_1", "heatmap_subclones"), lbl[["MPNST_1"]]) /
        native_panel(nf("MPNST_2", "heatmap_subclones"), lbl[["MPNST_2"]]) /
        native_panel(nf("MPNST_4", "heatmap_subclones"), lbl[["MPNST_4"]])) +
  plot_layout(heights = c(1.35, 1, 1.35)) +
  plot_annotation(
    title = "Figure 20. Native SCEVAN CNA heatmaps, subclone-annotated - the three reliable patients",
    subtitle = wrap_sub(paste0(
      "SCEVAN output embedded verbatim (primary run, norm_cell = NULL). Rows = cells grouped by CNA-defined subclone; ",
      "columns = genes in genomic order, chr1-22.\n",
      "CNA architecture is patient-specific: pairwise Jaccard of broad events among these three patients is only 0.20-0.23, ",
      "with a small shared core (chr18 loss, chr2 and chr7 gain).\n",
      "MPNST_3 is excluded: its SCEVAN partition failed the immune sanity gate and its malignant promotions were disabled (figure 25)."), 150),
    caption = M35A_CAPTION,
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 9, colour = "grey25"),
                  plot.caption = element_text(size = 7.5, colour = "grey35", hjust = 0)))
out20 <- file.path(M35A_FIG, "20_scevan_cna_reliable_patients.pdf")
pdf(out20, width = 13.0, height = 16.0); print(p20); dev.off()
m35a_msg("  [fig] %s", out20)
m35a_msg("== figures 17-20 done ==")
