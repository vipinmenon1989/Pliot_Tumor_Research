# =============================================================================
# Phase 4 shared plotting helpers
#
# Two lessons carried forward from Phase 3 are encoded here:
#   * ggplot2 (>= 4.0) REJECTS dpi = NA in ggsave. Always pass a number.
#   * A figure that clips its legend or overplots its points is not a
#     deliverable. save_fig() takes explicit page dimensions per figure.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(scales) })

P4_DPI_PDF <- 200
P4_DPI_PNG <- 220

save_fig <- function(plot, path, width, height, png = TRUE) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  ggsave(path, plot, width = width, height = height, units = "in",
         device = "pdf", dpi = P4_DPI_PDF, limitsize = FALSE)
  if (png) {
    p2 <- sub("\\.pdf$", ".png", path)
    ggsave(p2, plot, width = width, height = height, units = "in",
           device = "png", dpi = P4_DPI_PNG, bg = "white", limitsize = FALSE)
  }
  cat(sprintf("  [fig] %s (%.1f x %.1f in)\n", path, width, height))
  invisible(path)
}

p4_theme <- function(base = 10) {
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey94", colour = "grey70"),
          strip.text = element_text(face = "bold", size = base - 1),
          plot.title = element_text(face = "bold", size = base + 2),
          plot.subtitle = element_text(size = base - 1, colour = "grey25"),
          plot.caption = element_text(size = base - 2, colour = "grey35", hjust = 0),
          legend.key.size = unit(0.42, "cm"))
}

# Stable, readable palettes. Colour carries meaning consistently across every
# Phase 4 figure so a reader can move between panels without relearning.
P4_MAL <- c(Malignant = "#B2182B", `Non-malignant` = "#2166AC",
            Ambiguous = "#F0A202", `Excluded-low-quality` = "grey65")
P4_SCEVAN <- c(malignant = "#B2182B", `non-malignant` = "#2166AC",
               `not-assessed` = "grey75")
P4_CONF <- c(High = "#08519C", Moderate = "#6BAED6", Low = "#C6DBEF")
P4_SAMPLE <- c(MPNST_1 = "#1B9E77", MPNST_2 = "#D95F02",
               MPNST_3 = "#7570B3", MPNST_4 = "#E7298A")

# UMAP panels: draw the highlighted class last so it is never buried, and keep
# point size/alpha tied to n so 19,716 cells do not become a solid blob.
umap_layer <- function(df, x, y, colour_col, palette, title, subtitle = NULL,
                       caption = NULL, highlight = NULL, legend_title = colour_col) {
  n <- nrow(df)
  ps <- if (n > 15000) 0.30 else if (n > 5000) 0.45 else 0.8
  al <- if (n > 15000) 0.55 else 0.7
  if (!is.null(highlight)) {
    df$.ord <- ifelse(as.character(df[[colour_col]]) %in% highlight, 2L, 1L)
    df <- df[order(df$.ord), , drop = FALSE]
  }
  ggplot(df, aes(x = .data[[x]], y = .data[[y]], colour = .data[[colour_col]])) +
    geom_point(size = ps, alpha = al, stroke = 0) +
    scale_colour_manual(values = palette, name = legend_title, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 2.6, alpha = 1))) +
    labs(title = title, subtitle = subtitle, caption = caption,
         x = "UMAP 1 (Phase 2 Harmony)", y = "UMAP 2 (Phase 2 Harmony)") +
    coord_equal() + p4_theme() +
    theme(axis.text = element_blank(), axis.ticks = element_blank())
}

P4_CAPTION_CNV <- paste(
  "SCEVAN infers copy number from gene-expression patterns; this is not DNA sequencing.",
  "Inference is reliable for broad chromosomal, arm-level and large-segment events, not single genes.",
  "A SCEVAN non-malignant call is not proof of non-malignancy: some MPNST cells may be copy-number quiet.",
  sep = "\n")
