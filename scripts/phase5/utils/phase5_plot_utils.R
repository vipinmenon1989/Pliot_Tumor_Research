# =============================================================================
# Phase 5 / Phase 6 shared plotting helpers.
# Carries forward the Phase 3/4 lessons: ggplot2 >= 4 rejects dpi = NA, and a
# figure that clips its legend is not a deliverable, so every figure declares
# its own page size.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(scales); library(patchwork) })

P5_DPI_PDF <- 200; P5_DPI_PNG <- 220

save_fig5 <- function(plot, path, width, height, png = TRUE) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  ggsave(path, plot, width = width, height = height, units = "in",
         device = "pdf", dpi = P5_DPI_PDF, limitsize = FALSE)
  if (png) ggsave(sub("\\.pdf$", ".png", path), plot, width = width,
                  height = height, units = "in", device = "png",
                  dpi = P5_DPI_PNG, bg = "white", limitsize = FALSE)
  cat(sprintf("  [fig] %s (%.1f x %.1f in)\n", path, width, height))
  invisible(path)
}

p5_theme <- function(base = 10) {
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey94", colour = "grey70"),
          strip.text = element_text(face = "bold", size = base - 1),
          plot.title = element_text(face = "bold", size = base + 2),
          plot.subtitle = element_text(size = base - 1, colour = "grey25"),
          plot.caption = element_text(size = base - 2, colour = "grey35", hjust = 0),
          legend.key.size = unit(0.42, "cm"))
}

P5_SAMPLE <- c(MPNST_1 = "#1B9E77", MPNST_2 = "#D95F02",
               MPNST_3 = "#7570B3", MPNST_4 = "#E7298A")
P5_REC <- c(recurrent = "#1A9850", `shared-limited` = "#66BD63",
            `patient-private` = "#D73027", uncertain = "grey65")
P5_MAL <- c(Malignant = "#B2182B", `Non-malignant` = "#2166AC",
            Ambiguous = "#F0A202", `Excluded-low-quality` = "grey65")

wrap5 <- function(x, width = 150) {
  paste(vapply(strsplit(x, "\n", fixed = TRUE)[[1]],
               function(l) paste(strwrap(l, width = width), collapse = "\n"),
               character(1)), collapse = "\n")
}

P5_CAPTION <- paste(
  "n = 4 patients and sample_id = patient = dataset, so cells are NOT biological replicates:",
  "conclusions rest on patient-level direction, recurrence, effect size and robustness, not on cell-level p-values.",
  "A cNMF program is a continuous pattern of co-regulated expression - not a cell type, a lineage, a state or a clone.",
  sep = "\n")
