#!/usr/bin/env Rscript
# =============================================================================
# Phase 4 · M30 figures (§27)
#   01_phase2_vs_scevan_malignancy      05_malignancy_by_annotation
#   02_scevan_malignancy_umap           06_malignancy_by_patient
#   03_refined_malignancy_umap          07_ambiguous_population_malignancy
#   04_malignancy_by_phase2_cluster     08_cnv_profiles_by_annotation
# PDF primary, PNG alongside.
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(Matrix)
})
source("scripts/phase4/utils/phase4_plot_utils.R")

MAL <- "results/phase4/malignancy"
FIG <- "results/phase4/figures"
TAB <- "results/phase4/tables"
SC  <- "results/phase4/scevan/by_sample"
SAMPLES <- c("MPNST_1","MPNST_2","MPNST_3","MPNST_4")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

md <- readRDS(file.path(MAL, "phase4_malignancy_metadata.rds"))
cat("cells:", nrow(md), "\n")
UX <- "postint_umap_harmony_1"; UY <- "postint_umap_harmony_2"
stopifnot(all(c(UX, UY) %in% colnames(md)))

md$malignancy_refined <- factor(md$malignancy_refined,
  levels = c("Malignant","Non-malignant","Ambiguous","Excluded-low-quality"))
md$malignancy_scevan <- factor(md$malignancy_scevan,
  levels = c("malignant","non-malignant","not-assessed"))
md$malignancy_confidence <- factor(md$malignancy_confidence,
  levels = c("High","Moderate","Low"))
ann_ord <- md |> count(annotation_ccc) |> arrange(desc(n)) |> pull(annotation_ccc)
md$annotation_ccc <- factor(md$annotation_ccc, levels = ann_ord)

# ---- 01 Phase 2 vs SCEVAN ------------------------------------------------
d1 <- md |>
  count(annotation_ccc, malignancy_scevan) |>
  group_by(annotation_ccc) |> mutate(frac = n / sum(n), tot = sum(n)) |> ungroup()
lab <- md |> count(annotation_ccc, name = "tot") |>
  mutate(l = sprintf("%s  (n=%s)", annotation_ccc, format(tot, big.mark = ",")))
p1 <- ggplot(d1, aes(y = annotation_ccc, x = frac, fill = malignancy_scevan)) +
  geom_col(width = 0.78) +
  geom_text(data = lab, aes(y = annotation_ccc, x = 1.005, label = format(tot, big.mark=",")),
            inherit.aes = FALSE, hjust = 0, size = 2.7, colour = "grey30") +
  scale_fill_manual(values = P4_SCEVAN, name = "SCEVAN call\n(primary run)") +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.13))) +
  scale_y_discrete(limits = rev(ann_ord)) +
  labs(title = "Phase 2 annotation versus independent SCEVAN malignancy call",
       subtitle = paste("Primary SCEVAN run: norm_cell = NULL, so confident-normal detection saw",
                        "no Phase 2 label.\nRow labels are the frozen Phase 2 annotation_ccc; n at right."),
       x = "fraction of cells in the population", y = NULL,
       caption = P4_CAPTION_CNV) + p4_theme()
save_fig(p1, file.path(FIG, "01_phase2_vs_scevan_malignancy.pdf"), 10.5, 6.4)

# ---- 02 SCEVAN malignancy UMAP -------------------------------------------
p2 <- umap_layer(md, UX, UY, "malignancy_scevan", P4_SCEVAN,
  "SCEVAN malignancy call on the Phase 2 Harmony UMAP",
  paste0("The UMAP is for VISUALISATION ONLY. Classification comes from raw counts and inferred ",
         "copy number,\nnot from these coordinates. Malignant cells drawn last."),
  P4_CAPTION_CNV, highlight = "malignant", legend_title = "SCEVAN call")
save_fig(p2, file.path(FIG, "02_scevan_malignancy_umap.pdf"), 8.2, 7.4)

# ---- 03 refined malignancy UMAP (+ confidence facet) ---------------------
p3a <- umap_layer(md, UX, UY, "malignancy_refined", P4_MAL,
  "Evidence-integrated refined malignancy",
  "malignancy_refined, from the a priori rules in MALIGNANCY_DECISION_RULES.md",
  NULL, highlight = "Malignant", legend_title = "malignancy_refined")
p3b <- md |> filter(malignancy_refined == "Malignant") |>
  umap_layer(UX, UY, "malignancy_confidence", P4_CONF,
    "Confidence of the Malignant calls", "Malignant cells only",
    NULL, legend_title = "confidence")
p3 <- {
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    (p3a | p3b) + plot_annotation(caption = P4_CAPTION_CNV,
      theme = theme(plot.caption = element_text(size = 8, colour = "grey35", hjust = 0)))
  } else p3a
}
save_fig(p3, file.path(FIG, "03_refined_malignancy_umap.pdf"), 14.5, 7.2)

# ---- 04 malignancy by Phase 2 cluster ------------------------------------
cl <- md |>
  mutate(cluster = factor(postint_harmony_primary_cluster,
                          levels = sort(unique(as.integer(as.character(postint_harmony_primary_cluster)))))) |>
  count(cluster, malignancy_refined) |>
  group_by(cluster) |> mutate(frac = n/sum(n), tot = sum(n)) |> ungroup()
domann <- md |> group_by(cluster = factor(postint_harmony_primary_cluster,
                  levels = levels(cl$cluster))) |>
  summarise(dom = names(sort(table(annotation_ccc), decreasing = TRUE))[1],
            tot = n(), .groups = "drop") |>
  mutate(l = sprintf("C%s · %s (n=%s)", cluster, dom, format(tot, big.mark=",")))
cl <- cl |> left_join(domann |> select(cluster, l), by = "cluster")
cl$l <- factor(cl$l, levels = domann$l)
p4 <- ggplot(cl, aes(y = l, x = frac, fill = malignancy_refined)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = P4_MAL, name = "malignancy_refined", drop = FALSE) +
  scale_x_continuous(labels = percent_format(accuracy = 1), expand = expansion(0)) +
  scale_y_discrete(limits = rev(domann$l)) +
  labs(title = "Refined malignancy by frozen Phase 2 cluster",
       subtitle = "Phase 2 Harmony resolution 1.0 clusters, labelled with their dominant annotation_ccc",
       x = "fraction of cluster", y = NULL, caption = P4_CAPTION_CNV) +
  p4_theme() + theme(axis.text.y = element_text(size = 7.5))
save_fig(p4, file.path(FIG, "04_malignancy_by_phase2_cluster.pdf"), 10.5, 8.6)

# ---- 05 malignancy by annotation, with confidence -----------------------
d5 <- md |> count(annotation_ccc, malignancy_refined, malignancy_confidence) |>
  group_by(annotation_ccc) |> mutate(frac = n/sum(n)) |> ungroup() |>
  mutate(cell = paste(malignancy_refined, malignancy_confidence, sep = " / "))
p5 <- ggplot(d5, aes(y = annotation_ccc, x = frac,
                     fill = malignancy_refined, alpha = malignancy_confidence)) +
  geom_col(width = 0.78) +
  scale_fill_manual(values = P4_MAL, name = "malignancy_refined", drop = FALSE) +
  scale_alpha_manual(values = c(High = 1, Moderate = 0.62, Low = 0.32),
                     name = "confidence", drop = FALSE) +
  scale_x_continuous(labels = percent_format(accuracy = 1), expand = expansion(0)) +
  scale_y_discrete(limits = rev(ann_ord)) +
  labs(title = "Refined malignancy and confidence, by Phase 2 annotation",
       subtitle = "Opacity encodes confidence; a pale bar is a Low-confidence call and should be read as such",
       x = "fraction of population", y = NULL, caption = P4_CAPTION_CNV) + p4_theme()
save_fig(p5, file.path(FIG, "05_malignancy_by_annotation.pdf"), 10.5, 6.4)

# ---- 06 malignancy by patient -------------------------------------------
d6 <- md |> count(sample_id, malignancy_refined) |>
  group_by(sample_id) |> mutate(frac = n/sum(n), tot = sum(n)) |> ungroup()
p2f <- md |> group_by(sample_id) |>
  summarise(phase2 = mean(malignancy_phase2 == "malignant"),
            scevan = mean(malignancy_scevan == "malignant"),
            refined = mean(malignancy_refined == "Malignant"), .groups = "drop") |>
  pivot_longer(-sample_id, names_to = "stage", values_to = "frac") |>
  mutate(stage = factor(stage, levels = c("phase2","scevan","refined"),
    labels = c("Phase 2\nconservative","SCEVAN\n(primary)","Phase 4\nrefined")))
p6a <- ggplot(d6, aes(x = sample_id, y = frac, fill = malignancy_refined)) +
  geom_col(width = 0.72) +
  geom_text(data = distinct(d6, sample_id, tot),
            aes(x = sample_id, y = 1.02, label = format(tot, big.mark=",")),
            inherit.aes = FALSE, size = 2.8, colour = "grey30") +
  scale_fill_manual(values = P4_MAL, name = "malignancy_refined", drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.07))) +
  labs(title = "Refined malignancy composition per patient", x = NULL,
       y = "fraction of the patient's cells",
       subtitle = "sample_id = patient = dataset, so patient and batch cannot be separated") +
  p4_theme()
p6b <- ggplot(p2f, aes(x = stage, y = frac, group = sample_id, colour = sample_id)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.4) +
  scale_colour_manual(values = P4_SAMPLE, name = "patient") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(title = "Malignant fraction: Phase 2 → SCEVAN → refined",
       subtitle = "Each line is one patient", x = NULL, y = "malignant fraction") +
  p4_theme()
p6 <- if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork); (p6a | p6b) + plot_annotation(caption = P4_CAPTION_CNV,
    theme = theme(plot.caption = element_text(size = 8, colour = "grey35", hjust = 0))) } else p6a
save_fig(p6, file.path(FIG, "06_malignancy_by_patient.pdf"), 13.5, 5.6)

# ---- 07 ambiguous / disputed populations --------------------------------
AUD <- c("Fibroblast","Candidate-Malignant-Unresolved","Pericyte-VSMC",
         "MPNST-Tumor","Uncertain")
d7 <- md |> filter(annotation_ccc %in% AUD) |>
  count(annotation_ccc, sample_id, malignancy_scevan) |>
  group_by(annotation_ccc, sample_id) |> mutate(frac = n/sum(n), tot = sum(n)) |> ungroup() |>
  mutate(annotation_ccc = factor(as.character(annotation_ccc), levels = AUD))
p7a <- ggplot(d7, aes(x = sample_id, y = frac, fill = malignancy_scevan)) +
  geom_col(width = 0.75) +
  geom_text(data = distinct(d7, annotation_ccc, sample_id, tot),
            aes(x = sample_id, y = 1.03, label = tot), inherit.aes = FALSE,
            size = 2.4, colour = "grey30") +
  facet_wrap(~ annotation_ccc, nrow = 1) +
  scale_fill_manual(values = P4_SCEVAN, name = "SCEVAN call") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.09))) +
  labs(title = "The populations Phase 2 could not resolve — SCEVAN call per patient",
       subtitle = paste("Fibroblast and Pericyte-VSMC were NEVER used as normal references, and",
                        "Candidate-Malignant-Unresolved\nand Uncertain were excluded from Phase 3.",
                        "n above each bar."),
       x = NULL, y = "fraction of population in that patient") +
  p4_theme() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
d7b <- md |> filter(annotation_ccc %in% AUD) |>
  mutate(annotation_ccc = factor(as.character(annotation_ccc), levels = AUD)) |>
  filter(!is.na(cnv_burden))
p7b <- ggplot(d7b, aes(x = annotation_ccc, y = cnv_burden, fill = malignancy_scevan)) +
  geom_boxplot(outlier.size = 0.25, outlier.alpha = 0.3, linewidth = 0.32,
               position = position_dodge(preserve = "single")) +
  scale_fill_manual(values = P4_SCEVAN, name = "SCEVAN call") +
  labs(title = "Inferred CNA burden in the same populations",
       subtitle = paste0("Burden = fraction of assessed genes with |relative CNA| > 0.10 ",
                         "(gene-density-weighted, not base-pair-weighted)"),
       x = NULL, y = "CNA burden") +
  p4_theme() + theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 7.5))
p7 <- if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork); (p7a / p7b) + plot_layout(heights = c(1, 0.95)) +
    plot_annotation(caption = P4_CAPTION_CNV,
      theme = theme(plot.caption = element_text(size = 8, colour = "grey35", hjust = 0))) } else p7a
save_fig(p7, file.path(FIG, "07_ambiguous_population_malignancy.pdf"), 12.5, 9.6)

# ---- 08 CNA profiles along the genome, by annotation ---------------------
cat("\nbuilding genome-ordered mean CNA profiles\n")
prof <- list()
for (s in SAMPLES) {
  fm <- file.path(SC, s, "primary", "output", paste0(s, "_primary_CNAmtx.RData"))
  fa <- file.path(SC, s, "primary", "output", paste0(s, "_primary_count_mtx_annot.RData"))
  if (!file.exists(fm) || !file.exists(fa)) { cat("  skip", s, "\n"); next }
  e <- new.env(); load(fm, envir = e)
  nm <- ls(e)[vapply(ls(e), function(n) is.matrix(get(n, e)), logical(1))][1]
  M <- get(nm, e)
  a <- new.env(); load(fa, envir = a); ann <- get("count_mtx_annot", a)
  ann <- ann[!duplicated(ann$gene_name), ]
  rownames(ann) <- ann$gene_name
  g <- intersect(rownames(M), rownames(ann))
  M <- M[g, , drop = FALSE]; ann <- ann[g, , drop = FALSE]
  o <- order(as.numeric(as.character(ann$seqnames)), as.numeric(as.character(ann$start)))
  M <- M[o, , drop = FALSE]; ann <- ann[o, , drop = FALSE]
  ii <- match(colnames(M), md$cell_id)
  grp <- as.character(md$annotation_ccc)[ii]
  keep <- !is.na(grp) & grp %in% AUD
  if (!any(keep)) { rm(M); next }
  agg <- vapply(split(which(keep), grp[keep]),
                function(idx) rowMeans(M[, idx, drop = FALSE]), numeric(nrow(M)))
  prof[[s]] <- data.frame(sample_id = s, chr = as.numeric(as.character(ann$seqnames)),
                          idx = seq_len(nrow(M)), agg, check.names = FALSE)
  rm(M, e, a); invisible(gc(FALSE))
}
if (length(prof)) {
  pl <- bind_rows(lapply(prof, function(d)
    pivot_longer(d, cols = -c(sample_id, chr, idx),
                 names_to = "annotation_ccc", values_to = "cna"))) |>
    mutate(annotation_ccc = factor(annotation_ccc, levels = AUD))
  # chromosome boundaries, computed per sample because gene sets differ
  bnd <- pl |> distinct(sample_id, chr, idx) |> group_by(sample_id, chr) |>
    summarise(mid = mean(idx), end = max(idx), .groups = "drop")
  p8 <- ggplot(pl, aes(x = idx, y = cna, colour = annotation_ccc)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_vline(data = bnd, aes(xintercept = end), colour = "grey88", linewidth = 0.22) +
    geom_line(linewidth = 0.34) +
    facet_grid(annotation_ccc ~ sample_id, scales = "free_x") +
    scale_colour_brewer(palette = "Dark2", guide = "none") +
    scale_x_continuous(breaks = bnd$mid[bnd$sample_id == bnd$sample_id[1]],
                       labels = bnd$chr[bnd$sample_id == bnd$sample_id[1]],
                       expand = expansion(0)) +
    labs(title = "Mean inferred CNA profile along chromosomes 1–22, by population and patient",
         subtitle = paste("Genes ordered by genomic position. Deviation from zero indicates inferred",
                          "gain (up) or loss (down).\nChromosomes X and Y are absent: SCEVAN's",
                          "annotation covers chr1–22 only."),
         x = "chromosome", y = "mean relative CNA", caption = P4_CAPTION_CNV) +
    p4_theme(9) + theme(axis.text.x = element_text(size = 5.2),
                        panel.spacing.x = unit(0.12, "lines"))
  save_fig(p8, file.path(FIG, "08_cnv_profiles_by_annotation.pdf"), 15.5, 9.8)
  write.table(pl, gzfile(file.path(TAB, "CNA_PROFILES_BY_ANNOTATION.tsv.gz")),
              sep = "\t", quote = FALSE, row.names = FALSE)
} else cat("  no CNA matrices available; figure 08 skipped\n")

cat("\nM30 FIGURES COMPLETE\n")
