#!/usr/bin/env Rscript
# =============================================================================
# Phase 4 - M35A - FIGURE 5
#   25_MPNST3_scevan_failure_qc.pdf
#
# MPNST_3 is shown failing, not quietly dropped. A reader is entitled to see
# why one of four patients was barred from contributing malignant promotions.
# =============================================================================
source("scripts/phase4/figures/m35a_common.R")
suppressPackageStartupMessages(library(jsonlite))

md <- m35a_load_meta()
m35a_msg("== M35A figure 25: MPNST_3 SCEVAN failure QC ==")
m35a_validate(md)

# ---- A. primary vs sensitivity agreement, from the frozen M29 summaries -----
agr <- do.call(rbind, lapply(M35A_SAMPLES, function(s) {
  j <- jsonlite::fromJSON(file.path(M35A_SCEVAN, s,
        sprintf("m29_%s_summary.json", s)))
  data.frame(sample_id = s, agreement = j$concordance$agreement,
             cells_common = j$concordance$cells_common)
}))
agr$sample_id <- factor(agr$sample_id, levels = M35A_SAMPLES)
agr$reliable <- agr$sample_id %in% M35A_RELIABLE
m35a_msg("-- primary vs sensitivity agreement --"); print(agr)
stopifnot(all(abs(agr$agreement - M35A_FROZEN$agreement[as.character(agr$sample_id)]) < 1e-6))

pA <- ggplot(agr, aes(x = sample_id, y = agreement, fill = reliable)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = sprintf("%.4f", agreement)), vjust = -0.45, size = 3.2,
            fontface = "bold") +
  geom_hline(yintercept = 0.9, linetype = 2, linewidth = 0.35, colour = "grey40") +
  scale_fill_manual(values = c(`TRUE` = "#2166AC", `FALSE` = "#B2182B"),
                    labels = c(`TRUE` = "reliable", `FALSE` = "FAILED sanity gate"),
                    name = NULL) +
  scale_y_continuous(limits = c(0, 1.09), expand = c(0, 0),
                     labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "agreement", title = "A. Primary vs sensitivity run agreement",
       subtitle = wrap_sub(paste0("Primary: norm_cell = NULL. Sensitivity: immune cells as norm_cell, FIXED_NORMAL_CELLS = FALSE.\n",
                         "MPNST_3's two runs are not merely noisy - they are INVERTED (0.0864)."), 62)) +
  p4_theme(10) + theme(legend.position = "bottom")

# ---- B. MPNST_3 clone composition ------------------------------------------
c3 <- md |> filter(sample_id == "MPNST_3", !is.na(scevan_clone), scevan_clone != "")
sz3 <- c3 |> count(scevan_clone, name = "clone_size")
comp3 <- c3 |> count(scevan_clone, annotation_ccc, name = "n") |>
  left_join(sz3, by = "scevan_clone") |> mutate(frac = n / clone_size)
lv3 <- sz3 |> arrange(desc(clone_size)) |> pull(scevan_clone)
comp3$scevan_clone <- factor(comp3$scevan_clone, levels = lv3)
lab3 <- setNames(sprintf("%s\n(n=%s)", sub("^MPNST_3_clone", "clone ", lv3),
                         format(sz3$clone_size[match(lv3, sz3$scevan_clone)],
                                big.mark = ",")), lv3)
ann_lv <- comp3 |> group_by(annotation_ccc) |> summarise(n = sum(n)) |>
  arrange(desc(n)) |> pull(annotation_ccc)
comp3$annotation_ccc <- factor(comp3$annotation_ccc, levels = ann_lv)
pal3 <- setNames(scales::hue_pal()(length(ann_lv)), ann_lv)
for (k in intersect(names(M35A_ANN_COL), ann_lv)) pal3[k] <- M35A_ANN_COL[[k]]

pB <- ggplot(comp3, aes(x = scevan_clone, y = frac, fill = annotation_ccc)) +
  geom_col(width = 0.78, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = pal3, name = "Phase 2 annotation_ccc") +
  scale_x_discrete(labels = lab3) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(x = NULL, y = "fraction of clone",
       title = "B. MPNST_3 'clones' are canonical immune lineages",
       subtitle = wrap_sub(paste0("clone 3 predominantly T/NK; clone 2 plasma / pDC / B; clone 1 plasma-dominated.\n",
                         "A CNA-defined tumour partition should not reproduce the immune compartment."), 78)) +
  p4_theme(10)

# ---- C. malignant fraction of canonical immune populations ------------------
pf <- read.delim(file.path(M35A_TAB, "SCEVAN_POPULATION_MALIGNANT_FRACTION.tsv"),
                 check.names = FALSE)
imm <- pf |> filter(annotation_ccc %in% M35A_OTHER_IMMUNE, n_assessed >= 20) |>
  mutate(sample_id = factor(sample_id, levels = M35A_SAMPLES),
         annotation_ccc = factor(annotation_ccc, levels = M35A_OTHER_IMMUNE))
pC <- ggplot(imm, aes(x = annotation_ccc, y = pop_frac, fill = sample_id)) +
  geom_col(position = position_dodge2(preserve = "single", padding = 0.12),
           width = 0.85) +
  geom_hline(yintercept = 0.25, linetype = 2, linewidth = 0.4, colour = "#B2182B") +
  annotate("text", x = 0.6, y = 0.28, label = "a priori 25% sanity gate",
           hjust = 0, size = 2.7, colour = "#B2182B") +
  scale_fill_manual(values = P4_SAMPLE, name = "patient") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.03),
                     expand = c(0, 0)) +
  labs(x = NULL, y = "SCEVAN malignant fraction",
       title = "C. SCEVAN malignant fraction of canonical immune populations (primary run)",
       subtitle = wrap_sub(paste0("Populations with >= 20 assessed cells. MPNST_3 calls B cells, CD4-T, CD8-T, NK, T-other and pDC ",
                         "essentially 100% malignant - 6 of 9 populations fail on rate AND on breadth. ",
                         "Plasma cells are excluded from the sanity denominator by Amendment A2 (immunoglobulin-locus artefact); ",
                         "they are shown separately in the handoff, not here."), 150)) +
  p4_theme(10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# ---- D. depth context, explicitly NOT offered as the cause ------------------
dep <- md |> group_by(sample_id) |>
  summarise(median_genes = median(nFeature_RNA),
            median_umi = median(nCount_RNA), n = n(), .groups = "drop")
m35a_msg("-- depth --"); print(as.data.frame(dep))
pD <- ggplot(md, aes(x = sample_id, y = nFeature_RNA, fill = sample_id)) +
  geom_violin(scale = "width", linewidth = 0.2, alpha = 0.75, colour = "grey30") +
  geom_boxplot(width = 0.13, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  geom_text(data = dep, aes(x = sample_id, y = max(md$nFeature_RNA),
            label = sprintf("median %s", format(median_genes, big.mark = ","))),
            inherit.aes = FALSE, size = 2.8, colour = "grey25", vjust = -0.6) +
  scale_fill_manual(values = P4_SAMPLE, guide = "none") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0.02, 0.09))) +
  labs(x = NULL, y = "detected genes per cell (nFeature_RNA)",
       title = "D. Sequencing depth - CONTEXT ONLY",
       subtitle = wrap_sub(paste0("MPNST_3 has the lowest median depth of the four, and only 25 confident normal cells were found in it. ",
                         "Low depth alone is NOT claimed to have caused the failure: depth is a plausible contributor, not a demonstrated cause."), 150)) +
  p4_theme(10)

p <- (pA | pB) / pC / pD + plot_layout(heights = c(1.15, 1, 1)) +
  plot_annotation(
    title = "Figure 25. MPNST_3 - a failed SCEVAN partition, shown rather than silently omitted",
    subtitle = wrap_sub(paste0(
      "MPNST_3 malignant promotions were disabled because SCEVAN's inferred tumour partitions corresponded predominantly to ",
      "canonical immune lineages and the primary and sensitivity classifications were incompatible.\n",
      "Cost stated plainly: any malignant fibroblast-like population that exists in MPNST_3 cannot be detected by this analysis. ",
      "The gate is deliberately asymmetric - an unreliable run may still support a Non-malignant call, which always needs concordant lineage evidence too."), 170),
    caption = M35A_CAPTION,
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 9.5, colour = "grey20"),
                  plot.caption = element_text(size = 7.5, colour = "grey35", hjust = 0)))
save_fig(p, file.path(M35A_FIG, "25_MPNST3_scevan_failure_qc.pdf"), 15.8, 15.5)
m35a_msg("== figure 25 done ==")
