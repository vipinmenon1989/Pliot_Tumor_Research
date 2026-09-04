#!/usr/bin/env Rscript
# =============================================================================
# Phase 4 - M35A - MAIN INTEGRATED FIGURE
#   28_phase4_scevan_evidence_summary.pdf
#
# One page, one result:
#   CNA-based analysis reveals that a substantial fibroblast-like transcriptional
#   compartment belongs to CNA-defined malignant MPNST populations, while a
#   smaller fibroblast population remains non-malignant and discordant cells are
#   explicitly retained as Ambiguous.
# SCEVAN provides inferred copy number. It does NOT provide DNA-level proof.
# =============================================================================
source("scripts/phase4/figures/m35a_common.R")
suppressPackageStartupMessages({ library(png); library(grid); library(ggalluvial) })

md <- m35a_load_meta()
m35a_msg("== M35A figure 28: integrated evidence summary ==")
m35a_validate(md)

UX <- "postint_umap_harmony_1"; UY <- "postint_umap_harmony_2"
md$ann_grp <- m35a_group_ann(md$annotation_ccc)

# ---- A. Phase 2 annotation UMAP --------------------------------------------
dA <- md[order(md$ann_grp %in% c("Fibroblast", "Candidate-Malignant-Unresolved",
                                 "MPNST-Tumor")), ]
pA <- ggplot(dA, aes(.data[[UX]], .data[[UY]], colour = ann_grp)) +
  geom_point(size = 0.30, alpha = 0.6, stroke = 0) +
  scale_colour_manual(values = M35A_ANN_COL, name = NULL, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 2.4, alpha = 1), ncol = 2)) +
  coord_equal() +
  labs(title = "A. Phase 2 marker annotation", x = "UMAP 1", y = "UMAP 2",
       subtitle = wrap_sub("Frozen annotation_ccc, grouped. Marker evidence split the disputed cells three ways.", 90)) +
  p4_theme(9) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom", legend.text = element_text(size = 7))

# ---- B. refined malignancy UMAP --------------------------------------------
dB <- md[order(md$malignancy_refined != "Malignant"), ]
pB <- ggplot(dB, aes(.data[[UX]], .data[[UY]], colour = malignancy_refined)) +
  geom_point(size = 0.30, alpha = 0.6, stroke = 0) +
  scale_colour_manual(values = P4_MAL, name = NULL, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 2.4, alpha = 1), ncol = 2)) +
  coord_equal() +
  labs(title = "B. Phase 4 copy-number-refined malignancy",
       x = "UMAP 1", y = "UMAP 2",
       subtitle = wrap_sub("6,434 Malignant | 9,078 Non-malignant | 3,766 Ambiguous | 438 Excluded. The UMAP is for display only; classification came from raw counts.", 90)) +
  p4_theme(9) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom", legend.text = element_text(size = 7))

# ---- H. the disputed compartment on its own --------------------------------
dH <- md |> filter(annotation_ccc == "Fibroblast")
dH <- dH[order(dH$malignancy_refined != "Malignant"), ]
pH <- ggplot(md, aes(.data[[UX]], .data[[UY]])) +
  geom_point(colour = "grey88", size = 0.22, alpha = 0.45, stroke = 0) +
  geom_point(data = dH, aes(colour = malignancy_refined), size = 0.34,
             alpha = 0.75, stroke = 0) +
  scale_colour_manual(values = P4_MAL, name = NULL, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 2.4, alpha = 1), ncol = 2)) +
  coord_equal() +
  labs(title = "C. The 5,064 Phase-2 Fibroblast cells alone",
       x = "UMAP 1", y = "UMAP 2",
       subtitle = wrap_sub("4,036 Malignant, 908 Non-malignant, 120 Ambiguous. All other cells in grey.", 90)) +
  p4_theme(9) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom", legend.text = element_text(size = 7))

# ---- C. native SCEVAN CNA / subclone architecture --------------------------
read_native <- function(path, pad = 6L) {
  x <- png::readPNG(path)
  if (length(dim(x)) == 2L) x <- array(x, c(dim(x), 1L))
  ink <- apply(x[, , 1:min(3L, dim(x)[3]), drop = FALSE], c(1, 2),
               function(v) any(v < 0.985))
  rr <- range(which(rowSums(ink) > 0)); cc <- range(which(colSums(ink) > 0))
  x[max(1L, rr[1] - pad):min(dim(x)[1], rr[2] + pad),
    max(1L, cc[1] - pad):min(dim(x)[2], cc[2] + pad), , drop = FALSE]
}
hm <- file.path(M35A_SCEVAN, "MPNST_2", "primary", "output",
                "MPNST_2_primaryheatmap_subclones.png")
pC <- wrap_elements(full = rasterGrob(read_native(hm), interpolate = TRUE,
                                      width = unit(1, "npc"), height = unit(1, "npc"))) +
  labs(title = "D. Native SCEVAN CNA / subclone architecture (MPNST_2, primary run)",
       subtitle = wrap_sub(paste0("SCEVAN output embedded verbatim. Rows = cells grouped by CNA-defined subclone (left bar); columns = genes, chr1-22. ",
                         "All four MPNST_2 subclones are fibroblast-dominated. Full set: figures 17-20."), 165)) +
  theme(plot.title = element_text(face = "bold", size = 10, margin = margin(b = 3)),
        plot.subtitle = element_text(size = 8, colour = "grey25", margin = margin(b = 8)),
        plot.margin = margin(t = 12, r = 4, b = 4, l = 4))

# ---- D. clone composition by Phase 2 identity ------------------------------
cl <- md |> filter(!is.na(scevan_clone), scevan_clone != "",
                   sample_id %in% M35A_RELIABLE) |>
  mutate(ann = m35a_group_ann(annotation_ccc))
sz <- cl |> count(sample_id, scevan_clone, name = "clone_size")
compD <- cl |> count(sample_id, scevan_clone, ann, name = "n") |>
  left_join(sz, by = c("sample_id", "scevan_clone")) |>
  mutate(frac = n / clone_size,
         clone_label = sub("^MPNST_[0-9]+_clone", "clone ", scevan_clone))
ordD <- sz |> arrange(sample_id, desc(clone_size)) |>
  mutate(clone_label = sub("^MPNST_[0-9]+_clone", "clone ", scevan_clone))
compD$clone_label <- factor(compD$clone_label, levels = unique(ordD$clone_label))
compD$sample_id <- factor(compD$sample_id, levels = M35A_RELIABLE)
fibdom <- cl |> count(sample_id, scevan_clone, annotation_ccc) |>
  group_by(sample_id, scevan_clone) |> slice_max(n, n = 1, with_ties = FALSE) |>
  ungroup() |> filter(annotation_ccc == "Fibroblast") |>
  mutate(clone_label = factor(sub("^MPNST_[0-9]+_clone", "clone ", scevan_clone),
                              levels = levels(compD$clone_label)),
         sample_id = factor(sample_id, levels = M35A_RELIABLE))

pD <- ggplot(compD, aes(clone_label, frac, fill = ann)) +
  geom_col(width = 0.82, colour = "white", linewidth = 0.15) +
  geom_point(data = fibdom, aes(x = clone_label, y = -0.05), inherit.aes = FALSE,
             shape = 17, size = 1.9, colour = "#D95F02") +
  scale_fill_manual(values = M35A_ANN_COL, name = NULL, drop = FALSE) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(-0.08, 1),
                     expand = c(0, 0)) +
  facet_grid(~ sample_id, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "fraction of clone",
       title = "E. Composition of each CNA-defined clone by Phase 2 identity - the most direct evidence",
       subtitle = wrap_sub(paste0("SCEVAN built these clones from copy number alone and never saw a Phase 2 label. ",
                         "Orange triangle = fibroblast-dominated. MPNST_2: 4 of 4. MPNST_4: 7 of 8. MPNST_1: all 7 clones mix Candidate-Malignant-Unresolved with MPNST-Tumor, 6 of 7 also carry Fibroblast."), 165)) +
  p4_theme(9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "bottom", legend.text = element_text(size = 7))

# ---- E. fibroblast malignant vs non-malignant CNA burden -------------------
gE <- md |>
  mutate(group = case_when(
    annotation_ccc == "Fibroblast"  & malignancy_refined == "Malignant"     ~ "Fib -> Mal",
    annotation_ccc == "Fibroblast"  & malignancy_refined == "Non-malignant" ~ "Fib -> Non-mal",
    annotation_ccc == "MPNST-Tumor" & malignancy_refined == "Malignant"     ~ "MPNST-Tumor -> Mal",
    annotation_ccc %in% M35A_OTHER_IMMUNE                                   ~ "immune ctrl",
    TRUE ~ NA_character_)) |>
  filter(!is.na(group), !is.na(cnv_burden), sample_id %in% M35A_RELIABLE) |>
  mutate(group = factor(group, levels = c("Fib -> Mal", "Fib -> Non-mal",
                                          "MPNST-Tumor -> Mal", "immune ctrl")),
         sample_id = factor(sample_id, levels = M35A_RELIABLE))
nE <- gE |> count(sample_id, group)
pE <- ggplot(gE, aes(group, cnv_burden, fill = group)) +
  geom_violin(scale = "width", linewidth = 0.18, alpha = 0.75, colour = "grey35") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.28) +
  geom_text(data = nE, aes(x = group, y = -0.03, label = n), inherit.aes = FALSE,
            size = 2.1, colour = "grey35") +
  scale_fill_manual(values = c("Fib -> Mal" = "#D95F02", "Fib -> Non-mal" = "#2166AC",
                               "MPNST-Tumor -> Mal" = "#B2182B", "immune ctrl" = "grey60"),
                    guide = "none") +
  facet_wrap(~ sample_id, nrow = 1) +
  labs(x = NULL, y = "cnv_burden",
       title = "F. Fibroblast CNA burden, malignant vs non-malignant",
       subtitle = wrap_sub(paste0("Per patient, never pooled. Separation is clear in MPNST_1 and MPNST_2; ",
                         "MPNST_4 retains 1 non-malignant fibroblast, so its contrast is NOT EVALUABLE. n below each violin."), 105)) +
  p4_theme(9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7))

# ---- F. Phase 2 -> Phase 4 transition --------------------------------------
DISPUTED <- c("Fibroblast", "Candidate-Malignant-Unresolved", "MPNST-Tumor",
              "Pericyte-VSMC", "Endothelial", "Uncertain")
dF <- md |> filter(annotation_ccc %in% DISPUTED) |>
  count(annotation_ccc, malignancy_refined, name = "n") |>
  mutate(annotation_ccc = factor(annotation_ccc, levels = DISPUTED))
annF <- M35A_ANN_COL[c("Fibroblast", "Candidate-Malignant-Unresolved", "MPNST-Tumor",
                       "Pericyte-VSMC", "Endothelial")]; annF["Uncertain"] <- "grey60"
stratum_mid <- function(df, key) {
  tot <- df |> group_by(.data[[key]]) |> summarise(n = sum(n), .groups = "drop") |>
    arrange(.data[[key]])
  top <- sum(tot$n) - c(0, cumsum(tot$n)[-nrow(tot)])
  data.frame(stratum = as.character(tot[[key]]), n = tot$n, y = top - tot$n / 2)
}
labF <- rbind(stratum_mid(dF, "annotation_ccc") |> mutate(x = 1 - 0.17, hj = 1),
              stratum_mid(dF, "malignancy_refined") |> mutate(x = 2 + 0.17, hj = 0)) |>
  filter(n >= 1000) |>   # the small strata are labelled in full in figure 24
  mutate(txt = sprintf("%s\nn = %s", stratum, format(n, big.mark = ",")))

pF <- ggplot(dF, aes(y = n, axis1 = annotation_ccc, axis2 = malignancy_refined)) +
  geom_alluvium(aes(fill = annotation_ccc), width = 0.22, alpha = 0.72,
                knot.pos = 0.32, curve_type = "sigmoid") +
  geom_stratum(width = 0.22, fill = "grey96", colour = "grey35", linewidth = 0.25) +
  geom_text(data = labF, inherit.aes = FALSE,
            aes(x = x, y = y, label = txt, hjust = hj), size = 2.5,
            lineheight = 0.95, colour = "grey10") +
  geom_text(data = dF |> filter(n >= 300), stat = "alluvium",
            aes(label = format(n, big.mark = ",")), size = 2.4, colour = "grey10",
            fontface = "bold") +
  scale_fill_manual(values = annF, guide = "none") +
  scale_x_continuous(breaks = c(1, 2), labels = c("Phase 2", "Phase 4"),
                     limits = c(0.28, 2.72)) +
  scale_y_continuous(labels = comma) +
  labs(x = NULL, y = "cells",
       title = "G. Phase 2 identity to Phase 4 malignancy (disputed compartments)",
       subtitle = wrap_sub(paste0("Immune populations excluded so they cannot bury the transition. ",
         "Strata of >= 1,000 cells and flows of >= 300 cells are labelled; full detail in figure 24. Ambiguous is terminal and never becomes tumour."), 105)) +
  p4_theme(9) + theme(panel.grid.major.x = element_blank())

# ---- G. threshold robustness -----------------------------------------------
ts <- read.delim(file.path(M35A_TAB, "MALIGNANCY_THRESHOLD_SENSITIVITY.tsv"),
                 check.names = FALSE)
dG <- ts |> select(pop_frac_low,
                   `total refined Malignant` = Malignant,
                   `retained MPNST-Tumor` = mpnst_tumor_still_malignant,
                   `Fibroblast -> Malignant` = fibroblast_malignant) |>
  pivot_longer(-pop_frac_low, names_to = "series", values_to = "cells")
pG <- ggplot(dG, aes(pop_frac_low, cells, colour = series)) +
  geom_vline(xintercept = 0.25, linetype = 2, linewidth = 0.35, colour = "grey45") +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  geom_text(aes(label = format(cells, big.mark = ",")), vjust = -1.0, size = 2.3,
            show.legend = FALSE) +
  scale_colour_manual(values = c("total refined Malignant" = "#444444",
                                 "retained MPNST-Tumor" = "#B2182B",
                                 "Fibroblast -> Malignant" = "#D95F02"), name = NULL) +
  scale_x_continuous(breaks = ts$pop_frac_low) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0.1, 0.2))) +
  labs(x = "pop_frac_low", y = "cells",
       title = "H. Threshold robustness - what is solid and what is not",
       subtitle = wrap_sub(paste0("Fibroblast -> Malignant is flat at 4,036 across 0.15-0.40. The retained MPNST-Tumor count, and with it the ",
                         "32.63% cohort fraction, is threshold-sensitive and is not patient-robust."), 160)) +
  p4_theme(9) + theme(legend.position = "bottom", legend.text = element_text(size = 7))

# ---- assemble ---------------------------------------------------------------
design <- "
AABBCC
DDDDDD
EEEEEE
FFFGGG
HHHHHH
"
# patchwork assigns plots to design areas in ALPHABETICAL order of the area
# letters, NOT by variable name, so the chain below is written in the order the
# areas A..H appear on the page: three UMAPs, native CNA, clone composition,
# burden, transition, threshold.
p <- pA + pB + pH + pC + pD + pE + pF + pG +
  plot_layout(design = design, heights = c(1.35, 1.35, 1.05, 1.30, 0.95)) +
  plot_annotation(
    title = "Figure 28. Phase 4 SCEVAN evidence summary - a fibroblast-like compartment inside the CNA-defined malignant population",
    subtitle = wrap_sub(paste0(
      "Copy-number-based analysis places a substantial fibroblast-like transcriptional compartment inside CNA-defined malignant MPNST populations: ",
      "4,036 of 5,064 Phase-2 Fibroblast cells are refined Malignant, majority-malignant in 3 of 4 patients and threshold-independent.\n",
      "A smaller fibroblast population (908 cells) remains Non-malignant, and cells where lineage and copy-number evidence disagree are retained as ",
      "Ambiguous (3,766) rather than forced either way.\n",
      "Phase 2 annotation is not overturned - it found the Schwann-like and cycling malignant cells and could not see the mesenchymal/ECM ones, ",
      "because an ECM programme is what a fibroblast looks like."), 175),
    caption = paste0(M35A_CAPTION,
      "\nMPNST_3 is excluded from malignant promotions throughout (failed immune sanity gate; figure 25). ",
      "The 32.63% cohort malignant fraction is NOT patient-robust; the fibroblast conclusion is a separate and robust claim."),
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 9.5, colour = "grey20"),
                  plot.caption = element_text(size = 7.5, colour = "grey35", hjust = 0)))

save_fig(p, file.path(M35A_FIG, "28_phase4_scevan_evidence_summary.pdf"), 17.0, 24.0)
m35a_msg("== figure 28 done ==")
