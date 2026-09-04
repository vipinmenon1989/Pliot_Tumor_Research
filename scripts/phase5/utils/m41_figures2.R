#!/usr/bin/env Rscript
# Phase 5 final figures 08, 09, 12 - the malignant-ECM comparison and the model.
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr)
                                 library(patchwork); library(scales) })
source("scripts/phase5/utils/phase5_common.R")
source("scripts/phase5/utils/phase5_plot_utils.R")
set.seed(42)
args <- commandArgs(trailingOnly = TRUE); KSTAR <- as.integer(args[1])
rd <- function(p) read.delim(p, check.names = FALSE)

prog <- rd(file.path(P5_TAB, "MALIGNANT_PROGRAMS.tsv"))
PROGS <- prog$program_id
LABS <- setNames(sprintf("%s %s", prog$program_id, prog$program_label), PROGS)
ev   <- rd(file.path(P5_ECM, "M39_EVALUABILITY_BY_PATIENT.tsv"))
full <- rd(file.path(P5_TAB, "MALIGNANT_ECM_VS_TRUE_FIBROBLAST.tsv"))
sig  <- rd(file.path(P5_TAB, "MALIGNANT_ECM_SIGNATURE.tsv"))
usg  <- rd(file.path(P5_ECM, "M39_PROGRAM_USAGE_A_VS_B.tsv"))
cnv  <- rd(file.path(P5_ECM, "M39_CNV_METRICS_BY_GROUP.tsv"))
amb  <- rd(file.path(P5_TAB, "AMBIGUOUS_FIBROBLAST_PROGRAM_PROJECTION.tsv"))
rob  <- rd(file.path(P5_TAB, "PHASE5_PROGRAM_ROBUSTNESS.tsv"))
EVAL <- as.character(ev$sample_id[ev$evaluable])

# ---- 08 malignant ECM vs true fibroblast ------------------------------------
p8a <- ev |>
  pivot_longer(c(A_malignant_ECM, B_true_fibroblast), names_to = "group",
               values_to = "n") |>
  ggplot(aes(sample_id, n, fill = group)) +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8) +
  geom_text(aes(label = n), position = position_dodge2(width = 0.8),
            vjust = -0.3, size = 2.7) +
  geom_hline(yintercept = 30, linetype = 2, linewidth = 0.4, colour = "#B2182B") +
  geom_text(data = ev |> filter(!evaluable),
            aes(x = sample_id, y = 1, label = "NOT EVALUABLE"), inherit.aes = FALSE,
            angle = 90, hjust = 0, size = 2.8, colour = "grey25") +
  scale_fill_manual(values = c(A_malignant_ECM = "#D95F02",
                               B_true_fibroblast = "#2166AC"),
                    labels = c("A: malignant ECM-like", "B: true fibroblast"),
                    name = NULL) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 10, 30, 100, 300, 1000, 3000)) +
  labs(x = NULL, y = "cells (log1p scale)",
       title = "A. evaluability, checked before any test was run",
       subtitle = wrap5(paste0(
         "Red line = the declared minimum of 30 cells on BOTH sides. Only ", paste(EVAL, collapse = " and "),
         " qualify: MPNST_4 retains a single non-malignant fibroblast and MPNST_3 contributes no malignant ",
         "fibroblasts because Phase 4 disabled its promotions. A four-patient paired test is therefore not constructed."), 78)) +
  p5_theme(10)

lf <- full |> select(gene, all_of(EVAL)) |>
  rename_with(~paste0("lfc_", .x), all_of(EVAL))
au <- full |> select(gene, all_of(paste0("auc_", EVAL)))
# `is_sig` rather than `sig`: naming the new column `sig` would shadow the
# signature data frame inside the same mutate() call.
d8 <- lf |> left_join(au, by = "gene") |>
  mutate(is_sig = gene %in% sig$gene,
         dir = sig$direction[match(gene, sig$gene)])
p8b <- ggplot(d8, aes(.data[[paste0("lfc_", EVAL[1])]],
                      .data[[paste0("lfc_", EVAL[2])]])) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_point(data = ~subset(.x, !is_sig), colour = "grey82", size = 0.5, stroke = 0) +
  geom_point(data = ~subset(.x, is_sig), aes(colour = dir), size = 1.1, stroke = 0) +
  scale_colour_manual(values = c(up_in_malignant_ECM = "#D95F02",
                                 up_in_true_fibroblast = "#2166AC"),
                      na.value = "grey70", name = "signature direction") +
  labs(x = sprintf("pseudobulk log2FC (A vs B), %s", EVAL[1]),
       y = sprintf("pseudobulk log2FC (A vs B), %s", EVAL[2]),
       title = "B. the two evaluable patients agree on direction",
       subtitle = wrap5(paste0("Each point is a gene. Signature genes must move the same way, by >= 1 log2 unit, ",
         "with cell-level AUC >= 0.65, in EVERY evaluable patient."), 78)) +
  p5_theme(10)

pu <- usg |> filter(evaluable) |>
  mutate(plab = factor(LABS[program], levels = LABS[PROGS]))
p8c <- ggplot(pu, aes(plab, delta_median, fill = sample_id)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8) +
  scale_fill_manual(values = P5_SAMPLE, name = "patient") +
  labs(x = NULL, y = "median program usage\n(malignant ECM-like - true fibroblast)",
       title = "C. programs separate the two groups, per patient") +
  p5_theme(10) + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7.5))

p8d <- cnv |> filter(group != "C_ambiguous") |>
  ggplot(aes(sample_id, cnv_burden, fill = group)) +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8) +
  geom_text(aes(label = n), position = position_dodge2(width = 0.8),
            vjust = -0.3, size = 2.5) +
  scale_fill_manual(values = c(A_malignant_ECM = "#D95F02",
                               B_true_fibroblast = "#2166AC"), name = NULL) +
  labs(x = NULL, y = "median Phase 4 cnv_burden",
       title = "D. Phase 4 CNA burden, shown as context",
       subtitle = wrap5(paste0("These are the frozen Phase 4 metrics that defined the groups. They are context, ",
         "NOT independent validation of the split they produced."), 78),
       caption = P5_CAPTION) + p5_theme(10)
save_fig5((p8a | p8b) / (p8c | p8d),
          file.path(P5_FIG, "08_malignant_ecm_vs_true_fibroblast.pdf"), 15.0, 11.0)

# ---- 09 signature heatmap ----------------------------------------------------
top_up <- sig |> filter(direction == "up_in_malignant_ECM") |>
  arrange(desc(mean_lfc)) |> head(30)
top_dn <- sig |> filter(direction == "up_in_true_fibroblast") |>
  arrange(mean_lfc) |> head(30)
sg <- bind_rows(top_up, top_dn)
if (nrow(sg)) {
  h <- sg |> select(gene, direction, all_of(EVAL)) |>
    pivot_longer(all_of(EVAL), names_to = "sample_id", values_to = "lfc") |>
    mutate(gene = factor(gene, levels = rev(sg$gene)))
  p9 <- ggplot(h, aes(sample_id, gene, fill = lfc)) +
    geom_tile(colour = "white", linewidth = 0.15) +
    facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = "#2166AC", mid = "grey96", high = "#D95F02",
                         midpoint = 0, name = "pseudobulk\nlog2FC\n(A vs B)") +
    labs(x = NULL, y = NULL,
         title = "Figure 09. Transparent malignant-ECM signature",
         subtitle = wrap5(paste0(
           "Top 30 genes per direction. Declared criteria, applied in EVERY evaluable patient: |pseudobulk log2FC| >= 1, ",
           "cell-level AUC >= 0.65 (or <= 0.35), detected in >= 10% of the higher group. ",
           "No classifier was trained and no random cell-level train/test split was used, so no pseudo-independent ",
           "accuracy is reported. n = ", sum(sig$direction == "up_in_malignant_ECM"), " genes up in malignant ECM-like and ",
           sum(sig$direction == "up_in_true_fibroblast"), " up in true fibroblast."), 140),
         caption = P5_CAPTION) +
    p5_theme(9) + theme(axis.text.y = element_text(size = 6.5))
  save_fig5(p9, file.path(P5_FIG, "09_malignant_ecm_signature_heatmap.pdf"),
            7.0, max(9, 0.16 * nrow(sg) + 4))
} else cat("[warn] no signature genes met the declared criteria - figure 09 skipped\n")

# ---- 12 the Phase 5 model ----------------------------------------------------
m1 <- prog |> mutate(plab = factor(LABS[program_id], levels = LABS[PROGS])) |>
  ggplot(aes(plab, n_cells_dominant, fill = recurrence_status)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = n_cells_dominant), vjust = -0.4, size = 2.6) +
  scale_fill_manual(values = P5_REC, name = "recurrence", drop = FALSE) +
  labs(x = NULL, y = "malignant cells where\nthis is the top program",
       title = "A. programme size and recurrence") +
  p5_theme(9) + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7.5))
m2 <- rob |> mutate(plab = factor(LABS[program], levels = LABS[PROGS])) |>
  ggplot(aes(plab, loo_cosine_without_dominant_patient,
             fill = robust_overall)) +
  geom_col(width = 0.75) +
  geom_hline(yintercept = 0.60, linetype = 2, colour = "#B2182B", linewidth = 0.4) +
  scale_fill_manual(values = c(`TRUE` = "#A6DBA0", `FALSE` = "#F4A582"),
                    name = "robust overall") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "spectra cosine when the\ndominant patient is withheld",
       title = "B. the decisive robustness test") +
  p5_theme(9) + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7.5))
m3 <- ev |> mutate(lab = ifelse(evaluable, "evaluable", "NOT EVALUABLE")) |>
  ggplot(aes(sample_id, A_malignant_ECM, fill = lab)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("A=%d\nB=%d", A_malignant_ECM, B_true_fibroblast)),
            vjust = -0.2, size = 2.5) +
  scale_fill_manual(values = c(evaluable = "#A6DBA0",
                               `NOT EVALUABLE` = "grey80"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = NULL, y = "malignant ECM-like cells",
       title = "C. malignant ECM-like vs true fibroblast: where it can be tested") +
  p5_theme(9)
m4 <- amb |> mutate(plab = factor(LABS[program], levels = LABS[PROGS])) |>
  ggplot(aes(plab, median_usage, fill = sample_id)) +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8) +
  scale_fill_manual(values = P5_SAMPLE, name = "patient") +
  labs(x = NULL, y = "median program usage",
       title = "D. the 120 Ambiguous fibroblasts, projected - descriptive only",
       subtitle = wrap5("They are NOT reclassified: malignancy_refined stays Ambiguous for every one of them.", 110)) +
  p5_theme(9) + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7.5))

nrec <- sum(prog$recurrence_status == "recurrent")
p12 <- (m1 | m2) / (m3 | m4) +
  plot_annotation(
    title = "Figure 12. Phase 5 malignant-program model",
    subtitle = wrap5(paste0(
      "Phase 4 reported 0 of 8 DISCRETE malignant states recurrent across >= 3 patients, with 97.2% of malignant cells in ",
      "patient-private states. Phase 5 asked the continuous version of that question on the same 6,434 cells with the same ",
      "recurrence rule shape, and found ", nrec, " of ", nrow(prog), " programs recurrent. ",
      "A program is a continuous pattern of co-regulated expression: it is not a cell type, a lineage, a state or a clone, ",
      "and the two results are complementary rather than contradictory."), 165),
    caption = P5_CAPTION,
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 9.5, colour = "grey20"),
                  plot.caption = element_text(size = 7.5, colour = "grey35", hjust = 0)))
save_fig5(p12, file.path(P5_FIG, "12_phase5_malignant_program_model.pdf"), 15.0, 11.5)
cat("figures 08, 09, 12 written\n")
