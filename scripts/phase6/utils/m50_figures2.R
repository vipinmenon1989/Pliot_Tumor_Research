#!/usr/bin/env Rscript
# Phase 6 final figures 10-13.
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr)
                                 library(patchwork); library(scales) })
source("scripts/phase6/utils/phase6_common.R")
source("scripts/phase5/utils/phase5_plot_utils.R")
set.seed(42)
rd <- function(p) read.delim(p, check.names = FALSE)
lab   <- rd("results/phase5/tables/final/MALIGNANT_PROGRAMS.tsv")
LAB   <- setNames(sprintf("%s %s", lab$program_id, lab$program_label), lab$program_id)
PID   <- lab$program_id
port  <- rd(file.path(P6_TAB, "PATIENT_TUMOR_ARCHITECTURE.tsv"))
integ <- rd(file.path(P6_TAB, "PHASE6_INTEGRATED_EVIDENCE.tsv"))
rob   <- rd(file.path(P6_TAB, "PHASE6_ROBUSTNESS.tsv"))
plas  <- rd(file.path(P6_TAB, "PLASTICITY_METRICS.tsv"))
assoc <- rd(file.path(P6_TAB, "CLONE_PROGRAM_ASSOCIATION.tsv"))
tf    <- rd(file.path(P6_TAB, "PROGRAM_TF_ACTIVITY.tsv"))
pw    <- rd(file.path(P6_TAB, "PROGRAM_PATHWAY_ACTIVITY.tsv"))
md    <- readRDS(file.path(P6_CLONE, "phase6_cell_metadata.rds"))
arch  <- jsonlite::fromJSON(file.path(P6_VAL, "m48_integrated_facts.json"),
                            simplifyVector = FALSE)$architecture

# ---- 10 per-patient tumour architectures ------------------------------------
mal <- md |> filter(malignancy_refined == "Malignant")
comp <- mal |> count(sample_id, dominant_malignant_program) |>
  group_by(sample_id) |> mutate(frac = n / sum(n)) |> ungroup() |>
  mutate(plab = factor(LAB[dominant_malignant_program], levels = LAB[PID]))
p10a <- ggplot(comp, aes(sample_id, frac, fill = plab)) +
  geom_col(width = 0.75, colour = "white", linewidth = 0.15) +
  scale_fill_brewer(palette = "Set3", name = "dominant program") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(x = NULL, y = "fraction of the patient's malignant cells",
       title = "A. malignant program composition, all four patients") +
  p5_theme(10)
p10b <- port |> mutate(lab = ifelse(clone_structure_reliable,
                                    "clone structure usable", "clone structure EXCLUDED")) |>
  ggplot(aes(sample_id, n_malignant, fill = lab)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%d cells\n%s clones", n_malignant,
            ifelse(is.na(n_clones_evaluable), paste0(n_clones_total, " (unreliable)"),
                   paste0(n_clones_evaluable, " evaluable")))),
            vjust = -0.15, size = 2.6) +
  scale_fill_manual(values = c(`clone structure usable` = "#A6DBA0",
                               `clone structure EXCLUDED` = "#F4A582"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.28))) +
  labs(x = NULL, y = "refined malignant cells",
       title = "B. which patients support clone-based conclusions",
       subtitle = wrap5(CLONE_EXCLUSION_REASON, 120)) + p5_theme(10)
p10c <- plas |> filter(n_cells >= CLONE_SIZE_PRIMARY) |>
  ggplot(aes(sample_id, hard_effective_n, fill = sample_id)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.65) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.85) +
  geom_hline(yintercept = 2, linetype = 2, colour = "#B2182B", linewidth = 0.4) +
  scale_fill_manual(values = P5_SAMPLE, guide = "none") +
  labs(x = NULL, y = "effective programs per clone",
       title = "C. within-clone program diversity per patient") + p5_theme(10)
p10d <- assoc |> mutate(plab = factor(LAB[program], levels = LAB[PID])) |>
  ggplot(aes(sample_id, eta2, fill = sample_id)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.65) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.85) +
  scale_fill_manual(values = P5_SAMPLE, guide = "none") +
  labs(x = NULL, y = expression(eta^2 * " per program"),
       title = "D. clone-program coupling per patient",
       caption = P6_CAPTION) + p5_theme(10)
save_fig5((p10a | p10b) / (p10c | p10d) +
  plot_annotation(title = "Figure 10. Per-patient tumour architectures",
    subtitle = wrap5(paste0("With n = 4 patients - three of them usable for clone inference - deep per-patient portraits ",
      "are the honest deliverable, not a population-level statistical claim."), 160),
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 9.5, colour = "grey20"))),
  file.path(P6_FIG, "10_patient_tumor_architectures.pdf"), 14.5, 10.5)

# ---- 11 malignant-ECM regulatory profile ------------------------------------
ecm <- lab$program_id[grepl("ECM|Mesenchymal", lab$program_label, ignore.case = TRUE)]
if (!length(ecm)) ecm <- integ$program[which.max(integ$n_clones_dominated)]
p11a <- tf |> filter(program %in% ecm) |> group_by(program) |>
  slice_max(abs(pooled_rho), n = 15) |> ungroup() |>
  mutate(plab = factor(LAB[program], levels = LAB[PID])) |>
  ggplot(aes(pooled_rho, reorder(tf, pooled_rho), fill = recurrent_direction)) +
  geom_col(width = 0.7) + geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
  facet_wrap(~ plab, scales = "free_y") +
  scale_fill_manual(values = c(`TRUE` = "#1A9850", `FALSE` = "grey72"), name = NULL,
                    labels = c(`TRUE` = "concordant >= 3 patients", `FALSE` = "pooled only")) +
  labs(x = "pooled Spearman rho", y = NULL, title = "A. TF activity (CollecTRI)")+
  p5_theme(9) + theme(axis.text.y = element_text(size = 7))
p11b <- pw |> filter(program %in% ecm, layer == "PROGENy") |>
  mutate(plab = factor(LAB[program], levels = LAB[PID])) |>
  ggplot(aes(pooled_rho, reorder(pathway, pooled_rho), fill = recurrent_direction)) +
  geom_col(width = 0.7) + geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
  facet_wrap(~ plab, scales = "free_y") +
  scale_fill_manual(values = c(`TRUE` = "#1A9850", `FALSE` = "grey72"), guide = "none") +
  labs(x = "pooled Spearman rho", y = NULL, title = "B. PROGENy pathway activity") +
  p5_theme(9)
sigf <- "results/phase5/tables/final/MALIGNANT_ECM_SIGNATURE.tsv"
p11c <- if (file.exists(sigf) && nrow(rd(sigf))) {
  s <- rd(sigf) |> group_by(direction) |> slice_max(abs(mean_lfc), n = 20) |> ungroup()
  ggplot(s, aes(mean_lfc, reorder(gene, mean_lfc), fill = direction)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c(up_in_malignant_ECM = "#D95F02",
                                 up_in_true_fibroblast = "#2166AC"), name = NULL) +
    labs(x = "mean pseudobulk log2FC across evaluable patients", y = NULL,
         title = "C. the Phase 5 malignant-ECM signature it regulates",
         caption = P6_CAPTION) +
    p5_theme(9) + theme(axis.text.y = element_text(size = 6.5))
} else {
  ggplot() + annotate("text", 0, 0, size = 4,
    label = "No malignant-ECM signature gene met the declared criteria.") +
    theme_void() + labs(caption = P6_CAPTION)
}
save_fig5((p11a / p11b) | p11c +
  plot_annotation(title = "Figure 11. Regulatory profile of the malignant ECM-like compartment",
    subtitle = wrap5(paste0("Programs shown: ", paste(LAB[ecm], collapse = ", "), ". ",
      "Only green bars are treated as recurrent (concordant direction in >= 3 evaluable patients and |pooled rho| >= 0.20)."), 160),
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 9.5, colour = "grey20"))),
  file.path(P6_FIG, "11_malignant_ecm_regulatory_profile.pdf"), 15.0, 10.0)

# ---- 12 robustness -----------------------------------------------------------
p12 <- ggplot(rob, aes(reorder(paste(test, setting, sep = " | "), value), value,
                       fill = test)) +
  geom_col(width = 0.72) + coord_flip() +
  geom_text(aes(label = sprintf("%.2f  (n=%s)", value, n)), hjust = -0.05, size = 2.6) +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.30))) +
  labs(x = NULL, y = "value (metric named in the table)",
       title = "Figure 12. Phase 6 robustness",
       subtitle = wrap5(paste0(
         "Clone-size, dominant-program-margin and high-confidence sensitivity, plus the patient-level direction check. ",
         "Leave-one-patient-out is deliberately NOT applied to the clone analyses: they are inherently within-patient, so ",
         "withholding a patient deletes the analysis instead of testing it. Clones failing a size threshold are reported ",
         "NOT EVALUABLE, never hidden."), 150),
       caption = P6_CAPTION) +
  p5_theme(9) + theme(axis.text.y = element_text(size = 7), legend.position = "bottom")
save_fig5(p12, file.path(P6_FIG, "12_phase6_robustness_summary.pdf"), 13.0, 9.5)

# ---- 13 final architecture ----------------------------------------------------
ev <- plas |> filter(n_cells >= CLONE_SIZE_PRIMARY)
f1 <- ggplot(ev, aes(dominant_program_share_hard, hard_effective_n,
                     colour = sample_id, size = n_cells)) +
  geom_point(alpha = 0.85) +
  geom_hline(yintercept = 2, linetype = 2, colour = "#B2182B", linewidth = 0.4) +
  geom_vline(xintercept = 0.8, linetype = 2, colour = "#2166AC", linewidth = 0.4) +
  scale_colour_manual(values = P5_SAMPLE, name = "patient") +
  scale_size_continuous(range = c(1.5, 6), name = "clone size") +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  annotate("text", x = 0.83, y = max(ev$hard_effective_n, na.rm = TRUE),
           hjust = 0, vjust = 1, size = 2.7, colour = "#2166AC",
           label = "Model A region:\nclone-constrained") +
  annotate("text", x = min(ev$dominant_program_share_hard, na.rm = TRUE), y = 2.1,
           hjust = 0, vjust = 0, size = 2.7, colour = "#B2182B",
           label = "Model B region: clone spans\nmultiple programs") +
  labs(x = "share of the clone's cells in its top program",
       y = "effective number of programs per clone",
       title = "A. every evaluable clone placed against the two competing architectures") +
  p5_theme(10)
# Each row is scaled to its own maximum: the TF counts run into the hundreds
# while the others are single digits, so a shared scale would flatten every row
# except one. The printed number is always the raw value.
f2 <- integ |> mutate(plab = factor(LAB[program], levels = LAB[PID])) |>
  select(plab, median_eta2, n_patients_clone_associated, n_clones_dominated,
         n_recurrent_tfs, n_recurrent_progeny) |>
  pivot_longer(-plab) |>
  mutate(name = factor(name, levels = c("median_eta2", "n_patients_clone_associated",
    "n_clones_dominated", "n_recurrent_tfs", "n_recurrent_progeny"),
    labels = c("median eta^2", "patients where clone-associated", "clones dominated",
               "cross-patient-concordant TFs", "cross-patient-concordant PROGENy"))) |>
  group_by(name) |> mutate(rel = if (max(value) > 0) value / max(value) else 0) |>
  ungroup()
f2 <- ggplot(f2, aes(plab, name, fill = rel)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(value == round(value), sprintf("%.0f", value),
                               sprintf("%.2f", value)),
                colour = rel > 0.6), size = 2.7) +
  scale_fill_gradient(low = "grey96", high = "#2C7FB8", guide = "none") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"), guide = "none") +
  labs(x = NULL, y = NULL, title = "B. the three evidence layers per program, unmerged",
       subtitle = wrap5(paste0("Colour is scaled within each row because the counts live on different scales; ",
         "the printed number is the raw value. The layers are deliberately not combined into one score."), 105)) +
  p5_theme(10) + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8))

f3 <- port |> mutate(lab = ifelse(clone_structure_reliable, "usable", "EXCLUDED")) |>
  ggplot(aes(sample_id, ifelse(is.na(median_eta2), 0, median_eta2), fill = lab)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = ifelse(is.na(median_eta2), "clone structure\nEXCLUDED",
                               sprintf("eta^2 %.3f\n%.1f programs/clone", median_eta2,
                                       median_clone_effective_programs))),
            vjust = -0.15, size = 2.5) +
  scale_fill_manual(values = c(usable = "#A6DBA0", EXCLUDED = "#F4A582"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
  labs(x = NULL, y = expression("median " * eta^2), title = "C. per patient") +
  p5_theme(10)

p13 <- (f1 | f2) / f3 + plot_layout(heights = c(1.6, 1)) +
  plot_annotation(
    title = "Figure 13. Integrated MPNST malignant architecture",
    subtitle = wrap5(paste0("SELECTED MODEL: ", arch$selected, ". ", arch$verdict,
      "  Evidence: ", arch$evidence$n_multi_program, " of ", arch$evidence$n_evaluable_clones,
      " evaluable clones span >= 2 programs, ", arch$evidence$n_concentrated,
      " are concentrated in one, median eta^2 = ", sprintf("%.3f", arch$evidence$median_eta2),
      " (so roughly ", sprintf("%.0f%%", 100 * (1 - arch$evidence$median_eta2)),
      " of program variance sits WITHIN clones), and ", arch$evidence$n_clone_associated_pairs,
      " of ", arch$evidence$n_pairs, " program x patient pairs are clone-associated."), 165),
    caption = P6_CAPTION,
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 9.5, colour = "grey20"),
                  plot.caption = element_text(size = 7.5, colour = "grey35", hjust = 0)))
save_fig5(p13, file.path(P6_FIG, "13_final_mpnst_malignant_architecture.pdf"), 16.0, 12.0)
cat("figures 10-13 written\n")
