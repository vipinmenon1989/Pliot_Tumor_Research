#!/usr/bin/env Rscript
# Phase 6 final figures 01-13.
# Every clone figure distinguishes the three clone-reliable patients from
# MPNST_3, whose SCEVAN clone structure is excluded; MPNST_3 clones are never
# pooled with the reliable ones without a visible warning.
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr)
                                 library(patchwork); library(scales); library(ggalluvial) })
source("scripts/phase6/utils/phase6_common.R")
source("scripts/phase5/utils/phase5_plot_utils.R")
set.seed(42)
rd <- function(p) read.delim(p, check.names = FALSE)

lab   <- rd("results/phase5/tables/final/MALIGNANT_PROGRAMS.tsv")
LAB   <- setNames(sprintf("%s %s", lab$program_id, lab$program_label), lab$program_id)
PID   <- lab$program_id
cpa   <- rd(file.path(P6_TAB, "CLONE_PROGRAM_ACTIVITY.tsv"))
assoc <- rd(file.path(P6_TAB, "CLONE_PROGRAM_ASSOCIATION.tsv"))
div   <- rd(file.path(P6_TAB, "CLONE_PROGRAM_DIVERSITY.tsv"))
plas  <- rd(file.path(P6_TAB, "PLASTICITY_METRICS.tsv"))
btw   <- rd(file.path(P6_PLAS, "M44_BETWEEN_VS_WITHIN_CLONE.tsv"))
conf  <- rd(file.path(P6_PLAS, "M44_DIVERSITY_CONFOUNDERS.tsv"))
tf    <- rd(file.path(P6_TAB, "PROGRAM_TF_ACTIVITY.tsv"))
pw    <- rd(file.path(P6_TAB, "PROGRAM_PATHWAY_ACTIVITY.tsv"))
evd   <- rd(file.path(P6_TAB, "PROGRAM_REGULATORY_EVIDENCE.tsv"))
cx    <- rd(file.path(P6_TAB, "BROAD_CNA_EXPRESSION_EFFECTS.tsv"))
cpr   <- if (file.exists(file.path(P6_TAB, "CNA_PROGRAM_ASSOCIATION.tsv")))
  rd(file.path(P6_TAB, "CNA_PROGRAM_ASSOCIATION.tsv")) else NULL
port  <- rd(file.path(P6_TAB, "PATIENT_TUMOR_ARCHITECTURE.tsv"))
integ <- rd(file.path(P6_TAB, "PHASE6_INTEGRATED_EVIDENCE.tsv"))
rob   <- rd(file.path(P6_TAB, "PHASE6_ROBUSTNESS.tsv"))
md    <- readRDS(file.path(P6_CLONE, "phase6_cell_metadata.rds"))
shortclone <- function(x) sub("^MPNST_[0-9]+_clone", "c", x)

# ---- 01 clone x program heatmap ---------------------------------------------
d1 <- cpa |> mutate(plab = factor(LAB[program], levels = LAB[PID]),
                    cl = shortclone(clone))
p1 <- ggplot(d1, aes(plab, cl, fill = median_usage)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", median_usage)), size = 2.3) +
  facet_grid(sample_id ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient(low = "grey97", high = "#762A83", name = "median\nrelative usage") +
  labs(x = NULL, y = "SCEVAN CNA-defined clone",
       title = "Figure 01. Malignant program activity in each CNA-defined clone",
       subtitle = wrap5(paste0(
         "Only the three clone-reliable patients are shown. MPNST_3 is absent by design: ", CLONE_EXCLUSION_REASON,
         " Clone labels are patient-scoped - MPNST_1_clone1 and MPNST_4_clone1 are unrelated names and are never ",
         "compared as if they were the same clone."), 150),
       caption = P6_CAPTION) +
  p5_theme(10) + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8))
save_fig5(p1, file.path(P6_FIG, "01_clone_program_heatmap.pdf"), 12.0, 11.0)

# ---- 02 clone -> program alluvial -------------------------------------------
mal <- md |> filter(malignancy_refined == "Malignant", !is.na(tumor_clone_phase4),
                    sample_id %in% CLONE_RELIABLE)
d2 <- mal |> count(sample_id, tumor_clone_phase4, dominant_malignant_program) |>
  mutate(clone = shortclone(tumor_clone_phase4),
         plab = factor(LAB[dominant_malignant_program], levels = LAB[PID]))
p2 <- ggplot(d2, aes(y = n, axis1 = clone, axis2 = plab)) +
  geom_alluvium(aes(fill = plab), width = 0.26, alpha = 0.75, knot.pos = 0.3,
                curve_type = "sigmoid") +
  geom_stratum(width = 0.26, fill = "grey96", colour = "grey40", linewidth = 0.25) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.3) +
  scale_fill_brewer(palette = "Set3", name = "dominant program") +
  scale_x_discrete(limits = c("CNA clone", "dominant program"), expand = c(0.2, 0.2)) +
  facet_wrap(~ sample_id, scales = "free_y", nrow = 1) +
  labs(y = "malignant cells", x = NULL,
       title = "Figure 02. Where each clone's cells land in program space",
       subtitle = wrap5(paste0(
         "A clone whose ribbons fan out to several programs contains cells spanning multiple malignant transcriptional ",
         "programs. That is DIVERSITY, not observed switching: cross-sectional scRNA-seq cannot show a transition, so no ",
         "transition rate, direction or trajectory is claimed."), 150),
       caption = P6_CAPTION) +
  p5_theme(10) + theme(panel.grid.major.x = element_blank())
save_fig5(p2, file.path(P6_FIG, "02_clone_program_alluvial.pdf"), 15.0, 8.5)

# ---- 03 diversity by clone ---------------------------------------------------
d3 <- plas |> mutate(cl = shortclone(clone),
                     ev = n_cells >= CLONE_SIZE_PRIMARY)
p3a <- ggplot(d3, aes(reorder(cl, -n_cells), hard_effective_n, fill = sample_id)) +
  geom_col(width = 0.75) +
  geom_hline(yintercept = 2, linetype = 2, linewidth = 0.4, colour = "#B2182B") +
  geom_text(aes(label = n_cells), vjust = -0.35, size = 2.4) +
  geom_text(data = d3 |> filter(!ev), aes(y = 0.15, label = "NOT EVALUABLE"),
            angle = 90, hjust = 0, size = 2.3, colour = "grey30") +
  facet_grid(~ sample_id, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = P5_SAMPLE, guide = "none") +
  labs(x = NULL, y = "effective number of programs\n(hard assignment)",
       title = "Figure 03. Within-clone transcriptional-program diversity",
       subtitle = wrap5(paste0(
         "Effective number of programs = exp(Shannon entropy) of the clone's dominant-program composition, using the ",
         "declared margin rule (a cell keeps its top program only when it leads the runner-up by >= 0.10 relative usage). ",
         "Red line = 2 programs. Clone size is printed above each bar; clones below ", CLONE_SIZE_PRIMARY,
         " cells are marked NOT EVALUABLE rather than hidden."), 150)) +
  p5_theme(10) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
p3b <- conf |> mutate(diversity_metric = sub("_", " ", diversity_metric)) |>
  ggplot(aes(covariate, diversity_metric, fill = spearman_rho)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(is.na(spearman_rho), "NA", sprintf("%.2f", spearman_rho))),
            size = 2.8) +
  scale_fill_gradient2(low = "#2166AC", mid = "grey96", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), na.value = "grey88",
                       name = "Spearman rho") +
  labs(x = NULL, y = NULL,
       title = "B. technical confounders, assessed BEFORE any plasticity wording",
       subtitle = wrap5(paste0("Diversity that tracks sequencing depth, clone size, proliferation or malignancy ",
         "confidence is technical, not biological. |rho| >= 0.7 would be disqualifying."), 130),
       caption = P6_CAPTION) +
  p5_theme(10) + theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 8))
save_fig5(p3a / p3b + plot_layout(heights = c(1.5, 1)),
          file.path(P6_FIG, "03_program_diversity_by_clone.pdf"), 13.0, 11.0)

# ---- 04 within vs between clone ---------------------------------------------
p4a <- assoc |> mutate(plab = factor(LAB[program], levels = LAB[PID])) |>
  ggplot(aes(plab, eta2, fill = sample_id)) +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.85) +
  geom_point(aes(y = eta2_null_p95, group = sample_id),
             position = position_dodge2(width = 0.85), shape = 95, size = 3,
             colour = "grey20") +
  scale_fill_manual(values = P5_SAMPLE, name = "patient") +
  labs(x = NULL, y = expression(eta^2 * " (between-clone share of variance)"),
       title = "Figure 04. How much program variance sits BETWEEN clones rather than within them",
       subtitle = wrap5(paste0(
         "Bars are eta^2 within one patient; the dash is the 95th percentile of a null built by shuffling clone labels ",
         "INSIDE that patient (1,000 permutations), which preserves clone sizes and the score distribution. ",
         "Even a clearly non-null eta^2 of 0.10 leaves 90% of the variance within clones."), 150)) +
  p5_theme(10) + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8))
p4b <- btw |> pivot_longer(c(mean_between_clone_jsd, mean_within_clone_dispersion),
                           names_to = "k", values_to = "v") |>
  ggplot(aes(sample_id, v, fill = k)) +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.7) +
  scale_fill_manual(values = c(mean_between_clone_jsd = "#762A83",
                               mean_within_clone_dispersion = "#A6DBA0"),
                    labels = c("between-clone divergence (JSD)",
                               "within-clone dispersion"), name = NULL) +
  labs(x = NULL, y = "value", title = "B. between-clone divergence vs within-clone spread",
       caption = P6_CAPTION) + p5_theme(10)
save_fig5(p4a / p4b + plot_layout(heights = c(1.6, 1)),
          file.path(P6_FIG, "04_within_vs_between_clone_variation.pdf"), 12.0, 10.0)

# ---- 05 TF activity ----------------------------------------------------------
t5 <- tf |> group_by(program) |> slice_max(abs(pooled_rho), n = 12) |> ungroup() |>
  mutate(plab = factor(LAB[program], levels = LAB[PID]))
p5f <- ggplot(t5, aes(pooled_rho, reorder(tf, pooled_rho), fill = recurrent_direction)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
  facet_wrap(~ plab, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(`TRUE` = "#1A9850", `FALSE` = "grey72"),
                    labels = c(`TRUE` = "concordant in >= 3 patients",
                               `FALSE` = "pooled only"), name = NULL) +
  labs(x = "pooled Spearman rho (TF activity vs program usage)", y = NULL,
       title = "Figure 05. Transcription-factor activity distinguishing the malignant programs",
       subtitle = wrap5(paste0(
         "decoupleR run_ulm over CollecTRI on the frozen RNA log-normalised layer. Green bars are the only ones treated as ",
         "recurrent: direction concordant in >= 3 evaluable patients AND |pooled rho| >= 0.20. A strong pooled association ",
         "alone is never sufficient - with n = 4 and sample_id = patient = dataset it can be one patient's signal."), 150),
       caption = P6_CAPTION) +
  p5_theme(9) + theme(axis.text.y = element_text(size = 6.5), legend.position = "bottom")
save_fig5(p5f, file.path(P6_FIG, "05_tf_activity_by_program.pdf"),
          16.0, max(8, 2.6 * ceiling(length(PID) / 4) + 3))

# ---- 06 pathway activity -----------------------------------------------------
p6a <- pw |> filter(layer == "PROGENy") |>
  mutate(plab = factor(LAB[program], levels = LAB[PID])) |>
  ggplot(aes(plab, pathway, fill = pooled_rho)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(recurrent_direction, sprintf("%.2f*", pooled_rho),
                               sprintf("%.2f", pooled_rho))), size = 2.4) +
  scale_fill_gradient2(low = "#2166AC", mid = "grey96", high = "#B2182B",
                       midpoint = 0, name = "pooled rho") +
  labs(x = NULL, y = NULL, title = "A. PROGENy pathway activity vs program usage",
       subtitle = wrap5("* = direction concordant in >= 3 evaluable patients and |pooled rho| >= 0.20.", 130)) +
  p5_theme(10) + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8))
h6 <- pw |> filter(layer == "Hallmark") |> group_by(program) |>
  slice_max(abs(pooled_rho), n = 6) |> ungroup() |>
  mutate(plab = factor(LAB[program], levels = LAB[PID]),
         pathway = sub("^HALLMARK_", "", pathway))
p6b <- ggplot(h6, aes(pooled_rho, reorder(pathway, pooled_rho),
                      fill = recurrent_direction)) +
  geom_col(width = 0.7) + geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
  facet_wrap(~ plab, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(`TRUE` = "#1A9850", `FALSE` = "grey72"), name = NULL,
                    labels = c(`TRUE` = "concordant in >= 3 patients", `FALSE` = "pooled only")) +
  labs(x = "pooled Spearman rho", y = NULL, title = "B. Hallmark gene-set score vs program usage",
       subtitle = wrap5(paste0("PROGENy and Hallmark are kept as two independent layers and are never merged into a ",
         "composite score or a single numeric rank."), 150),
       caption = P6_CAPTION) +
  p5_theme(9) + theme(axis.text.y = element_text(size = 6.2), legend.position = "bottom")
save_fig5(p6a / p6b + plot_layout(heights = c(1, 1.5)),
          file.path(P6_FIG, "06_pathway_activity_by_program.pdf"),
          16.0, max(12, 2.4 * ceiling(length(PID) / 4) + 7))

# ---- 07 program - TF - pathway evidence map ---------------------------------
e7 <- evd |> mutate(plab = factor(LAB[program], levels = LAB[PID]))
p7 <- if (nrow(e7)) {
  ggplot(e7, aes(pooled_rho, reorder(paste(feature, layer), pooled_rho),
                 fill = layer)) +
    geom_col(width = 0.7) + geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    facet_wrap(~ plab, scales = "free_y", ncol = 3) +
    scale_fill_brewer(palette = "Dark2", name = "evidence layer") +
    labs(x = "pooled Spearman rho", y = NULL,
         title = "Figure 07. The three evidence layers per program, side by side",
         subtitle = wrap5(paste0(
           "Only cross-patient-concordant associations are shown (>= 3 evaluable patients agreeing in direction and ",
           "|pooled rho| >= 0.20). Transcriptional program, TF activity and pathway activity are INDEPENDENT layers; ",
           "no weighted composite is computed and no single numeric rank is imposed."), 150),
         caption = P6_CAPTION) +
    p5_theme(9) + theme(axis.text.y = element_text(size = 6.2))
} else {
  ggplot() + annotate("text", 0, 0, size = 5, label = paste(
    "No TF or pathway association was concordant across >= 3 patients\n",
    "at |pooled rho| >= 0.20. Reported as a negative result.")) +
    theme_void() + labs(title = "Figure 07. Program - TF - pathway evidence",
                        caption = P6_CAPTION)
}
save_fig5(p7, file.path(P6_FIG, "07_program_tf_pathway_network.pdf"),
          15.0, max(8, 2.6 * ceiling(length(PID) / 3) + 3))

# ---- 08 broad CNA expression effects ----------------------------------------
c8 <- cx |> filter(evaluable) |>
  mutate(seg = paste0(sample_id, " ", segment))
p8 <- if (nrow(c8)) {
  ggplot(c8, aes(cohens_d, reorder(seg, cohens_d), fill = event)) +
    geom_col(width = 0.7) + geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    scale_fill_manual(values = c(gain = "#B2182B", loss = "#2166AC"),
                      name = "inferred broad event") +
    labs(x = "Cohen's d of segment-level mean expression (carrier clones - non-carrier clones)",
         y = NULL,
         title = "Figure 08. Transcriptional consequences of broad inferred CNA events",
         subtitle = wrap5(paste0(
           "Within each patient, cells in clones carrying an inferred broad segment are compared with cells in clones of the ",
           "SAME patient that do not carry it, using the mean expression of all genes in that interval. Segment level only - ",
           "no single-gene deletion or amplification is asserted anywhere. ",
           "NOT VALIDATION: SCEVAN inferred these events from expression in the first place, so agreement is expected and is ",
           "an internal consistency check, not independent evidence for the copy-number call."), 150),
         caption = P6_CAPTION) +
    p5_theme(9) + theme(axis.text.y = element_text(size = 6.5))
} else {
  ggplot() + annotate("text", 0, 0, size = 4.5, label = paste(
    "No broad inferred CNA event was evaluable:\n",
    "an event needs carrier AND non-carrier clones in the same patient,\n",
    "each above the declared cell minimum, and >= 30 genes in the interval.")) +
    theme_void() + labs(title = "Figure 08. Broad CNA expression effects", caption = P6_CAPTION)
}
save_fig5(p8, file.path(P6_FIG, "08_broad_cna_expression_effects.pdf"),
          13.0, max(7, 0.22 * max(1, nrow(c8)) + 5))

# ---- 09 CNA - program association -------------------------------------------
p9 <- if (!is.null(cpr) && nrow(cpr)) {
  d9 <- cpr |> mutate(plab = factor(LAB[program], levels = LAB[PID]),
                      seg = paste0(sample_id, " ", segment))
  ggplot(d9, aes(plab, seg, fill = cohens_d)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    scale_fill_gradient2(low = "#2166AC", mid = "grey96", high = "#B2182B",
                         midpoint = 0, name = "Cohen's d") +
    labs(x = NULL, y = NULL,
         title = "Figure 09. Broad CNA events versus malignant program activity, within patient",
         subtitle = wrap5(paste0(
           "Each row is one inferred broad event in one patient; colour is the standardised difference in program usage ",
           "between the clones carrying it and the clones of the same patient that do not. CNA architecture is patient-specific ",
           "and is never pooled across patients."), 150),
         caption = P6_CAPTION) +
    p5_theme(9) + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8),
                        axis.text.y = element_text(size = 6))
} else {
  ggplot() + annotate("text", 0, 0, size = 4.5,
    label = "No evaluable CNA - program comparison.") + theme_void() +
    labs(title = "Figure 09. CNA - program associations", caption = P6_CAPTION)
}
save_fig5(p9, file.path(P6_FIG, "09_cna_program_associations.pdf"), 13.0, 10.0)
cat("figures 01-09 written\n")
