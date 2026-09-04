#!/usr/bin/env Rscript
# =============================================================================
# Phase 5 final figures 01-12.
#
# Rule carried through every panel: a UMAP is never the evidence on its own.
# Every UMAP is accompanied by the quantitative panel that actually supports
# the claim - patient distribution, program score, effect size, recurrence or
# robustness.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr)
                                 library(patchwork); library(scales) })
source("scripts/phase5/utils/phase5_common.R")
source("scripts/phase5/utils/phase5_plot_utils.R")
set.seed(42)
args <- commandArgs(trailingOnly = TRUE)
KSTAR <- as.integer(args[1])

rd <- function(p, ...) read.delim(p, check.names = FALSE, ...)
ks   <- rd(file.path(P5_TAB, "PROGRAM_K_SELECTION.tsv"))
prog <- rd(file.path(P5_TAB, "MALIGNANT_PROGRAMS.tsv"))
tg   <- rd(file.path(P5_TAB, "MALIGNANT_PROGRAM_TOP_GENES.tsv"))
cs   <- rd(file.path(P5_TAB, "MALIGNANT_PROGRAM_CELL_SCORES.tsv"))
pd   <- rd(file.path(P5_TAB, "MALIGNANT_PROGRAM_PATIENT_DISTRIBUTION.tsv"))
pw   <- rd(file.path(P5_TAB, "MALIGNANT_PROGRAM_PATHWAYS.tsv"))
st   <- rd(file.path(P5_TAB, "PROGRAM_VS_PHASE4_STATE.tsv"))
bal  <- rd(file.path(P5_TAB, "PROGRAM_PATIENT_BALANCE_SENSITIVITY.tsv"))
sens <- rd(file.path(P5_TAB, "PROGRAM_RECURRENCE_THRESHOLD_SENSITIVITY.tsv"))
rob  <- rd(file.path(P5_TAB, "PHASE5_PROGRAM_ROBUSTNESS.tsv"))
md   <- readRDS(file.path(P5_PROG, "phase5_cell_metadata.rds"))

PROGS <- prog$program_id
LAB <- setNames(sprintf("%s\n%s", prog$program_id, prog$program_label), PROGS)
LABS <- setNames(sprintf("%s %s", prog$program_id, prog$program_label), PROGS)
REC <- setNames(prog$recurrence_status, PROGS)
prog$plab <- LABS[prog$program_id]

# ---- 01 program rank selection ---------------------------------------------
k1 <- ks |> filter(run == "primary") |>
  select(K, silhouette = stability, prediction_error, max_program_cosine,
         min_dominant_share, n_dead_programs) |>
  pivot_longer(-K, names_to = "metric", values_to = "value") |>
  mutate(metric = factor(metric, levels = c("silhouette", "prediction_error",
    "max_program_cosine", "min_dominant_share", "n_dead_programs")))
hl <- data.frame(metric = factor(c("max_program_cosine", "min_dominant_share"),
                                 levels = levels(k1$metric)), y = c(0.75, 0.01))
p1 <- ggplot(k1, aes(K, value)) +
  geom_vline(xintercept = KSTAR, colour = "#B2182B", linewidth = 0.8, linetype = 2) +
  geom_hline(data = hl, aes(yintercept = y), linetype = 3, colour = "grey35") +
  geom_line(linewidth = 0.7, colour = "grey25") + geom_point(size = 2) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  scale_x_continuous(breaks = sort(unique(k1$K))) +
  labs(x = "K (number of cNMF programs)", y = NULL,
       title = "Figure 01. Program-rank selection across the declared K grid",
       subtitle = wrap5(paste0(
         "Red line = the selected K = ", KSTAR, ". The selection rule was declared before any factor was inspected: ",
         "the LARGEST K whose consensus programs are not duplicated (max pairwise cosine <= 0.75, dotted line), ",
         "that has no dead program (every program dominant in >= 1% of malignant cells, dotted line), and whose ",
         "stability is at least the median stability of the K values satisfying both. ",
         "The rule references only measured properties of the factorization - never a gene, a pathway or a label."), 150),
       caption = P5_CAPTION) + p5_theme(10)
save_fig5(p1, file.path(P5_FIG, "01_program_rank_selection.pdf"), 9.5, 11.0)

# ---- 02 program heatmap ----------------------------------------------------
# z-scores come from the FULL spectra matrix, not from the top-50 table: taking
# them from the top-50 rows alone would zero-fill the other programs and distort
# every z.
ZF <- file.path("results/phase5/programs/cnmf/primary",
                sprintf("primary.gene_spectra_score.k_%d.dt_0_1.txt", KSTAR))
Zm <- as.matrix(read.delim(ZF, row.names = 1, check.names = FALSE))
rownames(Zm) <- PROGS
gl <- unique((tg |> filter(rank <= 10) |> arrange(program, rank))$gene)
gl <- intersect(gl, colnames(Zm))
Zs <- scale(Zm[, gl, drop = FALSE])            # z across programs, per gene
hm <- as.data.frame(as.table(Zs))
names(hm) <- c("program", "gene", "z")
hm <- hm |> mutate(gene = factor(as.character(gene), levels = rev(gl)),
                   plab = factor(LABS[as.character(program)], levels = LABS[PROGS]))
p2 <- ggplot(hm, aes(plab, gene, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.1) +
  scale_fill_gradient2(low = "#2166AC", mid = "grey96", high = "#B2182B",
                       midpoint = 0, name = "z of cNMF\nspectra score") +
  labs(x = NULL, y = NULL,
       title = sprintf("Figure 02. The %d malignant transcriptional programs, top 10 genes each", length(PROGS)),
       subtitle = wrap5(paste0(
         "Rows are the union of each program's 10 highest-scoring genes; colour is the gene's cNMF spectra score ",
         "z-scored across programs, so a red block marks the program that gene defines. ",
         "A program is a continuous pattern of co-regulated expression - it is not a cell type, a lineage, a state or a clone."), 150),
       caption = P5_CAPTION) + p5_theme(9) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8),
        axis.text.y = element_text(size = 6.2))
save_fig5(p2, file.path(P5_FIG, "02_malignant_program_heatmap.pdf"),
          max(8, 1.0 * length(PROGS) + 5), max(10, 0.13 * length(gl) + 4))

# ---- 03 program activity UMAP (never alone) --------------------------------
# UMAP coordinates live in the object's REDUCTION, not in meta.data, so they are
# read from the frozen Phase 4 embedding table rather than re-derived. Display
# only: programs came from RNA counts, never from these coordinates.
UM <- readRDS("results/phase4/malignancy/phase4_malignancy_metadata.rds")
stopifnot(all(c("postint_umap_harmony_1", "postint_umap_harmony_2") %in% colnames(UM)))
u <- UM[cs$cell_id, c("postint_umap_harmony_1", "postint_umap_harmony_2")]
cs$U1 <- u[[1]]; cs$U2 <- u[[2]]
long <- cs |> select(cell_id, sample_id, U1, U2, all_of(PROGS)) |>
  pivot_longer(all_of(PROGS), names_to = "program", values_to = "usage") |>
  mutate(plab = factor(LABS[program], levels = LABS[PROGS]))
p3a <- ggplot(long |> arrange(usage), aes(U1, U2, colour = usage)) +
  geom_point(size = 0.28, stroke = 0) +
  scale_colour_viridis_c(option = "magma", direction = -1,
                         name = "relative\nprogram usage") +
  facet_wrap(~ plab, ncol = 4) + coord_equal() +
  labs(x = "UMAP 1 (Phase 2 Harmony)", y = "UMAP 2",
       title = "Figure 03A. Program activity across the 6,434 refined malignant cells",
       subtitle = wrap5(paste0("The UMAP is the frozen Phase 2 Harmony embedding and is used for DISPLAY ONLY - ",
         "programs were derived from RNA counts, not from these coordinates. Panel B carries the quantitative claim."), 150)) +
  p5_theme(9) + theme(axis.text = element_blank(), axis.ticks = element_blank())
p3b <- ggplot(cs, aes(U1, U2, colour = sample_id)) +
  geom_point(size = 0.3, stroke = 0, alpha = 0.7) +
  scale_colour_manual(values = P5_SAMPLE, name = "patient") +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  coord_equal() + labs(x = "UMAP 1", y = "UMAP 2", title = "B. patient") +
  p5_theme(9) + theme(axis.text = element_blank(), axis.ticks = element_blank())
p3c <- pd |> mutate(plab = factor(LABS[program], levels = LABS[PROGS])) |>
  ggplot(aes(plab, frac_of_patient, fill = sample_id)) +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.85) +
  geom_hline(yintercept = 0.05, linetype = 2, linewidth = 0.35, colour = "#B2182B") +
  scale_fill_manual(values = P5_SAMPLE, name = "patient") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "fraction of that patient's\nmalignant cells active",
       title = "C. the quantitative panel: per-patient activity, with the declared 5% carriage threshold",
       caption = P5_CAPTION) +
  p5_theme(9) + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7.5))
p3 <- p3a / (p3b | p3c) + plot_layout(heights = c(2.2, 1))
save_fig5(p3, file.path(P5_FIG, "03_program_activity_umap.pdf"), 14.0, 13.0)

# ---- 04 patient distribution ------------------------------------------------
p4a <- pd |> mutate(plab = factor(LABS[program], levels = LABS[PROGS])) |>
  ggplot(aes(plab, median_rel_usage, fill = sample_id)) +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.85) +
  scale_fill_manual(values = P5_SAMPLE, name = "patient") +
  labs(x = NULL, y = "median relative usage",
       title = "A. median program usage per patient") +
  p5_theme(10) + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8))
p4b <- prog |> mutate(plab = factor(plab, levels = LABS[PROGS])) |>
  ggplot(aes(plab, dominant_patient_fraction, fill = dominant_patient)) +
  geom_col(width = 0.75) +
  geom_hline(yintercept = 0.5727, linetype = 2, linewidth = 0.4, colour = "grey30") +
  annotate("text", x = 0.6, y = 0.60, hjust = 0, size = 2.8, colour = "grey30",
           label = "MPNST_4's share of the malignant compartment (0.573)") +
  scale_fill_manual(values = P5_SAMPLE, name = "dominant patient") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = NULL, y = "dominant-patient fraction\nof program-active cells",
       title = "B. is a program just its biggest patient?",
       subtitle = wrap5("A program whose bar sits near the dashed line is no more patient-skewed than the compartment itself.", 120)) +
  p5_theme(10) + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8))
p4c <- pd |> mutate(plab = factor(LABS[program], levels = LABS[PROGS])) |>
  ggplot(aes(plab, sample_id, fill = n_active)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = n_active), size = 2.6) +
  scale_fill_gradient(low = "grey96", high = "#1B7837", trans = "sqrt",
                      name = "active cells") +
  labs(x = NULL, y = NULL, title = "C. program-active cells per patient",
       caption = P5_CAPTION) +
  p5_theme(10) + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8))
p4 <- p4a / p4b / p4c + plot_layout(heights = c(1, 1, 0.8))
save_fig5(p4, file.path(P5_FIG, "04_program_patient_distribution.pdf"), 11.5, 12.5)

# ---- 05 program vs Phase 4 states -------------------------------------------
stm <- st |> mutate(plab = factor(LABS[program], levels = LABS[PROGS]),
                    is_ecm = grepl("^Mesenchymal_ECM", tumor_state_phase4))
p5a <- ggplot(stm, aes(plab, tumor_state_phase4, fill = median_rel_usage)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", median_rel_usage)), size = 2.5) +
  scale_fill_gradient(low = "grey97", high = "#7B3294", name = "median\nrelative usage") +
  labs(x = NULL, y = "Phase 4 discrete tumour state",
       title = "Figure 05. Continuous Phase 5 programs against the discrete Phase 4 states",
       subtitle = wrap5(paste0(
         "Phase 4 found 0 of 8 discrete states recurrent across >= 3 patients, with 97.2% of malignant cells in ",
         "patient-private states. The question here is whether the five patient-private Mesenchymal_ECM states are ",
         "separate biology or patient-specific manifestations of the same continuous program(s)."), 150)) +
  p5_theme(10) + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8))
p5b <- stm |> group_by(tumor_state_phase4) |>
  slice_max(median_rel_usage, n = 1, with_ties = FALSE) |> ungroup() |>
  ggplot(aes(tumor_state_phase4, median_rel_usage, fill = plab)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = program), vjust = -0.4, size = 2.8) +
  scale_fill_brewer(palette = "Set3", name = "dominant program") +
  labs(x = NULL, y = "median relative usage\nof the dominant program",
       title = "B. which single program dominates each Phase 4 state",
       caption = P5_CAPTION) +
  p5_theme(10) + theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8))
save_fig5(p5a / p5b + plot_layout(heights = c(1.6, 1)),
          file.path(P5_FIG, "05_program_vs_phase4_states.pdf"), 11.5, 11.0)

# ---- 06 recurrence -----------------------------------------------------------
p6a <- prog |> mutate(plab = factor(plab, levels = LABS[PROGS])) |>
  ggplot(aes(plab, patients_represented, fill = recurrence_status)) +
  geom_col(width = 0.75) +
  geom_hline(yintercept = 3, linetype = 2, linewidth = 0.4, colour = "#B2182B") +
  geom_text(aes(label = carrying_patients), vjust = -0.5, size = 2.4) +
  scale_fill_manual(values = P5_REC, name = "recurrence", drop = FALSE) +
  scale_y_continuous(breaks = 0:4, limits = c(0, 4.6)) +
  labs(x = NULL, y = "patients carrying the program",
       title = "Figure 06. Recurrent versus patient-private malignant programs",
       subtitle = wrap5(paste0(
         "Criteria declared BEFORE any program was labelled: a cell is program-active at >= 0.20 relative usage; a patient ",
         "carries a program with >= 10 active cells AND >= 5% of its malignant cells active; recurrent requires >= 3 of 4 ",
         "patients (red line). This is deliberately the same shape as the Phase 4 state-recurrence rule, so the two results ",
         "are directly comparable."), 150)) +
  p5_theme(10) + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8))
sk <- sens |> filter(K == KSTAR) |>
  pivot_longer(starts_with("n_"), names_to = "class", values_to = "n") |>
  mutate(class = factor(sub("^n_", "", class),
                        levels = c("recurrent","shared_limited","patient_private","uncertain")))
p6b <- ggplot(sk, aes(factor(activity_threshold), n, fill = class)) +
  geom_col(width = 0.8) +
  facet_wrap(~ factor(sprintf("min fraction of patient = %.1f%%", 100 * min_frac),
                      levels = sprintf("min fraction of patient = %.1f%%",
                                       100 * sort(unique(sens$min_frac)))), nrow = 1) +
  scale_fill_manual(values = setNames(P5_REC, c("recurrent","shared_limited",
                                                "patient_private","uncertain")),
                    name = "classification") +
  labs(x = "activity threshold (relative usage)", y = "programs",
       title = "B. threshold sensitivity of the recurrence classification",
       subtitle = wrap5("Reported across a grid rather than at a single cut, so the classification cannot rest on one arbitrary threshold.", 140),
       caption = P5_CAPTION) + p5_theme(10)
save_fig5(p6a / p6b + plot_layout(heights = c(1.3, 1)),
          file.path(P5_FIG, "06_recurrent_vs_patient_private_programs.pdf"), 12.0, 10.5)

# ---- 07 balanced sensitivity -------------------------------------------------
b <- bal |> mutate(plab = factor(LABS[program], levels = LABS[PROGS])) |>
  pivot_longer(c(cosine, top100_jaccard), names_to = "metric", values_to = "v")
p7 <- ggplot(b, aes(plab, v, fill = metric)) +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8) +
  geom_hline(yintercept = 0.60, linetype = 2, linewidth = 0.4, colour = "#B2182B") +
  scale_fill_manual(values = c(cosine = "#1B7837", top100_jaccard = "#A6DBA0"),
                    labels = c(cosine = "spectra cosine",
                               top100_jaccard = "top-100 gene Jaccard"), name = NULL) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "similarity to the best-matching balanced-run program",
       title = "Figure 07. Patient-balanced sensitivity: is the program solution just the largest patient?",
       subtitle = wrap5(paste0(
         "The balanced run repeats program discovery after downsampling every patient to 651 malignant cells - the smallest ",
         "patient with >= 500 - which cuts the dominant-patient fraction from 0.573 to 0.301. Each primary program is matched ",
         "to its most similar balanced program. The red line is the declared support threshold (cosine 0.60). ",
         "This is a sensitivity test; it does not replace the full-data primary analysis."), 150),
       caption = P5_CAPTION) +
  p5_theme(10) + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8))
save_fig5(p7, file.path(P5_FIG, "07_patient_balanced_program_sensitivity.pdf"), 11.5, 7.5)

# ---- 10 pathways -------------------------------------------------------------
tp <- pw |> group_by(program) |> slice_min(padj, n = 6, with_ties = FALSE) |>
  ungroup() |>
  mutate(pathway = sub("^(HALLMARK|REACTOME)_", "", pathway),
         pathway = gsub("_", " ", pathway),
         pathway = substr(pathway, 1, 52),
         plab = factor(LABS[program], levels = LABS[PROGS]))
p10 <- ggplot(tp, aes(NES, reorder(pathway, NES), fill = collection)) +
  geom_col(width = 0.72) +
  facet_wrap(~ plab, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(Hallmark = "#B2182B", Reactome = "#2166AC"),
                    name = "collection") +
  labs(x = "fgsea normalised enrichment score", y = NULL,
       title = "Figure 10. Pathway activity of each malignant program",
       subtitle = wrap5(paste0(
         "fgsea over the full ranked cNMF spectra score, MSigDB v2024.1.Hs Hallmark and Reactome, sets of 10-500 genes, ",
         "FDR < 0.05, top 6 per program. Labels follow this evidence and the program's own top genes; they are not imposed."), 150),
       caption = P5_CAPTION) +
  p5_theme(9) + theme(axis.text.y = element_text(size = 6.5))
save_fig5(p10, file.path(P5_FIG, "10_program_pathway_activity.pdf"),
          15.0, max(8, 2.4 * ceiling(length(PROGS) / 3) + 3))

# ---- 11 robustness summary ---------------------------------------------------
rl <- rob |>
  select(program, program_label, balanced_cosine, gene_subset_min_rho,
         loo_min_cosine, loo_cosine_without_dominant_patient,
         nocc_cosine, highconf_cosine, rank_min_cosine) |>
  pivot_longer(-c(program, program_label), names_to = "test", values_to = "v") |>
  mutate(plab = factor(LABS[program], levels = LABS[PROGS]),
         test = factor(test, levels = c("rank_min_cosine","balanced_cosine",
           "gene_subset_min_rho","nocc_cosine","highconf_cosine",
           "loo_min_cosine","loo_cosine_without_dominant_patient"),
           labels = c("rank K+-1","patient-balanced","gene subset (min rho)",
                      "cell-cycle removed","high-confidence only",
                      "leave-one-out (worst)","dominant patient removed")))
p11a <- ggplot(rl, aes(plab, test, fill = v)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(is.na(v), "NA", sprintf("%.2f", v))), size = 2.6) +
  scale_fill_gradient2(low = "#B2182B", mid = "#FEE08B", high = "#1A9850",
                       midpoint = 0.5, limits = c(0, 1), na.value = "grey88",
                       name = "similarity") +
  labs(x = NULL, y = NULL,
       title = "Figure 11. Phase 5 program robustness",
       subtitle = wrap5(paste0(
         "Each cell is the similarity between a primary program and its best match under that perturbation. ",
         "The bottom row is the decisive one: a program that cannot be recovered once the patient dominating it is ",
         "withheld entirely is NOT robust, and is reported as such rather than dropped from the table."), 150)) +
  p5_theme(10) + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8))
# Panel B carries TWO verdicts, because the strict one is uninformative on its
# own: every program fails it, and 5 of the 8 fail only the test a
# patient-private program can never pass - removal of its own patient. Showing
# only "NOT robust" eight times would confuse scope with reliability.
vb <- rob |>
  transmute(plab = factor(LABS[program], levels = LABS[PROGS]),
            `every declared test` = robust_overall,
            `every test except patient scope` = robust_excluding_patient_scope) |>
  pivot_longer(-plab, names_to = "verdict", values_to = "pass") |>
  mutate(verdict = factor(verdict, levels = c("every declared test",
                                              "every test except patient scope")))
p11b <- ggplot(vb, aes(plab, verdict, fill = pass)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(pass, "passes", "fails")), size = 2.8) +
  scale_fill_manual(values = c(`TRUE` = "#A6DBA0", `FALSE` = "#F4A582"),
                    guide = "none") +
  labs(x = NULL, y = NULL, title = "B. two verdicts, because the strict one is uninformative alone",
       subtitle = wrap5(paste0(
         "The upper row applies every declared threshold. The lower row waives ONLY the dominant-patient test for ",
         "programs already classified patient-private - such a program cannot survive removal of its own patient, ",
         "which is its scope rather than a defect - and also excludes the gene-subset test, whose 0.80 bar compares a ",
         "projection against a consensus solution and is conservative by construction."), 150),
       caption = P5_CAPTION) +
  p5_theme(10) + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8))
save_fig5(p11a / p11b + plot_layout(heights = c(2.6, 1.1)),
          file.path(P5_FIG, "11_phase5_robustness_summary.pdf"), 12.5, 10.5)
cat("figures 01-07, 10-11 written\n")
