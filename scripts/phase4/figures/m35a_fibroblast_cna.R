#!/usr/bin/env Rscript
# =============================================================================
# Phase 4 - M35A - FIGURE 3
#   22_fibroblast_cna_burden_by_patient.pdf
#   23_fibroblast_malignant_vs_nonmalignant_cna_profile.pdf
#
# Question: do Phase-2 `Fibroblast` cells that Phase 4 called Malignant carry
# CNA architecture closer to the malignant compartment than to fibroblasts that
# stayed Non-malignant?
#
# Reported as effect size and patient consistency, NOT as a pooled cell-level
# p-value: cells within a patient are not independent, so a pooled test on
# 5,064 cells would manufacture significance from n = 4 patients. Where a
# patient has too few cells on one side the comparison is printed as NOT
# EVALUABLE rather than pooled away.
# =============================================================================
source("scripts/phase4/figures/m35a_common.R")

md <- m35a_load_meta()
m35a_msg("== M35A figures 22-23: fibroblast CNA evidence ==")
m35a_validate(md)

MIN_N <- 20L   # fixed before looking: a group below this is NOT EVALUABLE

IMMUNE <- M35A_OTHER_IMMUNE
grp <- md |>
  mutate(group = case_when(
    annotation_ccc == "Fibroblast"  & malignancy_refined == "Malignant"      ~ "Fibroblast -> Malignant",
    annotation_ccc == "Fibroblast"  & malignancy_refined == "Non-malignant"  ~ "Fibroblast -> Non-malignant",
    annotation_ccc == "Fibroblast"  & malignancy_refined == "Ambiguous"      ~ "Fibroblast -> Ambiguous",
    annotation_ccc == "MPNST-Tumor" & malignancy_refined == "Malignant"      ~ "MPNST-Tumor -> Malignant",
    annotation_ccc %in% IMMUNE                                               ~ "Canonical immune (control)",
    TRUE ~ NA_character_)) |>
  filter(!is.na(group)) |>
  mutate(group = factor(group, levels = c(
    "Fibroblast -> Malignant", "Fibroblast -> Ambiguous", "Fibroblast -> Non-malignant",
    "MPNST-Tumor -> Malignant", "Canonical immune (control)")))

GCOL <- c("Fibroblast -> Malignant"     = "#D95F02",
          "Fibroblast -> Ambiguous"     = "#F0A202",
          "Fibroblast -> Non-malignant" = "#2166AC",
          "MPNST-Tumor -> Malignant"    = "#B2182B",
          "Canonical immune (control)"  = "grey55")

METRICS <- c(cnv_burden = "cnv_burden", cnv_frac_gain = "cnv_frac_gain",
             cnv_frac_loss = "cnv_frac_loss", cnv_mean_abs = "cnv_mean_abs")

long <- grp |>
  select(cell_id, sample_id, group, all_of(unname(METRICS))) |>
  pivot_longer(all_of(unname(METRICS)), names_to = "metric", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(metric = factor(metric, levels = unname(METRICS)))

# ---- effect sizes: Cliff's delta, malignant vs non-malignant fibroblasts -----
cliffs <- function(a, b) {           # +1 = a strictly greater than b
  if (length(a) < 1 || length(b) < 1) return(NA_real_)
  r <- rank(c(a, b)); n1 <- length(a); n2 <- length(b)
  U <- sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2
  2 * U / (n1 * n2) - 1
}
eff <- long |>
  filter(group %in% c("Fibroblast -> Malignant", "Fibroblast -> Non-malignant")) |>
  group_by(sample_id, metric) |>
  summarise(n_mal = sum(group == "Fibroblast -> Malignant"),
            n_non = sum(group == "Fibroblast -> Non-malignant"),
            med_mal = median(value[group == "Fibroblast -> Malignant"]),
            med_non = median(value[group == "Fibroblast -> Non-malignant"]),
            delta = cliffs(value[group == "Fibroblast -> Malignant"],
                           value[group == "Fibroblast -> Non-malignant"]),
            .groups = "drop") |>
  mutate(evaluable = n_mal >= MIN_N & n_non >= MIN_N,
         reliable  = sample_id %in% M35A_RELIABLE,
         delta_plot = ifelse(evaluable, delta, NA_real_))
m35a_msg("-- Cliff's delta, Fibroblast Malignant vs Non-malignant --")
print(as.data.frame(eff))
write.table(eff, file.path(M35A_TABF, "M35A_FIBROBLAST_CNA_EFFECT_SIZES.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

nlab <- long |> count(sample_id, group, metric) |>
  filter(metric == "cnv_burden") |>
  mutate(lab = format(n, big.mark = ","))

unrel <- data.frame(sample_id = factor("MPNST_3", levels = M35A_SAMPLES))

nlab <- long |> count(sample_id, group, metric, name = "n") |>
  left_join(long |> group_by(metric) |> summarise(ymin = min(value), .groups = "drop"),
            by = "metric")

pA <- ggplot(long, aes(x = group, y = value, fill = group)) +
  geom_rect(data = unrel, inherit.aes = FALSE, xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, fill = "#FDECEC", alpha = 0.9) +
  geom_violin(scale = "width", width = 0.85, linewidth = 0.2, colour = "grey30",
              alpha = 0.75, na.rm = TRUE) +
  geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.3,
               fill = "white", alpha = 0.9, na.rm = TRUE) +
  geom_text(data = nlab, inherit.aes = FALSE,
            aes(x = group, y = ymin, label = format(n, big.mark = ",")),
            size = 2.2, colour = "grey30", vjust = 1.6) +
  scale_fill_manual(values = GCOL, name = NULL, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.05))) +
  facet_grid(metric ~ sample_id, scales = "free_y", switch = "y") +
  labs(x = NULL, y = "Phase 4 cell-level CNA metric",
       title = "Cell-level CNA metrics of Phase-2 fibroblast-labelled cells, stratified by patient",
       subtitle = wrap_sub(paste0(
         "Groups are Phase 2 annotation_ccc -> Phase 4 malignancy_refined. Distributions are shown per patient, never pooled; ",
         "cell counts are printed below each violin.\n",
         "MPNST_3 (pink panel) is UNRELIABLE: its SCEVAN partition failed the immune sanity gate and malignant promotions were disabled, ",
         "so it contributes no Fibroblast -> Malignant cells.\n",
         "MPNST_4 retains only 1 Non-malignant fibroblast, so its malignant/non-malignant fibroblast contrast is NOT EVALUABLE."), 155)) +
  p4_theme(10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom", strip.placement = "outside")

pB <- ggplot(eff, aes(x = sample_id, y = delta_plot, fill = sample_id)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey30") +
  geom_col(width = 0.75, na.rm = TRUE) +
  geom_text(data = eff |> filter(evaluable),
            aes(x = sample_id, y = delta_plot, label = sprintf("%.2f", delta_plot)),
            vjust = -0.45, size = 2.6, colour = "grey15") +
  geom_text(data = eff |> filter(!evaluable),
            aes(x = sample_id, y = 0.03,
                label = sprintf("NOT EVALUABLE (n=%d vs %d)", n_mal, n_non)),
            angle = 90, hjust = 0, size = 2.3, colour = "grey35") +
  scale_fill_manual(values = P4_SAMPLE, guide = "none", drop = FALSE) +
  scale_y_continuous(limits = c(0, 1.06), expand = c(0, 0)) +
  facet_grid(~ metric) +
  labs(x = NULL, y = "Cliff's delta\n(Fibroblast Malignant vs Non-malignant)",
       title = "Effect size, per patient - not a pooled cell-level test",
       subtitle = wrap_sub(paste0(
         "delta = +1 means every malignant-called fibroblast exceeds every non-malignant one. Bars are drawn only where both groups have >= ",
         MIN_N, " cells, so the contrast exists in MPNST_1 and MPNST_2 only. ",
         "Where it exists it is large (0.78-0.92) and consistent across all four metrics and both patients."), 155),
       caption = M35A_CAPTION) +
  p4_theme(10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8))

p <- pA / pB + plot_layout(heights = c(2.9, 1))
save_fig(p, file.path(M35A_FIG, "22_fibroblast_cna_burden_by_patient.pdf"), 13.5, 12.0)

# =============================================================================
# Figure 23 - genome-wide mean inferred CNA profile per group, per patient.
# Built from the STORED native SCEVAN CNA matrices. SCEVAN is not re-run.
# =============================================================================
prof_groups <- c("Fibroblast -> Malignant", "Fibroblast -> Non-malignant",
                 "MPNST-Tumor -> Malignant", "Canonical immune (control)")
gmap <- setNames(as.character(grp$group), grp$cell_id)

prof <- list(); corr <- list()
for (s in M35A_RELIABLE) {
  fcna <- file.path(M35A_SCEVAN, s, "primary", "output",
                    sprintf("%s_primary_CNAmtx.RData", s))
  fann <- file.path(M35A_SCEVAN, s, "primary", "output",
                    sprintf("%s_primary_count_mtx_annot.RData", s))
  if (!file.exists(fcna) || !file.exists(fann)) {
    m35a_msg("  [skip] %s: stored CNA matrix missing", s); next
  }
  e <- new.env(); load(fcna, envir = e); load(fann, envir = e)
  M   <- get("CNA_mtx_relat", e)
  ann <- get("count_mtx_annot", e)
  stopifnot(nrow(M) == nrow(ann))
  keep <- ann$seqnames %in% as.character(1:22)
  M <- M[keep, , drop = FALSE]; ann <- ann[keep, , drop = FALSE]
  ann$chr <- factor(ann$seqnames, levels = as.character(1:22))
  o <- order(ann$chr, ann$start); M <- M[o, , drop = FALSE]; ann <- ann[o, ]
  # a shared genome coordinate, so patients that retain different numbers of
  # genes are still drawn on the same axis and remain directly comparable
  ann$idx <- M35A_CHROFF[as.character(ann$chr)] + ann$start

  g <- gmap[colnames(M)]
  means <- list()
  for (gg in prof_groups) {
    cid <- which(!is.na(g) & g == gg)
    if (length(cid) < MIN_N) { m35a_msg("  [%s] %-30s n=%d  NOT EVALUABLE", s, gg, length(cid)); next }
    means[[gg]] <- rowMeans(M[, cid, drop = FALSE])
    prof[[length(prof) + 1L]] <- data.frame(
      sample_id = s, group = gg, idx = ann$idx, chr = ann$chr,
      cna = means[[gg]], n_cells = length(cid))
    m35a_msg("  [%s] %-30s n=%d", s, gg, length(cid))
  }
  # does the malignant-called fibroblast profile sit closer to the malignant
  # compartment or to the non-malignant fibroblasts?
  if (!is.null(means[["Fibroblast -> Malignant"]])) {
    for (ref in c("MPNST-Tumor -> Malignant", "Fibroblast -> Non-malignant",
                  "Canonical immune (control)")) {
      if (is.null(means[[ref]])) next
      corr[[length(corr) + 1L]] <- data.frame(
        sample_id = s, reference = ref,
        pearson = cor(means[["Fibroblast -> Malignant"]], means[[ref]]))
    }
  }
  rm(M, e); invisible(gc())
}
prof <- do.call(rbind, prof); corr <- do.call(rbind, corr)

if (is.null(prof)) {
  m35a_msg("  figure 23 NOT produced: no stored CNA matrix yielded an evaluable group")
} else {
  bnd <- data.frame(chr = factor(names(M35A_CHRLEN), levels = names(M35A_CHRLEN)),
                    mx = M35A_CHROFF + M35A_CHRLEN, mid = M35A_CHRMID)
  m35a_msg("-- profile correlation with Fibroblast -> Malignant --")
  print(as.data.frame(corr))
  write.table(corr, file.path(M35A_TABF, "M35A_FIBROBLAST_CNA_PROFILE_CORRELATION.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  cl_lab <- corr |>
    mutate(t = sprintf("r(Fib->Mal, %s) = %.3f", reference, pearson)) |>
    group_by(sample_id) |>
    summarise(lab = paste(t, collapse = "\n"), .groups = "drop") |>
    mutate(idx = 0, cna = Inf)

  p23 <- ggplot(prof, aes(x = idx, y = cna, colour = group)) +
    geom_vline(data = bnd, aes(xintercept = mx), colour = "grey85",
               linewidth = 0.25, inherit.aes = FALSE) +
    geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.3) +
    geom_line(linewidth = 0.45, alpha = 0.9) +
    geom_text(data = cl_lab, aes(x = idx, y = cna, label = lab),
              inherit.aes = FALSE, hjust = 0, vjust = 1.1, size = 2.5,
              colour = "grey25", lineheight = 1.05) +
    scale_colour_manual(values = GCOL, name = NULL, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
    scale_x_continuous(breaks = M35A_CHRMID, labels = names(M35A_CHRLEN),
                       limits = c(0, M35A_GTOT), expand = c(0.004, 0)) +
    facet_wrap(~ sample_id, ncol = 1) +
    labs(x = "chromosome (hg38 layout, chr1-22)",
         y = "mean inferred relative CNA",
         title = "Genome-wide inferred CNA profile: fibroblast-labelled cells called Malignant vs Non-malignant",
         subtitle = wrap_sub(paste0(
           "Mean per-gene relative CNA from the STORED native SCEVAN CNA matrices (primary run). SCEVAN was not re-run.\n",
           "Positive = inferred relative gain, negative = inferred relative loss. Broad segments only; no single-gene CNV is asserted.\n",
           "A group is drawn only where the patient has >= ", MIN_N, " such cells - MPNST_4 has 1 non-malignant fibroblast, so that ",
           "contrast exists only in MPNST_1 and MPNST_2."), 150),
         caption = M35A_CAPTION) +
    p4_theme(10) + theme(legend.position = "bottom",
                         axis.text.x = element_text(size = 7))
  save_fig(p23, file.path(M35A_FIG,
    "23_fibroblast_malignant_vs_nonmalignant_cna_profile.pdf"), 13.0, 9.5, png = FALSE)
}
m35a_msg("== figures 22-23 done ==")
