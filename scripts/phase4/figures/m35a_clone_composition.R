#!/usr/bin/env Rscript
# =============================================================================
# Phase 4 - M35A - FIGURE 2
#   21_scevan_clone_composition_phase2.pdf
#   results/phase4/tables/final/SCEVAN_CLONE_COMPOSITION_VISUALIZATION.tsv
#
# The central figure of M35A. SCEVAN builds its subclones from inferred CNA
# profiles WITHOUT seeing any Phase 2 label, so the Phase 2 composition of a
# CNA-defined clone is the most direct malignancy evidence Phase 4 holds.
# Clone assignments are read from the frozen metadata; nothing is recomputed.
# =============================================================================
source("scripts/phase4/figures/m35a_common.R")

md <- m35a_load_meta()
m35a_msg("== M35A figure 21: clone composition by Phase 2 annotation ==")
m35a_validate(md)

cl <- md |>
  filter(!is.na(scevan_clone), scevan_clone != "") |>
  mutate(ann = m35a_group_ann(annotation_ccc))

# ---- clone-level attributes -------------------------------------------------
clone_stats <- cl |>
  group_by(sample_id, scevan_clone) |>
  summarise(clone_size = n(),
            n_refined_malignant = sum(malignancy_refined == "Malignant"),
            frac_refined_malignant = mean(malignancy_refined == "Malignant"),
            n_fibroblast = sum(annotation_ccc == "Fibroblast"),
            frac_fibroblast = mean(annotation_ccc == "Fibroblast"),
            scevan_reliable = unique(as.character(scevan_sample_reliable))[1],
            .groups = "drop") |>
  group_by(sample_id) |>
  mutate(dominant_annotation = NA_character_) |>
  ungroup()

# dominant Phase 2 identity, on the UNGROUPED annotation_ccc so "fibroblast-
# dominated" means what the frozen SCEVAN_CLONES.tsv means by it
dom <- cl |>
  count(sample_id, scevan_clone, annotation_ccc) |>
  group_by(sample_id, scevan_clone) |>
  slice_max(n, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(sample_id, scevan_clone, dominant_annotation = annotation_ccc,
         dominant_n = n)
clone_stats <- clone_stats |>
  select(-dominant_annotation) |>
  left_join(dom, by = c("sample_id", "scevan_clone")) |>
  mutate(fibroblast_dominated = dominant_annotation == "Fibroblast",
         disputed_triple = NA)

# does the clone carry all three disputed Phase 2 identities?
trip <- cl |>
  group_by(sample_id, scevan_clone) |>
  summarise(n_fib = sum(annotation_ccc == "Fibroblast"),
            n_cmu = sum(annotation_ccc == "Candidate-Malignant-Unresolved"),
            n_mt  = sum(annotation_ccc == "MPNST-Tumor"), .groups = "drop") |>
  mutate(disputed_triple = n_fib > 0 & n_cmu > 0 & n_mt > 0)
clone_stats <- clone_stats |> select(-disputed_triple) |>
  left_join(trip, by = c("sample_id", "scevan_clone"))

# ---- plotting / publication table -------------------------------------------
comp <- cl |>
  count(sample_id, scevan_clone, ann, name = "n_cells") |>
  complete(nesting(sample_id, scevan_clone), ann, fill = list(n_cells = 0L)) |>
  left_join(clone_stats, by = c("sample_id", "scevan_clone")) |>
  mutate(frac_of_clone = n_cells / clone_size,
         clone_label = sub("^MPNST_[0-9]+_clone", "clone ", scevan_clone)) |>
  arrange(sample_id, desc(clone_size), ann)

out_tsv <- file.path(M35A_TABF, "SCEVAN_CLONE_COMPOSITION_VISUALIZATION.tsv")
dir.create(M35A_TABF, showWarnings = FALSE, recursive = TRUE)
write.table(comp |>
  select(sample_id, scevan_clone, clone_label, clone_size, annotation_group = ann,
         n_cells, frac_of_clone, dominant_annotation, dominant_n,
         fibroblast_dominated, disputed_triple, n_fib, n_cmu, n_mt,
         n_refined_malignant, frac_refined_malignant, scevan_reliable),
  out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
m35a_msg("  [tsv] %s (%d rows)", out_tsv, nrow(comp))

# ---- console assertions on the claims the figure must make ------------------
fd <- clone_stats |> group_by(sample_id) |>
  summarise(n_clones = n(), n_fib_dom = sum(fibroblast_dominated),
            n_triple = sum(disputed_triple), .groups = "drop")
print(as.data.frame(fd))
stopifnot(fd$n_fib_dom[fd$sample_id == "MPNST_2"] == 4,
          fd$n_clones[fd$sample_id == "MPNST_2"] == 4,
          fd$n_fib_dom[fd$sample_id == "MPNST_4"] == 7,
          fd$n_clones[fd$sample_id == "MPNST_4"] == 8)

# ---- figure -----------------------------------------------------------------
pd <- comp |> filter(sample_id %in% M35A_RELIABLE) |>
  mutate(sample_id = factor(sample_id, levels = M35A_RELIABLE))
ord <- pd |> distinct(sample_id, scevan_clone, clone_label, clone_size) |>
  arrange(sample_id, desc(clone_size))
pd$clone_label <- factor(pd$clone_label, levels = unique(ord$clone_label))

hdr <- pd |> distinct(sample_id, clone_label, clone_size, frac_refined_malignant,
                      fibroblast_dominated, disputed_triple)

pA <- ggplot(pd, aes(x = clone_label, y = frac_of_clone, fill = ann)) +
  geom_col(width = 0.82, colour = "white", linewidth = 0.18) +
  geom_text(data = hdr, inherit.aes = FALSE,
            aes(x = clone_label, y = 1.015,
                label = sprintf("n=%s", format(clone_size, big.mark = ","))),
            size = 2.6, colour = "grey25", vjust = 0) +
  geom_point(data = hdr |> filter(fibroblast_dominated), inherit.aes = FALSE,
             aes(x = clone_label, y = -0.040), shape = 17, size = 2.3,
             colour = "#D95F02") +
  geom_point(data = hdr |> filter(disputed_triple), inherit.aes = FALSE,
             aes(x = clone_label, y = -0.088), shape = 15, size = 1.9,
             colour = "grey35") +
  scale_fill_manual(values = M35A_ANN_COL, name = "Phase 2 annotation_ccc",
                    drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     breaks = seq(0, 1, 0.25),
                     limits = c(-0.115, 1.075), expand = c(0, 0)) +
  facet_grid(~ sample_id, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "fraction of cells in the CNA-defined clone",
       title = "SCEVAN CNA-defined clones, composed by the frozen Phase 2 marker annotation",
       subtitle = wrap_sub(paste0(
         "SCEVAN built these clones from inferred copy number alone - it never saw a Phase 2 label. ",
         "Orange triangle = fibroblast-dominated clone; grey square = clone carrying all three disputed identities.\n",
         "MPNST_2: 4 of 4 clones fibroblast-dominated. MPNST_4: 7 of 8 (clone 8 is endothelial). ",
         "MPNST_1: all 7 clones mix Candidate-Malignant-Unresolved with MPNST-Tumor and 6 of 7 also carry Fibroblast cells ",
         "(clone 6 has none, clone 7 has one).\n",
         "MPNST_3 is omitted: its partition failed the sanity gate (figure 25)."), 148)) +
  p4_theme(10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "right")

pB <- ggplot(hdr, aes(x = clone_label, y = frac_refined_malignant)) +
  geom_col(width = 0.82, fill = "#B2182B", alpha = 0.85) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * frac_refined_malignant)),
            vjust = -0.35, size = 2.5, colour = "#7A1119") +
  geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.3, colour = "grey40") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.18),
                     breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  facet_grid(~ sample_id, scales = "free_x", space = "free_x") +
  labs(x = "SCEVAN clone", y = "fraction refined\nMalignant",
       caption = M35A_CAPTION) +
  p4_theme(10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_blank(), strip.background = element_blank())

p <- pA / pB + plot_layout(heights = c(3.4, 1), guides = "collect") &
  theme(legend.position = "right")
save_fig(p, file.path(M35A_FIG, "21_scevan_clone_composition_phase2.pdf"), 14.0, 8.8)
m35a_msg("== figure 21 done ==")
