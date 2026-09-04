#!/usr/bin/env Rscript
# =============================================================================
# Phase 4 - M35A - FIGURE 4
#   24_phase2_to_phase4_malignancy_transition.pdf
#
# Only the disputed compartments. The canonical immune populations are large
# and their calls were never in question, so including them would bury the
# transition this figure exists to show.
# =============================================================================
source("scripts/phase4/figures/m35a_common.R")
suppressPackageStartupMessages(library(ggalluvial))

md <- m35a_load_meta()
m35a_msg("== M35A figure 24: Phase 2 -> Phase 4 identity transition ==")
m35a_validate(md)

DISPUTED <- c("Fibroblast", "Candidate-Malignant-Unresolved", "MPNST-Tumor",
              "Pericyte-VSMC", "Endothelial", "Uncertain")

d <- md |>
  filter(annotation_ccc %in% DISPUTED) |>
  count(annotation_ccc, malignancy_refined, name = "n") |>
  mutate(annotation_ccc = factor(annotation_ccc, levels = DISPUTED))

tot <- d |> group_by(annotation_ccc) |> summarise(n = sum(n), .groups = "drop")
m35a_msg("-- flows plotted --"); print(as.data.frame(d))

# the numbers this figure must make explicit
chk <- function(a, b, want) {
  got <- sum(d$n[d$annotation_ccc == a & d$malignancy_refined == b])
  stopifnot(got == want); m35a_msg("  [OK] %-32s -> %-14s %d", a, b, got)
}
chk("Fibroblast", "Malignant", 4036); chk("Fibroblast", "Non-malignant", 908)
chk("Fibroblast", "Ambiguous", 120)
chk("Candidate-Malignant-Unresolved", "Malignant", 836)
chk("Candidate-Malignant-Unresolved", "Ambiguous", 395)
chk("MPNST-Tumor", "Malignant", 1405); chk("MPNST-Tumor", "Ambiguous", 2015)

ann_col <- M35A_ANN_COL[c("Fibroblast", "Candidate-Malignant-Unresolved",
                          "MPNST-Tumor", "Pericyte-VSMC", "Endothelial")]
ann_col["Uncertain"] <- "grey60"

# Stratum labels are placed OUTSIDE the strata so the flow counts inside stay
# readable. ggalluvial stacks strata with the first factor level on top, so the
# midpoints are computed directly and asserted against the plotted totals.
LAB_MIN <- 100L
stratum_mid <- function(df, key) {
  tot <- df |> group_by(.data[[key]]) |> summarise(n = sum(n), .groups = "drop") |>
    arrange(.data[[key]])
  top <- sum(tot$n) - c(0, cumsum(tot$n)[-nrow(tot)])
  data.frame(stratum = as.character(tot[[key]]), n = tot$n,
             y = top - tot$n / 2)
}
lab1 <- stratum_mid(d, "annotation_ccc") |> mutate(x = 1 - 0.17, hj = 1)
lab2 <- stratum_mid(d, "malignancy_refined") |> mutate(x = 2 + 0.17, hj = 0)
labs_all <- rbind(lab1, lab2) |>
  mutate(txt = sprintf("%s\nn = %s", stratum, format(n, big.mark = ",")))
stopifnot(sum(lab1$n) == sum(d$n), sum(lab2$n) == sum(d$n))

p <- ggplot(d, aes(y = n, axis1 = annotation_ccc, axis2 = malignancy_refined)) +
  geom_alluvium(aes(fill = annotation_ccc), width = 0.24, alpha = 0.72,
                knot.pos = 0.32, curve_type = "sigmoid") +
  geom_stratum(width = 0.24, fill = "grey96", colour = "grey35",
               linewidth = 0.3) +
  geom_text(data = labs_all, inherit.aes = FALSE,
            aes(x = x, y = y, label = txt, hjust = hj),
            size = 3.0, lineheight = 0.95, colour = "grey10") +
  geom_text(data = d |> filter(n >= LAB_MIN), stat = "alluvium",
            aes(label = format(n, big.mark = ",")), size = 2.8,
            colour = "grey10", fontface = "bold") +
  scale_fill_manual(values = ann_col, name = "Phase 2 annotation_ccc") +
  scale_x_continuous(breaks = c(1, 2),
                     labels = c("Phase 2 annotation_ccc", "Phase 4 malignancy_refined"),
                     limits = c(0.35, 2.65)) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0.02, 0.04))) +
  labs(y = "cells", x = NULL,
       title = "Phase 2 marker identity to Phase 4 copy-number-refined malignancy - disputed compartments only",
       subtitle = wrap_sub(paste0(
         "Restricted to the six populations whose malignant status marker expression could not settle ",
         "(", format(sum(d$n), big.mark = ","), " of 19,716 cells). ",
         "Canonical immune populations are deliberately excluded: their calls were never in dispute and their size would bury this transition.\n",
         "Fibroblast: 4,036 Malignant, 908 Non-malignant, 120 Ambiguous. ",
         "Candidate-Malignant-Unresolved: 836 Malignant, 395 Ambiguous, 0 Non-malignant. ",
         "MPNST-Tumor: 1,405 Malignant, 2,015 Ambiguous.\n",
         "Flow counts are printed for flows of at least ", LAB_MIN, " cells; every flow is tabulated in PHASE4_MALIGNANCY_CALLS.tsv. ",
         "Ambiguous is a terminal outcome: lineage and copy-number evidence disagree and neither is allowed to win."), 148),
       caption = M35A_CAPTION) +
  p4_theme(11) +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(face = "bold", size = 10))

save_fig(p, file.path(M35A_FIG, "24_phase2_to_phase4_malignancy_transition.pdf"),
         13.5, 9.5)
m35a_msg("== figure 24 done ==")
