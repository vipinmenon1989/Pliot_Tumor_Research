#!/usr/bin/env Rscript
# =============================================================================
# Phase 4 - M35A - FIGURE 6
#   26_fibroblast_malignancy_threshold_robustness.pdf
#
# Separates the robust claim from the fragile one. Fibroblast -> Malignant is
# 4,036 at every tested pop_frac_low; the retained MPNST-Tumor count, and with
# it the 32.63% cohort malignant fraction, is not.
# =============================================================================
source("scripts/phase4/figures/m35a_common.R")

m35a_msg("== M35A figure 26: threshold robustness ==")
ts <- read.delim(file.path(M35A_TAB, "MALIGNANCY_THRESHOLD_SENSITIVITY.tsv"),
                 check.names = FALSE)
print(ts)

stopifnot(all(ts$fibroblast_malignant == 4036),
          all(ts$candidate_malignant == 836),
          ts$Malignant[ts$is_a_priori == "TRUE" | ts$is_a_priori == TRUE] == 6434)
m35a_msg("  [OK] Fibroblast -> Malignant is 4036 at every tested pop_frac_low")

apx <- ts$pop_frac_low[as.character(ts$is_a_priori) %in% c("TRUE", "True")][1]

d <- ts |>
  select(pop_frac_low,
         `total refined Malignant`        = Malignant,
         `retained original MPNST-Tumor`  = mpnst_tumor_still_malignant,
         `Fibroblast -> Malignant`        = fibroblast_malignant) |>
  pivot_longer(-pop_frac_low, names_to = "series", values_to = "cells") |>
  mutate(series = factor(series, levels = c("total refined Malignant",
    "retained original MPNST-Tumor", "Fibroblast -> Malignant")))

SCOL <- c("total refined Malignant" = "#444444",
          "retained original MPNST-Tumor" = "#B2182B",
          "Fibroblast -> Malignant" = "#D95F02")

rng <- d |> group_by(series) |>
  summarise(lo = min(cells), hi = max(cells),
            stable = lo == hi, .groups = "drop") |>
  mutate(lab = ifelse(stable,
    sprintf("%s: %s at every threshold - THRESHOLD-STABLE", series, format(lo, big.mark = ",")),
    sprintf("%s: %s to %s - THRESHOLD-SENSITIVE", series,
            format(lo, big.mark = ","), format(hi, big.mark = ","))))
m35a_msg("%s", paste(rng$lab, collapse = "\n"))

pA <- ggplot(d, aes(x = pop_frac_low, y = cells, colour = series)) +
  geom_vline(xintercept = apx, linetype = 2, linewidth = 0.4, colour = "grey45") +
  annotate("text", x = apx, y = Inf, label = "  a priori 0.25", hjust = 0,
           vjust = 1.6, size = 3, colour = "grey35") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.6) +
  geom_text(aes(label = format(cells, big.mark = ",")), vjust = -1.1,
            size = 2.8, show.legend = FALSE) +
  scale_colour_manual(values = SCOL, name = NULL) +
  scale_x_continuous(breaks = ts$pop_frac_low, limits = c(0.13, 0.42)) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0.08, 0.16))) +
  labs(x = "pop_frac_low (population malignant-fraction threshold)", y = "cells",
       title = "A. What moves with the threshold, and what does not",
       subtitle = wrap_sub(paste0(
         "The a priori value 0.25 was fixed in MALIGNANCY_DECISION_RULES.md before any SCEVAN result was inspected. ",
         "Fibroblast -> Malignant is a flat line at 4,036 across 0.15-0.40. The retained MPNST-Tumor count steps from 3,266 to 1,405 ",
         "because MPNST_1's pop_frac is 0.222, just under 0.25."), 125)) +
  p4_theme(11) + theme(legend.position = "bottom")

frac <- ts |>
  mutate(`refined malignant fraction of 19,716` = Malignant / 19716)
pB <- ggplot(frac, aes(x = pop_frac_low, y = `refined malignant fraction of 19,716`)) +
  geom_vline(xintercept = apx, linetype = 2, linewidth = 0.4, colour = "grey45") +
  geom_line(linewidth = 0.9, colour = "#444444") +
  geom_point(size = 2.6, colour = "#444444") +
  geom_text(aes(label = sprintf("%.2f%%", 100 * Malignant / 19716)),
            vjust = -1.2, size = 2.8, colour = "grey25") +
  geom_hline(yintercept = 0.1735, linetype = 3, linewidth = 0.4, colour = "#2166AC") +
  annotate("text", x = 0.40, y = 0.1735, label = "Phase 2 conservative 17.35%",
           hjust = 1, vjust = -0.6, size = 2.8, colour = "#2166AC") +
  scale_x_continuous(breaks = ts$pop_frac_low, limits = c(0.13, 0.42)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0.15, 0.47)) +
  labs(x = "pop_frac_low", y = "refined malignant fraction",
       title = "B. The cohort malignant fraction is the fragile number",
       subtitle = wrap_sub(paste0(
         "32.63% at the a priori threshold, 42.07% at 0.15-0.20. It is also not patient-robust: dropping MPNST_4 gives 0.214 ",
         "against a Phase 2 fraction of 0.215 (fold change 1.00). ",
         "The fibroblast conclusion and the cohort fraction are separate claims. Only the first is robust."), 125),
       caption = M35A_CAPTION) +
  p4_theme(11)

p <- pA / pB + plot_layout(heights = c(1.35, 1)) +
  plot_annotation(
    title = "Figure 26. Threshold robustness of the fibroblast conclusion",
    subtitle = wrap_sub(paste0("4,036 Phase-2 fibroblast-labelled cells are refined Malignant at every tested value of pop_frac_low ",
                      "(0.15, 0.20, 0.25, 0.30, 0.40). Candidate-Malignant-Unresolved is likewise flat at 836."), 130),
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 10, colour = "grey20")))
save_fig(p, file.path(M35A_FIG, "26_fibroblast_malignancy_threshold_robustness.pdf"),
         11.5, 11.0)
m35a_msg("== figure 26 done ==")
