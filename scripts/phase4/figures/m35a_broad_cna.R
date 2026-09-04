#!/usr/bin/env Rscript
# =============================================================================
# Phase 4 - M35A - OPTIONAL FIGURE
#   27_broad_cna_recurrence.pdf
#
# Secondary evidence. Broad recurrent patterns only, drawn from the stored
# clonal segment calls; no gene-level CNV claim is made anywhere on this page.
# Where the NF1 / NF2 loci are marked, the wording is "broad segment containing
# the <gene> locus" - this analysis cannot distinguish a focal deletion from a
# large-segment loss that happens to include the locus, and does not try.
# =============================================================================
source("scripts/phase4/figures/m35a_common.R")

m35a_msg("== M35A figure 27: broad CNA recurrence ==")
seg <- read.delim(file.path(M35A_TAB, "SCEVAN_CNV_SUMMARY.tsv"), check.names = FALSE)
seg <- seg |> filter(level == "clonal", broad %in% c(TRUE, "TRUE"), event != "neutral",
                     sample_id %in% M35A_RELIABLE, Chr %in% 1:22)
m35a_msg("  broad clonal non-neutral segments in the reliable patients: %d", nrow(seg))

CHRLEN <- M35A_CHRLEN; off <- M35A_CHROFF
mid <- M35A_CHRMID; gtot <- M35A_GTOT

seg <- seg |>
  mutate(chr = as.character(Chr),
         gs = off[chr] + Pos, ge = off[chr] + End,
         sample_id = factor(sample_id, levels = rev(M35A_RELIABLE)),
         event = factor(event, levels = c("gain", "loss")))

# recurrent core: present in all three reliable patients
rec <- read.delim(file.path(M35A_TAB, "SCEVAN_RECURRENT_BROAD_EVENTS.tsv"),
                  check.names = FALSE) |>
  filter(n_patients_reliable3 == 3)
m35a_msg("-- recurrent in 3/3 reliable patients --"); print(rec)
rec <- rec |> mutate(chr = as.character(Chr),
                     xmin = off[chr], xmax = off[chr] + CHRLEN[chr],
                     lab = sprintf("chr%s %s\n3/3 reliable patients", Chr, event))

# loci marked as CONTAINED IN a broad segment - never as a focal deletion
loci <- data.frame(
  gene = c("NF1", "NF2"), chr = c("17", "22"),
  pos  = c(31094927, 29603556))
loci <- loci |> mutate(g = off[chr] + pos,
                       lab = sprintf("broad segment containing\nthe %s locus", gene),
                       hj = 1, ylab = c(3.84, 3.42))  # staggered: the two loci are only 130 Mb apart

ECOL <- c(gain = "#B2182B", loss = "#2166AC")

p <- ggplot() +
  geom_rect(data = rec, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "#FFF3C4", alpha = 0.75) +
  geom_vline(xintercept = c(off, gtot), colour = "grey86", linewidth = 0.25) +
  geom_segment(data = seg,
               aes(x = gs, xend = ge, y = sample_id, yend = sample_id,
                   colour = event), linewidth = 5.4, lineend = "butt") +
  geom_vline(data = loci, aes(xintercept = g), linetype = 2, linewidth = 0.4,
             colour = "grey25") +
  geom_text(data = loci, aes(x = g, y = ylab, label = lab, hjust = hj), size = 2.5,
            colour = "grey20", lineheight = 0.95, nudge_x = -8e6) +
  geom_text(data = rec, aes(x = (xmin + xmax) / 2, y = 0.36, label = lab),
            size = 2.5, colour = "#8A6D00", lineheight = 0.95) +
  scale_colour_manual(values = ECOL, name = "inferred broad event",
                      labels = c(gain = "relative gain", loss = "relative loss")) +
  scale_x_continuous(breaks = mid, labels = names(CHRLEN),
                     limits = c(0, gtot), expand = c(0.004, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.75, 1.35))) +
  labs(x = "chromosome (hg38 layout, chr1-22)", y = NULL,
       title = "Figure 27. Broad recurrent inferred copy-number architecture in the three reliable patients",
       subtitle = wrap_sub(paste0(
         "Clonal segments >= 10 Mb and non-neutral only (", nrow(seg), " segments). Shaded bands mark the shared core: ",
         "chr18 loss, chr2 gain and chr7 gain, each present in 3 of 3 reliable patients.\n",
         "The architecture is otherwise patient-specific - pairwise Jaccard of broad-event sets is only 0.20-0.23 and just 1 of 59 events is shared by all four patients.\n",
         "NO gene-level CNV claim is made. The NF1 and NF2 marks indicate only that a broad segment spans the locus; distinguishing a focal deletion from ",
         "large-segment loss requires DNA sequencing.\n",
         "MPNST_3 is excluded (failed sanity gate). chrX and chrY are not assessable by this method."), 175),
       caption = M35A_CAPTION) +
  p4_theme(11) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.text.x = element_text(size = 7.5),
        legend.position = "bottom")

save_fig(p, file.path(M35A_FIG, "27_broad_cna_recurrence.pdf"), 14.0, 6.2, png = FALSE)
m35a_msg("== figure 27 done ==")
