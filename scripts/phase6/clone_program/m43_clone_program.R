#!/usr/bin/env Rscript
# =============================================================================
# Phase 6 - M43 - clone <-> program coupling, WITHIN each reliable patient.
#
# Clone labels are patient-scoped: MPNST_1_clone1 and MPNST_4_clone1 are
# unrelated names. Every association is therefore computed inside one patient
# and never pooled across patients as if clones matched.
#
# Association strength is eta^2 (the fraction of a program's variance that
# sits between clones rather than within them), with a permutation null built
# by shuffling clone labels INSIDE the patient - which preserves patient
# structure, clone sizes and the program score distribution.
# =============================================================================
source("scripts/phase6/utils/phase6_common.R")
set.seed(42)
NPERM <- 1000L
facts <- list(milestone = "M43", n_perm = NPERM,
              generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))

md <- readRDS(file.path(P6_CLONE, "phase6_cell_metadata.rds"))
PROGS <- grep("^program_P[0-9]+_score$", colnames(md), value = TRUE)
PID <- sub("^program_(P[0-9]+)_score$", "\\1", PROGS)
lab <- read.delim("results/phase5/tables/final/MALIGNANT_PROGRAMS.tsv")
LABEL <- setNames(lab$program_label, lab$program_id)
p6_msg("programs: %s", paste(sprintf("%s=%s", PID, LABEL[PID]), collapse = "  "))

mal <- md[md$malignancy_refined == "Malignant" & !is.na(md$tumor_clone_phase4), ]
p6_msg("malignant cells with a clone label: %d", nrow(mal))
mal_rel <- mal[mal$sample_id %in% CLONE_RELIABLE, ]
p6_msg("of which in clone-reliable patients (%s): %d",
       paste(CLONE_RELIABLE, collapse = ", "), nrow(mal_rel))
p6_msg("EXCLUDED from every clone analysis: %s (%d cells) - %s",
       CLONE_EXCLUDED, sum(mal$sample_id == CLONE_EXCLUDED), CLONE_EXCLUSION_REASON)

p6_sec("1. Program activity by clone (reliable patients only)")
act <- mal_rel |>
  select(cell_id, sample_id, clone = tumor_clone_phase4, all_of(PROGS)) |>
  pivot_longer(all_of(PROGS), names_to = "program", values_to = "usage") |>
  mutate(program = sub("^program_(P[0-9]+)_score$", "\\1", program))
cpa <- act |> group_by(sample_id, clone, program) |>
  summarise(n_cells = n(), median_usage = median(usage), mean_usage = mean(usage),
            sd_usage = sd(usage), .groups = "drop") |>
  mutate(program_label = unname(LABEL[program]))
p6_tsv(cpa, file.path(P6_TAB, "CLONE_PROGRAM_ACTIVITY.tsv"))

dom <- mal_rel |> rename(clone = tumor_clone_phase4) |>
  group_by(sample_id, clone) |>
  summarise(n_cells = n(),
            dominant_program = names(sort(table(dominant_malignant_program),
                                          decreasing = TRUE))[1],
            dominant_program_share = max(table(dominant_malignant_program)) / n(),
            n_programs_present = length(unique(dominant_malignant_program)),
            .groups = "drop") |>
  mutate(dominant_program_label = unname(LABEL[dominant_program]),
         evaluable_ge20 = n_cells >= CLONE_SIZE_PRIMARY)
print(as.data.frame(dom), row.names = FALSE, digits = 3)
p6_tsv(dom, file.path(P6_CLONE, "M43_CLONE_DOMINANT_PROGRAM.tsv"))

p6_sec("2. Clone-program association within each patient (eta^2 + permutation)")
eta2 <- function(y, g) {
  gm <- tapply(y, g, mean); n <- table(g)
  ssb <- sum(n * (gm - mean(y))^2); sst <- sum((y - mean(y))^2)
  if (sst <= 0) return(NA_real_)
  ssb / sst
}
assoc <- list()
for (s in CLONE_RELIABLE) {
  d <- mal_rel[mal_rel$sample_id == s, ]
  # clones below the declared minimum are reported, not silently dropped
  keep <- names(which(table(d$tumor_clone_phase4) >= CLONE_SIZE_PRIMARY))
  small <- setdiff(unique(d$tumor_clone_phase4), keep)
  if (length(small))
    p6_msg("  %s: clones below %d cells reported NOT EVALUABLE: %s", s,
           CLONE_SIZE_PRIMARY, paste(small, collapse = ", "))
  d <- d[d$tumor_clone_phase4 %in% keep, ]
  g <- factor(d$tumor_clone_phase4)
  for (i in seq_along(PROGS)) {
    y <- d[[PROGS[i]]]
    e <- eta2(y, g)
    null <- replicate(NPERM, eta2(y, sample(g)))
    assoc[[length(assoc) + 1L]] <- data.frame(
      sample_id = s, program = PID[i], program_label = unname(LABEL[PID[i]]),
      n_cells = nrow(d), n_clones = nlevels(g), eta2 = e,
      eta2_null_mean = mean(null, na.rm = TRUE),
      eta2_null_p95 = quantile(null, 0.95, na.rm = TRUE, names = FALSE),
      perm_p = (1 + sum(null >= e, na.rm = TRUE)) / (NPERM + 1),
      excess_over_null = e - mean(null, na.rm = TRUE))
  }
  p6_msg("  %s: %d clones >= %d cells, %d cells used", s, nlevels(g),
         CLONE_SIZE_PRIMARY, nrow(d))
}
assoc <- do.call(rbind, assoc)
assoc$clone_associated <- assoc$eta2 >= 0.10 & assoc$perm_p < 0.01
p6_tsv(assoc, file.path(P6_TAB, "CLONE_PROGRAM_ASSOCIATION.tsv"))
print(as.data.frame(assoc |> select(sample_id, program, program_label, eta2,
                                    eta2_null_p95, perm_p, clone_associated)),
      row.names = FALSE, digits = 3)

p6_sec("3. How much program variance is between clones at all?")
summ <- assoc |> group_by(sample_id) |>
  summarise(n_programs = n(), median_eta2 = median(eta2, na.rm = TRUE),
            max_eta2 = max(eta2, na.rm = TRUE),
            program_at_max = program[which.max(eta2)],
            n_clone_associated = sum(clone_associated, na.rm = TRUE),
            .groups = "drop")
print(as.data.frame(summ), row.names = FALSE, digits = 3)
p6_tsv(summ, file.path(P6_CLONE, "M43_ASSOCIATION_SUMMARY_BY_PATIENT.tsv"))
p6_msg("Interpretation guard: eta^2 is the share of a program's variance sitting BETWEEN clones.")
p6_msg("Even a strongly significant eta^2 of 0.10 leaves 90%% of the variance WITHIN clones.")

p6_sec("4. Does one program appear across several distinct clones?")
spread <- cpa |> filter(n_cells >= CLONE_SIZE_PRIMARY) |>
  group_by(sample_id, program) |>
  summarise(n_clones_evaluable = n(),
            n_clones_active = sum(median_usage >= 0.20),
            max_median = max(median_usage), min_median = min(median_usage),
            range_median = max(median_usage) - min(median_usage),
            .groups = "drop") |>
  mutate(program_label = unname(LABEL[program]))
p6_tsv(spread, file.path(P6_CLONE, "M43_PROGRAM_SPREAD_ACROSS_CLONES.tsv"))
print(as.data.frame(spread), row.names = FALSE, digits = 3)

facts$association <- assoc; facts$summary <- summ
facts$clone_dominant <- dom
facts$excluded <- list(patient = CLONE_EXCLUDED, reason = CLONE_EXCLUSION_REASON)
facts$eta2_note <- "eta^2 = between-clone share of program variance; a permutation null was built by shuffling clone labels WITHIN each patient, preserving clone sizes and the score distribution"
p6_json(facts, file.path(P6_VAL, "m43_clone_program_facts.json"))
p6_sec("M43 complete")
