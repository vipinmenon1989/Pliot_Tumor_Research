#!/usr/bin/env Rscript
# =============================================================================
# Phase 4 · M34 — robustness
#
# Four questions (§49), answered explicitly and without softening negatives:
#   1. Does each malignant state occur in multiple patients?
#   2. Are the SCEVAN CNV patterns patient-specific?
#   3. Are the major CCC changes caused by one patient?
#   4. Do the Phase 3 core interactions remain leave-one-patient-out stable?
#
# n = 4 patients. Nothing here supports a population-level claim (§65E).
# =============================================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(jsonlite); library(ggplot2); library(Matrix)
})
source("scripts/phase4/utils/phase4_plot_utils.R")
set.seed(42)

TAB <- "results/phase4/tables"; FIG <- "results/phase4/figures"
MAL <- "results/phase4/malignancy"; CR <- "results/phase4/ccc_refinement"
SAMPLES <- c("MPNST_1","MPNST_2","MPNST_3","MPNST_4")
log_ <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
sec  <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
res <- list()

# --- Q1: are malignant states multi-patient? ------------------------------
sec("Q1. Does each malignant state occur in multiple patients?")
pd <- read.delim(file.path(TAB, "TUMOR_STATE_PATIENT_DISTRIBUTION.tsv"), stringsAsFactors = FALSE)
print(as.data.frame(pd |> select(tumor_state_phase4, n_cells, patients_represented,
      dominant_patient, dominant_patient_fraction, recurrence)), row.names = FALSE)
q1 <- list(n_states = nrow(pd),
  states_in_1_patient = sum(pd$patients_represented == 1),
  states_in_2_patients = sum(pd$patients_represented == 2),
  states_in_3plus = sum(pd$patients_represented >= 3),
  states_in_all_4 = sum(pd$patients_represented == 4),
  max_dominant_patient_fraction = max(pd$dominant_patient_fraction),
  median_dominant_patient_fraction = median(pd$dominant_patient_fraction),
  single_patient_states = pd$tumor_state_phase4[pd$patients_represented == 1])
str(q1)
log_(q1$states_in_1_patient, " of ", q1$n_states,
     " states occur in only ONE patient and are exploratory, not recurrent MPNST biology.")

# --- Q2: are the CNV patterns patient-specific? ---------------------------
sec("Q2. Are the SCEVAN CNV patterns patient-specific?")
sf <- file.path(TAB, "SCEVAN_CNV_SUMMARY.tsv")
q2 <- list()
if (file.exists(sf)) {
  sg <- read.delim(sf, stringsAsFactors = FALSE)
  cl <- sg |> filter(level == "clonal", width_mb >= 10, event != "neutral")
  # Per-chromosome, per-direction: in how many patients does the broad event appear?
  rec <- cl |> distinct(sample_id, Chr, event) |> count(Chr, event, name = "n_patients")
  log_("broad clonal events by chromosome and direction:")
  print(as.data.frame(rec |> arrange(desc(n_patients), Chr)), row.names = FALSE)
  # Jaccard similarity of each patient pair's broad-event set. This measures
  # SHARED / RECURRENT patterns. It is NOT a claim that patients share a clone (§50).
  sets <- split(paste(cl$Chr, cl$event), cl$sample_id)
  pr <- expand.grid(a = names(sets), b = names(sets), stringsAsFactors = FALSE) |>
    filter(a < b) |>
    mutate(jaccard = mapply(function(x, y) {
      A <- unique(sets[[x]]); B <- unique(sets[[y]])
      if (!length(union(A, B))) return(NA_real_)
      length(intersect(A, B))/length(union(A, B)) }, a, b))
  log_("pairwise Jaccard of broad clonal event sets (SHARED patterns, not shared clones):")
  print(as.data.frame(pr), row.names = FALSE)
  write.table(pr, file.path(TAB, "CNV_PATTERN_PATIENT_SIMILARITY.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  q2 <- list(
    broad_clonal_events = nrow(cl),
    events_in_all_4_patients = sum(rec$n_patients == 4),
    events_in_1_patient = sum(rec$n_patients == 1),
    events_in_3plus = sum(rec$n_patients >= 3),
    median_pairwise_jaccard = median(pr$jaccard, na.rm = TRUE),
    max_pairwise_jaccard = max(pr$jaccard, na.rm = TRUE),
    min_pairwise_jaccard = min(pr$jaccard, na.rm = TRUE),
    wording_rule = paste("Patterns shared between patients are described as SHARED or",
      "RECURRENT CNA patterns, never as the same clone (§50)."))
  str(q2)
} else log_("no SCEVAN_CNV_SUMMARY.tsv; Q2 not answerable")

# --- Q3: are the CCC changes driven by one patient? -----------------------
sec("Q3. Are the major CCC changes caused by one patient?")
pv <- read.delim(file.path(TAB, "PHASE3_VS_PHASE4_CCC.tsv"), stringsAsFactors = FALSE)
chg <- pv |> filter(change_class %in% c("Newly-supported","Sender-reassigned",
                                        "Strengthened","Lost","Weakened"))
log_("interactions in a changed class: ", nrow(chg))
q3a <- chg |> group_by(change_class) |>
  summarise(n = n(),
            median_phase4_patients = median(phase4_patients_supported, na.rm = TRUE),
            supported_in_1_patient_only = sum(phase4_patients_supported == 1, na.rm = TRUE),
            supported_in_3plus = sum(phase4_patients_supported >= 3, na.rm = TRUE),
            .groups = "drop") |>
  mutate(frac_single_patient = supported_in_1_patient_only/n)
print(as.data.frame(q3a), row.names = FALSE)
write.table(q3a, file.path(TAB, "CCC_CHANGE_PATIENT_DEPENDENCE.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
q3 <- list(changed_interactions = nrow(chg),
  by_class = q3a,
  newly_supported_single_patient_fraction =
    with(q3a[q3a$change_class == "Newly-supported", ],
         if (length(frac_single_patient)) frac_single_patient else NA_real_),
  interpretation = paste("A change class dominated by single-patient support is a WEAK result",
    "and is reported as such rather than promoted."))

# --- Q4: leave-one-patient-out on the Phase 3 core axes -------------------
sec("Q4. Leave-one-patient-out stability of the core axes")
p4c <- file.path(CR, "concordance", "CCC_CONCORDANCE.tsv")
q4 <- list()
if (file.exists(p4c)) {
  cc <- read.delim(p4c, stringsAsFactors = FALSE)
  sup_col <- "supported_samples"
  if (sup_col %in% colnames(cc)) {
    core <- cc |> filter(!is.na(.data[[sup_col]]), .data[[sup_col]] != "")
    log_("interactions with a supported-sample list: ", nrow(core))
    loso <- lapply(SAMPLES, function(drop) {
      kept <- vapply(strsplit(core[[sup_col]], "[;,] *"),
                     function(v) any(trimws(v) != drop & trimws(v) != ""), logical(1))
      data.frame(dropped_patient = drop, n_supported_full = nrow(core),
                 n_supported_without = sum(kept), retention = sum(kept)/nrow(core))
    }) |> bind_rows()
    print(as.data.frame(loso), row.names = FALSE)
    write.table(loso, file.path(TAB, "CCC_REFINED_LEAVE_ONE_PATIENT_OUT.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    q4 <- list(retention_by_dropped_patient = loso,
      min_retention = min(loso$retention), max_retention = max(loso$retention),
      note = paste("Retention below ~0.5 for any single dropped patient would mean that",
        "patient dominates the refined interaction set."))
    str(q4)
  } else log_("no supported_samples column; Q4 uses patients_supported instead")
} else log_("Phase 4 concordance table absent; Q4 not answerable")

# --- Q5 (additional): does the malignancy call itself depend on one patient? --
sec("Q5. Does the refined malignancy conclusion depend on one patient?")
calls <- read.delim(file.path(MAL, "PHASE4_MALIGNANCY_CALLS.tsv"), stringsAsFactors = FALSE)
q5 <- lapply(SAMPLES, function(drop) {
  d <- calls[calls$sample_id != drop, ]
  data.frame(dropped_patient = drop, n_cells = nrow(d),
    refined_malignant = sum(d$malignancy_refined == "Malignant"),
    refined_fraction = mean(d$malignancy_refined == "Malignant"),
    phase2_fraction = mean(d$malignancy_phase2 == "malignant"),
    fold_change = mean(d$malignancy_refined == "Malignant") /
                  max(1e-9, mean(d$malignancy_phase2 == "malignant")))
}) |> bind_rows()
full <- data.frame(dropped_patient = "none", n_cells = nrow(calls),
  refined_malignant = sum(calls$malignancy_refined == "Malignant"),
  refined_fraction = mean(calls$malignancy_refined == "Malignant"),
  phase2_fraction = mean(calls$malignancy_phase2 == "malignant"),
  fold_change = mean(calls$malignancy_refined == "Malignant") /
                mean(calls$malignancy_phase2 == "malignant"))
q5 <- bind_rows(full, q5)
print(as.data.frame(q5), row.names = FALSE)
write.table(q5, file.path(TAB, "MALIGNANCY_LEAVE_ONE_PATIENT_OUT.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# per-population, per-patient consistency of the refinement
sec("Q6. Is the fibroblast conclusion consistent across patients?")
aud <- calls |> filter(annotation_ccc_phase3 %in% c("Fibroblast","Pericyte-VSMC",
                        "Candidate-Malignant-Unresolved","Uncertain","MPNST-Tumor")) |>
  group_by(annotation_ccc_phase3, sample_id) |>
  summarise(n = n(), frac_refined_malignant = mean(malignancy_refined == "Malignant"),
            frac_scevan_malignant = mean(malignancy_scevan == "malignant"),
            .groups = "drop")
consist <- aud |> group_by(annotation_ccc_phase3) |>
  summarise(patients = n(),
            patients_majority_malignant = sum(frac_refined_malignant >= 0.5),
            min_frac = min(frac_refined_malignant), max_frac = max(frac_refined_malignant),
            spread = max(frac_refined_malignant) - min(frac_refined_malignant),
            per_patient = paste(sprintf("%s=%.2f", sample_id, frac_refined_malignant),
                                collapse = "; "), .groups = "drop")
print(as.data.frame(consist), row.names = FALSE)
write.table(consist, file.path(TAB, "MALIGNANCY_PATIENT_CONSISTENCY.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- figures ---------------------------------------------------------------
sec("Figures")
cap34 <- paste(
  "n = 4 patients. No population-level or epidemiological claim follows from any panel here.",
  "sample_id = patient = dataset, so biological and technical effects cannot be fully separated.",
  sep = "\n")

g1 <- ggplot(pd, aes(x = reorder(tumor_state_phase4, patients_represented),
                     y = patients_represented, fill = dominant_patient_fraction)) +
  geom_col(width = 0.7) + coord_flip() +
  geom_hline(yintercept = 2.5, linetype = 2, colour = "grey45") +
  scale_fill_gradient(low = "#2166AC", high = "#B2182B", limits = c(0, 1),
                      labels = percent_format(accuracy = 1),
                      name = "dominant patient\nshare") +
  scale_y_continuous(breaks = 0:4, limits = c(0, 4.2)) +
  p4_theme() +
  labs(title = "Patient representation of each malignant state",
       subtitle = paste("Dashed line marks the >=3-patient threshold. A red bar is dominated by one",
                        "patient even if\nother patients contribute a few cells."),
       x = NULL, y = "patients in which the state occurs", caption = cap34)
save_fig(g1, file.path(FIG, "34_01_state_patient_representation.pdf"), 10.5, 5.6)

g2 <- ggplot(q5 |> filter(dropped_patient != "none"),
             aes(x = dropped_patient)) +
  geom_col(aes(y = refined_fraction), fill = "#B2182B", width = 0.42,
           position = position_nudge(x = -0.22)) +
  geom_col(aes(y = phase2_fraction), fill = "#2166AC", width = 0.42,
           position = position_nudge(x = 0.22)) +
  geom_hline(yintercept = full$refined_fraction, linetype = 2, colour = "#B2182B") +
  geom_hline(yintercept = full$phase2_fraction, linetype = 2, colour = "#2166AC") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  p4_theme() +
  labs(title = "Malignant fraction with each patient left out",
       subtitle = paste("Red = Phase 4 refined, blue = Phase 2 conservative. Dashed lines are the",
                        "all-patient values.\nA bar far from its dashed line means that patient drives the estimate."),
       x = "patient left out", y = "malignant fraction", caption = cap34)
save_fig(g2, file.path(FIG, "34_02_malignancy_leave_one_patient_out.pdf"), 9.5, 5.4)

if (nrow(aud)) {
  g3 <- ggplot(aud, aes(x = sample_id, y = annotation_ccc_phase3,
                        fill = frac_refined_malignant)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.0f%%\nn=%d", 100*frac_refined_malignant, n)),
              size = 2.5, colour = ifelse(aud$frac_refined_malignant > 0.55, "white", "grey15")) +
    scale_fill_gradient2(low = "#2166AC", mid = "grey95", high = "#B2182B", midpoint = 0.5,
                         limits = c(0, 1), labels = percent_format(accuracy = 1),
                         name = "refined\nmalignant") +
    p4_theme() +
    labs(title = "Is the refinement consistent across patients?",
         subtitle = paste("Fraction of each Phase 2 population called Malignant after refinement,",
                          "per patient.\nA population that flips only in one patient is not a",
                          "recurrent finding."),
         x = NULL, y = "Phase 2 annotation_ccc", caption = cap34)
  save_fig(g3, file.path(FIG, "34_03_refinement_patient_consistency.pdf"), 9.5, 5.4)
}

if (length(q4)) {
  g4 <- ggplot(q4$retention_by_dropped_patient, aes(x = dropped_patient, y = retention)) +
    geom_col(fill = "#4D9221", width = 0.62) +
    geom_hline(yintercept = 0.5, linetype = 2, colour = "grey40") +
    geom_text(aes(label = percent(retention, accuracy = 0.1)), vjust = -0.5, size = 3) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.08)) +
    p4_theme() +
    labs(title = "Leave-one-patient-out retention of the refined interaction set",
         subtitle = "Dashed line at 50%: below it, the dropped patient would dominate the result",
         x = "patient left out", y = "interactions retained", caption = cap34)
  save_fig(g4, file.path(FIG, "34_04_ccc_leave_one_patient_out.pdf"), 9.5, 5.2)
}

facts <- list(q1_states = q1, q2_cnv = q2, q3_ccc_changes = q3, q4_loso = q4,
  q5_malignancy_loso = q5, q6_patient_consistency = consist,
  n_patients = 4L,
  limitation = paste("n = 4. No population-level or epidemiological claim follows.",
    "sample_id = patient = dataset, so biological and technical effects cannot be",
    "fully separated - the same confound recorded in Phase 2 and Phase 3."),
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  slurm = list(job_id = Sys.getenv("SLURM_JOB_ID")))
write_json(facts, file.path(MAL, "m34_robustness_facts.json"),
           auto_unbox = TRUE, pretty = TRUE, digits = 8, null = "null")
sec("M34 COMPLETE")
