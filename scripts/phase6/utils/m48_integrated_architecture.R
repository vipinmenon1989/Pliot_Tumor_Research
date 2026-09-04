#!/usr/bin/env Rscript
# =============================================================================
# Phase 6 - M48 - integrated malignant architecture and per-patient portraits.
#
# Three candidate architectures are evaluated against measured quantities and
# whichever the data support is reported. The conclusion is NOT chosen first:
#
#   Model A  clone-constrained  - a clone's cells concentrate in one program
#   Model B  program-diverse    - clones span several programs
#   Model C  mixed              - some programs clone-associated, others not
#
# With n = 4 patients (3 with usable clone structure) the honest deliverable is
# a deep per-patient portrait, not a population-level statistical claim.
# =============================================================================
source("scripts/phase6/utils/phase6_common.R")
set.seed(42)
facts <- list(milestone = "M48", generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))

md   <- readRDS(file.path(P6_CLONE, "phase6_cell_metadata.rds"))
assoc<- read.delim(file.path(P6_TAB, "CLONE_PROGRAM_ASSOCIATION.tsv"))
div  <- read.delim(file.path(P6_TAB, "CLONE_PROGRAM_DIVERSITY.tsv"))
cpa  <- read.delim(file.path(P6_TAB, "CLONE_PROGRAM_ACTIVITY.tsv"))
tf   <- read.delim(file.path(P6_TAB, "PROGRAM_TF_ACTIVITY.tsv"))
pw   <- read.delim(file.path(P6_TAB, "PROGRAM_PATHWAY_ACTIVITY.tsv"))
lab  <- read.delim("results/phase5/tables/final/MALIGNANT_PROGRAMS.tsv")
LABEL <- setNames(lab$program_label, lab$program_id)
cx   <- if (file.exists(file.path(P6_TAB, "BROAD_CNA_EXPRESSION_EFFECTS.tsv")))
  read.delim(file.path(P6_TAB, "BROAD_CNA_EXPRESSION_EFFECTS.tsv")) else NULL
btw  <- read.delim(file.path(P6_PLAS, "M44_BETWEEN_VS_WITHIN_CLONE.tsv"))

p6_sec("1. Evidence for each candidate architecture")
ev <- div[div$n_cells >= CLONE_SIZE_PRIMARY, ]
n_ev <- nrow(ev)
n_multi <- sum(ev$hard_effective_n >= 2, na.rm = TRUE)
n_conc <- sum(ev$dominant_program_share_hard >= 0.80, na.rm = TRUE)
med_eta2 <- median(assoc$eta2, na.rm = TRUE)
max_eta2 <- max(assoc$eta2, na.rm = TRUE)
n_assoc <- sum(assoc$clone_associated, na.rm = TRUE)
n_pairs <- nrow(assoc)
p6_msg("evaluable clones (>= %d cells): %d", CLONE_SIZE_PRIMARY, n_ev)
p6_msg("  spanning >= 2 programs (hard effective n >= 2):  %d (%.0f%%)",
       n_multi, 100 * n_multi / max(1, n_ev))
p6_msg("  concentrated (>= 80%% of cells in one program):   %d (%.0f%%)",
       n_conc, 100 * n_conc / max(1, n_ev))
p6_msg("clone-program eta^2: median %.3f, max %.3f", med_eta2, max_eta2)
p6_msg("  clone-associated program x patient pairs: %d of %d (%.0f%%)",
       n_assoc, n_pairs, 100 * n_assoc / max(1, n_pairs))
p6_msg("  => median eta^2 %.3f means about %.0f%% of program variance sits WITHIN clones",
       med_eta2, 100 * (1 - med_eta2))

# The decisive quantity, and the one that separates "clone-constrained" from
# everything else: how large is the divergence BETWEEN a patient's clones
# compared with the spread of cells WITHIN one clone?
ratio <- median(btw$ratio, na.rm = TRUE)
p6_msg("between-clone divergence / within-clone dispersion: %s (median %.4f)",
       paste(sprintf("%s=%.4f", btw$sample_id, btw$ratio), collapse = "  "), ratio)
p6_msg("  a ratio far below 1 means a patient's clones are transcriptionally near-interchangeable")

modelA <- n_conc / max(1, n_ev) >= 0.70 && med_eta2 >= 0.30 && ratio >= 0.50
modelB <- n_multi / max(1, n_ev) >= 0.70 && n_assoc / max(1, n_pairs) <= 0.25
selected <- if (modelA) "A_clone_constrained" else
            if (modelB) "B_within_clone_program_diversity" else "C_mixed"
p6_msg("selected architecture: %s", selected)
verdict <- switch(selected,
  A_clone_constrained = "Malignant transcriptional phenotype is strongly constrained by CNA-defined clone structure.",
  B_within_clone_program_diversity = paste(
    "Individual CNA-defined clones occupy multiple transcriptional programs, indicating substantial",
    "within-clone transcriptional-program diversity. This is CONSISTENT WITH phenotypic plasticity;",
    "cross-sectional data cannot demonstrate an actual transition."),
  C_mixed = paste(
    "MPNST shows a mixed architecture: a minority of transcriptional programs associate detectably with",
    "specific CNA-defined clones while the great majority of program variance sits WITHIN clones.",
    sprintf("Median eta^2 is %.3f, so roughly %.0f%% of each program's variance is within-clone,", med_eta2, 100 * (1 - med_eta2)),
    sprintf("and between-clone divergence is only %.1f%% of within-clone dispersion (median across patients),", 100 * ratio),
    "meaning a patient's CNA-defined clones are transcriptionally near-interchangeable.",
    "How much program diversity a clone contains is itself patient-specific, so this is a mixed",
    "architecture rather than a single rule that holds across the cohort."))
p6_msg("%s", verdict)
facts$architecture <- list(
  evaluated = c("A_clone_constrained", "B_within_clone_program_diversity", "C_mixed"),
  selected = selected, verdict = verdict,
  criteria = list(
    A = "clones concentrated (>= 80% of cells in one program) in >= 70% of evaluable clones AND median eta^2 >= 0.30 AND between-clone divergence >= 50% of within-clone dispersion",
    B = "clones spanning >= 2 programs in >= 70% of evaluable clones AND <= 25% of program x patient pairs clone-associated",
    C = "otherwise"),
  evidence = list(n_evaluable_clones = n_ev, n_multi_program = n_multi,
                  n_concentrated = n_conc, median_eta2 = med_eta2,
                  max_eta2 = max_eta2, n_clone_associated_pairs = n_assoc,
                  n_pairs = n_pairs,
                  between_within_ratio_median = ratio,
                  between_within_by_patient = btw))

p6_sec("2. Per-patient tumour portraits")
port <- list()
for (s in SAMPLES) {
  reliable <- s %in% CLONE_RELIABLE
  m <- md[md$sample_id == s & md$malignancy_refined == "Malignant", ]
  progcomp <- sort(table(m$dominant_malignant_program), decreasing = TRUE)
  pc <- paste(sprintf("%s(%s)=%d", names(progcomp), LABEL[names(progcomp)],
                      as.integer(progcomp)), collapse = "; ")
  if (reliable) {
    d <- div[div$sample_id == s, ]
    de <- d[d$n_cells >= CLONE_SIZE_PRIMARY, ]
    a <- assoc[assoc$sample_id == s, ]
    tfs <- tf[tf$recurrent_direction & !is.na(tf[[paste0("rho_", s)]]) &
                abs(tf[[paste0("rho_", s)]]) >= 0.20, ]
    tfs <- tfs[order(-abs(tfs[[paste0("rho_", s)]])), ]
    pws <- pw[pw$layer == "PROGENy" & !is.na(pw[[paste0("rho_", s)]]), ]
    pws <- pws[order(-abs(pws[[paste0("rho_", s)]])), ]
    cxs <- if (!is.null(cx)) cx[cx$sample_id == s & cx$evaluable, ] else NULL
    port[[s]] <- data.frame(
      sample_id = s, clone_structure_reliable = TRUE,
      n_malignant = nrow(m),
      n_clones_total = length(unique(na.omit(m$tumor_clone_phase4))),
      n_clones_evaluable = nrow(de),
      malignant_program_composition = pc,
      dominant_program = names(progcomp)[1],
      dominant_program_label = unname(LABEL[names(progcomp)[1]]),
      median_clone_effective_programs = median(de$hard_effective_n, na.rm = TRUE),
      median_clone_dominant_share = median(de$dominant_program_share_hard, na.rm = TRUE),
      n_clones_multi_program = sum(de$hard_effective_n >= 2, na.rm = TRUE),
      median_eta2 = median(a$eta2, na.rm = TRUE),
      n_clone_associated_programs = sum(a$clone_associated, na.rm = TRUE),
      broad_cna_events_evaluable = if (!is.null(cxs)) nrow(cxs) else NA_integer_,
      broad_cna_direction_match = if (!is.null(cxs) && nrow(cxs))
        sprintf("%d/%d", sum(cxs$direction_matches_event), nrow(cxs)) else NA_character_,
      top_tf_activities = paste(head(tfs$tf, 6), collapse = ","),
      top_pathways = paste(head(pws$pathway, 5), collapse = ","),
      note = "clone-based conclusions are supported for this patient")
  } else {
    port[[s]] <- data.frame(
      sample_id = s, clone_structure_reliable = FALSE,
      n_malignant = nrow(m),
      n_clones_total = length(unique(na.omit(m$tumor_clone_phase4))),
      n_clones_evaluable = NA_integer_,
      malignant_program_composition = pc,
      dominant_program = names(progcomp)[1],
      dominant_program_label = unname(LABEL[names(progcomp)[1]]),
      median_clone_effective_programs = NA_real_,
      median_clone_dominant_share = NA_real_, n_clones_multi_program = NA_integer_,
      median_eta2 = NA_real_, n_clone_associated_programs = NA_integer_,
      broad_cna_events_evaluable = NA_integer_, broad_cna_direction_match = NA_character_,
      top_tf_activities = NA_character_, top_pathways = NA_character_,
      note = paste("Phase 5 transcriptional programs are described for this patient,",
                   "but its SCEVAN clone architecture is UNRELIABLE and excluded.",
                   CLONE_EXCLUSION_REASON))
  }
}
port <- do.call(rbind, port)
p6_tsv(port, file.path(P6_TAB, "PATIENT_TUMOR_ARCHITECTURE.tsv"))
print(as.data.frame(port |> select(sample_id, clone_structure_reliable, n_malignant,
  n_clones_total, n_clones_evaluable, dominant_program_label,
  median_clone_effective_programs, median_eta2)), row.names = FALSE, digits = 3)

p6_sec("3. Integrated evidence table")
integ <- lab |> select(program = program_id, program_label, recurrence_status,
                       dominant_patient, dominant_patient_fraction) |>
  left_join(assoc |> group_by(program) |>
              summarise(median_eta2 = median(eta2, na.rm = TRUE),
                        max_eta2 = max(eta2, na.rm = TRUE),
                        n_patients_clone_associated = sum(clone_associated, na.rm = TRUE),
                        .groups = "drop"), by = "program") |>
  left_join(tf |> filter(recurrent_direction) |> group_by(program) |>
              summarise(n_recurrent_tfs = n(),
                        top_tfs = paste(head(tf[order(-abs(pooled_rho))], 5), collapse = ","),
                        .groups = "drop"), by = "program") |>
  left_join(pw |> filter(recurrent_direction, layer == "PROGENy") |>
              group_by(program) |>
              summarise(n_recurrent_progeny = n(),
                        top_progeny = paste(head(pathway[order(-abs(pooled_rho))], 4),
                                            collapse = ","), .groups = "drop"),
            by = "program") |>
  left_join(div |> filter(n_cells >= CLONE_SIZE_PRIMARY) |>
              group_by(dominant_program) |>
              summarise(n_clones_dominated = n(), .groups = "drop") |>
              rename(program = dominant_program), by = "program") |>
  mutate(across(c(n_recurrent_tfs, n_recurrent_progeny, n_clones_dominated),
                ~tidyr::replace_na(.x, 0L)),
         layers_note = "program, TF activity and pathway activity are independent layers - not merged")
p6_tsv(integ, file.path(P6_TAB, "PHASE6_INTEGRATED_EVIDENCE.tsv"))
print(as.data.frame(integ |> select(program, program_label, recurrence_status,
  median_eta2, n_patients_clone_associated, n_recurrent_tfs, n_recurrent_progeny,
  n_clones_dominated)), row.names = FALSE, digits = 3)

facts$patient_portraits <- port
facts$integrated <- integ
p6_json(facts, file.path(P6_VAL, "m48_integrated_facts.json"))
p6_sec("M48 complete")
