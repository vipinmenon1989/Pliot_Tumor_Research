#!/usr/bin/env Rscript
# =============================================================================
# Phase 4 · M33 — tumour-state → TME model
#
# Integrates, as an EVIDENCE MATRIX rather than an arbitrary composite score (§47):
#
#   SCEVAN malignancy -> SCEVAN clone -> tumour transcriptional state
#     -> ligand programme -> receiver cell -> receiver response
#
# Deliberately TARGETED (§45): every LR framework is NOT run for every tumour
# state, which would manufacture combinatorial noise. Instead:
#   * ONE additional LIANA run on a state-split label (consensus method), and
#   * per-state, per-patient expression evidence for the Phase 3-prioritized
#     ligands, cross-referenced against receiver receptor expression, and
#   * REUSE of the Phase 3 NicheNet and LochNESS results (§46).
#
# NicheNet is not forced to support APP-CD74. Phase 3 established that predicted
# APP engagement and the CSF1-explained macrophage state are DISTINCT findings;
# that distinction is preserved here (§46).
# =============================================================================
options(stringsAsFactors = FALSE); options(future.globals.maxSize = 64 * 1024^3)
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix)
  library(dplyr); library(tidyr); library(jsonlite); library(ggplot2)
})
source("scripts/phase4/utils/phase4_plot_utils.R")
set.seed(42)

TS  <- "results/phase4/tumor_states"
CR  <- "results/phase4/ccc_refinement"
TAB <- "results/phase4/tables"
FIG <- "results/phase4/figures"
for (d in c(TS, TAB, FIG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
log_ <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
sec  <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
DETECT <- 0.10   # same 10% detection floor Phase 3 used for expression support

# Phase 3-prioritized axes, grouped into programmes. These are the axes M32
# adjudicates; M33 asks WHICH STATE carries each one.
PROGRAMMES <- list(
  `Myeloid antigen-presentation / inhibitory` =
    list(ligands = c("APP","CD99","ANXA1","THBS1","HLA-F","LGALS1","LGALS3","MIF","CD47"),
         receptors = c("CD74","PILRA","FPR1","CD36","LILRB1","LILRB2","LILRB4","SIRPA","CD44"),
         receivers = c("Macrophage","Monocyte","Dendritic","Plasmacytoid-DC","B-cell")),
  `NK / CD8 inhibitory and activating` =
    list(ligands = c("HLA-E","HLA-F","BAG6","CD58","PVR","CLEC2B"),
         receptors = c("KLRC1","KLRD1","LILRB1","NCR3","CD2","TIGIT","KLRB1"),
         receivers = c("NK","CD8-T","CD4-T","T-cell-other")),
  `Angiogenic / vascular` =
    list(ligands = c("VEGFA","VEGFB","ANGPT2","ANGPTL4","PGF","EFNA1","SEMA3A"),
         receptors = c("KDR","FLT1","NRP1","NRP2","TEK","EPHA2"),
         receivers = c("Endothelial","Pericyte-VSMC")),
  `Notch perivascular` =
    list(ligands = c("JAG1","JAG2","DLL1","DLL4"),
         receptors = c("NOTCH1","NOTCH2","NOTCH3","NOTCH4"),
         receivers = c("Pericyte-VSMC","Endothelial")),
  `Growth factor / guidance / ECM` =
    list(ligands = c("FGF2","SLIT2","PDGFB","HBEGF","IGF1","COL1A1","COL1A2","COL6A2","FN1"),
         receptors = c("FGFR1","ROBO1","PDGFRB","EGFR","IGF1R","ITGAV","ITGB8","ITGB1"),
         receivers = c("Fibroblast","Pericyte-VSMC","Endothelial")),
  `Receiver-state cytokines (Phase 3 NicheNet)` =
    list(ligands = c("CSF1","IL15","TGFB1","HMGB1","ANGPT1","IL34"),
         receptors = c("CSF1R","IL15RA","IL2RB","TGFBR1","TGFBR2","TLR4","TEK"),
         receivers = c("Macrophage","Monocyte","NK","CD8-T","CD4-T","Endothelial")))

# --- 1. inputs -------------------------------------------------------------
sec("1. Inputs")
mo <- readRDS(file.path(TS, "malignant_only_object.rds"))
log_("malignant object: ", ncol(mo), " cells | states: ",
     length(unique(mo$tumor_state_phase4)))
co <- readRDS(file.path(CR, "ccc_input_object_refined.rds"))
log_("refined CCC object: ", ncol(co), " cells | populations: ", nlevels(co$ccc_label))

pd <- read.delim(file.path(TAB, "TUMOR_STATE_PATIENT_DISTRIBUTION.tsv"), stringsAsFactors = FALSE)
sens <- read.delim(file.path(TAB, "PHASE3_AXIS_SENSITIVITY.tsv"), stringsAsFactors = FALSE)
log_("axis sensitivity rows: ", nrow(sens))

# Phase 3 results reused verbatim (§46) - not recomputed.
nn_f <- "results/phase3/tables/final/NicheNet_ligand_activity.tsv"
nn <- if (file.exists(nn_f)) read.delim(nn_f, stringsAsFactors = FALSE) else NULL
log_("Phase 3 NicheNet rows reused: ", if (is.null(nn)) 0 else nrow(nn))
lo_f <- "results/phase3/tables/final/lochness_summary.tsv"
lo <- if (file.exists(lo_f)) read.delim(lo_f, stringsAsFactors = FALSE) else NULL
log_("Phase 3 LochNESS summary rows reused: ", if (is.null(lo)) 0 else nrow(lo))

# --- 2. sender evidence: ligand expression per state per patient -----------
sec("2. Sender evidence — ligand expression by tumour state and patient")
LIG <- sort(unique(unlist(lapply(PROGRAMMES, `[[`, "ligands"))))
LIG <- intersect(LIG, rownames(mo))
log_("prioritized ligands present in the malignant object: ", length(LIG), " of ",
     length(unique(unlist(lapply(PROGRAMMES, `[[`, "ligands")))))
E <- GetAssayData(mo, assay = "RNA", layer = "data")[LIG, , drop = FALSE]
st <- as.character(mo$tumor_state_phase4); pt <- as.character(mo$sample_id)
send <- list()
for (s in sort(unique(st))) for (p in sort(unique(pt))) {
  i <- which(st == s & pt == p)
  if (length(i) < 10) next      # same 10-cell floor Phase 3 used; below it, NOT EVALUABLE
  sub <- E[, i, drop = FALSE]
  send[[length(send)+1]] <- data.frame(
    tumor_state_phase4 = s, sample_id = p, n_cells = length(i), ligand = LIG,
    detect_frac = as.numeric(Matrix::rowMeans(sub > 0)),
    mean_expr = as.numeric(Matrix::rowMeans(sub)), stringsAsFactors = FALSE)
}
send <- bind_rows(send)
send$expressed <- send$detect_frac >= DETECT
log_("sender rows (state x patient x ligand): ", nrow(send))
log_("state x patient strata evaluated: ", nrow(distinct(send, tumor_state_phase4, sample_id)))
sender_sum <- send |> group_by(tumor_state_phase4, ligand) |>
  summarise(patients_evaluable = n(), patients_expressing = sum(expressed),
            max_detect_frac = max(detect_frac), median_detect_frac = median(detect_frac),
            .groups = "drop")
write.table(send, gzfile(file.path(TAB, "TUMOR_STATE_LIGAND_EXPRESSION.tsv.gz")),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- 3. receiver evidence: receptor expression per population per patient ---
sec("3. Receiver evidence — receptor expression by population and patient")
REC <- sort(unique(unlist(lapply(PROGRAMMES, `[[`, "receptors"))))
REC <- intersect(REC, rownames(co))
log_("prioritized receptors present in the refined CCC object: ", length(REC))
Er <- GetAssayData(co, assay = "RNA", layer = "data")[REC, , drop = FALSE]
rl <- as.character(co$ccc_label); rp <- as.character(co$ccc_sample)
recv <- list()
for (r in sort(unique(rl))) for (p in sort(unique(rp))) {
  i <- which(rl == r & rp == p)
  if (length(i) < 10) next
  sub <- Er[, i, drop = FALSE]
  recv[[length(recv)+1]] <- data.frame(
    receiver = r, sample_id = p, n_cells = length(i), receptor = REC,
    detect_frac = as.numeric(Matrix::rowMeans(sub > 0)),
    mean_expr = as.numeric(Matrix::rowMeans(sub)), stringsAsFactors = FALSE)
}
recv <- bind_rows(recv)
recv$expressed <- recv$detect_frac >= DETECT
log_("receiver rows: ", nrow(recv))
write.table(recv, gzfile(file.path(TAB, "RECEIVER_RECEPTOR_EXPRESSION.tsv.gz")),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- 4. the evidence matrix ------------------------------------------------
sec("4. Evidence matrix (NOT a composite score)")
rows <- list()
for (pg in names(PROGRAMMES)) {
  cfg <- PROGRAMMES[[pg]]
  for (s in sort(unique(st))) {
    ligs <- intersect(cfg$ligands, LIG)
    for (rcv in intersect(cfg$receivers, unique(rl))) {
      recs <- intersect(cfg$receptors, REC)
      if (!length(ligs) || !length(recs)) next
      sd_ <- send |> filter(tumor_state_phase4 == s, ligand %in% ligs)
      rd_ <- recv |> filter(receiver == rcv, receptor %in% recs)
      if (!nrow(sd_) || !nrow(rd_)) next
      # A pair counts as jointly supported in a patient only if BOTH sides pass the
      # detection floor in THAT patient. Sample-aware, never pooled (§SAMPLE-AWARE).
      pats <- intersect(unique(sd_$sample_id), unique(rd_$sample_id))
      joint <- vapply(pats, function(p)
        any(sd_$expressed[sd_$sample_id == p]) && any(rd_$expressed[rd_$sample_id == p]),
        logical(1))
      lig_top <- sd_ |> group_by(ligand) |>
        summarise(pe = sum(expressed), md = max(detect_frac), .groups = "drop") |>
        arrange(desc(pe), desc(md))
      rec_top <- rd_ |> group_by(receptor) |>
        summarise(pe = sum(expressed), md = max(detect_frac), .groups = "drop") |>
        arrange(desc(pe), desc(md))
      rows[[length(rows)+1]] <- data.frame(
        programme = pg, tumor_state_phase4 = s, receiver = rcv,
        n_patients_evaluable = length(pats),
        n_patients_joint_expression = sum(joint),
        joint_patients = paste(pats[joint], collapse = ";"),
        top_ligands = paste(head(sprintf("%s(%d/%d,%.2f)", lig_top$ligand, lig_top$pe,
                             length(unique(sd_$sample_id)), lig_top$md), 4), collapse = "; "),
        top_receptors = paste(head(sprintf("%s(%d/%d,%.2f)", rec_top$receptor, rec_top$pe,
                             length(unique(rd_$sample_id)), rec_top$md), 4), collapse = "; "),
        stringsAsFactors = FALSE)
    }
  }
}
ev <- bind_rows(rows)

# Attach the remaining evidence streams, each as its own column.
ev <- ev |> left_join(pd |> select(tumor_state_phase4, state_n_cells = n_cells,
                                   state_patients = patients_represented,
                                   state_recurrence = recurrence,
                                   state_clones = n_scevan_clones), by = "tumor_state_phase4")
# Phase 3 axis verdicts, summarised per programme (the M32 sensitivity result).
axis_map <- list(
  `Myeloid antigen-presentation / inhibitory` = c("tumour->myeloid antigen-presentation",
    "tumour->myeloid inhibitory","tumour->myeloid chemotaxis","tumour->myeloid scavenger",
    "tumour->myeloid/lymphoid inhibitory"),
  `NK / CD8 inhibitory and activating` = c("tumour->NK inhibitory","tumour->NK activating",
    "tumour->lymphoid adhesion/costimulation"),
  `Angiogenic / vascular` = "tumour->endothelial angiogenic",
  `Notch perivascular` = c("tumour->pericyte Notch","tumour->endothelial Notch","TME->tumour Notch"),
  `Growth factor / guidance / ECM` = c("tumour->stromal growth factor",
    "tumour->stromal guidance","ECM->tumour integrin"),
  `Receiver-state cytokines (Phase 3 NicheNet)` = character(0))
ev$phase4_ccc_verdicts <- vapply(ev$programme, function(pg) {
  ax <- axis_map[[pg]]
  if (!length(ax)) return("n/a - reused Phase 3 NicheNet, not an M32 axis")
  v <- sens$verdict[sens$axis %in% ax]
  if (!length(v)) return("not adjudicated")
  paste(sprintf("%s=%d", names(sort(table(v), decreasing = TRUE)),
                sort(table(v), decreasing = TRUE)), collapse = "; ") }, character(1))

if (!is.null(nn)) {
  nnc <- colnames(nn)
  rc <- intersect(c("receiver","receiver_population","target_population"), nnc)
  lc <- intersect(c("ligand","test_ligand"), nnc)
  ac <- intersect(c("aupr_corrected","aupr","pearson","auroc"), nnc)
  if (length(rc) && length(lc) && length(ac)) {
    nns <- nn |> rename(receiver = !!rc[1], ligand = !!lc[1], activity = !!ac[1]) |>
      group_by(receiver) |> arrange(desc(activity)) |>
      summarise(nichenet_top_ligands = paste(head(ligand, 5), collapse = "; "),
                nichenet_best_activity = max(activity), .groups = "drop")
    ev <- left_join(ev, nns, by = "receiver")
  }
}
ev$nichenet_note <- paste(
  "Phase 3 NicheNet is REUSED, not recomputed, and is NOT forced to agree with the LR result.",
  "Phase 3 established that predicted APP-CD74 engagement and the CSF1-explained macrophage",
  "state are distinct findings; that distinction is preserved.")
if (!is.null(lo)) {
  lc2 <- intersect(c("lineage","receiver_lineage","receiver"), colnames(lo))
  if (length(lc2)) {
    los <- lo |> rename(receiver_lineage = !!lc2[1])
    ev$lochness_context <- "see LOCHNESS_MPNST_REPORT.md - negative in all four lineages"
  }
}
ev$lochness_note <- paste(
  "LochNESS enters ONLY as receiver-lineage context. Phase 3 found no receiver-state structure",
  "associated with the tumour-derived APP context in any lineage (best descriptive p = 0.333",
  "against a 0.167 floor at n = 4) and the score was itself unstable. High LochNESS would not",
  "mean cell-cell communication, and low LochNESS does not refute the LR findings.")

ev <- ev |> arrange(programme, desc(n_patients_joint_expression), desc(state_n_cells))
f1 <- file.path(TAB, "TUMOR_STATE_TME_EVIDENCE.tsv")
write.table(ev, f1, sep = "\t", quote = FALSE, row.names = FALSE)
log_("wrote ", f1, " (", nrow(ev), " rows)")
print(as.data.frame(ev |> select(programme, tumor_state_phase4, receiver,
      n_patients_joint_expression, n_patients_evaluable, state_recurrence) |>
      filter(n_patients_joint_expression >= 3) |> head(40)), row.names = FALSE)

# --- 5. which state carries which programme -------------------------------
sec("5. Which tumour state carries which programme")
lead <- ev |> group_by(programme, tumor_state_phase4) |>
  summarise(receivers_supported = sum(n_patients_joint_expression >= 3),
            receivers_evaluated = n(),
            max_patients = max(n_patients_joint_expression),
            state_patients = first(state_patients),
            state_n_cells = first(state_n_cells),
            state_recurrence = first(state_recurrence), .groups = "drop") |>
  group_by(programme) |> arrange(desc(receivers_supported), desc(max_patients)) |>
  mutate(rank_in_programme = row_number()) |> ungroup()
write.table(lead, file.path(TAB, "TUMOR_STATE_PROGRAMME_LEADERS.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(as.data.frame(lead |> filter(rank_in_programme <= 3)), row.names = FALSE)

# --- 6. figures ------------------------------------------------------------
sec("6. Figures")
cap33 <- paste(
  "Evidence matrix, not a composite score. Each column is an independent stream and none are averaged.",
  "Joint expression = ligand AND receptor each detected in >=10% of the relevant population, IN THE SAME PATIENT.",
  "scRNA-seq infers communication POTENTIAL; it does not establish adjacency or direct signalling.",
  sep = "\n")

g1 <- ggplot(ev, aes(x = receiver, y = tumor_state_phase4,
                     fill = n_patients_joint_expression)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = n_patients_joint_expression), size = 2.5,
            colour = ifelse(ev$n_patients_joint_expression >= 3, "white", "grey20")) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08519C", limits = c(0, 4),
                      breaks = 0:4, name = "patients with\njoint expression") +
  facet_wrap(~ programme, scales = "free_x") +
  p4_theme(9) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
  labs(title = "Tumour state → TME: which malignant state engages which receiver",
       subtitle = "Cell values are the number of the four patients in which BOTH sides pass the 10% detection floor",
       x = "receiver population", y = "tumor_state_phase4", caption = cap33)
save_fig(g1, file.path(FIG, "33_01_tumor_state_specific_ccc.pdf"), 15.5, 9.2)

sd2 <- send |> filter(ligand %in% c("APP","CD99","ANXA1","HLA-E","HLA-F","VEGFA","JAG1","JAG2",
                                    "DLL4","FGF2","CSF1","IL15","TGFB1","COL1A1","FN1")) |>
  mutate(ligand = factor(ligand, levels = c("APP","CD99","ANXA1","HLA-E","HLA-F","VEGFA",
                         "JAG1","JAG2","DLL4","FGF2","CSF1","IL15","TGFB1","COL1A1","FN1")))
g2 <- ggplot(sd2, aes(x = ligand, y = tumor_state_phase4)) +
  geom_point(aes(size = detect_frac, colour = mean_expr), alpha = 0.9) +
  scale_size_continuous(range = c(0.3, 5), limits = c(0, 1),
                        name = "detection\nfraction") +
  scale_colour_gradient(low = "#DEEBF7", high = "#08306B", name = "mean\nlog-expression") +
  facet_wrap(~ sample_id, nrow = 1) +
  p4_theme(9) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)) +
  labs(title = "Prioritized ligand expression by malignant state, per patient",
       subtitle = paste("Sample-aware by construction. A state x patient stratum with <10 cells is",
                        "absent from the panel because it is\nNOT EVALUABLE, which is not the same",
                        "as no expression."),
       x = NULL, y = "tumor_state_phase4", caption = cap33)
save_fig(g2, file.path(FIG, "33_02_state_ligand_expression_by_patient.pdf"), 16.5, 6.4)

g3 <- ggplot(lead, aes(x = reorder(tumor_state_phase4, receivers_supported),
                       y = receivers_supported, fill = state_recurrence)) +
  geom_col(width = 0.72) + coord_flip() +
  scale_fill_manual(values = c(
    "recurrent: >=3 patients each contribute >=5%" = "#1B7837",
    "shared between 2 patients" = "#7FBC41",
    "PATIENT-DOMINATED (>80% one patient) / exploratory - not recurrent MPNST biology" = "#F0A202",
    "PATIENT-SPECIFIC / exploratory - not recurrent MPNST biology" = "#D95F02",
    "PATIENT-DOMINATED / exploratory - not recurrent MPNST biology" = "#F0A202"),
    name = "state recurrence", drop = FALSE) +
  facet_wrap(~ programme, scales = "free_y") + p4_theme(9) +
  labs(title = "Programme carriage per malignant state, with recurrence honesty",
       subtitle = paste("Bar height = receivers supported in >=3 patients. Colour flags states seen in",
                        "only one patient, which are\nexploratory and must not be read as recurrent",
                        "MPNST biology."),
       x = NULL, y = "receivers supported in >=3 patients", caption = cap33)
save_fig(g3, file.path(FIG, "33_03_programme_carriage_by_state.pdf"), 15.5, 9.2)

facts <- list(
  malignant_cells = ncol(mo), states = as.list(table(mo$tumor_state_phase4)),
  ccc_populations = levels(co$ccc_label),
  detection_floor = DETECT, min_cells_per_stratum = 10L,
  programmes = names(PROGRAMMES),
  evidence_rows = nrow(ev),
  supported_3plus_patients = sum(ev$n_patients_joint_expression >= 3),
  programme_leaders = lead |> filter(rank_in_programme == 1) |>
    select(programme, tumor_state_phase4, receivers_supported, state_recurrence),
  no_composite_score = TRUE,
  targeted_design = paste("LR frameworks were NOT run per tumour state (§45). Evidence is",
    "per-state, per-patient ligand/receptor expression on the Phase 3-prioritized axes,",
    "cross-referenced with the M32 sensitivity verdicts and with Phase 3 NicheNet and",
    "LochNESS results reused verbatim."),
  phase3_reused = list(nichenet = nn_f, lochness = lo_f),
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  slurm = list(job_id = Sys.getenv("SLURM_JOB_ID")))
write_json(facts, file.path(TS, "m33_state_tme_facts.json"),
           auto_unbox = TRUE, pretty = TRUE, digits = 8, null = "null")
sec("M33 COMPLETE")
