#!/usr/bin/env Rscript
# =============================================================================
# Phase 5 - M39 - malignant ECM-like MPNST cells versus true fibroblasts.
#
#   Group A  historical Fibroblast + malignancy_refined == Malignant      4,036
#   Group B  historical Fibroblast + malignancy_refined == Non-malignant    908
#   Group C  historical Fibroblast + malignancy_refined == Ambiguous        120
#            (sensitivity / projection ONLY - never in the primary comparison)
#
# "historical" is the field M36 VERIFIED reproduces the frozen fibroblast split.
# annotation_ccc_refined is barred from this role: it already encodes the
# Phase 4 conclusion, so using it would make the comparison circular.
#
# Statistics: patient is the replicate. Per-patient pseudobulk fold change and
# per-patient cell-level effect size (AUC) are primary; a patient with too few
# cells on one side is printed NOT EVALUABLE rather than pooled away. No
# classifier is trained and no random cell-level train/test split is used.
# =============================================================================
suppressPackageStartupMessages({ library(edgeR); library(presto); library(Matrix) })
source("scripts/phase5/utils/phase5_common.R")
source("scripts/phase5/utils/program_projection.R")
set.seed(42)
args  <- commandArgs(trailingOnly = TRUE)
KSTAR <- as.integer(args[1]); DT <- if (length(args) > 1) args[2] else "0_1"
facts <- list(milestone = "M39", K = KSTAR,
              generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))

MIN_CELLS_EVAL <- 30L   # declared: a patient needs >= 30 cells on BOTH sides

p5_sec("1. Groups, from the VERIFIED historical annotation field")
md <- readRDS(file.path(P5_PROG, "phase5_cell_metadata.rds"))
cnts <- readRDS(file.path(P5_PROG, "phase5_working_counts.rds"))
HIST <- "annotation_ccc_phase3"
md$hist <- md[[HIST]]
p5_msg("historical field: %s (verified in M36)", HIST)

fib <- md[md$hist == "Fibroblast", ]
stopifnot(nrow(fib) == P4_FROZEN$fibroblast[["total"]])
grp <- setNames(rep(NA_character_, nrow(fib)), fib$cell_id)
grp[fib$malignancy_refined == "Malignant"]     <- "A_malignant_ECM"
grp[fib$malignancy_refined == "Non-malignant"] <- "B_true_fibroblast"
grp[fib$malignancy_refined == "Ambiguous"]     <- "C_ambiguous"
p5_msg("A malignant ECM-like %d | B true fibroblast %d | C ambiguous %d",
       sum(grp == "A_malignant_ECM"), sum(grp == "B_true_fibroblast"),
       sum(grp == "C_ambiguous"))
stopifnot(sum(grp == "A_malignant_ECM") == P4_FROZEN$fibroblast[["Malignant"]],
          sum(grp == "B_true_fibroblast") == P4_FROZEN$fibroblast[["Non-malignant"]],
          sum(grp == "C_ambiguous") == P4_FROZEN$fibroblast[["Ambiguous"]])
fib$group <- unname(grp[fib$cell_id])

p5_sec("2. Evaluability, per patient - checked BEFORE any test is run")
ev <- fib |> filter(group != "C_ambiguous") |>
  count(sample_id, group) |>
  pivot_wider(names_from = group, values_from = n, values_fill = 0) |>
  mutate(evaluable = A_malignant_ECM >= MIN_CELLS_EVAL &
           B_true_fibroblast >= MIN_CELLS_EVAL,
         reason = case_when(
           A_malignant_ECM < MIN_CELLS_EVAL & B_true_fibroblast < MIN_CELLS_EVAL ~
             "both groups below the declared minimum",
           A_malignant_ECM < MIN_CELLS_EVAL ~
             "too few malignant ECM-like cells (MPNST_3 promotions were disabled in Phase 4)",
           B_true_fibroblast < MIN_CELLS_EVAL ~
             "too few non-malignant fibroblasts to compare against",
           TRUE ~ "evaluable"))
print(as.data.frame(ev), row.names = FALSE)
p5_tsv(ev, file.path(P5_ECM, "M39_EVALUABILITY_BY_PATIENT.tsv"))
EVAL <- as.character(ev$sample_id[ev$evaluable])
p5_msg("evaluable patients: %s (of 4)", paste(EVAL, collapse = ", "))
p5_msg("a four-patient paired test is NOT constructed - the groups do not coexist in 2 of 4 patients")
facts$evaluability <- ev; facts$evaluable_patients <- EVAL
facts$min_cells_declared <- MIN_CELLS_EVAL

p5_sec("3. Pseudobulk fold change per patient (patient is the replicate)")
ab <- fib[fib$group %in% c("A_malignant_ECM", "B_true_fibroblast"), ]
key <- paste(ab$sample_id, ab$group, sep = "|")
pb <- t(rowsum(t(as.matrix(cnts[, ab$cell_id, drop = FALSE])), group = key))
p5_msg("pseudobulk matrix: %d genes x %d sample-group columns", nrow(pb), ncol(pb))
keep <- rowSums(pb) >= 10
pb <- pb[keep, , drop = FALSE]
cpm <- edgeR::cpm(pb, log = TRUE, prior.count = 3)
p5_msg("genes retained for pseudobulk: %d", nrow(pb))

lfc <- do.call(cbind, lapply(EVAL, function(s) {
  a <- paste0(s, "|A_malignant_ECM"); b <- paste0(s, "|B_true_fibroblast")
  cpm[, a] - cpm[, b]
}))
colnames(lfc) <- EVAL
p5_msg("per-patient pseudobulk log2FC (A vs B) computed for: %s",
       paste(EVAL, collapse = ", "))
if (length(EVAL) >= 2)
  p5_msg("cross-patient correlation of log2FC: %.3f",
         cor(lfc[, 1], lfc[, 2], method = "spearman"))

p5_sec("4. Cell-level effect size per patient (AUC; not a pooled p-value)")
auc <- lapply(EVAL, function(s) {
  ids <- ab$cell_id[ab$sample_id == s]
  m <- cnts[, ids, drop = FALSE]
  ln <- log1p(t(t(m) / pmax(Matrix::colSums(m), 1) * 1e4))
  w <- presto::wilcoxauc(as.matrix(ln), ab$group[match(ids, ab$cell_id)])
  w <- w[w$group == "A_malignant_ECM", c("feature", "auc", "logFC", "pct_in", "pct_out")]
  names(w) <- c("gene", paste0(c("auc_", "logFC_", "pct_A_", "pct_B_"), s))
  w
})
auc <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), auc)
p5_msg("AUC table: %d genes x %d patients", nrow(auc), length(EVAL))

p5_sec("5. Transparent malignant-ECM signature (declared criteria)")
# Declared before inspection: a signature gene must move in the SAME direction
# in EVERY evaluable patient, by >= 1 log2 unit in pseudobulk, with a cell-level
# AUC >= 0.65 (or <= 0.35) in every evaluable patient, and be detected in >= 10%
# of cells in the higher group. No classifier, no train/test split.
LFC_MIN <- 1.0; AUC_MIN <- 0.65; PCT_MIN <- 0.10
sig <- data.frame(gene = rownames(lfc), lfc, check.names = FALSE)
sig <- merge(sig, auc, by = "gene")
aucc <- as.matrix(sig[, paste0("auc_", EVAL), drop = FALSE])
lfcc <- as.matrix(sig[, EVAL, drop = FALSE])
pctA <- as.matrix(sig[, paste0("pct_A_", EVAL), drop = FALSE])
pctB <- as.matrix(sig[, paste0("pct_B_", EVAL), drop = FALSE])
up <- rowSums(lfcc >= LFC_MIN, na.rm = TRUE) == length(EVAL) &
      rowSums(aucc >= AUC_MIN, na.rm = TRUE) == length(EVAL) &
      rowSums(pctA >= PCT_MIN * 100, na.rm = TRUE) == length(EVAL)
dn <- rowSums(lfcc <= -LFC_MIN, na.rm = TRUE) == length(EVAL) &
      rowSums(aucc <= 1 - AUC_MIN, na.rm = TRUE) == length(EVAL) &
      rowSums(pctB >= PCT_MIN * 100, na.rm = TRUE) == length(EVAL)
sig$direction <- ifelse(up, "up_in_malignant_ECM",
                 ifelse(dn, "up_in_true_fibroblast", NA))
sig$mean_lfc <- rowMeans(lfcc, na.rm = TRUE)
sig$min_abs_auc_distance <- apply(abs(aucc - 0.5), 1, min)
sg <- sig[!is.na(sig$direction), ]
sg <- sg[order(sg$direction, -abs(sg$mean_lfc)), ]
p5_msg("signature genes: %d up in malignant ECM-like, %d up in true fibroblast",
       sum(sg$direction == "up_in_malignant_ECM"),
       sum(sg$direction == "up_in_true_fibroblast"))
p5_msg("  malignant ECM-like: %s", paste(utils::head(
  sg$gene[sg$direction == "up_in_malignant_ECM"], 25), collapse = " "))
p5_msg("  true fibroblast   : %s", paste(utils::head(
  sg$gene[sg$direction == "up_in_true_fibroblast"], 25), collapse = " "))
p5_tsv(sg, file.path(P5_TAB, "MALIGNANT_ECM_SIGNATURE.tsv"))
p5_tsv(sig[order(-abs(sig$mean_lfc)), ],
       file.path(P5_TAB, "MALIGNANT_ECM_VS_TRUE_FIBROBLAST.tsv"))
facts$signature <- list(criteria = list(min_abs_pseudobulk_log2FC = LFC_MIN,
                                        min_auc = AUC_MIN, min_pct_detected = PCT_MIN,
                                        must_hold_in = "every evaluable patient"),
                        n_up_malignant = sum(sg$direction == "up_in_malignant_ECM"),
                        n_up_fibroblast = sum(sg$direction == "up_in_true_fibroblast"),
                        evaluable_patients = EVAL,
                        no_classifier_trained = TRUE)

p5_sec("6. Programs, CNV metrics and curated scores across the groups")
CN <- file.path("results/phase5/programs/cnmf/primary")
Wt <- t(as.matrix(read.delim(file.path(CN, sprintf(
  "primary.gene_spectra_tpm.k_%d.dt_%s.txt", KSTAR, DT)),
  row.names = 1, check.names = FALSE)))                      # genes x K
colnames(Wt) <- paste0("P", seq_len(ncol(Wt)))
p5_msg("spectra for projection: %d genes x %d programs", nrow(Wt), ncol(Wt))
allf <- fib$cell_id
Xf <- tpm_on_genes(cnts[, allf, drop = FALSE], rownames(Wt))
Hf <- relative_usage(project_programs(Xf, Wt))
p5_msg("projected %d historical fibroblasts onto the %d programs", nrow(Hf), ncol(Hf))

proj <- data.frame(cell_id = allf, sample_id = fib$sample_id, group = fib$group,
                   malignancy_refined = fib$malignancy_refined,
                   malignancy_confidence = fib$malignancy_confidence,
                   cnv_burden = fib$cnv_burden, cnv_frac_gain = fib$cnv_frac_gain,
                   cnv_frac_loss = fib$cnv_frac_loss, cnv_mean_abs = fib$cnv_mean_abs,
                   as.data.frame(Hf), check.names = FALSE)
proj$dominant_program <- colnames(Hf)[max.col(Hf, ties.method = "first")]
p5_tsv(proj, file.path(P5_ECM, "M39_FIBROBLAST_PROGRAM_PROJECTION_ALL.tsv"))

cmp <- proj |> filter(group != "C_ambiguous") |>
  pivot_longer(all_of(colnames(Hf)), names_to = "program", values_to = "usage") |>
  group_by(sample_id, program, group) |>
  summarise(n = n(), median_usage = median(usage), mean_usage = mean(usage),
            .groups = "drop") |>
  pivot_wider(names_from = group, values_from = c(n, median_usage, mean_usage)) |>
  mutate(evaluable = sample_id %in% EVAL,
         delta_median = median_usage_A_malignant_ECM - median_usage_B_true_fibroblast)
p5_tsv(cmp, file.path(P5_ECM, "M39_PROGRAM_USAGE_A_VS_B.tsv"))
p5_msg("program usage difference (A - B), evaluable patients only:")
print(as.data.frame(cmp |> filter(evaluable) |>
        select(sample_id, program, delta_median) |>
        pivot_wider(names_from = sample_id, values_from = delta_median)),
      row.names = FALSE, digits = 3)

# CNV metrics, which Phase 4 already computed - reported as context, not as a
# new result, and never as independent validation of the malignancy call
cnv <- proj |> group_by(sample_id, group) |>
  summarise(n = n(), across(c(cnv_burden, cnv_frac_gain, cnv_frac_loss,
                              cnv_mean_abs), ~median(.x, na.rm = TRUE)),
            .groups = "drop")
p5_tsv(cnv, file.path(P5_ECM, "M39_CNV_METRICS_BY_GROUP.tsv"))
print(as.data.frame(cnv), row.names = FALSE, digits = 3)

p5_sec("7. Ambiguous fibroblast projection - descriptive only")
amb <- proj[proj$group == "C_ambiguous", ]
p5_msg("Ambiguous historical fibroblasts projected: %d", nrow(amb))
p5_msg("NOT reclassified: malignancy_refined stays Ambiguous for every one of them")
stopifnot(all(amb$malignancy_refined == "Ambiguous"))
ambs <- amb |> pivot_longer(all_of(colnames(Hf)), names_to = "program",
                            values_to = "usage") |>
  group_by(sample_id, program) |>
  summarise(n_cells = n(), median_usage = median(usage),
            mean_usage = mean(usage), .groups = "drop")
ref <- proj |> filter(group != "C_ambiguous") |>
  pivot_longer(all_of(colnames(Hf)), names_to = "program", values_to = "usage") |>
  group_by(group, program) |>
  summarise(median_usage = median(usage), .groups = "drop") |>
  pivot_wider(names_from = group, values_from = median_usage)
ambo <- ambs |> left_join(ref, by = "program") |>
  mutate(closer_to = ifelse(abs(median_usage - A_malignant_ECM) <=
                              abs(median_usage - B_true_fibroblast),
                            "A_malignant_ECM", "B_true_fibroblast"),
         note = "descriptive projection only; these cells remain Ambiguous")
p5_tsv(ambo, file.path(P5_TAB, "AMBIGUOUS_FIBROBLAST_PROGRAM_PROJECTION.tsv"))
p5_msg("Ambiguous cells per patient: %s",
       paste(sprintf("%s=%d", names(table(amb$sample_id)),
                     as.integer(table(amb$sample_id))), collapse = " "))

facts$program_projection <- list(
  method = "non-negative least squares by multiplicative update with the consensus spectra held FIXED; projected cells cannot influence the programs",
  n_projected = nrow(proj), programs = ncol(Hf))
p5_json(facts, file.path(P5_VAL, "m39_malignant_ecm_facts.json"))
p5_sec("M39 complete")
