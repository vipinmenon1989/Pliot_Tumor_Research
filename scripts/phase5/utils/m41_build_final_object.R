#!/usr/bin/env Rscript
# =============================================================================
# Phase 5 - M41 - build, validate and freeze the Phase 5 object.
#
# The Phase 4 object is opened READ-ONLY, Phase 5 fields are ADDED, and nothing
# pre-existing is touched. Preservation is asserted column by column against a
# copy taken before any field was added - the lesson carried from Phase 2 M13
# and Phase 4 M35 - not by a whole-frame digest that can hide a single change.
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject) })
source("scripts/phase5/utils/phase5_common.R")
source("scripts/phase5/utils/program_projection.R")
set.seed(42)
args <- commandArgs(trailingOnly = TRUE)
KSTAR <- as.integer(args[1]); DT <- if (length(args) > 1) args[2] else "0_1"
CNMF <- "results/phase5/programs/cnmf"
OUT  <- "results/phase5/phase5_final_object.rds"
facts <- list(milestone = "M41", K = KSTAR,
              generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))

p5_sec("1. Verify and load the frozen Phase 4 object")
stopifnot(p5_md5(P4_OBJECT) == P4_MD5, p5_sha256(P4_OBJECT) == P4_SHA256)
p5_msg("Phase 4 object verified by md5 and sha256")
obj <- readRDS(P4_OBJECT)
stopifnot(ncol(obj) == P4_FROZEN$total_cells)
md_before <- obj@meta.data                    # the preservation reference
red_before <- lapply(Reductions(obj), function(r) Embeddings(obj, r))
names(red_before) <- Reductions(obj)
cells_before <- colnames(obj)
p5_msg("cells %d | metadata columns %d | reductions %s", ncol(obj),
       ncol(md_before), paste(Reductions(obj), collapse = ", "))

p5_sec("2. Phase 5 program scores for EVERY cell")
# Malignant cells carry their own cNMF consensus usages. All other cells are
# PROJECTED onto the fixed spectra so the field is defined object-wide - a
# projection is descriptive and changes no Phase 4 call.
U <- as.matrix(read.delim(file.path(CNMF, "primary", sprintf(
  "primary.usages.k_%d.dt_%s.consensus.txt", KSTAR, DT)),
  row.names = 1, check.names = FALSE))
PROGS <- paste0("P", seq_len(ncol(U))); colnames(U) <- PROGS
REL <- U / rowSums(U)
p5_msg("cNMF consensus usages: %d malignant cells x %d programs", nrow(REL), ncol(REL))

Wt <- t(as.matrix(read.delim(file.path(CNMF, "primary", sprintf(
  "primary.gene_spectra_tpm.k_%d.dt_%s.txt", KSTAR, DT)),
  row.names = 1, check.names = FALSE)))
colnames(Wt) <- PROGS

rna <- JoinLayers(obj[["RNA"]])
cnts_all <- LayerData(rna, layer = "counts")
other <- setdiff(colnames(obj), rownames(REL))
p5_msg("projecting %d non-malignant / ambiguous / excluded cells onto the fixed spectra", length(other))
Xo <- tpm_on_genes(cnts_all[, other, drop = FALSE], rownames(Wt))
Ho <- relative_usage(project_programs(Xo, Wt))
rm(Xo, rna); invisible(gc())

ALL <- matrix(NA_real_, nrow = ncol(obj), ncol = length(PROGS),
              dimnames = list(colnames(obj), PROGS))
ALL[rownames(REL), ] <- REL
ALL[rownames(Ho), ]  <- Ho
stopifnot(!anyNA(ALL))

lab <- read.delim(file.path(P5_TAB, "MALIGNANT_PROGRAMS.tsv"))
lmap <- setNames(lab$program_label, lab$program_id)
p5_msg("program labels: %s", paste(sprintf("%s=%s", names(lmap), lmap), collapse = "  "))

p5_sec("3. Add Phase 5 fields - nothing overwritten")
newcols <- list()
for (p in PROGS) newcols[[paste0("program_", p, "_score")]] <- ALL[, p]
dom <- PROGS[max.col(ALL, ties.method = "first")]
newcols$dominant_malignant_program <- dom
newcols$dominant_malignant_program_label <- unname(lmap[dom])
newcols$dominant_program_score <- apply(ALL, 1, max)
# confidence: how far the top program is above the second, on the same cell
srt <- t(apply(ALL, 1, function(v) sort(v, decreasing = TRUE)[1:2]))
gap <- srt[, 1] - srt[, 2]
newcols$dominant_program_margin <- gap
newcols$malignant_program_confidence <- ifelse(
  gap >= 0.20, "High", ifelse(gap >= 0.10, "Moderate", "Low"))
# a cell that was never part of program discovery is flagged as projected
newcols$program_score_source <- ifelse(colnames(obj) %in% rownames(REL),
                                       "cnmf_consensus", "projected")
clash <- intersect(names(newcols), colnames(md_before))
if (length(clash)) stop("M41: Phase 5 field name collides with an existing field: ",
                        paste(clash, collapse = ", "))
for (n in names(newcols)) obj[[n]] <- unname(newcols[[n]])
p5_msg("added %d Phase 5 metadata fields", length(newcols))
print(table(obj$program_score_source))
print(table(obj$malignant_program_confidence))
print(table(obj$dominant_malignant_program_label,
            obj$malignancy_refined)[, "Malignant", drop = FALSE])

p5_sec("4. Preservation guards - asserted column by column")
md_after <- obj@meta.data
g <- list()
g$cells_unchanged <- identical(colnames(obj), cells_before)
g$cell_order_unchanged <- identical(rownames(md_after), rownames(md_before))
same <- vapply(colnames(md_before), function(cn)
  identical(md_before[[cn]], md_after[[cn]]), logical(1))
g$all_phase1_4_columns_identical <- all(same)
if (!all(same)) p5_msg("CHANGED COLUMNS: %s", paste(names(same)[!same], collapse = ", "))
g$only_expected_columns_added <- setequal(setdiff(colnames(md_after),
                                                  colnames(md_before)),
                                          names(newcols))
g$reductions_unchanged <- identical(Reductions(obj), names(red_before)) &&
  all(vapply(names(red_before), function(r)
    identical(Embeddings(obj, r), red_before[[r]]), logical(1)))
g$assays_unchanged <- identical(sort(Assays(obj)), sort(c("RNA", "SCT")))
g$malignancy_refined_unchanged <- identical(md_before$malignancy_refined,
                                            md_after$malignancy_refined)
g$scevan_clone_unchanged <- identical(md_before$scevan_clone, md_after$scevan_clone)
g$tumor_state_unchanged <- identical(md_before$tumor_state_phase4,
                                     md_after$tumor_state_phase4)
g$frozen_counts_hold <- identical(as.integer(table(md_after$malignancy_refined)[
  names(P4_FROZEN$refined)]), as.integer(P4_FROZEN$refined))
for (n in names(g)) p5_msg("  [%s] %s", if (isTRUE(g[[n]])) "OK  " else "FAIL", n)
if (!all(unlist(g))) stop("M41: a preservation guard failed - STOP")
facts$preservation <- g

p5_sec("5. Save, reload, validate")
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
saveRDS(obj, OUT)
sz <- file.info(OUT)$size
md5 <- p5_md5(OUT); sha <- p5_sha256(OUT)
p5_msg("saved %s | %s bytes | md5 %s", OUT, format(sz, big.mark = ","), md5)
p5_msg("sha256 %s", sha)
rm(obj); invisible(gc())
o2 <- readRDS(OUT)
v <- list(
  reload_cells = ncol(o2) == P4_FROZEN$total_cells,
  reload_cell_order = identical(colnames(o2), cells_before),
  reload_assays = identical(sort(Assays(o2)), sort(c("RNA", "SCT"))),
  reload_reductions = identical(Reductions(o2), names(red_before)),
  reload_phase4_columns = all(vapply(colnames(md_before), function(cn)
    identical(md_before[[cn]], o2@meta.data[[cn]]), logical(1))),
  reload_phase5_columns = all(names(newcols) %in% colnames(o2@meta.data)),
  reload_no_na_program_scores = !anyNA(o2@meta.data[[paste0("program_", PROGS[1], "_score")]]),
  reload_md5_stable = p5_md5(OUT) == md5)
for (n in names(v)) p5_msg("  [%s] %s", if (isTRUE(v[[n]])) "OK  " else "FAIL", n)
if (!all(unlist(v))) stop("M41: a reload validation failed - STOP")
facts$validation <- v
facts$phase5_object <- list(path = OUT, size_bytes = sz, md5 = md5, sha256 = sha,
                            cells = ncol(o2),
                            metadata_columns = ncol(o2@meta.data),
                            phase5_fields_added = names(newcols))

p5_sec("6. Phase 4 object still untouched")
stopifnot(p5_md5(P4_OBJECT) == P4_MD5)
p5_msg("Phase 4 md5 unchanged: %s", P4_MD5)
facts$phase4_unchanged <- TRUE
p5_json(facts, file.path(P5_VAL, "m41_final_object_validation.json"))
p5_sec("M41 object build complete")
