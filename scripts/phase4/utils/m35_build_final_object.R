#!/usr/bin/env Rscript
# =============================================================================
# Phase 4 · M35 — build, save and validate the final Phase 4 object
#
# Phase 4 ONLY ADDS. Every Phase 1 provenance field, Phase 2 reduction, Phase 2
# cluster, Phase 2 annotation and the Phase 3 annotation_ccc survive verbatim,
# and that is asserted rather than assumed (§54).
# =============================================================================
options(stringsAsFactors = FALSE); options(future.globals.maxSize = +Inf)
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix); library(jsonlite); library(digest)
})
source("scripts/R/provenance_utils.R")
set.seed(42)

PHASE2 <- "results/phase2/phase2_final_object.rds"
PHASE2_MD5 <- "153d5f6acc70f9c05aa48cabc4f4ac2d"
OUT    <- "results/phase4/phase4_final_object.rds"
MAL    <- "results/phase4/malignancy"
TS     <- "results/phase4/tumor_states"
log_ <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
sec  <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

sec("1. Load the frozen Phase 2 object")
md5 <- calculate_file_checksum(PHASE2)
log_("Phase 2 md5 ", md5)
stopifnot(identical(md5, PHASE2_MD5))
obj <- readRDS(PHASE2)
meta_before   <- colnames(obj@meta.data)
red_before    <- Reductions(obj)
assay_before  <- Assays(obj)
cells_before  <- colnames(obj)
log_("cells ", length(cells_before), " | metadata cols ", length(meta_before),
     " | assays ", paste(assay_before, collapse=","),
     " | reductions ", paste(red_before, collapse=","))

sec("2. Attach Phase 4 layers (adds only)")
calls <- read.delim(file.path(MAL, "PHASE4_MALIGNANCY_CALLS.tsv"), stringsAsFactors = FALSE)
rownames(calls) <- calls$cell_id
stopifnot(setequal(calls$cell_id, cells_before))
calls <- calls[cells_before, , drop = FALSE]
P4_FIELDS <- c("scevan_call","scevan_confident_normal","scevan_clone","scevan_subclone_raw",
               "scevan_sample_reliable","population_class",
               "cnv_burden","cnv_frac_gain","cnv_frac_loss","cnv_mean_abs",
               "cnv_ref_threshold","cnv_elevated",
               "malignancy_phase2","malignancy_scevan","malignancy_scevan_sens",
               "malignancy_refined","malignancy_confidence","malignancy_rule",
               "malignancy_reason","annotation_ccc_phase3","annotation_ccc_refined",
               "marker_evidence","score_malignant_schwann_nc","score_nonmalignant_max",
               "pop_frac","pop_n_assessed","runs_agree")
added <- character(0)
for (f in intersect(P4_FIELDS, colnames(calls))) {
  if (f %in% meta_before) stop("COLLISION: Phase 4 field '", f, "' already exists in Phase 2 metadata.")
  obj[[f]] <- calls[[f]]; added <- c(added, f)
}
# scevan_sample: which SCEVAN run produced this cell's call
obj$scevan_sample <- as.character(obj$sample_id); added <- c(added, "scevan_sample")
# A compact human-readable CNA summary per cell.
obj$cnv_summary <- sprintf("burden=%.3f gain=%.3f loss=%.3f elevated=%s",
  ifelse(is.na(obj$cnv_burden), NA_real_, obj$cnv_burden),
  ifelse(is.na(obj$cnv_frac_gain), NA_real_, obj$cnv_frac_gain),
  ifelse(is.na(obj$cnv_frac_loss), NA_real_, obj$cnv_frac_loss), obj$cnv_elevated)
added <- c(added, "cnv_summary")

# Tumour states exist only for malignant cells; everything else is NA by design.
asg <- read.delim("results/phase4/tables/TUMOR_STATE_ASSIGNMENTS.tsv", stringsAsFactors = FALSE)
rownames(asg) <- asg$cell_id
obj$tumor_state_phase4 <- asg$tumor_state_phase4[match(cells_before, asg$cell_id)]
obj$tumor_clone_phase4 <- asg$tumor_clone_phase4[match(cells_before, asg$cell_id)]
obj$tumor_cluster_phase4 <- asg$tumor_cluster_phase4[match(cells_before, asg$cell_id)]
added <- c(added, "tumor_state_phase4","tumor_clone_phase4","tumor_cluster_phase4")
log_("Phase 4 fields added: ", length(added))
log_("cells with a tumour state: ", sum(!is.na(obj$tumor_state_phase4)))

sec("3. Preservation assertions (§6, §54)")
bad <- character(0)
# Re-read the frozen object and compare column by column, so a preservation
# failure names the offending column rather than reporting a whole-frame digest.
# (Phase 2 M13 taught this: a blunt digest fired on a change that had not happened.)
p2 <- readRDS(PHASE2)
for (f in meta_before) {
  x <- p2@meta.data[[f]]; y <- obj@meta.data[cells_before, f]
  same <- if (is.factor(x)) identical(as.character(x), as.character(y)) else isTRUE(all.equal(x, y))
  if (!same) bad <- c(bad, f)
}
if (length(bad)) stop("PRESERVATION FAILURE: Phase 2/3 metadata columns changed: ",
                      paste(bad, collapse = ", "))
log_("GUARD PASSED: all ", length(meta_before), " pre-existing metadata columns are byte-identical")
stopifnot(identical(Reductions(obj), red_before))
log_("GUARD PASSED: reductions unchanged (", paste(red_before, collapse=", "), ")")
for (r in red_before) {
  stopifnot(isTRUE(all.equal(Embeddings(obj, r), Embeddings(p2, r))))
}
log_("GUARD PASSED: every Phase 2 embedding is numerically identical (PCA and Harmony included)")
stopifnot(identical(Assays(obj), assay_before), identical(colnames(obj), cells_before))
log_("GUARD PASSED: assays and cell order unchanged")
stopifnot(identical(as.character(obj$annotation_ccc), as.character(obj$annotation_ccc_phase3)))
log_("GUARD PASSED: annotation_ccc_phase3 is a verbatim copy of the frozen Phase 3 layer")
stopifnot(!any(obj$annotation_ccc_refined == "MPNST-Tumor" & obj$malignancy_refined != "Malignant"))
log_("GUARD PASSED: no non-Malignant cell carries the refined MPNST-Tumor label")
stopifnot(all(is.na(obj$tumor_state_phase4) | obj$malignancy_refined == "Malignant"))
log_("GUARD PASSED: tumour states exist only on Malignant cells")
rm(p2); invisible(gc(FALSE))

sec("4. Save")
saveRDS(obj, OUT)
sz <- file.info(OUT)$size
md5o <- calculate_file_checksum(OUT)
sha <- digest(OUT, file = TRUE, algo = "sha256")
log_(sprintf("wrote %s (%.2f GB)", OUT, sz/1024^3))
log_("md5    ", md5o)
log_("sha256 ", sha)

sec("5. Reload and validate (§55)")
rm(obj); invisible(gc(FALSE))
o2 <- readRDS(OUT)
checks <- list()
chk <- function(name, ok, detail = "") {
  checks[[length(checks)+1]] <<- list(check = name, pass = isTRUE(ok), detail = detail)
  log_(sprintf("  [%s] %-46s %s", if (isTRUE(ok)) "PASS" else "FAIL", name, detail)) }
chk("cell count == 19716", ncol(o2) == 19716, paste("got", ncol(o2)))
chk("feature count (RNA) == 31764", nrow(o2[["RNA"]]) == 31764, paste("got", nrow(o2[["RNA"]])))
chk("assays preserved", identical(Assays(o2), assay_before), paste(Assays(o2), collapse=","))
chk("RNA layers present", length(Layers(o2[["RNA"]])) > 0, paste(length(Layers(o2[["RNA"]])), "layers"))
chk("SCT layers present", length(Layers(o2[["SCT"]])) > 0, paste(length(Layers(o2[["SCT"]])), "layers"))
chk("reductions preserved", identical(Reductions(o2), red_before), paste(Reductions(o2), collapse=","))
chk("no duplicate cell IDs", !any(duplicated(colnames(o2))), paste(sum(duplicated(colnames(o2))), "dups"))
chk("no duplicate feature names", !any(duplicated(rownames(o2[["RNA"]]))), "")
chk("all Phase 2 metadata columns retained", all(meta_before %in% colnames(o2@meta.data)),
    paste(length(setdiff(meta_before, colnames(o2@meta.data))), "missing"))
chk("Phase 4 fields present", all(added %in% colnames(o2@meta.data)),
    paste(length(setdiff(added, colnames(o2@meta.data))), "missing"))
chk("no missing malignancy_refined", !any(is.na(o2$malignancy_refined)),
    paste(sum(is.na(o2$malignancy_refined)), "NA"))
chk("no missing malignancy_confidence", !any(is.na(o2$malignancy_confidence)),
    paste(sum(is.na(o2$malignancy_confidence)), "NA"))
chk("no missing annotation_ccc_refined", !any(is.na(o2$annotation_ccc_refined)),
    paste(sum(is.na(o2$annotation_ccc_refined)), "NA"))
chk("malignancy_refined levels valid",
    all(o2$malignancy_refined %in% c("Malignant","Non-malignant","Ambiguous","Excluded-low-quality")),
    paste(sort(unique(o2$malignancy_refined)), collapse="/"))
for (r in Reductions(o2)) {
  e <- Embeddings(o2, r)
  chk(paste0("embedding finite: ", r), all(is.finite(e)),
      paste(sum(!is.finite(e)), "non-finite of", length(e)))
}
chk("annotation_ccc == annotation_ccc_phase3",
    identical(as.character(o2$annotation_ccc), as.character(o2$annotation_ccc_phase3)), "")
chk("refined MPNST-Tumor implies Malignant",
    !any(o2$annotation_ccc_refined == "MPNST-Tumor" & o2$malignancy_refined != "Malignant"), "")
chk("Ambiguous cells are never MPNST-Tumor",
    !any(o2$malignancy_refined == "Ambiguous" & o2$annotation_ccc_refined == "MPNST-Tumor"), "")
chk("tumour states only on Malignant cells",
    all(is.na(o2$tumor_state_phase4) | o2$malignancy_refined == "Malignant"), "")
chk("save/load md5 stable", identical(calculate_file_checksum(OUT), md5o), md5o)
npass <- sum(vapply(checks, function(x) x$pass, logical(1)))
log_("checks: ", npass, "/", length(checks), " passed")
if (npass != length(checks)) stop("M35 VALIDATION FAILED: ", length(checks)-npass, " check(s) failed.")

summary_out <- list(
  path = OUT, size_bytes = sz, md5 = md5o, sha256 = sha,
  cells = ncol(o2), features_rna = nrow(o2[["RNA"]]), features_sct = nrow(o2[["SCT"]]),
  assays = Assays(o2), default_assay = DefaultAssay(o2), reductions = Reductions(o2),
  rna_layers = Layers(o2[["RNA"]]), sct_layers = Layers(o2[["SCT"]]),
  n_metadata_columns = ncol(o2@meta.data),
  phase4_fields_added = added,
  phase2_metadata_columns_preserved = length(meta_before),
  validation_checks = checks, validation_passed = npass, validation_total = length(checks),
  malignancy_refined = as.list(table(o2$malignancy_refined)),
  malignancy_confidence = as.list(table(o2$malignancy_confidence)),
  annotation_ccc_refined = as.list(table(o2$annotation_ccc_refined)),
  tumor_states = as.list(table(o2$tumor_state_phase4)),
  tumor_clones = length(unique(na.omit(o2$tumor_clone_phase4))),
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  slurm = list(job_id = Sys.getenv("SLURM_JOB_ID")))
write_json(summary_out, "results/phase4/phase4_final_object_validation.json",
           auto_unbox = TRUE, pretty = TRUE, digits = 8, null = "null")
sec("M35 OBJECT BUILD AND VALIDATION COMPLETE")
