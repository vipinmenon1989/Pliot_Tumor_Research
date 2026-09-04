#!/usr/bin/env Rscript
# =============================================================================
# Phase 4 · M32 — refined CCC input preparation
#
# Builds a CCC input object keyed on `annotation_ccc_refined` instead of the
# frozen Phase 3 `annotation_ccc`, honouring EXACTLY the Phase 3 contract:
#
#     ccc_label   factor   the population identity every runner reads
#     ccc_sample  chr      sample_id
#     RNA assay, JoinLayers + LogNormalize(1e4), SCT dropped
#
# Because the contract is identical, scripts/phase3/ccc/run_liana.R,
# run_cellchat.R and run_cellphonedb.py are reused **completely unmodified**,
# at the same versions, thresholds, resources and expression basis. That is what
# makes the Phase 3 -> Phase 4 comparison a label-sensitivity analysis rather
# than a confound of tool drift (§39, §40).
#
# NOT re-run: the whole Phase 3 interaction space. Only what is needed for a
# fair tumour-label comparison (§37).
# =============================================================================
options(stringsAsFactors = FALSE); options(future.globals.maxSize = +Inf)
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix); library(jsonlite); library(digest)
})
source("scripts/R/provenance_utils.R")
set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
getopt <- function(f, d = NULL) { i <- which(args == paste0("--", f))
  if (!length(i)) return(d); args[i + 1] }
PHASE2   <- getopt("input", "results/phase2/phase2_final_object.rds")
CALLS    <- getopt("calls", "results/phase4/malignancy/PHASE4_MALIGNANCY_CALLS.tsv")
OUT      <- getopt("out-dir", "results/phase4/ccc_refinement")
CPDB     <- getopt("cpdb-dir", "results/phase4/ccc_refinement/cellphonedb_inputs")
TAB      <- getopt("tab-dir", "results/phase4/tables")
MIN_CELLS <- as.integer(getopt("min-cells", "10"))   # identical to Phase 3
for (d in c(OUT, CPDB, TAB)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

log_ <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
sec  <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

sec("1. Inputs")
md5_in <- digest(PHASE2, file = TRUE, algo = "md5")
log_(PHASE2, " md5 ", md5_in)
if (!identical(md5_in, "153d5f6acc70f9c05aa48cabc4f4ac2d"))
  stop("Phase 2 object checksum changed; refusing to proceed.")
calls <- read.delim(CALLS, stringsAsFactors = FALSE)
rownames(calls) <- calls$cell_id
log_("malignancy calls: ", nrow(calls), " cells")
log_("annotation_ccc_refined composition:")
print(sort(table(calls$annotation_ccc_refined), decreasing = TRUE))

obj <- readRDS(PHASE2)
stopifnot(setequal(colnames(obj), calls$cell_id))
calls <- calls[colnames(obj), , drop = FALSE]
log_("loaded Phase 2 object: ", ncol(obj), " cells; call table aligned")

# Attach the Phase 4 layers. Nothing existing is overwritten.
for (f in c("annotation_ccc_refined","annotation_ccc_phase3","malignancy_refined",
            "malignancy_confidence","malignancy_rule","scevan_call","scevan_clone",
            "population_class","cnv_burden")) {
  if (f %in% colnames(calls)) obj[[f]] <- calls[[f]]
}
stopifnot(identical(as.character(obj$annotation_ccc), as.character(obj$annotation_ccc_phase3)))
log_("GUARD PASSED: frozen annotation_ccc is unchanged and equals annotation_ccc_phase3")

# ---------------------------------------------------------------------------
sec("2. CCC-readiness of the refined labels")
# Same policy shape as Phase 3: populations that are not interpretable as a
# communicating cell identity are excluded outright, not silently included.
NOT_READY <- c("Ambiguous-unresolved", "Low-quality-excluded", "Uncertain",
               "Candidate-Malignant-Unresolved")
lab <- as.character(obj$annotation_ccc_refined)
smp <- as.character(obj$sample_id)
ready_pop <- sort(setdiff(unique(lab), NOT_READY))
log_("populations EXCLUDED as not CCC-ready: ",
     paste(sort(intersect(unique(lab), NOT_READY)), collapse = ", "))
log_("populations retained (", length(ready_pop), "): ", paste(ready_pop, collapse = ", "))

lv <- sort(unique(smp))
ct <- table(factor(lab, levels = sort(unique(lab))), factor(smp, levels = lv))
counts_df <- do.call(rbind, lapply(rownames(ct), function(p) do.call(rbind, lapply(lv, function(s) {
  n <- as.integer(ct[p, s])
  data.frame(annotation_ccc_refined = p, sample_id = s, patient = s, n_cells = n,
    phase4_ccc_ready = p %in% ready_pop,
    evaluable_in_sample = (p %in% ready_pop) && n >= MIN_CELLS,
    reason = if (!(p %in% ready_pop)) "excluded: not CCC-ready under the refined layer"
             else if (n < MIN_CELLS) sprintf("not evaluable: %d < %d cells in this sample", n, MIN_CELLS)
             else "evaluable", stringsAsFactors = FALSE) }))))
write.table(counts_df, file.path(TAB, "CCC_REFINED_SAMPLE_CELL_COUNTS.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(counts_df[counts_df$phase4_ccc_ready, ], row.names = FALSE)
log_("A population below the minimum is NOT EVALUABLE, never 'no signalling'.")

# ---------------------------------------------------------------------------
sec("3. Expression basis: RNA -> JoinLayers -> LogNormalize (identical to Phase 3)")
DefaultAssay(obj) <- "RNA"
log_("RNA layers before join: ", paste(Layers(obj[["RNA"]]), collapse = ";"))
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
log_("RNA layers after join : ", paste(Layers(obj[["RNA"]]), collapse = ";"))
stopifnot("counts" %in% Layers(obj[["RNA"]]))
obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize",
                     scale.factor = 1e4, verbose = FALSE)
log_("LogNormalize complete. Harmony coordinates are NEVER used as expression.")

# ---------------------------------------------------------------------------
sec("4. Build the refined CCC object")
keep_cells <- colnames(obj)[lab %in% ready_pop]
co <- subset(obj, cells = keep_cells)
co$ccc_label  <- factor(as.character(co$annotation_ccc_refined), levels = ready_pop)
co$ccc_sample <- as.character(co$sample_id)
DefaultAssay(co) <- "RNA"
if ("SCT" %in% Assays(co)) co[["SCT"]] <- NULL
for (r in setdiff(Reductions(co), c("postint_harmony","postint_umap_harmony"))) co[[r]] <- NULL
log_("refined CCC object: ", nrow(co), " features x ", ncol(co), " cells | ",
     nlevels(co$ccc_label), " populations")
print(table(co$ccc_label, co$ccc_sample))
rds <- file.path(OUT, "ccc_input_object_refined.rds")
saveRDS(co, rds)
md5_out <- digest(rds, file = TRUE, algo = "md5")
log_("saved ", rds, sprintf(" (%.2f GB, md5 %s)", file.info(rds)$size/1024^3, md5_out))

# ---------------------------------------------------------------------------
sec("5. Per-sample CellPhoneDB inputs (same format Phase 3 used)")
analysable <- lv
for (s in analysable) {
  sd <- file.path(CPDB, s); dir.create(sd, recursive = TRUE, showWarnings = FALSE)
  cells_s <- colnames(co)[co$ccc_sample == s]
  lab_s <- as.character(co$ccc_label[match(cells_s, colnames(co))])
  okp <- names(which(table(lab_s) >= MIN_CELLS))
  cells_s <- cells_s[lab_s %in% okp]; lab_s <- lab_s[lab_s %in% okp]
  m <- GetAssayData(co, assay = "RNA", layer = "data")[, cells_s, drop = FALSE]
  Matrix::writeMM(as(m, "dgCMatrix"), file.path(sd, "matrix.mtx"))
  writeLines(colnames(m), file.path(sd, "barcodes.tsv"))
  writeLines(rownames(m), file.path(sd, "features.tsv"))
  write.table(data.frame(Cell = cells_s, cell_type = lab_s),
              file.path(sd, "meta.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  log_(sprintf("  %s: %d cells x %d genes, %d populations -> %s",
               s, ncol(m), nrow(m), length(okp), sd))
}

# ---------------------------------------------------------------------------
sec("6. Phase 3 vs Phase 4 population comparison")
p3 <- as.character(obj$annotation_ccc_phase3)
p4 <- as.character(obj$annotation_ccc_refined)
cmp <- as.data.frame(table(phase3 = p3, phase4 = p4))
cmp <- cmp[cmp$Freq > 0, ]
cmp <- cmp[order(-cmp$Freq), ]
write.table(cmp, file.path(TAB, "CCC_LABEL_PHASE3_VS_PHASE4.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(cmp, row.names = FALSE)
moved <- sum(p3 != p4)
log_("cells whose CCC label changed: ", moved, sprintf(" (%.2f%%)", 100*moved/length(p3)))

facts <- list(
  phase2_input = list(path = PHASE2, md5 = md5_in),
  malignancy_calls = list(path = CALLS, md5 = digest(CALLS, file = TRUE, algo = "md5")),
  refined_ccc_object = list(path = rds, md5 = md5_out, cells = ncol(co),
                            populations = levels(co$ccc_label)),
  min_cells = MIN_CELLS,
  not_ready_excluded = NOT_READY,
  expression_basis = "RNA assay, JoinLayers, LogNormalize (scale.factor 1e4) - identical to Phase 3",
  reuse_statement = paste("scripts/phase3/ccc/run_liana.R, run_cellchat.R and",
    "run_cellphonedb.py are reused UNMODIFIED at Phase 3 versions; only the",
    "ccc_label definition differs. This is a label-sensitivity analysis."),
  cells_relabelled = moved,
  population_counts = as.list(table(co$ccc_label)),
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  slurm = list(job_id = Sys.getenv("SLURM_JOB_ID")))
write_json(facts, file.path(OUT, "m32_refined_ccc_input_facts.json"),
           auto_unbox = TRUE, pretty = TRUE, digits = 8, null = "null")
sec("REFINED CCC INPUT PREPARATION COMPLETE")
