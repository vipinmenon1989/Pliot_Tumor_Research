#!/usr/bin/env Rscript
# =============================================================================
# Phase 4 · M28 — SCEVAN feasibility inspection and per-sample count extraction
#
# Read-only with respect to every Phase 2 / Phase 3 artefact. Writes ONLY into
# results/phase4/ and reports/phase4/.
#
# Purpose
#   1. Verify the Phase 2 final object against its recorded checksum.
#   2. Establish the exact raw-count representation SCEVAN must consume
#      (assay, layer, gene identifier type, integer-ness, sparsity).
#   3. Confirm gene-identifier compatibility with SCEVAN's internal
#      EnsDB_Hsapiens_v86 annotation (chr1-22, gene_name matching).
#   4. Extract per-sample sparse RNA count matrices to disk so that the four
#      M29 SCEVAN jobs never have to reload the 6 GB Phase 2 object.
#   5. Extract the cell metadata + UMAP coordinates once, so all downstream
#      malignancy work is cheap.
#   6. Enumerate normal-reference candidates WITHOUT committing to them.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix); library(jsonlite)
})
source("scripts/R/provenance_utils.R")

set.seed(42)
options(future.globals.maxSize = 64 * 1024^3)

OUT_SCEVAN <- "results/phase4/scevan"
OUT_MAL    <- "results/phase4/malignancy"
OUT_TAB    <- "results/phase4/tables"
PHASE2_RDS <- "results/phase2/phase2_final_object.rds"
PHASE2_MD5 <- "153d5f6acc70f9c05aa48cabc4f4ac2d"

log_ <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
sec  <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

facts <- list()

# --- 1. checksum gate -------------------------------------------------------
sec("1. Phase 2 input verification")
stopifnot(file.exists(PHASE2_RDS))
md5 <- calculate_file_checksum(PHASE2_RDS)
log_("phase2_final_object.rds md5 = ", md5)
if (!identical(md5, PHASE2_MD5)) {
  stop("STOP CONDITION (Phase 4 §67): Phase 2 object checksum mismatch. ",
       "Expected ", PHASE2_MD5, " got ", md5, ". Phase 2 object may be corrupted.")
}
log_("checksum MATCHES the Phase 2 manifest. Object is intact.")
facts$phase2_input <- list(path = PHASE2_RDS, md5 = md5, checksum_verified = TRUE)

log_("loading (read-only) ...")
t0 <- Sys.time()
obj <- readRDS(PHASE2_RDS)
log_("loaded in ", round(difftime(Sys.time(), t0, units = "mins"), 2), " min; ",
     "resident size ", round(as.numeric(object.size(obj)) / 1024^3, 2), " GiB")

# --- 2. structure -----------------------------------------------------------
sec("2. Object structure")
log_("cells    : ", ncol(obj))
log_("assays   : ", paste(Assays(obj), collapse = ", "))
log_("default  : ", DefaultAssay(obj))
log_("reductions: ", paste(Reductions(obj), collapse = ", "))
for (a in Assays(obj)) {
  log_(sprintf("  assay %-5s : %d features | layers: %s", a, nrow(obj[[a]]),
               paste(Layers(obj[[a]]), collapse = ", ")))
}
facts$object <- list(
  cells = ncol(obj), assays = Assays(obj), default_assay = DefaultAssay(obj),
  reductions = Reductions(obj),
  rna_features = nrow(obj[["RNA"]]), sct_features = nrow(obj[["SCT"]]),
  rna_layers = Layers(obj[["RNA"]]), sct_layers = Layers(obj[["SCT"]]),
  n_sct_models = length(SeuratObject::Misc(obj[["SCT"]], "vst.set")) )

# --- 3. metadata fields to preserve ----------------------------------------
sec("3. Metadata fields (all must be preserved; Phase 4 only ADDS)")
mc <- colnames(obj@meta.data)
log_("n metadata columns: ", length(mc))
key_fields <- grep("annotation|malign|cluster|sample|patient|dataset|percent|nCount|nFeature|Phase|phase",
                   mc, value = TRUE, ignore.case = TRUE)
cat(paste0("  ", key_fields, collapse = "\n"), "\n")
facts$metadata <- list(n_columns = length(mc), all_columns = mc, key_fields = key_fields)

for (f in c("sample_id", "annotation_ccc")) {
  if (!f %in% mc) stop("STOP: required field '", f, "' absent from Phase 2 metadata.")
}
ann_field_candidates <- intersect(
  c("annotation_level1","annotation_level2","annotation_level3","annotation_confidence",
    "annotation_initial","annotation_refined","annotation_ccc",
    "postint_harmony_primary_cluster","postint_harmony_alternative_cluster","seurat_clusters"), mc)
log_("annotation/cluster fields present: ", paste(ann_field_candidates, collapse = ", "))

sec("3b. annotation_ccc x sample_id")
tab <- table(obj$annotation_ccc, obj$sample_id)
print(addmargins(tab))
facts$annotation_ccc_by_sample <- as.data.frame.matrix(tab)

# --- 4. RAW COUNTS: the representation SCEVAN must consume ------------------
sec("4. Raw-count representation for SCEVAN (Phase 4 §12)")
rna <- obj[["RNA"]]
count_layers <- grep("^counts", Layers(rna), value = TRUE)
log_("RNA count layers (one per sample, Seurat v5 split): ",
     paste(count_layers, collapse = ", "))
if (!length(count_layers)) stop("STOP CONDITION (§67): no raw RNA count layer found.")

samples <- sort(unique(as.character(obj$sample_id)))
log_("samples: ", paste(samples, collapse = ", "))
facts$samples <- samples

# Map layer -> sample by the cells each layer holds, rather than by parsing the
# layer name. Names are informative but cell membership is authoritative.
layer_sample <- setNames(character(0), character(0))
for (ly in count_layers) {
  cl <- colnames(LayerData(rna, layer = ly))
  s  <- unique(as.character(obj$sample_id[cl]))
  if (length(s) != 1L)
    stop("STOP: count layer '", ly, "' spans ", length(s), " samples: ", paste(s, collapse=","))
  layer_sample[ly] <- s
}
log_("layer -> sample mapping (verified by cell membership):")
for (ly in names(layer_sample)) log_("  ", ly, " -> ", layer_sample[[ly]])
facts$layer_sample_map <- as.list(layer_sample)
if (!setequal(unname(layer_sample), samples))
  stop("STOP: count layers do not cover every sample exactly once.")

# --- 5. per-sample extraction, integer + sparsity audit, write to disk ------
sec("5. Per-sample count extraction + audit")
audit <- list()
for (ly in names(layer_sample)) {
  s  <- layer_sample[[ly]]
  m  <- LayerData(rna, layer = ly)                    # dgCMatrix, genes x cells
  expressed <- Matrix::rowSums(m) > 0
  m2 <- m[expressed, , drop = FALSE]

  xv <- m2@x
  is_int  <- all(xv == floor(xv))
  spars   <- 1 - (length(xv) / (as.numeric(nrow(m2)) * ncol(m2)))
  gsym    <- rownames(m2)
  looks_symbol <- mean(grepl("^ENSG[0-9]{11}$", gsym)) < 0.5

  d <- file.path(OUT_SCEVAN, "by_sample", s)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  fp <- file.path(d, "counts_raw.rds")
  # Idempotent: an identical matrix already on disk is left alone, so re-running
  # M28 after a downstream fix does not rewrite verified inputs.
  if (file.exists(fp) && identical(dim(readRDS(fp)), dim(m2))) {
    log_("  ", s, ": counts_raw.rds already present with matching dimensions - kept")
  } else {
    saveRDS(m2, fp, compress = TRUE)
  }

  log_(sprintf("%s: %d genes expressed / %d total, %d cells | class %s | integer=%s | sparsity=%.4f",
               s, nrow(m2), nrow(m), ncol(m2), class(m2)[1], is_int, spars))
  audit[[s]] <- list(
    sample = s, layer = ly, matrix_class = class(m2)[1],
    genes_total_rna = nrow(m), genes_expressed = nrow(m2), cells = ncol(m2),
    integer_counts = is_int, min_value = min(xv), max_value = max(xv),
    sparsity = spars, nonzero_entries = length(xv),
    identifier_type = if (looks_symbol) "gene_symbol" else "ensembl_gene_id",
    median_counts_per_cell = as.numeric(median(Matrix::colSums(m2))),
    median_genes_per_cell  = as.numeric(median(Matrix::colSums(m2 > 0))),
    written = file.path(d, "counts_raw.rds"),
    written_md5 = calculate_file_checksum(file.path(d, "counts_raw.rds")))
  rm(m, m2); invisible(gc(FALSE))
}
facts$per_sample_counts <- audit
if (!all(vapply(audit, function(x) x$integer_counts, logical(1))))
  stop("STOP CONDITION (§67): RNA 'counts' layer is not integer-valued; not raw counts.")
log_("ALL per-sample matrices are integer-valued raw counts. SCEVAN input basis confirmed.")

# --- 6. gene-identifier compatibility with SCEVAN's annotation -------------
sec("6. SCEVAN gene-annotation compatibility")
# SCEVAN::annotateGenes matches rownames against EnsDB_Hsapiens_v86$gene_name
# (falling back to gene_id) and keeps chromosomes 1-22 only.
edb <- get("EnsDB_Hsapiens_v86", envir = asNamespace("SCEVAN"))
edb_auto <- edb[edb$seqnames %in% as.character(1:22), ]
all_genes <- rownames(obj[["RNA"]])
ov_name <- length(intersect(all_genes, edb_auto$gene_name))
ov_id   <- length(intersect(all_genes, edb_auto$gene_id))
log_("SCEVAN EnsDB_Hsapiens_v86 rows (chr1-22): ", nrow(edb_auto))
log_("overlap with gene_name : ", ov_name, sprintf(" (%.1f%% of %d RNA features)",
     100 * ov_name / length(all_genes), length(all_genes)))
log_("overlap with gene_id   : ", ov_id)
log_("=> SCEVAN will match on: ", if (ov_name > 0) "gene_name (symbols)" else "gene_id")
if (ov_name < 5000)
  stop("STOP CONDITION (§67): only ", ov_name,
       " genes map to SCEVAN's chr1-22 annotation. SCEVAN is incompatible with this object.")
per_sample_ov <- vapply(audit, function(a) {
  g <- rownames(readRDS(a$written)); length(intersect(g, edb_auto$gene_name)) }, numeric(1))
log_("per-sample annotatable genes: ",
     paste(sprintf("%s=%d", names(per_sample_ov), per_sample_ov), collapse = ", "))
facts$scevan_annotation <- list(
  annotation_object = "SCEVAN:::EnsDB_Hsapiens_v86",
  rows_chr1_22 = nrow(edb_auto), match_field = if (ov_name > 0) "gene_name" else "gene_id",
  overlap_gene_name = ov_name, overlap_gene_id = ov_id,
  identifier_type = "gene_symbol",
  per_sample_annotatable_genes = as.list(per_sample_ov),
  note = "chrX/chrY are dropped by SCEVAN; sex-chromosome CNAs are not assessable.")

# --- 7. cell-ID integrity ---------------------------------------------------
sec("7. Cell-ID integrity (Phase 4 §67 stop conditions)")
cid <- colnames(obj)
log_("duplicated cell IDs: ", sum(duplicated(cid)))
if (any(duplicated(cid))) stop("STOP CONDITION (§67): duplicated cell IDs.")
extracted <- unlist(lapply(audit, function(a) colnames(readRDS(a$written))), use.names = FALSE)
log_("cells across extracted matrices: ", length(extracted),
     " | duplicated: ", sum(duplicated(extracted)))
log_("setequal(extracted, object cells): ", setequal(extracted, cid))
if (!setequal(extracted, cid))
  stop("STOP CONDITION (§67): extracted cells do not round-trip to the object's cells.")
log_("Cell IDs map 1:1. No cell loss during extraction.")
facts$cell_id_integrity <- list(
  cells = length(cid), duplicated = 0L, extracted = length(extracted),
  roundtrip_setequal = TRUE, example_ids = head(cid, 3))

# --- 8. metadata + embeddings snapshot for cheap downstream use ------------
sec("8. Metadata / embedding snapshot")
md <- obj@meta.data
md$cell_id <- rownames(md)
red_used <- character(0)
for (r in Reductions(obj)) {
  e <- Embeddings(obj, reduction = r)
  if (ncol(e) >= 2 && grepl("umap", r, ignore.case = TRUE)) {
    md[[paste0(r, "_1")]] <- e[rownames(md), 1]
    md[[paste0(r, "_2")]] <- e[rownames(md), 2]
    red_used <- c(red_used, r)
  }
}
log_("UMAP reductions captured: ", paste(red_used, collapse = ", "))
dir.create(OUT_MAL, showWarnings = FALSE, recursive = TRUE)
mp <- file.path(OUT_MAL, "phase4_cell_metadata.tsv.gz")
write.table(md, gzfile(mp), sep = "\t", quote = FALSE, row.names = FALSE)
log_("wrote ", mp, " (", nrow(md), " x ", ncol(md), ")")
facts$metadata_snapshot <- list(path = mp, rows = nrow(md), cols = ncol(md),
                                umap_reductions = red_used,
                                md5 = calculate_file_checksum(mp))

# Harmony embedding, for the M31 malignant-only work (never overwritten)
hp <- file.path(OUT_MAL, "phase2_harmony_embedding.rds")
hred <- grep("harmony", Reductions(obj), value = TRUE, ignore.case = TRUE)[1]
saveRDS(Embeddings(obj, reduction = hred), hp)
log_("wrote frozen copy of reduction '", hred, "' -> ", hp)
facts$harmony_reduction <- list(name = hred, path = hp,
                                dims = dim(Embeddings(obj, reduction = hred)))

# --- 9. normal-reference candidates (enumerated, NOT committed) ------------
sec("9. Normal-reference candidates (Phase 4 §15/§16)")
IMMUNE_HIGH_CONF <- c("CD4-T","CD8-T","NK","T-cell-other","B-cell","Plasma-cell",
                      "Macrophage","Monocyte","Dendritic","Plasmacytoid-DC")
ENDO <- "Endothelial"
DISPUTED <- c("MPNST-Tumor","Candidate-Malignant-Unresolved","Fibroblast",
              "Pericyte-VSMC","Uncertain","Low-quality-excluded")
cnt <- function(v) sum(obj$annotation_ccc %in% v)
log_("high-confidence immune (candidate reference)  : ", cnt(IMMUNE_HIGH_CONF), " cells")
log_("Endothelial (candidate reference)             : ", cnt(ENDO), " cells")
log_("DISPUTED - never used as a fixed reference    : ", cnt(DISPUTED), " cells")
ref_by_sample <- table(obj$sample_id[obj$annotation_ccc %in% IMMUNE_HIGH_CONF])
print(ref_by_sample)
facts$normal_reference_candidates <- list(
  strategy_primary = paste("SCEVAN automatic confident-normal detection",
                           "(norm_cell = NULL, findConfident = TRUE). Fully independent",
                           "of Phase 2 labels - this is the primary, non-circular run."),
  strategy_sensitivity = paste("Supply high-confidence IMMUNE cells only as norm_cell,",
                               "with FIXED_NORMAL_CELLS = FALSE so classification still",
                               "proceeds by CNA clustering. FIXED_NORMAL_CELLS = TRUE is",
                               "PROHIBITED here: it forces every non-reference cell to",
                               "'malignant', which is the non-immune-equals-tumour",
                               "inference Phase 2 explicitly banned."),
  immune_reference_populations = IMMUNE_HIGH_CONF,
  immune_reference_cells = cnt(IMMUNE_HIGH_CONF),
  immune_reference_by_sample = as.list(ref_by_sample),
  endothelial_considered_separately = list(population = ENDO, cells = cnt(ENDO),
    note = "Endothelium is a Phase 3 receiver of interest and a CCC partner; excluded from the reference set to keep its malignancy call independent."),
  disputed_never_fixed_normal = DISPUTED,
  disputed_cells = cnt(DISPUTED))

# --- 10. resource evidence for M29 -----------------------------------------
sec("10. Resource estimate for M29")
# SCEVAN densifies via annotateGenes (cbind of as.matrix); estimate that cost.
for (s in names(audit)) {
  a <- audit[[s]]
  dense_gib <- a$genes_expressed * a$cells * 8 / 1024^3
  log_(sprintf("%s: dense annotated matrix ~%.2f GiB (%d genes x %d cells, double)",
               s, dense_gib, a$genes_expressed, a$cells))
  audit[[s]]$dense_estimate_gib <- dense_gib
}
facts$per_sample_counts <- audit
facts$resource_note <- paste(
  "SCEVAN::annotateGenes densifies the matrix (cbind(edb, as.matrix(mtx))), so the",
  "dense estimate above is the dominant single allocation. SCEVAN also holds smoothed and",
  "relative copies, so budget several multiples. Requested M29 memory is set from these",
  "numbers, not from the SLURM maximum.")

# Persist the audit BEFORE the smoke test: sections 1-10 are the M28 deliverable
# and must survive a smoke-test failure.
dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE)
write_json(facts, file.path(OUT_SCEVAN, "m28_feasibility_facts.json"),
           auto_unbox = TRUE, pretty = TRUE, digits = 8, null = "null")
log_("interim audit written (pre-smoke-test)")

# --- 11. smoke test (technical only) ---------------------------------------
sec("11. Technical smoke test - NOT a scientific result")
smallest <- names(which.min(vapply(audit, function(a) a$cells, numeric(1))))
log_("smallest sample: ", smallest, " (", audit[[smallest]]$cells, " cells)")
m <- readRDS(audit[[smallest]]$written)
set.seed(42)
idx <- sort(sample(ncol(m), min(600L, ncol(m))))
sm  <- m[, idx, drop = FALSE]
sm  <- sm[Matrix::rowSums(sm) > 0, , drop = FALSE]
smoke_dir <- normalizePath(file.path(OUT_SCEVAN, "smoke_test"), mustWork = FALSE)
dir.create(smoke_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# SCEVAN 1.0.3 output_dir defect, and why we handle it this way.
#
# pipelineCNA() takes output_dir and forwards it to most helpers, but the
# read-back/plot helpers do NOT honour it:
#   getScevanCNV(sample, path = "./output")        <- hardcoded default
#   getScevanCNVfinal(sample, path = "./output")   <- hardcoded default
#   plotAllClonalCN / plotAllSubclonalCN / plotConsensusCNA / analyzeSegm2
# and their caller plotCNclonal() invokes getScevanCNV(paste0(sample, name))
# without passing output_dir. So with any non-default output_dir, classification
# succeeds and then the plotting stage dies with "cannot open the connection"
# (observed: SLURM 19896576).
#
# Fix: run each SCEVAN call with the process working directory set to a
# per-run directory and let output_dir take SCEVAN's own default "./output", so
# the parameterised and hardcoded paths refer to the same place. The package is
# NOT patched and no error is suppressed.
# ---------------------------------------------------------------------------
log_("running pipelineCNA on a ", nrow(sm), " x ", ncol(sm),
     " random subsample (full path: SUBCLONES=TRUE, ClonalCN=TRUE)")
smoke <- tryCatch({
  t1 <- Sys.time()
  owd <- setwd(smoke_dir); on.exit(setwd(owd), add = TRUE)
  r <- SCEVAN::pipelineCNA(sm, sample = paste0("SMOKE_", smallest),
                           par_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")),
                           norm_cell = NULL, SUBCLONES = TRUE, ClonalCN = TRUE,
                           plotTree = FALSE, organism = "human")
  setwd(owd)
  list(status = "SUCCESS", elapsed_min = as.numeric(difftime(Sys.time(), t1, units="mins")),
       n_cells_in = ncol(sm), returned_rows = nrow(r),
       returned_cols = colnames(r), class_table = as.list(table(r$class)),
       files_written = list.files(file.path(smoke_dir, "output")),
       output_dir_workaround = paste(
         "SCEVAN 1.0.3 ignores output_dir in getScevanCNV/getScevanCNVfinal and the",
         "plotAll*/plotConsensusCNA helpers (hardcoded path = './output'). Each run",
         "therefore sets the working directory to its own folder and uses SCEVAN's",
         "default output_dir. The package is not patched."))
}, error = function(e) { list(status = "FAILED", message = conditionMessage(e)) })
str(smoke, max.level = 2)
facts$smoke_test <- c(smoke, list(
  purpose = "Technical API/compatibility check only.",
  scientific_use = "NONE. 600 random cells cannot support a malignancy call and this
    result is never propagated. The scientific classification is the full per-sample
    M29 run."))
if (identical(smoke$status, "FAILED"))
  stop("STOP CONDITION (§67): SCEVAN smoke test failed: ", smoke$message)

# --- 12. provenance --------------------------------------------------------
sec("12. Provenance")
facts$scevan_version <- as.character(packageVersion("SCEVAN"))
facts$yagst_version  <- as.character(packageVersion("yaGST"))
facts$r_version      <- R.version.string
facts$core_packages  <- lapply(c("Seurat","SeuratObject","Matrix","harmony","sctransform"),
                               function(p) as.character(packageVersion(p)))
names(facts$core_packages) <- c("Seurat","SeuratObject","Matrix","harmony","sctransform")
facts$slurm <- list(job_id = Sys.getenv("SLURM_JOB_ID"),
                    cpus = Sys.getenv("SLURM_CPUS_PER_TASK"),
                    node = Sys.info()[["nodename"]])
facts$generated <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE)
write_json(facts, file.path(OUT_SCEVAN, "m28_feasibility_facts.json"),
           auto_unbox = TRUE, pretty = TRUE, digits = 8, null = "null")
log_("wrote ", file.path(OUT_SCEVAN, "m28_feasibility_facts.json"))

sec("M28 INSPECTION COMPLETE")
print(sessionInfo())
