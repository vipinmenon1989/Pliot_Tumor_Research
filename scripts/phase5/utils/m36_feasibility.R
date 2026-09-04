#!/usr/bin/env Rscript
# =============================================================================
# Phase 5 - M36 - reconstruction and feasibility.
#
#   * verifies the frozen Phase 4 object by md5 AND sha256 before opening it
#   * asserts every frozen Phase 4 count
#   * VERIFIES (does not assume) which frozen field carries the pre-Phase-4
#     biological annotation
#   * audits assays, layers, cell order, patient balance, depth and confidence
#   * exports the small, durable Phase 5 working artefacts so the 6 GB object
#     never has to be re-opened by a downstream milestone
#
# The Phase 4 object is opened READ-ONLY. Nothing is written back to it.
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject) })
source("scripts/phase5/utils/phase5_common.R")
set.seed(42)
facts <- list(milestone = "M36",
              generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))

# ---- 1. verify the frozen object BEFORE loading it --------------------------
p5_sec("1. Phase 4 object integrity")
stopifnot(file.exists(P4_OBJECT))
sz  <- file.info(P4_OBJECT)$size
md5 <- p5_md5(P4_OBJECT); sha <- p5_sha256(P4_OBJECT)
p5_msg("size   %s bytes", format(sz, big.mark = ","))
p5_msg("md5    %s  (expected %s)", md5, P4_MD5)
p5_msg("sha256 %s", sha)
p5_msg("       %s (expected)", P4_SHA256)
stopifnot(md5 == P4_MD5, sha == P4_SHA256)
p5_msg("  [OK] Phase 4 object verified by md5 AND sha256")
facts$phase4_object <- list(path = P4_OBJECT, size_bytes = sz, md5 = md5,
                            sha256 = sha, verified = TRUE)

# ---- 2. load, read-only -----------------------------------------------------
p5_sec("2. Load")
t0 <- Sys.time()
obj <- readRDS(P4_OBJECT)
p5_msg("loaded in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
p5_msg("cells %d | assays %s | default %s", ncol(obj),
       paste(Assays(obj), collapse = ", "), DefaultAssay(obj))
stopifnot(ncol(obj) == P4_FROZEN$total_cells)
md <- obj@meta.data
md$cell_id <- rownames(md)
facts$object <- list(cells = ncol(obj), assays = Assays(obj),
                     default_assay = DefaultAssay(obj),
                     reductions = Reductions(obj),
                     meta_columns = ncol(md) - 1L,
                     features_RNA = nrow(obj[["RNA"]]),
                     features_SCT = tryCatch(nrow(obj[["SCT"]]), error = function(e) NA))

# ---- 3. assert every frozen Phase 4 count -----------------------------------
p5_sec("3. Frozen Phase 4 counts")
chk <- list(); ok <- TRUE
add <- function(nm, got, want) {
  good <- identical(as.integer(got), as.integer(want))
  chk[[nm]] <<- list(check = nm, expected = as.integer(want),
                     observed = as.integer(got), pass = good)
  p5_msg("  [%s] %-40s expected %-6s observed %s",
         if (good) "OK  " else "FAIL", nm, want, got)
  ok <<- ok && good
}
add("total_cells", ncol(obj), P4_FROZEN$total_cells)
ref <- table(md$malignancy_refined)
for (k in names(P4_FROZEN$refined)) add(paste0("refined_", k), ref[[k]],
                                        P4_FROZEN$refined[[k]])

# ---- 4. VERIFY the historical annotation field ------------------------------
# The pre-Phase-4 biological annotation must be identified, not assumed. It is
# the field that reproduces the frozen Fibroblast 5,064 -> 4,036/908/120 split.
# annotation_ccc_refined is explicitly NOT eligible: it already encodes the
# Phase 4 conclusion, so using it would make the M39 comparison circular.
p5_sec("4. Historical annotation field - verified, not assumed")
cands <- grep("^annotation", colnames(md), value = TRUE)
p5_msg("candidate annotation columns: %s", paste(cands, collapse = ", "))
audit <- lapply(cands, function(cn) {
  v <- as.character(md[[cn]])
  fib <- v == "Fibroblast"
  t <- table(md$malignancy_refined[fib])
  g <- function(k) if (k %in% names(t)) as.integer(t[[k]]) else 0L
  data.frame(column = cn, n_levels = length(unique(v)),
             n_fibroblast = sum(fib), fib_Malignant = g("Malignant"),
             fib_Non_malignant = g("Non-malignant"), fib_Ambiguous = g("Ambiguous"),
             matches_frozen = sum(fib) == P4_FROZEN$fibroblast[["total"]] &&
               g("Malignant") == P4_FROZEN$fibroblast[["Malignant"]] &&
               g("Non-malignant") == P4_FROZEN$fibroblast[["Non-malignant"]] &&
               g("Ambiguous") == P4_FROZEN$fibroblast[["Ambiguous"]],
             eligible = cn != "annotation_ccc_refined")
})
audit <- do.call(rbind, audit)
print(audit, row.names = FALSE)
p5_tsv(audit, file.path(P5_VAL, "M36_HISTORICAL_ANNOTATION_FIELD_AUDIT.tsv"))
hits <- audit$column[audit$matches_frozen & audit$eligible]
p5_msg("fields reproducing the frozen fibroblast split (and eligible): %s",
       paste(hits, collapse = ", "))
if (!length(hits)) stop("M36: no eligible frozen field reproduces the Phase 4 ",
                        "fibroblast split - STOP and diagnose")
# prefer annotation_ccc_phase3 when it qualifies, but only because it qualified
HIST <- if ("annotation_ccc_phase3" %in% hits) "annotation_ccc_phase3" else hits[1]
p5_msg("  [OK] historical annotation field = %s%s", HIST,
       if (HIST == "annotation_ccc_phase3") " (the expected field, verified)" else
         " (NOT the expected field - documented)")
facts$historical_annotation_field <- list(
  selected = HIST, expected = "annotation_ccc_phase3",
  matched_expectation = HIST == "annotation_ccc_phase3",
  eligible_matches = hits,
  excluded_by_design = "annotation_ccc_refined - encodes the Phase 4 conclusion, so using it would make the M39 fibroblast comparison circular",
  audit_table = file.path(P5_VAL, "M36_HISTORICAL_ANNOTATION_FIELD_AUDIT.tsv"))
md$hist_annotation <- as.character(md[[HIST]])

fibt <- table(md$malignancy_refined[md$hist_annotation == "Fibroblast"])
add("fibroblast_total", sum(md$hist_annotation == "Fibroblast"),
    P4_FROZEN$fibroblast[["total"]])
for (k in c("Malignant", "Non-malignant", "Ambiguous"))
  add(paste0("fibroblast_", k), fibt[[k]], P4_FROZEN$fibroblast[[k]])
for (s in SAMPLES)
  add(paste0("clones_", s),
      length(unique(na.omit(md$scevan_clone[md$sample_id == s]))),
      P4_FROZEN$n_clones[[s]])
if (!ok) stop("M36: a frozen Phase 4 count did not reproduce - STOP")
facts$frozen_checks <- unname(chk)
p5_msg("  all %d frozen checks passed", length(chk))

# ---- 5. required fields present and populated -------------------------------
p5_sec("5. Required Phase 4 fields")
REQ <- c("malignancy_refined", "malignancy_confidence", "malignancy_reason",
         HIST, "annotation_ccc_refined", "scevan_clone", "tumor_clone_phase4",
         "tumor_state_phase4", "sample_id", "cnv_burden", "cnv_frac_gain",
         "cnv_frac_loss", "cnv_mean_abs", "scevan_sample_reliable",
         "nCount_RNA", "nFeature_RNA")
fld <- do.call(rbind, lapply(REQ, function(f) data.frame(
  field = f, present = f %in% colnames(md),
  class = if (f %in% colnames(md)) class(md[[f]])[1] else NA_character_,
  n_non_na = if (f %in% colnames(md)) sum(!is.na(md[[f]])) else NA_integer_,
  n_unique = if (f %in% colnames(md)) length(unique(md[[f]])) else NA_integer_)))
print(fld, row.names = FALSE)
p5_tsv(fld, file.path(P5_VAL, "M36_REQUIRED_FIELD_AUDIT.tsv"))
if (any(!fld$present)) stop("M36: missing required field(s): ",
                            paste(fld$field[!fld$present], collapse = ", "))
facts$required_fields <- fld
# tumor_clone_phase4 exists only on malignant cells by Phase 4 construction
p5_msg("tumor_clone_phase4 non-NA: %d (malignant cells: %d)",
       sum(!is.na(md$tumor_clone_phase4)), sum(md$malignancy_refined == "Malignant"))
p5_msg("tumor_state_phase4 non-NA: %d", sum(!is.na(md$tumor_state_phase4)))

# ---- 6. assay / layer / cell-order audit ------------------------------------
p5_sec("6. Assays, layers, cell order")
rna <- obj[["RNA"]]
lay <- Layers(rna)
p5_msg("RNA layers: %s", paste(lay, collapse = ", "))
p5_msg("RNA features: %d", nrow(rna))
cnt_lay <- grep("^counts", lay, value = TRUE)
dat_lay <- grep("^data",   lay, value = TRUE)
p5_msg("counts layers: %d | data layers: %d", length(cnt_lay), length(dat_lay))
stopifnot(length(cnt_lay) >= 1)
same_order <- identical(colnames(obj), rownames(obj@meta.data))
p5_msg("colnames(obj) == rownames(meta.data): %s", same_order)
stopifnot(same_order)
facts$rna_assay <- list(layers = lay, counts_layers = cnt_lay,
                        data_layers = dat_lay, features = nrow(rna),
                        cell_order_consistent = same_order)

# ---- 7. patient distribution audit ------------------------------------------
p5_sec("7. Patient distribution of the malignant compartment")
mal <- md[md$malignancy_refined == "Malignant", ]
pat <- mal |> group_by(sample_id) |>
  summarise(n_malignant = n(),
            frac_of_malignant = n() / nrow(mal),
            median_nCount = median(nCount_RNA),
            median_nFeature = median(nFeature_RNA),
            conf_High = sum(malignancy_confidence == "High"),
            conf_Moderate = sum(malignancy_confidence == "Moderate"),
            conf_Low = sum(malignancy_confidence == "Low"),
            n_clones = length(unique(na.omit(scevan_clone))),
            clone_reliable = unique(as.character(scevan_sample_reliable))[1],
            .groups = "drop") |>
  arrange(desc(n_malignant))
print(as.data.frame(pat), row.names = FALSE, digits = 4)
p5_tsv(pat, file.path(P5_TAB, "M36_MALIGNANT_PATIENT_DISTRIBUTION.tsv"))
dom <- max(pat$frac_of_malignant)
p5_msg("dominant-patient fraction of the malignant compartment: %.4f (%s)",
       dom, pat$sample_id[which.max(pat$frac_of_malignant)])
p5_msg("smallest malignant patient: %s with %d cells",
       pat$sample_id[which.min(pat$n_malignant)], min(pat$n_malignant))
facts$patient_distribution <- pat
facts$patient_imbalance <- list(
  dominant_patient = as.character(pat$sample_id[which.max(pat$frac_of_malignant)]),
  dominant_fraction = dom,
  smallest_patient = as.character(pat$sample_id[which.min(pat$n_malignant)]),
  smallest_n = min(pat$n_malignant),
  imbalance_ratio = max(pat$n_malignant) / max(1, min(pat$n_malignant)),
  note = "patient imbalance is severe enough that a patient-balanced sensitivity analysis is mandatory, not optional")

# confidence composition overall
p5_msg("malignancy_confidence over malignant cells:")
print(table(mal$malignancy_confidence))
facts$malignant_confidence <- as.list(table(mal$malignancy_confidence))

# Phase 4 discrete states, for the M38 comparison
p5_msg("Phase 4 tumour states over malignant cells:")
st <- as.data.frame(table(mal$tumor_state_phase4, mal$sample_id))
names(st) <- c("tumor_state_phase4", "sample_id", "n")
p5_tsv(st |> pivot_wider(names_from = sample_id, values_from = n),
       file.path(P5_TAB, "M36_PHASE4_STATE_BY_PATIENT.tsv"))
print(table(mal$tumor_state_phase4))

# ---- 8. export the durable Phase 5 working artefacts ------------------------
# Downstream milestones must not have to re-open a 6 GB object.
p5_sec("8. Export Phase 5 working artefacts")
dir.create(P5_PROG, showWarnings = FALSE, recursive = TRUE)

# 8a. full 19,716-cell metadata (small) - carries every Phase 1-4 field
saveRDS(md, file.path(P5_PROG, "phase5_cell_metadata.rds"))
p5_msg("  [rds] %s (%d x %d)", file.path(P5_PROG, "phase5_cell_metadata.rds"),
       nrow(md), ncol(md))

# 8b. working cell set = malignant UNION historical fibroblast.
#     Malignant drives program discovery; the historical fibroblasts (including
#     the 908 Non-malignant and 120 Ambiguous) are needed for M39/M40.
keep <- md$cell_id[md$malignancy_refined == "Malignant" |
                     md$hist_annotation == "Fibroblast"]
p5_msg("working cells: %d (malignant %d, historical fibroblast %d, overlap %d)",
       length(keep), sum(md$malignancy_refined == "Malignant"),
       sum(md$hist_annotation == "Fibroblast"),
       sum(md$malignancy_refined == "Malignant" & md$hist_annotation == "Fibroblast"))

# 8c. counts and the frozen normalized data, joined across the per-sample layers
rna_j <- JoinLayers(rna)
cnts <- LayerData(rna_j, layer = "counts")[, keep, drop = FALSE]
p5_msg("counts matrix: %d genes x %d cells, %.1f%% non-zero",
       nrow(cnts), ncol(cnts), 100 * length(cnts@x) / (as.numeric(nrow(cnts)) * ncol(cnts)))
stopifnot(all(cnts@x == floor(cnts@x)))    # integers, as Phase 4 asserted
saveRDS(cnts, file.path(P5_PROG, "phase5_working_counts.rds"))
p5_msg("  [rds] phase5_working_counts.rds")

dat <- LayerData(rna_j, layer = "data")[, keep, drop = FALSE]
saveRDS(dat, file.path(P5_PROG, "phase5_working_lognorm.rds"))
p5_msg("  [rds] phase5_working_lognorm.rds (frozen Phase 2/3/4 RNA data layer)")
rm(rna_j); invisible(gc())

# 8d. per-cell gene detection over the malignant compartment, for gene filtering
malc <- cnts[, md$cell_id[md$malignancy_refined == "Malignant"], drop = FALSE]
gd <- data.frame(gene = rownames(malc),
                 n_cells_detected = Matrix::rowSums(malc > 0),
                 total_counts = Matrix::rowSums(malc))
gd$frac_cells <- gd$n_cells_detected / ncol(malc)
p5_msg("genes detected in >=1 malignant cell: %d / %d",
       sum(gd$n_cells_detected > 0), nrow(gd))
for (th in c(0.001, 0.005, 0.01, 0.02, 0.05))
  p5_msg("  genes in >= %.1f%% of malignant cells: %d", 100 * th,
         sum(gd$frac_cells >= th))
saveRDS(gd, file.path(P5_PROG, "phase5_malignant_gene_detection.rds"))
facts$gene_detection <- list(
  genes_total = nrow(gd), genes_detected = sum(gd$n_cells_detected > 0),
  genes_ge_0.5pct = sum(gd$frac_cells >= 0.005),
  genes_ge_1pct = sum(gd$frac_cells >= 0.01))

# 8e. cell-cycle gene list actually present, for the declared sensitivity run
cc <- intersect(p5_cc_genes(), rownames(cnts))
p5_msg("canonical cell-cycle genes present: %d", length(cc))
writeLines(cc, file.path(P5_PROG, "phase5_cellcycle_genes.txt"))
facts$cellcycle_genes_present <- length(cc)

facts$exports <- list(
  cell_metadata = file.path(P5_PROG, "phase5_cell_metadata.rds"),
  working_counts = file.path(P5_PROG, "phase5_working_counts.rds"),
  working_lognorm = file.path(P5_PROG, "phase5_working_lognorm.rds"),
  gene_detection = file.path(P5_PROG, "phase5_malignant_gene_detection.rds"),
  cellcycle_genes = file.path(P5_PROG, "phase5_cellcycle_genes.txt"),
  working_cells = length(keep),
  malignant_cells = sum(md$malignancy_refined == "Malignant"),
  historical_fibroblasts = sum(md$hist_annotation == "Fibroblast"))

# ---- 9. assert nothing was written back -------------------------------------
p5_sec("9. Phase 4 object unchanged")
md5b <- p5_md5(P4_OBJECT)
p5_msg("md5 after M36: %s (%s)", md5b,
       if (md5b == P4_MD5) "UNCHANGED" else "CHANGED - STOP")
stopifnot(md5b == P4_MD5)
facts$phase4_object_unchanged_after_m36 <- TRUE

p5_json(facts, file.path(P5_VAL, "m36_feasibility_facts.json"))
p5_sec("M36 complete")
