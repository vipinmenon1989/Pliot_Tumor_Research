#!/usr/bin/env Rscript
# =============================================================================
# Phase 6 - M50 - build, validate and freeze the Phase 6 object.
#
# Every Phase 1-5 field is preserved and asserted column by column against a
# copy taken before anything was added.
#
# What is added, and what deliberately is NOT: the full TF-activity matrix
# (hundreds of regulons) and the full Hallmark matrix stay as separate RDS
# artefacts referenced by the manifest. Pushing hundreds of columns into the
# object metadata would bloat it for no analytical gain. What goes IN is the
# compact, interpretable layer: the 14 PROGENy pathway activities, the top
# cross-patient-concordant TF activities, and the clone/program metrics.
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject) })
source("scripts/phase6/utils/phase6_common.R")
set.seed(42)
OUT <- "results/phase6/phase6_final_object.rds"
MAXTF <- 25L
facts <- list(milestone = "M50", generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))

p6_sec("1. Verify and load the Phase 5 object")
mf5 <- fromJSON("results/phase5/phase5_manifest.json", simplifyVector = FALSE)
stopifnot(p6_md5(P5_OBJECT) == mf5$phase5_object$md5,
          p6_sha256(P5_OBJECT) == mf5$phase5_object$sha256)
p6_msg("Phase 5 object verified: md5 %s", mf5$phase5_object$md5)
obj <- readRDS(P5_OBJECT)
md_before <- obj@meta.data
red_before <- lapply(Reductions(obj), function(r) Embeddings(obj, r))
names(red_before) <- Reductions(obj)
cells_before <- colnames(obj)
p6_msg("cells %d | metadata columns %d", ncol(obj), ncol(md_before))

p6_sec("2. Assemble the compact Phase 6 layer")
A <- readRDS(file.path(P6_REG, "m45_tf_activity_matrix.rds"))     # cells x TFs
P <- readRDS(file.path(P6_PATH, "m46_progeny_activity_matrix.rds"))
tf <- read.delim(file.path(P6_TAB, "PROGRAM_TF_ACTIVITY.tsv"))
top_tf <- tf |> filter(recurrent_direction) |>
  group_by(tf) |> summarise(m = max(abs(pooled_rho)), .groups = "drop") |>
  arrange(desc(m)) |> head(MAXTF) |> pull(tf)
if (!length(top_tf)) {
  top_tf <- tf |> group_by(tf) |> summarise(m = max(abs(pooled_rho)), .groups = "drop") |>
    arrange(desc(m)) |> head(MAXTF) |> pull(tf)
  p6_msg("no TF was cross-patient concordant; carrying the %d strongest pooled TFs instead, labelled as such", length(top_tf))
} else p6_msg("carrying %d cross-patient-concordant TF activities into the object", length(top_tf))
p6_msg("TFs: %s", paste(top_tf, collapse = ", "))

div <- read.delim(file.path(P6_TAB, "CLONE_PROGRAM_DIVERSITY.tsv"))
assoc <- read.delim(file.path(P6_TAB, "CLONE_PROGRAM_ASSOCIATION.tsv"))
nm <- length(cells_before)
blank <- function() setNames(rep(NA_real_, nm), cells_before)
newcols <- list()
for (p in colnames(P)) {
  v <- blank(); v[rownames(P)] <- P[, p]
  newcols[[paste0("progeny_", p)]] <- v
}
for (t in top_tf) {
  v <- blank(); v[rownames(A)] <- A[, t]
  newcols[[paste0("tfact_", make.names(t))]] <- v
}
# clone-level metrics broadcast to the cells of that clone
key <- paste(md_before$sample_id, md_before$tumor_clone_phase4, sep = "|")
dkey <- paste(div$sample_id, div$clone, sep = "|")
for (m in c("hard_effective_n", "continuous_effective_n",
            "dominant_program_share_hard", "program_dispersion")) {
  v <- blank(); i <- match(key, dkey)
  v[!is.na(i)] <- div[[m]][i[!is.na(i)]]
  newcols[[paste0("clone_", m)]] <- v
}
v <- rep(NA_character_, nm); names(v) <- cells_before
i <- match(key, dkey); v[!is.na(i)] <- div$dominant_program_label[i[!is.na(i)]]
newcols$clone_dominant_program_label <- v
v <- rep(NA, nm); names(v) <- cells_before
v[!is.na(i)] <- div$n_cells[i[!is.na(i)]] >= CLONE_SIZE_PRIMARY
newcols$clone_evaluable_for_diversity <- v
v <- rep(NA, nm); names(v) <- cells_before
v[] <- md_before$sample_id %in% CLONE_RELIABLE
newcols$clone_structure_reliable_phase6 <- v

clash <- intersect(names(newcols), colnames(md_before))
if (length(clash)) stop("M50: Phase 6 field name collides with an existing field: ",
                        paste(clash, collapse = ", "))
for (n in names(newcols)) obj[[n]] <- unname(newcols[[n]])
p6_msg("added %d Phase 6 metadata fields", length(newcols))
p6_msg("PROGENy %d | TF %d | clone metrics %d", ncol(P), length(top_tf),
       length(newcols) - ncol(P) - length(top_tf))

p6_sec("3. Preservation guards, column by column")
md_after <- obj@meta.data
g <- list()
g$cells_unchanged <- identical(colnames(obj), cells_before)
same <- vapply(colnames(md_before), function(cn)
  identical(md_before[[cn]], md_after[[cn]]), logical(1))
g$all_phase1_5_columns_identical <- all(same)
if (!all(same)) p6_msg("CHANGED: %s", paste(names(same)[!same], collapse = ", "))
g$only_expected_added <- setequal(setdiff(colnames(md_after), colnames(md_before)),
                                  names(newcols))
g$reductions_unchanged <- identical(Reductions(obj), names(red_before)) &&
  all(vapply(names(red_before), function(r)
    identical(Embeddings(obj, r), red_before[[r]]), logical(1)))
g$assays_unchanged <- identical(sort(Assays(obj)), sort(c("RNA", "SCT")))
g$malignancy_refined_unchanged <- identical(md_before$malignancy_refined,
                                            md_after$malignancy_refined)
g$phase5_program_scores_unchanged <- all(vapply(
  grep("^program_P[0-9]+_score$", colnames(md_before), value = TRUE),
  function(cn) identical(md_before[[cn]], md_after[[cn]]), logical(1)))
g$frozen_counts_hold <- identical(
  as.integer(table(md_after$malignancy_refined)[c("Malignant", "Non-malignant",
    "Ambiguous", "Excluded-low-quality")]), c(6434L, 9078L, 3766L, 438L))
for (n in names(g)) p6_msg("  [%s] %s", if (isTRUE(g[[n]])) "OK  " else "FAIL", n)
if (!all(unlist(g))) stop("M50: a preservation guard failed - STOP")
facts$preservation <- g

p6_sec("4. Save, reload, validate")
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
saveRDS(obj, OUT)
sz <- file.info(OUT)$size; md5 <- p6_md5(OUT); sha <- p6_sha256(OUT)
p6_msg("saved %s | %s bytes", OUT, format(sz, big.mark = ","))
p6_msg("md5    %s", md5); p6_msg("sha256 %s", sha)
rm(obj); invisible(gc())
o2 <- readRDS(OUT)
v <- list(
  reload_cells = ncol(o2) == 19716L,
  reload_cell_order = identical(colnames(o2), cells_before),
  reload_assays = identical(sort(Assays(o2)), sort(c("RNA", "SCT"))),
  reload_reductions = identical(Reductions(o2), names(red_before)),
  reload_phase1_5_columns = all(vapply(colnames(md_before), function(cn)
    identical(md_before[[cn]], o2@meta.data[[cn]]), logical(1))),
  reload_phase6_columns = all(names(newcols) %in% colnames(o2@meta.data)),
  reload_progeny_on_malignant = all(!is.na(o2@meta.data[[paste0("progeny_",
    colnames(P)[1])]][o2$malignancy_refined == "Malignant"])),
  reload_md5_stable = p6_md5(OUT) == md5)
for (n in names(v)) p6_msg("  [%s] %s", if (isTRUE(v[[n]])) "OK  " else "FAIL", n)
if (!all(unlist(v))) stop("M50: a reload validation failed - STOP")
facts$validation <- v
facts$phase6_object <- list(path = OUT, size_bytes = sz, md5 = md5, sha256 = sha,
                            cells = ncol(o2), metadata_columns = ncol(o2@meta.data),
                            phase6_fields_added = names(newcols),
                            not_embedded = list(
                              tf_activity_matrix = file.path(P6_REG, "m45_tf_activity_matrix.rds"),
                              hallmark_matrix = file.path(P6_PATH, "m46_hallmark_score_matrix.rds"),
                              progeny_matrix = file.path(P6_PATH, "m46_progeny_activity_matrix.rds"),
                              reason = "full activity matrices are kept as separate artefacts rather than bloating the object with hundreds of metadata columns"))

p6_sec("5. Upstream objects untouched")
stopifnot(p6_md5(P4_OBJECT) == P4_MD5,
          p6_md5(P5_OBJECT) == mf5$phase5_object$md5)
p6_msg("Phase 4 md5 unchanged: %s", P4_MD5)
p6_msg("Phase 5 md5 unchanged: %s", mf5$phase5_object$md5)
facts$upstream_unchanged <- TRUE
p6_json(facts, file.path(P6_VAL, "m50_final_object_validation.json"))
p6_sec("M50 object build complete")
