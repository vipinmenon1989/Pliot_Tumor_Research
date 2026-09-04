#!/usr/bin/env Rscript
# =============================================================================
# Phase 5 - M37a - prepare the cNMF expression input.
#
# Declared BEFORE any factorization is run, so no choice below can be made
# after seeing a program:
#   * cells      : malignancy_refined == "Malignant" ONLY (6,434)
#   * expression : RNA raw integer counts (cNMF normalizes internally to TPM).
#                  NOT Harmony, NOT UMAP, NOT PCA scores, NOT SCT residuals,
#                  NOT the SCEVAN CNA matrix.
#   * genes      : detected in >= 0.5% of malignant cells, mitochondrial
#                  (^MT-) genes removed. Ribosomal, cell-cycle, ECM, HLA,
#                  interferon, Schwann and neural-crest genes are all RETAINED
#                  in the primary run - removing them would delete the signal
#                  the phase exists to find.
#   * balanced   : equal-cell downsample capped at the smallest patient with
#                  >= 500 malignant cells; patients below the cap contribute
#                  all their cells. Seed 42.
#   * no-CC      : declared sensitivity gene universe with the 95 canonical
#                  Tirosh S/G2M genes removed.
# =============================================================================
source("scripts/phase5/utils/phase5_common.R")
set.seed(42)
facts <- list(milestone = "M37a", generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))

OUT <- file.path(P5_PROG, "cnmf_input")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

p5_sec("1. Load Phase 5 working artefacts")
md <- readRDS(file.path(P5_PROG, "phase5_cell_metadata.rds"))
cnts <- readRDS(file.path(P5_PROG, "phase5_working_counts.rds"))
p5_msg("metadata %d x %d | counts %d genes x %d cells",
       nrow(md), ncol(md), nrow(cnts), ncol(cnts))
stopifnot(nrow(md) == P4_FROZEN$total_cells)

mal_ids <- md$cell_id[md$malignancy_refined == "Malignant"]
stopifnot(length(mal_ids) == P4_FROZEN$refined[["Malignant"]],
          all(mal_ids %in% colnames(cnts)))
M <- cnts[, mal_ids, drop = FALSE]
p5_msg("malignant matrix: %d genes x %d cells", nrow(M), ncol(M))
mmd <- md[mal_ids, ]

p5_sec("2. Gene universe (declared in advance)")
det <- Matrix::rowSums(M > 0)
frac <- det / ncol(M)
MIN_FRAC <- 0.005
is_mt <- grepl("^MT-", rownames(M))
keep_gene <- frac >= MIN_FRAC & !is_mt
p5_msg("genes total                      %d", nrow(M))
p5_msg("  detected in >= %.1f%% of cells  %d", 100 * MIN_FRAC, sum(frac >= MIN_FRAC))
p5_msg("  mitochondrial (^MT-) removed   %d", sum(is_mt & frac >= MIN_FRAC))
p5_msg("  PRIMARY gene universe          %d", sum(keep_gene))
G <- M[keep_gene, , drop = FALSE]

cc <- readLines(file.path(P5_PROG, "phase5_cellcycle_genes.txt"))
cc_in <- intersect(cc, rownames(G))
p5_msg("canonical cell-cycle genes inside the primary universe: %d", length(cc_in))
p5_msg("  RETAINED in the primary run; removed only in the declared sensitivity run")
# sanity: the biology this phase exists to find must still be in the universe
fams <- list(ECM = c("COL1A1","COL1A2","COL3A1","FN1","POSTN"),
             HLA = c("HLA-A","HLA-B","HLA-E","B2M"),
             Interferon = c("ISG15","IFIT1","MX1","STAT1"),
             Schwann = c("S100B","PLP1","SOX10","MPZ"),
             NeuralCrest = c("SOX9","TWIST1","NES","ERBB3"),
             Angiogenic = c("VEGFA","HIF1A","NDRG1"),
             CellCycle = c("MKI67","TOP2A","CDK1","CCNB1"))
famtab <- do.call(rbind, lapply(names(fams), function(n)
  data.frame(family = n, n_expected = length(fams[[n]]),
             n_present = length(intersect(fams[[n]], rownames(G))),
             missing = paste(setdiff(fams[[n]], rownames(G)), collapse = ","))))
print(famtab, row.names = FALSE)
p5_tsv(famtab, file.path(P5_VAL, "M37_GENE_FAMILY_RETENTION.tsv"))
facts$gene_universe <- list(
  genes_input = nrow(M), min_frac_malignant_cells = MIN_FRAC,
  mitochondrial_removed = sum(is_mt & frac >= MIN_FRAC),
  genes_primary = sum(keep_gene), cellcycle_genes_in_universe = length(cc_in),
  retained_families = famtab)

p5_sec("3. Patient-balanced subset (sensitivity design, declared in advance)")
pt <- table(mmd$sample_id)
print(pt)
MIN_FOR_CAP <- 500L
cap <- min(pt[pt >= MIN_FOR_CAP])
cap_pat <- names(pt)[pt == cap][1]
p5_msg("cap = %d cells (%s: the smallest patient with >= %d malignant cells)",
       cap, cap_pat, MIN_FOR_CAP)
p5_msg("patients below the cap contribute all their malignant cells")
bal <- unlist(lapply(names(pt), function(s) {
  ids <- mmd$cell_id[mmd$sample_id == s]
  if (length(ids) <= cap) ids else sort(sample(ids, cap))
}), use.names = FALSE)
bt <- table(mmd$sample_id[match(bal, mmd$cell_id)])
p5_msg("balanced subset: %d cells", length(bal)); print(bt)
p5_msg("dominant-patient fraction: full %.4f -> balanced %.4f",
       max(pt) / sum(pt), max(bt) / sum(bt))
facts$balanced_design <- list(
  strategy = "equal-cell downsampling per patient, capped at the smallest patient with >= 500 malignant cells; smaller patients contribute all their cells",
  cap = as.integer(cap), cap_patient = cap_pat, seed = 42L,
  n_cells = length(bal),
  per_patient_full = as.list(pt), per_patient_balanced = as.list(bt),
  dominant_fraction_full = unname(max(pt) / sum(pt)),
  dominant_fraction_balanced = unname(max(bt) / sum(bt)))

p5_sec("4. Write matrices for cNMF")
write_run <- function(mat, cells, genes, tag) {
  d <- file.path(OUT, tag); dir.create(d, showWarnings = FALSE, recursive = TRUE)
  sub <- mat[genes, cells, drop = FALSE]
  sub <- sub[Matrix::rowSums(sub) > 0, , drop = FALSE]   # cNMF rejects all-zero genes
  Matrix::writeMM(as(sub, "dgCMatrix"), file.path(d, "counts.mtx"))
  writeLines(rownames(sub), file.path(d, "genes.txt"))
  writeLines(colnames(sub), file.path(d, "cells.txt"))
  p5_msg("  [%s] %d genes x %d cells -> %s", tag, nrow(sub), ncol(sub), d)
  list(tag = tag, genes = nrow(sub), cells = ncol(sub), dir = d)
}
runs <- list(
  primary  = write_run(G, colnames(G), rownames(G), "primary"),
  balanced = write_run(G, bal, rownames(G), "balanced"),
  nocc     = write_run(G, colnames(G), setdiff(rownames(G), cc_in), "nocc"))
facts$runs <- runs

# per-cell covariates travel with the matrices so Python never guesses
cov <- mmd |> select(cell_id, sample_id, malignancy_confidence,
                     tumor_state_phase4, tumor_clone_phase4, scevan_clone,
                     nCount_RNA, nFeature_RNA, percent.mt,
                     cnv_burden, cnv_frac_gain, cnv_frac_loss, cnv_mean_abs) |>
  mutate(hist_annotation = mmd$annotation_ccc_phase3)
p5_tsv(cov, file.path(OUT, "malignant_cell_covariates.tsv"))
writeLines(bal, file.path(OUT, "balanced_cells.txt"))

p5_json(facts, file.path(P5_VAL, "m37a_prepare_facts.json"))
p5_sec("M37a complete")
