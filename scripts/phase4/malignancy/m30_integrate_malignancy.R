#!/usr/bin/env Rscript
# =============================================================================
# Phase 4 · M30 — malignancy integration
#
# Applies the rules in reports/phase4/MALIGNANCY_DECISION_RULES.md, which were
# written and committed BEFORE any SCEVAN result was inspected. Nothing here is
# tuned to the answer.
#
# SCEVAN carries substantial weight but not absolute authority (§23). Every call
# records the rule that produced it, and disagreement between Phase 2 and SCEVAN
# is preserved rather than resolved by fiat (§25).
# =============================================================================
suppressPackageStartupMessages({
  library(Matrix); library(jsonlite); library(Seurat); library(SeuratObject)
})
source("scripts/R/provenance_utils.R")
set.seed(42)

SC   <- "results/phase4/scevan/by_sample"
MAL  <- "results/phase4/malignancy"
TAB  <- "results/phase4/tables"
SAMPLES <- c("MPNST_1","MPNST_2","MPNST_3","MPNST_4")
dir.create(TAB, showWarnings = FALSE, recursive = TRUE)

log_ <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
sec  <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

# --- a priori thresholds, restated from the rules document ------------------
THR <- list(
  pop_frac_high      = 0.50,   # "majority of the stratum"
  pop_frac_low       = 0.25,   # "a quarter of the stratum"
  pop_min_assessed   = 20L,    # below this, pop_frac is unreliable
  cna_abs            = 0.10,   # |relative CNA| counted as altered
  cna_normal_pctile  = 0.95)   # cnv_elevated reference among confident normals
sec("A priori thresholds (from MALIGNANCY_DECISION_RULES.md §2)")
str(THR)

# --- population classes ----------------------------------------------------
CLS <- list(
  EXCLUDED         = "Low-quality-excluded",
  PHASE2_MALIGNANT = "MPNST-Tumor",
  CANONICAL_IMMUNE = c("CD4-T","CD8-T","NK","T-cell-other","B-cell","Plasma-cell",
                       "Macrophage","Monocyte","Dendritic","Plasmacytoid-DC"),
  CANONICAL_ENDO   = "Endothelial",
  DISPUTED_STROMAL = c("Fibroblast","Pericyte-VSMC"),
  PROVISIONAL      = c("Candidate-Malignant-Unresolved","Uncertain"))
pop_class <- function(a) {
  out <- rep(NA_character_, length(a))
  for (k in names(CLS)) out[a %in% CLS[[k]]] <- k
  out
}

# --- 1. metadata -----------------------------------------------------------
sec("1. Phase 2 / Phase 3 metadata (frozen, read-only)")
md <- read.delim(file.path(MAL, "phase4_cell_metadata.tsv.gz"), stringsAsFactors = FALSE)
rownames(md) <- md$cell_id
log_("cells: ", nrow(md), " | annotation_ccc levels: ", length(unique(md$annotation_ccc)))
md$population_class <- pop_class(md$annotation_ccc)
stopifnot(!any(is.na(md$population_class)))
log_("population_class assigned to every cell:")
print(table(md$population_class))

# --- 2. SCEVAN classifications --------------------------------------------
sec("2. SCEVAN classifications (primary + sensitivity)")
read_cls <- function(s, tag) {
  f <- file.path(SC, s, tag, sprintf("scevan_classification_%s_%s.rds", s, tag))
  if (!file.exists(f)) stop("missing SCEVAN result: ", f)
  x <- readRDS(f)
  data.frame(cell_id = x$cell_id, call = x$scevan_call,
             confident_normal = x$scevan_confident_normal,
             subclone = x$scevan_subclone, stringsAsFactors = FALSE)
}
prim <- do.call(rbind, lapply(SAMPLES, read_cls, tag = "primary"))
sens <- do.call(rbind, lapply(SAMPLES, read_cls, tag = "sensitivity"))
log_("primary rows: ", nrow(prim), " | sensitivity rows: ", nrow(sens))
log_("primary call table:");     print(table(prim$call))
log_("sensitivity call table:"); print(table(sens$call))

# SCEVAN drops cells during preprocessing; those become "not-assessed".
md$scevan_call_raw <- prim$call[match(md$cell_id, prim$cell_id)]
md$scevan_call_raw[is.na(md$scevan_call_raw)] <- "not-assessed"
md$scevan_call_raw[md$scevan_call_raw == "filtered"] <- "not-assessed"
md$scevan_sens_raw <- sens$call[match(md$cell_id, sens$cell_id)]
md$scevan_sens_raw[is.na(md$scevan_sens_raw)] <- "not-assessed"
md$scevan_sens_raw[md$scevan_sens_raw == "filtered"] <- "not-assessed"
md$scevan_confident_normal <- prim$confident_normal[match(md$cell_id, prim$cell_id)]
md$scevan_confident_normal[is.na(md$scevan_confident_normal)] <- "no"
md$scevan_subclone_raw <- prim$subclone[match(md$cell_id, prim$cell_id)]
md$scevan_sample <- md$sample_id

norm3 <- function(x) c(tumor = "malignant", normal = "non-malignant",
                       `not-assessed` = "not-assessed")[x]
md$malignancy_scevan      <- unname(norm3(md$scevan_call_raw))
md$malignancy_scevan_sens <- unname(norm3(md$scevan_sens_raw))
md$scevan_call <- md$malignancy_scevan
md$runs_agree  <- md$malignancy_scevan == md$malignancy_scevan_sens
log_("malignancy_scevan:");      print(table(md$malignancy_scevan))
log_("primary vs sensitivity agreement (assessed cells only): ",
     sprintf("%.4f", mean(md$runs_agree[md$malignancy_scevan != "not-assessed" &
                                        md$malignancy_scevan_sens != "not-assessed"])))
md$scevan_clone <- ifelse(is.na(md$scevan_subclone_raw), NA_character_,
                          paste0(md$sample_id, "_clone", md$scevan_subclone_raw))
log_("clone labels: ", length(unique(na.omit(md$scevan_clone))))

# --- 3. malignancy_phase2 (derived from the frozen annotation) -------------
sec("3. malignancy_phase2")
md$malignancy_phase2 <- with(md, ifelse(
  population_class == "PHASE2_MALIGNANT", "malignant",
  ifelse(population_class == "EXCLUDED", "excluded",
  ifelse(population_class == "PROVISIONAL", "ambiguous", "non-malignant"))))
print(table(md$malignancy_phase2))

# --- 3b. PER-SAMPLE IMMUNE SANITY CHECK -----------------------------------
# ADDED AFTER M29, in response to an observed failure, and documented as such in
# MALIGNANCY_DECISION_RULES.md Amendment A1. Not an a priori rule.
#
# WHY: in MPNST_3 the primary run called 89/89 B cells, 327/328 CD4-T, 135/136
# CD8-T, 32/32 NK, 415/418 plasma cells and 165/165 pDC "tumor" while calling
# macrophages "normal"; the sensitivity run inverted it (macrophages/DC/monocytes
# "tumor", lymphoid "normal"). Primary/sensitivity agreement was 0.0864. SCEVAN
# found only 25 confident normal cells in that sample and MPNST_3 has the lowest
# depth of the four (median 1,594 genes/cell). With a poorly estimated baseline,
# classifyCluster() partitions the top-level CNA clustering along a lymphoid /
# myeloid axis rather than a malignant / normal one.
#
# A run that calls a large share of unambiguous immune cells malignant has failed
# a basic sanity check, and must not be allowed to drive PROMOTIONS of disputed
# stromal or provisional cells to Malignant.
sec("3b. Per-sample SCEVAN immune sanity check")
IMMUNE_SANITY_MAX  <- 0.25   # >25% of a sample's canonical immune called malignant
IMMUNE_POPS_MAX    <- 3L     # ... or >=3 immune populations individually over that rate
SANITY_MIN_CELLS   <- 20L    # a population needs this many assessed cells to be judged
# Amendment A2: Plasma-cell is excluded from the sanity DENOMINATOR on mechanistic
# grounds. Plasma cells express immunoglobulin loci (IGH 14q32, IGK 2p11, IGL 22q11)
# at extraordinary levels, which expression-derived CNV inference reads as large
# segmental gains. This is a documented artefact class for inferCNV/CopyKAT/SCEVAN-
# style methods, it is not specific to this dataset, and plasma cells are not part
# of the malignancy question. Both rates are computed and reported.
SANITY_ARTEFACT_POPS <- "Plasma-cell"

per_pop_fail <- function(d) {
  pops <- setdiff(unique(d$annotation_ccc), SANITY_ARTEFACT_POPS)
  n_tested <- 0L; n_fail <- 0L; failing <- character(0)
  for (p in pops) {
    dp <- d[d$annotation_ccc == p, ]
    a <- sum(dp$malignancy_scevan != "not-assessed")
    if (a < SANITY_MIN_CELLS) next
    n_tested <- n_tested + 1L
    r <- sum(dp$malignancy_scevan == "malignant") / a
    if (r > IMMUNE_SANITY_MAX) { n_fail <- n_fail + 1L
      failing <- c(failing, sprintf("%s=%.2f", p, r)) }
  }
  list(n_tested = n_tested, n_fail = n_fail, failing = paste(failing, collapse = "; "))
}

sanity <- do.call(rbind, lapply(SAMPLES, function(s) {
  d  <- md[md$sample_id == s & md$population_class == "CANONICAL_IMMUNE", ]
  d2 <- d[!d$annotation_ccc %in% SANITY_ARTEFACT_POPS, ]
  a_all <- sum(d$malignancy_scevan  != "not-assessed")
  m_all <- sum(d$malignancy_scevan  == "malignant")
  a_ex  <- sum(d2$malignancy_scevan != "not-assessed")
  m_ex  <- sum(d2$malignancy_scevan == "malignant")
  a_sen <- sum(d2$malignancy_scevan_sens != "not-assessed")
  m_sen <- sum(d2$malignancy_scevan_sens == "malignant")
  pf <- per_pop_fail(d)
  dp <- d[d$annotation_ccc %in% SANITY_ARTEFACT_POPS, ]
  data.frame(sample_id = s,
    immune_cells = nrow(d),
    immune_assessed_all = a_all, immune_malignant_all = m_all,
    rate_all_immune = if (a_all) m_all/a_all else NA_real_,
    plasma_assessed = sum(dp$malignancy_scevan != "not-assessed"),
    plasma_malignant = sum(dp$malignancy_scevan == "malignant"),
    immune_assessed_excl_plasma = a_ex, immune_malignant_excl_plasma = m_ex,
    rate_excl_plasma = if (a_ex) m_ex/a_ex else NA_real_,
    rate_excl_plasma_sensitivity = if (a_sen) m_sen/a_sen else NA_real_,
    n_immune_pops_tested = pf$n_tested,
    n_immune_pops_failing = pf$n_fail,
    failing_populations = pf$failing,
    stringsAsFactors = FALSE) }))

sanity$fail_by_rate   <- !is.na(sanity$rate_excl_plasma) & sanity$rate_excl_plasma > IMMUNE_SANITY_MAX
sanity$fail_by_breadth <- sanity$n_immune_pops_failing >= IMMUNE_POPS_MAX
sanity$scevan_reliable <- !(sanity$fail_by_rate | sanity$fail_by_breadth)
sanity$flat_gate_would_have_failed <- !is.na(sanity$rate_all_immune) &
                                      sanity$rate_all_immune > IMMUNE_SANITY_MAX
sanity$verdict <- ifelse(sanity$scevan_reliable, "PASS",
                  ifelse(sanity$fail_by_rate & sanity$fail_by_breadth,
                         "FAIL (rate AND breadth) - promotions disabled",
                  ifelse(sanity$fail_by_rate, "FAIL (rate) - promotions disabled",
                         "FAIL (breadth) - promotions disabled")))
print(sanity[, c("sample_id","rate_all_immune","plasma_malignant","plasma_assessed",
                 "rate_excl_plasma","rate_excl_plasma_sensitivity",
                 "n_immune_pops_failing","n_immune_pops_tested",
                 "flat_gate_would_have_failed","scevan_reliable","verdict")],
      row.names = FALSE, digits = 4)
for (i in seq_len(nrow(sanity))) if (nzchar(sanity$failing_populations[i]))
  log_("  ", sanity$sample_id[i], " failing populations: ", sanity$failing_populations[i])
write.table(sanity, file.path(TAB, "SCEVAN_SAMPLE_RELIABILITY.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
md$scevan_sample_reliable <- sanity$scevan_reliable[match(md$sample_id, sanity$sample_id)]
for (i in seq_len(nrow(sanity))) if (!sanity$scevan_reliable[i])
  log_("*** ", sanity$sample_id[i], " FAILS the immune sanity check (",
       sprintf("%.1f%%", 100*sanity$rate_primary[i]),
       " of canonical immune cells called malignant). Promotions from this sample are disabled. ***")
log_("samples passing: ", sum(sanity$scevan_reliable), " of ", nrow(sanity))

# --- 4. CNA burden ---------------------------------------------------------
sec("4. CNA burden (E4) and cnv_elevated")
md$cnv_burden <- NA_real_; md$cnv_mean_abs <- NA_real_
md$cnv_frac_gain <- NA_real_; md$cnv_frac_loss <- NA_real_
cna_ref <- list()
for (s in SAMPLES) {
  f <- file.path(SC, s, "primary", "output", paste0(s, "_primary_CNAmtx.RData"))
  if (!file.exists(f)) { log_("WARNING: no CNA matrix for ", s, " (", f, ")"); next }
  e <- new.env(); load(f, envir = e)
  nm <- ls(e)[vapply(ls(e), function(n) is.matrix(get(n, e)), logical(1))][1]
  M <- get(nm, e)
  log_(s, ": CNA matrix '", nm, "' ", nrow(M), " genes x ", ncol(M), " cells")
  # Per cell: fraction of assessed genes whose |relative CNA| exceeds THR$cna_abs.
  # CNA values are piecewise-constant along segments, so this is a gene-density-
  # weighted genome fraction, not a base-pair fraction. Stated, not implied.
  A <- abs(M) > THR$cna_abs
  burden <- colMeans(A)
  gain <- colMeans(M >  THR$cna_abs)
  loss <- colMeans(M < -THR$cna_abs)
  mabs <- colMeans(abs(M))
  ii <- match(colnames(M), md$cell_id)
  md$cnv_burden[ii]    <- burden
  md$cnv_frac_gain[ii] <- gain
  md$cnv_frac_loss[ii] <- loss
  md$cnv_mean_abs[ii]  <- mabs
  cn <- md$cell_id[md$scevan_confident_normal == "yes" & md$sample_id == s]
  cn <- intersect(cn, colnames(M))
  ref <- if (length(cn) >= 10) as.numeric(quantile(burden[cn], THR$cna_normal_pctile))
         else NA_real_
  cna_ref[[s]] <- list(confident_normals = length(cn),
                       threshold_p95 = ref,
                       median_burden_all = as.numeric(median(burden)))
  log_(sprintf("  confident normals %d | p%.0f burden threshold %s | median burden all %.4f",
       length(cn), 100*THR$cna_normal_pctile,
       if (is.na(ref)) "NA (too few normals)" else sprintf("%.4f", ref), median(burden)))
  rm(M, A, e); invisible(gc(FALSE))
}
md$cnv_ref_threshold <- vapply(md$sample_id, function(s)
  if (is.null(cna_ref[[s]])) NA_real_ else cna_ref[[s]]$threshold_p95, numeric(1))
md$cnv_elevated <- !is.na(md$cnv_burden) & !is.na(md$cnv_ref_threshold) &
                   md$cnv_burden > md$cnv_ref_threshold
log_("cnv_elevated cells: ", sum(md$cnv_elevated), " / ", sum(!is.na(md$cnv_burden)),
     " assessed")

# --- 5. pop_frac (E5) ------------------------------------------------------
sec("5. pop_frac: SCEVAN malignant fraction per (annotation_ccc x sample_id)")
md$.assessed  <- md$malignancy_scevan != "not-assessed"
md$.malignant <- md$malignancy_scevan == "malignant"
agg <- aggregate(cbind(n_assessed = .assessed, n_malignant = .malignant) ~
                 annotation_ccc + sample_id, data = md, FUN = sum)
agg$pop_frac <- ifelse(agg$n_assessed > 0, agg$n_malignant / agg$n_assessed, NA_real_)
agg$pop_frac_reliable <- agg$n_assessed >= THR$pop_min_assessed
key <- paste(md$annotation_ccc, md$sample_id)
akey <- paste(agg$annotation_ccc, agg$sample_id)
md$pop_frac <- agg$pop_frac[match(key, akey)]
md$pop_n_assessed <- agg$n_assessed[match(key, akey)]
md$pop_frac_reliable <- agg$pop_frac_reliable[match(key, akey)]
md$pop_frac_reliable[is.na(md$pop_frac_reliable)] <- FALSE
print(agg[order(-agg$pop_frac), c("annotation_ccc","sample_id","n_assessed",
                                  "n_malignant","pop_frac","pop_frac_reliable")],
      row.names = FALSE)
write.table(agg, file.path(TAB, "SCEVAN_POPULATION_MALIGNANT_FRACTION.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- 6. lineage-marker module scores (E6) ---------------------------------
sec("6. Lineage-marker module scores")
# Built from the per-sample raw counts already extracted in M28, on the SAME
# expression basis Phase 3 used (RNA counts -> LogNormalize). The 6 GB Phase 2
# object is not reloaded: these are the identical count matrices.
mats <- lapply(SAMPLES, function(s) readRDS(file.path(SC, s, "counts_raw.rds")))
genes <- sort(Reduce(union, lapply(mats, rownames)))
pad <- function(m) { miss <- setdiff(genes, rownames(m))
  if (length(miss)) m <- rbind(m, Matrix(0, nrow = length(miss), ncol = ncol(m),
      sparse = TRUE, dimnames = list(miss, colnames(m))))
  m[genes, , drop = FALSE] }
cm <- do.call(cbind, lapply(mats, pad)); rm(mats); invisible(gc(FALSE))
log_("combined counts: ", nrow(cm), " genes x ", ncol(cm), " cells")
so <- CreateSeuratObject(counts = cm, project = "phase4_markers", min.cells = 0, min.features = 0)
so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)

PANELS <- list(
  # MPNST / Schwann-lineage and neural-crest. Phase 2 recorded that SOX10, CNP,
  # PMP22 and NGFR are DOWN-regulated in MPNST relative to neurofibroma, so the
  # panel is deliberately broad rather than classic-Schwann-only.
  malignant_schwann_nc = c("S100B","SOX10","PLP1","MPZ","PMP22","NGFR","L1CAM","GFRA3",
                           "CRYAB","ERBB3","SOX9","TWIST1","PAX3","ZIC1","ETS1","FOXD3",
                           "NES","CDH19","PLEKHB1","MIA","SEMA3B","ABCB5","CCN3","RSPO3"),
  # Deliberately EXCLUDES COL1A1/COL1A2/FN1/POSTN, which Mes-NC-like malignant
  # cells can express - using those would beg the question.
  nonmalignant_fibroblast = c("DCN","LUM","PI16","CFD","FBLN1","FBLN2","MMP2","GSN",
                              "CLU","APOD","SFRP2","PDGFRA"),
  nonmalignant_immune  = c("PTPRC","CD3E","CD3D","CD2","CD68","CD14","LYZ","AIF1",
                           "CD79A","MS4A1","NKG7","GZMB","JCHAIN","TYROBP","FCER1G"),
  nonmalignant_endo    = c("PECAM1","VWF","CDH5","CLDN5","ERG","EGFL7","RAMP2","PLVAP"),
  nonmalignant_mural   = c("RGS5","NOTCH3","ACTA2","MYH11","PDGFRB","TAGLN","MYL9","HIGD1B"))
present <- lapply(PANELS, function(g) intersect(g, rownames(so)))
for (n in names(present)) log_(sprintf("  %-24s %d/%d genes present", n,
                                       length(present[[n]]), length(PANELS[[n]])))
empty <- names(present)[lengths(present) == 0]
if (length(empty)) {
  log_("WARNING: dropping panels with no genes present: ", paste(empty, collapse = ", "))
  present <- present[lengths(present) > 0]
}
stopifnot(length(present) > 0)
# ctrl must not exceed the smallest expression bin; 100 is the package default and
# is what real data of this size supports.
nbin <- 24
ctrl <- min(100L, floor(nrow(so) / (nbin * 3)))
log_("AddModuleScore: nbin=", nbin, " ctrl=", ctrl)
so <- AddModuleScore(so, features = present, name = "PANEL_", nbin = nbin,
                     ctrl = ctrl, seed = 42, verbose = FALSE)
sc_cols <- paste0("PANEL_", seq_along(present))
sco <- so@meta.data[, sc_cols, drop = FALSE]
colnames(sco) <- paste0("score_", names(present))
sco$cell_id <- rownames(sco)
md <- merge(md, sco, by = "cell_id", all.x = TRUE, sort = FALSE)
rownames(md) <- md$cell_id
rm(so, cm); invisible(gc(FALSE))

nm_cols <- intersect(paste0("score_", grep("^nonmalignant", names(present), value = TRUE)),
                     colnames(md))
stopifnot(length(nm_cols) > 0)
log_("non-malignant panels contributing to the max: ", paste(nm_cols, collapse = ", "))
md$score_nonmalignant_max <- do.call(pmax, c(md[nm_cols], list(na.rm = TRUE)))
med_mal <- median(md$score_malignant_schwann_nc, na.rm = TRUE)
med_non <- median(md$score_nonmalignant_max, na.rm = TRUE)
log_(sprintf("object-wide medians: malignant panel %.4f | non-malignant max %.4f",
             med_mal, med_non))
# Population-level marker evidence, per the rules document (population, not cell).
pm <- aggregate(cbind(score_malignant_schwann_nc, score_nonmalignant_max) ~ annotation_ccc,
                data = md, FUN = median)
pm$marker_evidence <- with(pm, ifelse(
  score_malignant_schwann_nc >  med_mal & score_nonmalignant_max <= med_non, "malignant-consistent",
  ifelse(score_malignant_schwann_nc <= med_mal & score_nonmalignant_max >  med_non,
         "non-malignant-consistent", "uninformative")))
print(pm, row.names = FALSE)
md$marker_evidence <- pm$marker_evidence[match(md$annotation_ccc, pm$annotation_ccc)]
write.table(pm, file.path(TAB, "LINEAGE_MARKER_EVIDENCE_BY_POPULATION.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- 7. apply the rules ----------------------------------------------------
sec("7. Applying R0-R99")
pc <- md$population_class; sv <- md$malignancy_scevan
ra <- md$runs_agree;        pf <- md$pop_frac
pfr <- md$pop_frac_reliable & !is.na(pf)
ce <- md$cnv_elevated
hi <- pfr & pf >= THR$pop_frac_high
lo <- pfr & pf <  THR$pop_frac_low
mid_or_hi <- pfr & pf >= THR$pop_frac_low

refined <- rep(NA_character_, nrow(md)); conf <- rep(NA_character_, nrow(md))
rule    <- rep(NA_character_, nrow(md))
set_rule <- function(mask, id, val, cf) {
  m <- mask & is.na(refined)
  refined[m] <<- val; conf[m] <<- cf; rule[m] <<- id
  log_(sprintf("  %-5s %-20s %-9s n=%d", id, val, cf, sum(m)))
}
na  <- sv == "not-assessed"
rel <- md$scevan_sample_reliable   # Amendment A1: promotions require a reliable sample
log_("cells in SCEVAN-reliable samples: ", sum(rel), " of ", nrow(md))
set_rule(pc == "EXCLUDED",                          "R0",  "Excluded-low-quality","High")
set_rule(na & pc == "CANONICAL_IMMUNE",             "R1a", "Non-malignant","Moderate")
set_rule(na,                                        "R1b", "Ambiguous","Low")
set_rule(pc=="PHASE2_MALIGNANT" & sv=="malignant" &  ra &  rel, "R2","Malignant","High")
set_rule(pc=="PHASE2_MALIGNANT" & sv=="malignant" &  ra & !rel, "R2u","Malignant","Moderate")
set_rule(pc=="PHASE2_MALIGNANT" & sv=="malignant" & !ra, "R3","Malignant","Moderate")
set_rule(pc=="PHASE2_MALIGNANT" & sv=="non-malignant" &  mid_or_hi, "R4","Malignant","Moderate")
set_rule(pc=="PHASE2_MALIGNANT" & sv=="non-malignant",             "R5","Ambiguous","Low")
# Amendment A2: plasma cells called malignant are attributed to the documented
# immunoglobulin-locus artefact rather than left as Ambiguous, which would remove
# 1,709 cells from the refined CCC analysis for a reason we can name.
set_rule(md$annotation_ccc=="Plasma-cell" & sv=="malignant", "R7p","Non-malignant","Moderate")
set_rule(pc=="CANONICAL_IMMUNE" & sv=="non-malignant","R6","Non-malignant","High")
set_rule(pc=="CANONICAL_IMMUNE" & sv=="malignant" & !hi,"R7","Non-malignant","Moderate")
set_rule(pc=="CANONICAL_IMMUNE" & sv=="malignant",     "R8","Ambiguous","Low")
set_rule(pc=="CANONICAL_ENDO"   & sv=="non-malignant","R9","Non-malignant","High")
set_rule(pc=="CANONICAL_ENDO"   & sv=="malignant" & !hi,"R10","Non-malignant","Moderate")
set_rule(pc=="CANONICAL_ENDO"   & sv=="malignant",     "R11","Ambiguous","Low")
set_rule(pc=="DISPUTED_STROMAL" & sv=="non-malignant" & lo, "R12","Non-malignant","High")
set_rule(pc=="DISPUTED_STROMAL" & sv=="malignant" & !rel, "R13x","Ambiguous","Low")
set_rule(pc=="DISPUTED_STROMAL" & sv=="malignant" & ra & hi & ce, "R13","Malignant","High")
set_rule(pc=="DISPUTED_STROMAL" & sv=="malignant" & mid_or_hi,    "R14","Malignant","Moderate")
set_rule(pc=="DISPUTED_STROMAL" & sv=="malignant",                "R15","Ambiguous","Low")
set_rule(pc=="DISPUTED_STROMAL" & sv=="non-malignant",            "R16","Non-malignant","Moderate")
set_rule(pc=="PROVISIONAL" & sv=="malignant" & !rel,    "R17x","Ambiguous","Low")
set_rule(pc=="PROVISIONAL" & sv=="malignant" & ra & hi, "R17","Malignant","High")
set_rule(pc=="PROVISIONAL" & sv=="malignant",           "R18","Malignant","Moderate")
set_rule(pc=="PROVISIONAL" & sv=="non-malignant" & lo,  "R19","Non-malignant","Moderate")
set_rule(pc=="PROVISIONAL" & sv=="non-malignant",       "R20","Ambiguous","Low")
set_rule(rep(TRUE, nrow(md)),                           "R99","Ambiguous","Low")
if (sum(rule == "R99") > 0)
  log_("WARNING: ", sum(rule == "R99"), " cells reached the R99 fail-safe (rule-coverage gap)")

md$malignancy_refined <- refined
md$malignancy_confidence <- conf
md$malignancy_rule <- rule

# --- 8. marker-evidence confidence modifier -------------------------------
sec("8. Marker-evidence confidence modifier (downgrade only)")
dn_mal <- md$malignancy_refined == "Malignant" & md$malignancy_confidence == "High" &
          md$marker_evidence == "non-malignant-consistent"
dn_non <- md$malignancy_refined == "Non-malignant" & md$malignancy_confidence == "High" &
          md$marker_evidence == "malignant-consistent"
log_("downgraded Malignant/High -> Moderate: ", sum(dn_mal, na.rm = TRUE))
log_("downgraded Non-malignant/High -> Moderate: ", sum(dn_non, na.rm = TRUE))
md$marker_modifier <- ""
md$marker_modifier[which(dn_mal | dn_non)] <- "marker_evidence_discordant"
md$malignancy_confidence[which(dn_mal | dn_non)] <- "Moderate"

md$malignancy_reason <- sprintf(
  "%s | phase2=%s scevan=%s sens=%s runs_agree=%s pop_frac=%s(n=%d,reliable=%s) cnv_elevated=%s marker=%s sample_scevan_reliable=%s%s",
  md$malignancy_rule, md$malignancy_phase2, md$malignancy_scevan, md$malignancy_scevan_sens,
  md$runs_agree, ifelse(is.na(md$pop_frac), "NA", sprintf("%.3f", md$pop_frac)),
  ifelse(is.na(md$pop_n_assessed), 0L, md$pop_n_assessed), md$pop_frac_reliable,
  md$cnv_elevated, md$marker_evidence, md$scevan_sample_reliable,
  ifelse(md$marker_modifier == "", "", paste0(" | ", md$marker_modifier)))

sec("9. Refined malignancy outcome")
print(table(md$malignancy_refined, md$malignancy_confidence))
print(table(md$malignancy_refined, md$sample_id))

# --- 10. annotation_ccc_refined -------------------------------------------
sec("10. annotation_ccc_refined (§29)")
md$annotation_ccc_phase3 <- md$annotation_ccc      # frozen Phase 3 layer, verbatim
md$annotation_ccc_refined <- with(md, ifelse(
  malignancy_refined == "Malignant", "MPNST-Tumor",
  ifelse(malignancy_refined == "Excluded-low-quality", "Low-quality-excluded",
  ifelse(malignancy_refined == "Ambiguous", "Ambiguous-unresolved", annotation_ccc))))
print(table(md$annotation_ccc_refined, md$sample_id))
# Ambiguous cells must never become tumour.
stopifnot(!any(md$annotation_ccc_refined == "MPNST-Tumor" &
               md$malignancy_refined != "Malignant"))
log_("GUARD PASSED: no Ambiguous or Non-malignant cell is labelled MPNST-Tumor")
stopifnot(identical(md$annotation_ccc_phase3, md$annotation_ccc))
log_("GUARD PASSED: annotation_ccc_phase3 is a verbatim copy of the frozen Phase 3 layer")

# --- 11. outputs ----------------------------------------------------------
sec("11. Writing outputs")
ord <- c("cell_id","sample_id","patient","postint_harmony_primary_cluster",
         "annotation_level1","annotation_level2","annotation_level3","annotation_confidence",
         "annotation_ccc","annotation_ccc_phase3","annotation_ccc_refined",
         "population_class","scevan_sample","scevan_call","scevan_confident_normal",
         "scevan_clone","scevan_subclone_raw","scevan_sample_reliable",
         "cnv_burden","cnv_frac_gain","cnv_frac_loss","cnv_mean_abs",
         "cnv_ref_threshold","cnv_elevated",
         "pop_frac","pop_n_assessed","pop_frac_reliable","runs_agree",
         "score_malignant_schwann_nc","score_nonmalignant_max","marker_evidence",
         "malignancy_phase2","malignancy_scevan","malignancy_scevan_sens",
         "malignancy_refined","malignancy_confidence","malignancy_rule","malignancy_reason")
ord <- intersect(ord, colnames(md))
md$.assessed <- NULL; md$.malignant <- NULL
calls <- md[, ord, drop = FALSE]
names(calls)[names(calls) == "postint_harmony_primary_cluster"] <- "primary_cluster"
f1 <- file.path(MAL, "PHASE4_MALIGNANCY_CALLS.tsv")
write.table(calls, f1, sep = "\t", quote = FALSE, row.names = FALSE)
log_("wrote ", f1, " (", nrow(calls), " x ", ncol(calls), ")")
saveRDS(md, file.path(MAL, "phase4_malignancy_metadata.rds"))

# SCEVAN vs Phase 2 confusion, per population and per patient (§21)
xt <- do.call(rbind, lapply(split(md, list(md$annotation_ccc, md$sample_id), drop = TRUE),
  function(d) data.frame(
    annotation_ccc = d$annotation_ccc[1], sample_id = d$sample_id[1],
    n_total = nrow(d),
    n_scevan_malignant = sum(d$malignancy_scevan == "malignant"),
    n_scevan_normal    = sum(d$malignancy_scevan == "non-malignant"),
    n_scevan_notassessed = sum(d$malignancy_scevan == "not-assessed"),
    fraction_malignant = mean(d$malignancy_scevan == "malignant"),
    fraction_malignant_of_assessed = if (sum(d$malignancy_scevan != "not-assessed")) 
      sum(d$malignancy_scevan == "malignant")/sum(d$malignancy_scevan != "not-assessed") else NA_real_,
    n_confident_normal = sum(d$scevan_confident_normal == "yes"),
    median_cnv_burden = median(d$cnv_burden, na.rm = TRUE),
    n_refined_malignant = sum(d$malignancy_refined == "Malignant"),
    n_refined_nonmalignant = sum(d$malignancy_refined == "Non-malignant"),
    n_refined_ambiguous = sum(d$malignancy_refined == "Ambiguous"),
    stringsAsFactors = FALSE)))
xt <- xt[order(xt$annotation_ccc, xt$sample_id), ]
f2 <- file.path(TAB, "SCEVAN_VS_PHASE2_ANNOTATION.tsv")
write.table(xt, f2, sep = "\t", quote = FALSE, row.names = FALSE)
log_("wrote ", f2)
print(xt[, c("annotation_ccc","sample_id","n_total","n_scevan_malignant",
             "n_scevan_normal","fraction_malignant_of_assessed","n_refined_malignant")],
      row.names = FALSE)

# per-cell SCEVAN classification table (final deliverable name)
f3 <- file.path(TAB, "SCEVAN_CELL_CLASSIFICATION.tsv")
write.table(md[, c("cell_id","sample_id","annotation_ccc","scevan_call",
                   "scevan_confident_normal","scevan_clone","malignancy_scevan",
                   "malignancy_scevan_sens","runs_agree","cnv_burden","cnv_elevated")],
            f3, sep = "\t", quote = FALSE, row.names = FALSE)
log_("wrote ", f3)

# summaries by cluster and by sample
by_cl <- do.call(rbind, lapply(split(md, md$postint_harmony_primary_cluster), function(d)
  data.frame(primary_cluster = d$postint_harmony_primary_cluster[1],
    dominant_annotation_ccc = names(sort(table(d$annotation_ccc), decreasing = TRUE))[1],
    n = nrow(d),
    n_phase2_malignant = sum(d$malignancy_phase2 == "malignant"),
    n_scevan_malignant = sum(d$malignancy_scevan == "malignant"),
    n_refined_malignant = sum(d$malignancy_refined == "Malignant"),
    frac_refined_malignant = mean(d$malignancy_refined == "Malignant"),
    n_refined_ambiguous = sum(d$malignancy_refined == "Ambiguous"),
    median_cnv_burden = median(d$cnv_burden, na.rm = TRUE),
    n_patients = length(unique(d$sample_id)),
    dominant_patient_frac = max(table(d$sample_id))/nrow(d),
    stringsAsFactors = FALSE)))
write.table(by_cl, file.path(TAB, "MALIGNANCY_SUMMARY_BY_CLUSTER.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(by_cl, row.names = FALSE)

by_s <- do.call(rbind, lapply(split(md, md$sample_id), function(d)
  data.frame(sample_id = d$sample_id[1], patient = d$sample_id[1], n_cells = nrow(d),
    phase2_malignant = sum(d$malignancy_phase2 == "malignant"),
    phase2_fraction = mean(d$malignancy_phase2 == "malignant"),
    scevan_assessed = sum(d$malignancy_scevan != "not-assessed"),
    scevan_malignant = sum(d$malignancy_scevan == "malignant"),
    scevan_fraction_of_all = mean(d$malignancy_scevan == "malignant"),
    scevan_fraction_of_assessed = sum(d$malignancy_scevan=="malignant")/
                                  max(1,sum(d$malignancy_scevan!="not-assessed")),
    refined_malignant = sum(d$malignancy_refined == "Malignant"),
    refined_fraction = mean(d$malignancy_refined == "Malignant"),
    refined_ambiguous = sum(d$malignancy_refined == "Ambiguous"),
    refined_nonmalignant = sum(d$malignancy_refined == "Non-malignant"),
    n_clones = length(unique(na.omit(d$scevan_clone))),
    stringsAsFactors = FALSE)))
write.table(by_s, file.path(TAB, "MALIGNANCY_SUMMARY_BY_SAMPLE.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(by_s, row.names = FALSE)

# ambiguous-population audit (§22)
AUD <- c("Fibroblast","Candidate-Malignant-Unresolved","Pericyte-VSMC","MPNST-Tumor","Uncertain")
aud <- do.call(rbind, lapply(AUD, function(a) {
  d <- md[md$annotation_ccc == a, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  ps <- do.call(rbind, lapply(split(d, d$sample_id), function(x) data.frame(
    sample_id = x$sample_id[1], n = nrow(x),
    frac_scevan_malignant = mean(x$malignancy_scevan == "malignant"),
    frac_refined_malignant = mean(x$malignancy_refined == "Malignant"))))
  data.frame(annotation_ccc = a, n_total = nrow(d),
    n_scevan_malignant = sum(d$malignancy_scevan == "malignant"),
    n_scevan_normal = sum(d$malignancy_scevan == "non-malignant"),
    n_not_assessed = sum(d$malignancy_scevan == "not-assessed"),
    frac_scevan_malignant_of_assessed =
      sum(d$malignancy_scevan=="malignant")/max(1,sum(d$malignancy_scevan!="not-assessed")),
    patients_with_majority_malignant = sum(ps$frac_scevan_malignant >= 0.5),
    patients_n = nrow(ps),
    per_patient_frac = paste(sprintf("%s=%.3f", ps$sample_id, ps$frac_scevan_malignant),
                             collapse = "; "),
    n_clones = length(unique(na.omit(d$scevan_clone))),
    clone_membership = {
      tt <- head(sort(table(na.omit(d$scevan_clone)), decreasing = TRUE), 6)
      if (length(tt)) paste(sprintf("%s=%d", names(tt), as.integer(tt)), collapse = "; ") else "" },
    median_cnv_burden = median(d$cnv_burden, na.rm = TRUE),
    frac_cnv_elevated = mean(d$cnv_elevated, na.rm = TRUE),
    median_score_malignant = median(d$score_malignant_schwann_nc, na.rm = TRUE),
    median_score_nonmalignant = median(d$score_nonmalignant_max, na.rm = TRUE),
    marker_evidence = d$marker_evidence[1],
    n_refined_malignant = sum(d$malignancy_refined == "Malignant"),
    n_refined_nonmalignant = sum(d$malignancy_refined == "Non-malignant"),
    n_refined_ambiguous = sum(d$malignancy_refined == "Ambiguous"),
    rules_applied = paste(sprintf("%s=%d", names(sort(table(d$malignancy_rule),
                          decreasing = TRUE)), sort(table(d$malignancy_rule),
                          decreasing = TRUE)), collapse = "; "),
    stringsAsFactors = FALSE)
}))
f4 <- file.path(TAB, "AMBIGUOUS_POPULATION_MALIGNANCY_AUDIT.tsv")
write.table(aud, f4, sep = "\t", quote = FALSE, row.names = FALSE)
log_("wrote ", f4)
print(aud[, c("annotation_ccc","n_total","frac_scevan_malignant_of_assessed",
              "patients_with_majority_malignant","median_cnv_burden",
              "n_refined_malignant","n_refined_ambiguous")], row.names = FALSE)

# --- 12. headline fractions (§28) -----------------------------------------
sec("12. Malignant fraction: Phase 2 -> SCEVAN -> refined")
N <- nrow(md)
hl <- list(
  n_cells = N,
  phase2_malignant = sum(md$malignancy_phase2 == "malignant"),
  phase2_fraction = mean(md$malignancy_phase2 == "malignant"),
  scevan_assessed = sum(md$malignancy_scevan != "not-assessed"),
  scevan_malignant = sum(md$malignancy_scevan == "malignant"),
  scevan_fraction_of_all = mean(md$malignancy_scevan == "malignant"),
  scevan_fraction_of_assessed = sum(md$malignancy_scevan == "malignant") /
                                sum(md$malignancy_scevan != "not-assessed"),
  refined_malignant = sum(md$malignancy_refined == "Malignant"),
  refined_fraction = mean(md$malignancy_refined == "Malignant"),
  refined_nonmalignant = sum(md$malignancy_refined == "Non-malignant"),
  refined_ambiguous = sum(md$malignancy_refined == "Ambiguous"),
  refined_excluded = sum(md$malignancy_refined == "Excluded-low-quality"),
  refined_high_confidence_malignant = sum(md$malignancy_refined == "Malignant" &
                                          md$malignancy_confidence == "High"))
str(hl)

# --- 12b. threshold sensitivity ---------------------------------------------
sec("12b. Sensitivity of the headline numbers to pop_frac_low")
# The a priori pop_frac_low = 0.25 separates R4 (Malignant/Moderate) from R5
# (Ambiguous/Low) for Phase 2 malignant cells that SCEVAN calls normal, and R14
# from R15 for disputed stromal cells. Some strata land close to it, so the
# headline numbers are recomputed across a grid. The PRIMARY result remains the
# a priori 0.25; this table exists so a reader can see the knife-edges rather
# than discover them, and NO threshold was chosen to improve an answer.
grid <- c(0.15, 0.20, 0.25, 0.30, 0.40)
thr_sens <- do.call(rbind, lapply(grid, function(g) {
  hi2 <- pfr & pf >= THR$pop_frac_high
  lo2 <- pfr & pf <  g
  mid2 <- pfr & pf >= g
  r <- rep(NA_character_, nrow(md))
  setr <- function(mask, val) { m <- mask & is.na(r); r[m] <<- val }
  setr(pc == "EXCLUDED", "Excluded-low-quality")
  setr(na & pc == "CANONICAL_IMMUNE", "Non-malignant")
  setr(na, "Ambiguous")
  setr(pc=="PHASE2_MALIGNANT" & sv=="malignant" &  ra &  rel, "Malignant")
  setr(pc=="PHASE2_MALIGNANT" & sv=="malignant" &  ra & !rel, "Malignant")
  setr(pc=="PHASE2_MALIGNANT" & sv=="malignant" & !ra, "Malignant")
  setr(pc=="PHASE2_MALIGNANT" & sv=="non-malignant" &  mid2, "Malignant")
  setr(pc=="PHASE2_MALIGNANT" & sv=="non-malignant", "Ambiguous")
  setr(pc=="CANONICAL_IMMUNE" & sv=="non-malignant", "Non-malignant")
  setr(pc=="CANONICAL_IMMUNE" & sv=="malignant" & !hi2, "Non-malignant")
  setr(pc=="CANONICAL_IMMUNE" & sv=="malignant", "Ambiguous")
  setr(pc=="CANONICAL_ENDO"   & sv=="non-malignant", "Non-malignant")
  setr(pc=="CANONICAL_ENDO"   & sv=="malignant" & !hi2, "Non-malignant")
  setr(pc=="CANONICAL_ENDO"   & sv=="malignant", "Ambiguous")
  setr(pc=="DISPUTED_STROMAL" & sv=="non-malignant" & lo2, "Non-malignant")
  setr(pc=="DISPUTED_STROMAL" & sv=="malignant" & !rel, "Ambiguous")
  setr(pc=="DISPUTED_STROMAL" & sv=="malignant" & ra & hi2 & ce, "Malignant")
  setr(pc=="DISPUTED_STROMAL" & sv=="malignant" & mid2, "Malignant")
  setr(pc=="DISPUTED_STROMAL" & sv=="malignant", "Ambiguous")
  setr(pc=="DISPUTED_STROMAL" & sv=="non-malignant", "Non-malignant")
  setr(pc=="PROVISIONAL" & sv=="malignant" & !rel, "Ambiguous")
  setr(pc=="PROVISIONAL" & sv=="malignant" & ra & hi2, "Malignant")
  setr(pc=="PROVISIONAL" & sv=="malignant", "Malignant")
  setr(pc=="PROVISIONAL" & sv=="non-malignant" & lo2, "Non-malignant")
  setr(pc=="PROVISIONAL" & sv=="non-malignant", "Ambiguous")
  setr(rep(TRUE, nrow(md)), "Ambiguous")
  data.frame(pop_frac_low = g,
    is_a_priori = identical(g, THR$pop_frac_low),
    Malignant = sum(r == "Malignant"),
    malignant_fraction = mean(r == "Malignant"),
    Non_malignant = sum(r == "Non-malignant"),
    Ambiguous = sum(r == "Ambiguous"),
    fibroblast_malignant = sum(r == "Malignant" & md$annotation_ccc == "Fibroblast"),
    mpnst_tumor_still_malignant = sum(r == "Malignant" & md$annotation_ccc == "MPNST-Tumor"),
    candidate_malignant = sum(r == "Malignant" &
      md$annotation_ccc == "Candidate-Malignant-Unresolved"))
}))
print(thr_sens, row.names = FALSE, digits = 4)
write.table(thr_sens, file.path(TAB, "MALIGNANCY_THRESHOLD_SENSITIVITY.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
log_("The a priori 0.25 row is the reported result. The grid is disclosure, not tuning.")

facts <- list(thresholds = THR, immune_sanity_max = IMMUNE_SANITY_MAX,
              immune_pops_max = IMMUNE_POPS_MAX,
              sanity_artefact_pops = SANITY_ARTEFACT_POPS,
              threshold_sensitivity = thr_sens,
              sample_reliability = sanity, population_classes = CLS, headline = hl,
              cna_reference = cna_ref,
              marker_panels_present = present,
              marker_medians = list(malignant = med_mal, nonmalignant_max = med_non),
              rule_counts = as.list(table(md$malignancy_rule)),
              refined_by_confidence = as.data.frame.matrix(
                table(md$malignancy_refined, md$malignancy_confidence)),
              refined_by_sample = as.data.frame.matrix(
                table(md$malignancy_refined, md$sample_id)),
              ccc_refined_counts = as.list(table(md$annotation_ccc_refined)),
              scevan_version = as.character(packageVersion("SCEVAN")),
              slurm = list(job_id = Sys.getenv("SLURM_JOB_ID")),
              generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
write_json(facts, file.path(MAL, "m30_malignancy_facts.json"),
           auto_unbox = TRUE, pretty = TRUE, digits = 8, null = "null")
sec("M30 INTEGRATION COMPLETE")
