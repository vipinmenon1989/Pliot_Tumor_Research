#!/usr/bin/env Rscript
# =============================================================================
# Phase 6 - M44 - within-clone transcriptional-program diversity.
#
# WORDING RULE, enforced throughout: this measures DIVERSITY, not switching.
# Cross-sectional scRNA-seq cannot show a transition. The permitted claim is
# "clone X contains cells spanning multiple malignant transcriptional
# programs"; it may be read as "consistent with transcriptional plasticity"
# only after the technical confounders below have been assessed.
#
# Diversity is measured with several metrics because no single one is
# definitive, and on BOTH the continuous program scores and a cautiously
# defined dominant-program assignment whose threshold is printed.
# =============================================================================
source("scripts/phase6/utils/phase6_common.R")
set.seed(42)
facts <- list(milestone = "M44", generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))

md <- readRDS(file.path(P6_CLONE, "phase6_cell_metadata.rds"))
PROGS <- grep("^program_P[0-9]+_score$", colnames(md), value = TRUE)
PID <- sub("^program_(P[0-9]+)_score$", "\\1", PROGS)
lab <- read.delim("results/phase5/tables/final/MALIGNANT_PROGRAMS.tsv")
LABEL <- setNames(lab$program_label, lab$program_id)

mal <- md[md$malignancy_refined == "Malignant" & !is.na(md$tumor_clone_phase4), ]
mal_rel <- mal[mal$sample_id %in% CLONE_RELIABLE, ]
S <- as.matrix(mal_rel[, PROGS]); colnames(S) <- PID
p6_msg("clone-reliable malignant cells with a clone: %d across %d clones",
       nrow(mal_rel), length(unique(mal_rel$tumor_clone_phase4)))
p6_msg("MPNST_3 excluded: %s", CLONE_EXCLUSION_REASON)

# hard assignment, with the threshold printed rather than hidden
MARGIN <- 0.10
mal_rel$hard_program <- ifelse(mal_rel$dominant_program_margin >= MARGIN,
                               mal_rel$dominant_malignant_program, "Mixed")
p6_msg("hard assignment: a cell keeps its top program only when it leads the runner-up by >= %.2f relative usage", MARGIN)
print(table(mal_rel$hard_program))

p6_sec("1. Diversity metrics per clone")
div <- do.call(rbind, lapply(split(seq_len(nrow(mal_rel)),
                                   mal_rel$tumor_clone_phase4), function(i) {
  d <- mal_rel[i, ]; M <- S[i, , drop = FALSE]
  # continuous view: mean program profile of the clone, and its entropy
  prof <- colMeans(M); prof <- prof / sum(prof)
  # hard view: composition over dominant-program labels (excluding Mixed)
  hp <- d$hard_program[d$hard_program != "Mixed"]
  hcomp <- if (length(hp)) prop.table(table(hp)) else numeric(0)
  # dispersion: mean Euclidean distance of a cell from the clone centroid
  cen <- colMeans(M)
  disp <- mean(sqrt(rowSums((M - matrix(cen, nrow(M), ncol(M), byrow = TRUE))^2)))
  data.frame(
    sample_id = d$sample_id[1], clone = d$tumor_clone_phase4[1], n_cells = nrow(d),
    continuous_entropy = shannon(prof),
    continuous_effective_n = effective_n(prof),
    hard_entropy = if (length(hcomp)) shannon(as.numeric(hcomp)) else NA_real_,
    hard_effective_n = if (length(hcomp)) effective_n(as.numeric(hcomp)) else NA_real_,
    dominant_program = names(which.max(prof)),
    dominant_program_share_continuous = max(prof),
    dominant_program_share_hard = if (length(hcomp)) max(hcomp) else NA_real_,
    n_programs_hard = length(hcomp),
    frac_mixed_cells = mean(d$hard_program == "Mixed"),
    program_dispersion = disp,
    median_nCount = median(d$nCount_RNA), median_nFeature = median(d$nFeature_RNA),
    frac_high_confidence = mean(d$malignancy_confidence == "High"))
}))
div$dominant_program_label <- unname(LABEL[div$dominant_program])
div$max_effective_n <- length(PID)
for (th in CLONE_SIZE_GRID)
  div[[paste0("evaluable_ge", th)]] <- div$n_cells >= th
div <- div[order(div$sample_id, -div$n_cells), ]
print(as.data.frame(div |> select(sample_id, clone, n_cells, continuous_effective_n,
                                  hard_effective_n, dominant_program_share_hard,
                                  n_programs_hard, frac_mixed_cells)),
      row.names = FALSE, digits = 3)
p6_tsv(div, file.path(P6_TAB, "CLONE_PROGRAM_DIVERSITY.tsv"))

p6_sec("2. Between-clone divergence within each patient")
jsd <- function(p, q) {
  m <- (p + q) / 2
  0.5 * sum(p[p > 0] * log(p[p > 0] / m[p > 0])) +
    0.5 * sum(q[q > 0] * log(q[q > 0] / m[q > 0]))
}
btw <- do.call(rbind, lapply(CLONE_RELIABLE, function(s) {
  d <- mal_rel[mal_rel$sample_id == s, ]
  cl <- names(which(table(d$tumor_clone_phase4) >= CLONE_SIZE_PRIMARY))
  if (length(cl) < 2) return(data.frame(sample_id = s, n_clones = length(cl),
    mean_between_clone_jsd = NA_real_, mean_within_clone_dispersion = NA_real_,
    ratio = NA_real_))
  profs <- t(vapply(cl, function(c0) {
    m <- colMeans(S[d$tumor_clone_phase4 == c0, , drop = FALSE]); m / sum(m)
  }, numeric(ncol(S))))
  pr <- combn(seq_along(cl), 2)
  j <- mean(apply(pr, 2, function(k) jsd(profs[k[1], ], profs[k[2], ])))
  w <- mean(div$program_dispersion[div$sample_id == s & div$clone %in% cl])
  data.frame(sample_id = s, n_clones = length(cl), mean_between_clone_jsd = j,
             mean_within_clone_dispersion = w, ratio = j / w)
}))
print(as.data.frame(btw), row.names = FALSE, digits = 3)
p6_tsv(btw, file.path(P6_PLAS, "M44_BETWEEN_VS_WITHIN_CLONE.tsv"))

p6_sec("3. Technical confounders - assessed before any plasticity wording")
ev <- div[div$n_cells >= CLONE_SIZE_PRIMARY, ]
cc_prog <- lab$program_id[grepl("Cycling", lab$program_label, ignore.case = TRUE)]
if (length(cc_prog)) {
  cc_col <- paste0("program_", cc_prog[1], "_score")
  ev$median_cellcycle_program <- vapply(seq_len(nrow(ev)), function(i)
    median(mal_rel[[cc_col]][mal_rel$tumor_clone_phase4 == ev$clone[i]]), numeric(1))
  p6_msg("cell-cycle program used as the proliferation covariate: %s (%s)",
         cc_prog[1], LABEL[cc_prog[1]])
} else {
  ev$median_cellcycle_program <- NA_real_
  p6_msg("no program was labelled Cycling; the proliferation covariate is NA")
}
conf <- do.call(rbind, lapply(c("n_cells", "median_nCount", "median_nFeature",
                                "frac_high_confidence", "median_cellcycle_program"),
                              function(v) {
  x <- ev[[v]]
  do.call(rbind, lapply(c("continuous_effective_n", "hard_effective_n",
                          "program_dispersion"), function(y) {
    yy <- ev[[y]]
    ok <- is.finite(x) & is.finite(yy)
    data.frame(covariate = v, diversity_metric = y, n = sum(ok),
               spearman_rho = if (sum(ok) >= 4)
                 suppressWarnings(cor(x[ok], yy[ok], method = "spearman")) else NA_real_)
  }))
}))
print(as.data.frame(conf), row.names = FALSE, digits = 3)
p6_tsv(conf, file.path(P6_PLAS, "M44_DIVERSITY_CONFOUNDERS.tsv"))
strong <- conf[!is.na(conf$spearman_rho) & abs(conf$spearman_rho) >= 0.7, ]
if (nrow(strong)) {
  p6_msg("CAUTION: diversity tracks a technical covariate at |rho| >= 0.7:")
  print(as.data.frame(strong), row.names = FALSE, digits = 3)
} else p6_msg("no covariate reaches |rho| >= 0.7 against a diversity metric")

p6_sec("4. Clone-size sensitivity of the diversity estimates")
sens <- do.call(rbind, lapply(CLONE_SIZE_GRID, function(th) {
  e <- div[div$n_cells >= th, ]
  data.frame(min_clone_size = th, n_clones_evaluable = nrow(e),
             n_clones_not_evaluable = sum(div$n_cells < th),
             median_continuous_effective_n = median(e$continuous_effective_n),
             median_hard_effective_n = median(e$hard_effective_n, na.rm = TRUE),
             median_dominant_share_hard = median(e$dominant_program_share_hard, na.rm = TRUE))
}))
print(as.data.frame(sens), row.names = FALSE, digits = 3)
p6_tsv(sens, file.path(P6_PLAS, "M44_CLONE_SIZE_SENSITIVITY.tsv"))

p6_sec("5. PLASTICITY_METRICS.tsv")
pm <- div |>
  select(sample_id, clone, n_cells, dominant_program, dominant_program_label,
         dominant_program_share_continuous, dominant_program_share_hard,
         n_programs_hard, continuous_entropy, continuous_effective_n,
         hard_entropy, hard_effective_n, program_dispersion, frac_mixed_cells,
         median_nCount, median_nFeature, frac_high_confidence,
         starts_with("evaluable_")) |>
  mutate(interpretation = ifelse(
    n_cells < CLONE_SIZE_PRIMARY, "NOT EVALUABLE - clone below the declared minimum",
    ifelse(hard_effective_n >= 2,
      "contains cells spanning multiple malignant transcriptional programs",
      "cells concentrate in a single dominant program")),
    wording_rule = "diversity, not observed switching; no transition rate, direction or trajectory is claimed")
p6_tsv(pm, file.path(P6_TAB, "PLASTICITY_METRICS.tsv"))
nmulti <- sum(pm$n_cells >= CLONE_SIZE_PRIMARY & pm$hard_effective_n >= 2, na.rm = TRUE)
nev <- sum(pm$n_cells >= CLONE_SIZE_PRIMARY)
p6_msg("clones with >= %d cells: %d; of those, %d span multiple programs (effective n >= 2)",
       CLONE_SIZE_PRIMARY, nev, nmulti)

facts$diversity <- div; facts$between_within <- btw
facts$confounders <- conf; facts$clone_size_sensitivity <- sens
facts$hard_assignment_margin <- MARGIN
facts$n_clones_evaluable <- nev; facts$n_clones_multi_program <- nmulti
facts$wording <- paste("Within-clone program DIVERSITY is what is measured.",
  "Cross-sectional scRNA-seq cannot demonstrate a state transition, so no",
  "switching, transition rate, direction or trajectory is claimed anywhere.")
p6_json(facts, file.path(P6_VAL, "m44_diversity_facts.json"))
p6_sec("M44 complete")
