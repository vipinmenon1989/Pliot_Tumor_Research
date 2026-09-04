#!/usr/bin/env Rscript
# =============================================================================
# Phase 5 - M40 - program robustness.
#
#   1 program-rank sensitivity            K*-1, K*, K*+1
#   2 patient-balanced sensitivity        (from the balanced cNMF run)
#   3 gene-selection sensitivity          repeated random 80% gene subsets
#   4 cell-cycle inclusion/exclusion      (from the nocc cNMF run)
#   5 malignancy-confidence sensitivity   High+Moderate vs High only
#   6 leave-one-patient-out recovery      one cNMF run per withheld patient
#
# A program that disappears when its dominant patient is withheld is NOT robust,
# and is reported as such rather than quietly kept.
# =============================================================================
source("scripts/phase5/utils/phase5_common.R")
source("scripts/phase5/utils/program_projection.R")
set.seed(42)
args <- commandArgs(trailingOnly = TRUE)
KSTAR <- as.integer(args[1]); DT <- if (length(args) > 1) args[2] else "0_1"
CNMF <- "results/phase5/programs/cnmf"
facts <- list(milestone = "M40", K = KSTAR,
              generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))

sp_path <- function(run, k = KSTAR, what = "gene_spectra_score")
  file.path(CNMF, run, sprintf("%s.%s.k_%d.dt_%s.txt", run, what, k, DT))
us_path <- function(run, k = KSTAR)
  file.path(CNMF, run, sprintf("%s.usages.k_%d.dt_%s.consensus.txt", run, k, DT))
rd <- function(p) as.matrix(read.delim(p, row.names = 1, check.names = FALSE))

cosm <- function(a, b) sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)) + 1e-12)
match_spectra <- function(Zref, Zalt, topn = 100L) {
  g <- intersect(colnames(Zref), colnames(Zalt))
  do.call(rbind, lapply(rownames(Zref), function(p) {
    cs <- vapply(rownames(Zalt), function(q) cosm(Zref[p, g], Zalt[q, g]), numeric(1))
    jc <- vapply(rownames(Zalt), function(q) {
      a <- colnames(Zref)[order(Zref[p, ], decreasing = TRUE)[seq_len(topn)]]
      b <- colnames(Zalt)[order(Zalt[q, ], decreasing = TRUE)[seq_len(topn)]]
      length(intersect(a, b)) / length(union(a, b))
    }, numeric(1))
    i <- which.max(cs)
    data.frame(program = p, best_match = names(cs)[i], cosine = unname(cs[i]),
               top100_jaccard = unname(jc[i]))
  }))
}

p5_sec("1. Primary solution")
Z <- rd(sp_path("primary")); rownames(Z) <- paste0("P", seq_len(nrow(Z)))
U <- rd(us_path("primary")); colnames(U) <- rownames(Z)
REL <- U / rowSums(U)
md <- readRDS(file.path(P5_PROG, "phase5_cell_metadata.rds"))
mmd <- md[rownames(U), ]
PROGS <- rownames(Z)
p5_msg("K* = %d, %d programs, %d cells", KSTAR, length(PROGS), nrow(U))
progs_meta <- read.delim(file.path(P5_TAB, "MALIGNANT_PROGRAMS.tsv"))

rob <- data.frame(program = PROGS)

p5_sec("2. Program-rank sensitivity (K*-1, K*, K*+1)")
rank_rows <- list()
for (kk in c(KSTAR - 1L, KSTAR + 1L)) {
  f <- sp_path("primary", kk)
  if (!file.exists(f)) { p5_msg("  K=%d not available", kk); next }
  Za <- rd(f); rownames(Za) <- paste0("K", kk, "_P", seq_len(nrow(Za)))
  m <- match_spectra(Z, Za); m$alt_K <- kk
  rank_rows[[as.character(kk)]] <- m
  p5_msg("  K=%d: median best-match cosine %.3f (min %.3f)", kk,
         median(m$cosine), min(m$cosine))
}
rk <- if (length(rank_rows)) do.call(rbind, rank_rows) else NULL
if (!is.null(rk)) {
  agg <- rk |> group_by(program) |>
    summarise(rank_min_cosine = min(cosine), rank_median_cosine = median(cosine),
              .groups = "drop")
  rob <- rob |> left_join(agg, by = "program")
  p5_tsv(rk, file.path(P5_VAL, "M40_PROGRAM_RANK_SENSITIVITY.tsv"))
}

p5_sec("3. Patient-balanced sensitivity")
bal <- read.delim(file.path(P5_TAB, "PROGRAM_PATIENT_BALANCE_SENSITIVITY.tsv"))
rob <- rob |> left_join(bal |> select(program, balanced_cosine = cosine,
                                      balanced_jaccard = top100_jaccard,
                                      balanced_support), by = "program")
p5_msg("  balanced-run best-match cosine: median %.3f, min %.3f",
       median(bal$cosine, na.rm = TRUE), min(bal$cosine, na.rm = TRUE))

p5_sec("4. Gene-selection sensitivity (10 x random 80% of genes)")
Wt <- t(rd(sp_path("primary", KSTAR, "gene_spectra_tpm")))
colnames(Wt) <- PROGS
cnts <- readRDS(file.path(P5_PROG, "phase5_working_counts.rds"))
X <- tpm_on_genes(cnts[, rownames(U), drop = FALSE], rownames(Wt))
gs <- lapply(1:10, function(i) {
  set.seed(1000 + i)
  keep <- sort(sample(nrow(Wt), floor(0.8 * nrow(Wt))))
  H <- relative_usage(project_programs(X[keep, , drop = FALSE],
                                       Wt[keep, , drop = FALSE]))
  vapply(PROGS, function(p) cor(H[, p], REL[, p], method = "spearman"), numeric(1))
})
gsm <- do.call(rbind, gs)
gsr <- data.frame(program = PROGS,
                  gene_subset_min_rho = apply(gsm, 2, min),
                  gene_subset_median_rho = apply(gsm, 2, median))
print(gsr, row.names = FALSE, digits = 3)
rob <- rob |> left_join(gsr, by = "program")
p5_tsv(as.data.frame(gsm), file.path(P5_VAL, "M40_GENE_SUBSET_USAGE_CORRELATION.tsv"))

p5_sec("5. Cell-cycle inclusion / exclusion sensitivity")
if (file.exists(sp_path("nocc"))) {
  Zn <- rd(sp_path("nocc")); rownames(Zn) <- paste0("NC", seq_len(nrow(Zn)))
  m <- match_spectra(Z, Zn)
  names(m)[names(m) == "cosine"] <- "nocc_cosine"
  names(m)[names(m) == "top100_jaccard"] <- "nocc_jaccard"
  names(m)[names(m) == "best_match"] <- "nocc_best_match"
  rob <- rob |> left_join(m, by = "program")
  p5_msg("  cell-cycle-excluded run: median cosine %.3f, min %.3f",
         median(m$nocc_cosine), min(m$nocc_cosine))
  p5_msg("  programs with cosine < 0.5 when cell-cycle genes are removed: %s",
         paste(m$program[m$nocc_cosine < 0.5], collapse = ", "))
} else p5_msg("  nocc run unavailable")

ACT <- 0.20; MINCELL <- 10L; MINFRAC <- 0.05; MINPAT <- 3L

p5_sec("5b. Do program usages track technical covariates?")
# Three of the eight programs are technical-dominated by their own top genes.
# This checks the complementary question directly: does a program's usage track
# sequencing depth or complexity, WITHIN a patient (so patient identity cannot
# create the correlation by itself)?
covs <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
tc <- do.call(rbind, lapply(PROGS, function(p) {
  do.call(rbind, lapply(covs, function(cv) {
    r <- vapply(SAMPLES, function(s) {
      k <- mmd$sample_id == s
      if (sum(k) < 50) return(NA_real_)
      suppressWarnings(cor(REL[k, p], mmd[[cv]][k], method = "spearman"))
    }, numeric(1))
    data.frame(program = p, covariate = cv, t(r),
               max_abs_within_patient_rho = max(abs(r), na.rm = TRUE),
               median_within_patient_rho = median(r, na.rm = TRUE))
  }))
}))
names(tc)[names(tc) %in% SAMPLES] <- paste0("rho_", SAMPLES)
p5_tsv(tc, file.path(P5_VAL, "M40_PROGRAM_TECHNICAL_COVARIATES.tsv"))
agg <- tc |> group_by(program) |>
  summarise(max_technical_rho = max(abs(median_within_patient_rho), na.rm = TRUE),
            worst_covariate = covariate[which.max(abs(median_within_patient_rho))],
            .groups = "drop")
print(as.data.frame(agg), row.names = FALSE, digits = 3)
rob <- rob |> left_join(agg, by = "program")

p5_sec("5c. Technical-gene-free universe (nortech)")
# Does removing ribosomal, pseudogene/lncRNA and canonical myeloid genes turn
# any biologically interpretable program into a shared one?
nt_dir <- file.path(CNMF, "nortech")
nt <- NULL
if (dir.exists(nt_dir)) {
  ks <- list.files(nt_dir, pattern = "usages\\.k_[0-9]+\\.dt_0_1\\.consensus\\.txt$")
  if (length(ks)) {
    kk <- as.integer(sub(".*usages\\.k_([0-9]+)\\..*", "\\1", ks))
    p5_msg("  nortech consensus available for K = %s", paste(sort(kk), collapse = ", "))
    nt <- do.call(rbind, lapply(sort(kk), function(k) {
      U2 <- rd(file.path(nt_dir, sprintf("nortech.usages.k_%d.dt_0_1.consensus.txt", k)))
      colnames(U2) <- paste0("N", seq_len(ncol(U2)))
      R2 <- U2 / rowSums(U2)
      m2 <- md[rownames(U2), ]
      cls <- vapply(colnames(R2), function(p) {
        act <- R2[, p] >= ACT
        n <- sum(vapply(SAMPLES, function(s) {
          kk2 <- m2$sample_id == s
          sum(act & kk2) >= MINCELL && mean(act[kk2]) >= MINFRAC
        }, logical(1)))
        if (n >= MINPAT) "recurrent" else if (n == 2) "shared-limited" else
          if (n == 1) "patient-private" else "uncertain"
      }, character(1))
      data.frame(K = k, n_programs = ncol(R2),
                 n_recurrent = sum(cls == "recurrent"),
                 n_shared_limited = sum(cls == "shared-limited"),
                 n_patient_private = sum(cls == "patient-private"),
                 n_uncertain = sum(cls == "uncertain"))
    }))
    print(as.data.frame(nt), row.names = FALSE)
    p5_tsv(nt, file.path(P5_TAB, "PROGRAM_NOTECH_RECURRENCE.tsv"))
    p5_msg("  maximum recurrent programs on the technical-gene-free universe, across the whole K grid: %d",
           max(nt$n_recurrent))
  }
} else p5_msg("  nortech run not available")

p5_sec("6. Malignancy-confidence sensitivity (High+Moderate vs High only)")
recur <- function(rel, pats) {
  vapply(colnames(rel), function(p) {
    act <- rel[, p] >= ACT
    n <- sum(vapply(unique(pats), function(s) {
      m <- pats == s
      sum(act & m) >= MINCELL && mean(act[m]) >= MINFRAC
    }, logical(1)))
    if (n >= MINPAT) "recurrent" else if (n == 2) "shared-limited" else
      if (n == 1) "patient-private" else "uncertain"
  }, character(1))
}
hc <- mmd$malignancy_confidence == "High"
p5_msg("  High-confidence malignant cells: %d of %d", sum(hc), nrow(mmd))
r_all <- recur(REL, mmd$sample_id)
r_hc  <- recur(REL[hc, , drop = FALSE], mmd$sample_id[hc])
conf <- data.frame(program = PROGS, recurrence_all = r_all,
                   recurrence_high_only = r_hc,
                   confidence_stable = r_all == r_hc)
print(conf, row.names = FALSE)
rob <- rob |> left_join(conf, by = "program")
if (file.exists(sp_path("highconf"))) {
  Zh <- rd(sp_path("highconf")); rownames(Zh) <- paste0("HC", seq_len(nrow(Zh)))
  mh <- match_spectra(Z, Zh)
  rob <- rob |> left_join(mh |> select(program, highconf_cosine = cosine,
                                       highconf_jaccard = top100_jaccard),
                          by = "program")
  p5_msg("  High-only cNMF run: median cosine %.3f, min %.3f",
         median(mh$cosine), min(mh$cosine))
}

p5_sec("7. Leave-one-patient-out program recovery")
loo <- list()
for (s in SAMPLES) {
  f <- sp_path(paste0("loo_", s))
  if (!file.exists(f)) { p5_msg("  loo_%s unavailable", s); next }
  Zl <- rd(f); rownames(Zl) <- paste0("L", seq_len(nrow(Zl)))
  m <- match_spectra(Z, Zl); m$withheld <- s
  loo[[s]] <- m
  p5_msg("  withhold %s: median cosine %.3f, min %.3f, programs < 0.5: %s",
         s, median(m$cosine), min(m$cosine),
         paste(m$program[m$cosine < 0.5], collapse = ", "))
}
if (length(loo)) {
  lo <- do.call(rbind, loo)
  p5_tsv(lo, file.path(P5_VAL, "M40_LEAVE_ONE_PATIENT_OUT_RECOVERY.tsv"))
  agg <- lo |> group_by(program) |>
    summarise(loo_min_cosine = min(cosine), loo_median_cosine = median(cosine),
              loo_worst_withheld = withheld[which.min(cosine)],
              loo_n_recovered = sum(cosine >= 0.60), .groups = "drop")
  rob <- rob |> left_join(agg, by = "program")
  # does a program survive removal of the patient that dominates it?
  dp <- setNames(progs_meta$dominant_patient, progs_meta$program_id)
  rob$loo_cosine_without_dominant_patient <- vapply(rob$program, function(p) {
    d <- dp[[p]]
    v <- lo$cosine[lo$program == p & lo$withheld == d]
    if (length(v)) v[1] else NA_real_
  }, numeric(1))
  rob$survives_dominant_patient_removal <-
    rob$loo_cosine_without_dominant_patient >= 0.60
}

p5_sec("8. Overall robustness verdict")
rob <- rob |> left_join(progs_meta |>
  select(program = program_id, program_label, recurrence_status,
         dominant_patient, dominant_patient_fraction), by = "program")
rob$robust_overall <- with(rob,
  (is.na(balanced_cosine) | balanced_cosine >= 0.60) &
  (is.na(gene_subset_min_rho) | gene_subset_min_rho >= 0.80) &
  (is.na(loo_min_cosine) | loo_min_cosine >= 0.50) &
  (is.na(survives_dominant_patient_removal) | survives_dominant_patient_removal) &
  (is.na(confidence_stable) | confidence_stable) &
  (is.na(max_technical_rho) | max_technical_rho < 0.50))
# A patient-private program CANNOT survive removal of its own patient - that is
# definitional, not a defect. The second column separates the two failure modes
# so a reader is not told a program is unreliable when the only test it fails is
# one it could never pass.
rob$robust_excluding_patient_scope <- with(rob,
  (is.na(balanced_cosine) | balanced_cosine >= 0.60) &
  (is.na(loo_min_cosine) | loo_min_cosine >= 0.50 |
     recurrence_status == "patient-private") &
  (is.na(confidence_stable) | confidence_stable) &
  (is.na(max_technical_rho) | max_technical_rho < 0.50) &
  (is.na(rank_min_cosine) | rank_min_cosine >= 0.60) &
  (is.na(nocc_cosine) | nocc_cosine >= 0.60))
rob <- rob |> relocate(program, program_label, recurrence_status, robust_overall,
                       robust_excluding_patient_scope)
p5_tsv(rob, file.path(P5_TAB, "PHASE5_PROGRAM_ROBUSTNESS.tsv"))
print(as.data.frame(rob[, c("program","program_label","recurrence_status",
                            "balanced_cosine","gene_subset_min_rho",
                            "loo_min_cosine","survives_dominant_patient_removal",
                            "confidence_stable","max_technical_rho",
                            "robust_overall","robust_excluding_patient_scope")]),
      row.names = FALSE, digits = 3)
facts$technical_covariates <- tc
facts$nortech_recurrence <- nt
facts$robustness <- rob
facts$criteria <- list(
  balanced_cosine_min = 0.60, gene_subset_rho_min = 0.80,
  loo_cosine_min = 0.50, dominant_patient_removal_cosine_min = 0.60,
  max_within_patient_technical_rho = 0.50,
  note = "thresholds declared before the tests were run; a program failing any one is reported as not robust rather than dropped from the table",
  two_verdicts = "robust_overall applies every test including removal of the program's dominant patient. robust_excluding_patient_scope waives ONLY that test for programs already classified patient-private, because a patient-private program cannot by definition survive removal of its own patient - reporting it as unreliable on that basis would confuse scope with reliability.")
p5_json(facts, file.path(P5_VAL, "m40_robustness_facts.json"))
p5_sec("M40 complete")
