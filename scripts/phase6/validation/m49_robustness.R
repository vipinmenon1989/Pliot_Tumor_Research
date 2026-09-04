#!/usr/bin/env Rscript
# =============================================================================
# Phase 6 - M49 - robustness.
#
#   1 clone-size threshold sensitivity        (>= 20 / 50 / 100 cells)
#   2 high-confidence malignant only          (vs High + Moderate)
#   3 dominant-program margin sensitivity
#   4 TF and pathway conclusions under the high-confidence restriction
#   5 patient-level direction consistency
#
# Leave-one-patient-out is deliberately NOT applied to the clone analyses: they
# are inherently within-patient, so withholding a patient deletes the analysis
# rather than testing it. Patient-level direction consistency is tested instead.
# =============================================================================
source("scripts/phase6/utils/phase6_common.R")
set.seed(42)
facts <- list(milestone = "M49", generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))

md <- readRDS(file.path(P6_CLONE, "phase6_cell_metadata.rds"))
PROGS <- grep("^program_P[0-9]+_score$", colnames(md), value = TRUE)
PID <- sub("^program_(P[0-9]+)_score$", "\\1", PROGS)
lab <- read.delim("results/phase5/tables/final/MALIGNANT_PROGRAMS.tsv")
LABEL <- setNames(lab$program_label, lab$program_id)
mal <- md[md$malignancy_refined == "Malignant" & !is.na(md$tumor_clone_phase4) &
            md$sample_id %in% CLONE_RELIABLE, ]
p6_msg("clone-reliable malignant cells with a clone: %d", nrow(mal))

eta2 <- function(y, g) {
  gm <- tapply(y, g, mean); n <- table(g)
  sst <- sum((y - mean(y))^2); if (sst <= 0) return(NA_real_)
  sum(n * (gm - mean(y))^2) / sst
}
run_assoc <- function(d, minsize) {
  out <- lapply(CLONE_RELIABLE, function(s) {
    x <- d[d$sample_id == s, ]
    keep <- names(which(table(x$tumor_clone_phase4) >= minsize))
    x <- x[x$tumor_clone_phase4 %in% keep, ]
    if (length(keep) < 2 || nrow(x) < 20) return(NULL)
    g <- factor(x$tumor_clone_phase4)
    do.call(rbind, lapply(seq_along(PROGS), function(i)
      data.frame(sample_id = s, program = PID[i], n_clones = nlevels(g),
                 n_cells = nrow(x), eta2 = eta2(x[[PROGS[i]]], g))))
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (!length(out)) NULL else do.call(rbind, out)
}

p6_sec("1. Clone-size threshold sensitivity")
cs <- do.call(rbind, lapply(CLONE_SIZE_GRID, function(th) {
  a <- run_assoc(mal, th)
  ncl <- sum(vapply(CLONE_RELIABLE, function(s)
    sum(table(mal$tumor_clone_phase4[mal$sample_id == s]) >= th), integer(1)))
  nne <- sum(vapply(CLONE_RELIABLE, function(s)
    sum(table(mal$tumor_clone_phase4[mal$sample_id == s]) < th), integer(1)))
  if (is.null(a)) return(data.frame(min_clone_size = th, n_pairs = 0L,
    median_eta2 = NA_real_, max_eta2 = NA_real_, n_clones_used = ncl,
    n_clones_not_evaluable = nne))
  data.frame(min_clone_size = th, n_pairs = nrow(a),
             median_eta2 = median(a$eta2, na.rm = TRUE),
             max_eta2 = max(a$eta2, na.rm = TRUE),
             n_clones_used = ncl, n_clones_not_evaluable = nne)
}))
print(as.data.frame(cs), row.names = FALSE, digits = 3)

p6_sec("2. High-confidence malignant only")
hc <- mal[mal$malignancy_confidence == "High", ]
p6_msg("High-confidence cells: %d of %d (%.0f%%)", nrow(hc), nrow(mal),
       100 * nrow(hc) / nrow(mal))
print(table(hc$sample_id))
a_all <- run_assoc(mal, CLONE_SIZE_PRIMARY)
a_hc  <- run_assoc(hc,  CLONE_SIZE_PRIMARY)
conf <- NULL
if (!is.null(a_hc) && !is.null(a_all)) {
  conf <- merge(a_all[, c("sample_id", "program", "eta2")],
                a_hc[, c("sample_id", "program", "eta2")],
                by = c("sample_id", "program"), suffixes = c("_all", "_high"))
  conf$program_label <- unname(LABEL[conf$program])
  conf$delta <- conf$eta2_high - conf$eta2_all
  print(as.data.frame(conf), row.names = FALSE, digits = 3)
  p6_msg("spearman(eta^2 all, eta^2 High-only) over %d patient x program pairs: %.3f",
         nrow(conf), suppressWarnings(cor(conf$eta2_all, conf$eta2_high,
                                          use = "complete.obs", method = "spearman")))
} else p6_msg("the High-only subset leaves too few evaluable clones in some patient")

p6_sec("3. Dominant-program margin sensitivity")
dt <- do.call(rbind, lapply(c(0.00, 0.05, 0.10, 0.15, 0.20), function(mg) {
  d <- mal
  d$hard <- ifelse(d$dominant_program_margin >= mg,
                   d$dominant_malignant_program, "Mixed")
  per <- do.call(rbind, lapply(split(seq_len(nrow(d)), d$tumor_clone_phase4),
    function(i) {
      x <- d[i, ]; hp <- x$hard[x$hard != "Mixed"]
      data.frame(clone = x$tumor_clone_phase4[1], n_cells = nrow(x),
                 eff_n = if (length(hp)) effective_n(as.numeric(prop.table(table(hp)))) else NA_real_,
                 frac_mixed = mean(x$hard == "Mixed"))
    }))
  per <- per[per$n_cells >= CLONE_SIZE_PRIMARY, ]
  data.frame(margin_threshold = mg, n_clones = nrow(per),
             median_effective_n = median(per$eff_n, na.rm = TRUE),
             frac_clones_multi = mean(per$eff_n >= 2, na.rm = TRUE),
             median_frac_mixed = median(per$frac_mixed))
}))
print(as.data.frame(dt), row.names = FALSE, digits = 3)

p6_sec("4. TF and pathway conclusions under the high-confidence restriction")
tf <- read.delim(file.path(P6_TAB, "PROGRAM_TF_ACTIVITY.tsv"))
pw <- read.delim(file.path(P6_TAB, "PROGRAM_PATHWAY_ACTIVITY.tsv"))
A <- readRDS(file.path(P6_REG, "m45_tf_activity_matrix.rds"))
P <- readRDS(file.path(P6_PATH, "m46_progeny_activity_matrix.rds"))
mm <- md[rownames(A), ]
hcm <- mm$malignancy_confidence == "High"
p6_msg("High-confidence cells in the activity matrices: %d of %d", sum(hcm), nrow(A))
recompute <- function(M, mask) {
  do.call(rbind, lapply(seq_along(PROGS), function(i) {
    y <- mm[[PROGS[i]]]
    data.frame(program = PID[i], feature = colnames(M),
               rho_all = suppressWarnings(cor(M, y, method = "spearman")[, 1]),
               rho_high = suppressWarnings(cor(M[mask, , drop = FALSE], y[mask],
                                               method = "spearman")[, 1]))
  }))
}
tfr <- recompute(A, hcm); tfr$same_sign <- sign(tfr$rho_all) == sign(tfr$rho_high)
pwr <- recompute(P, hcm); pwr$same_sign <- sign(pwr$rho_all) == sign(pwr$rho_high)
sumr <- rbind(
  data.frame(layer = "TF (CollecTRI)", n_features = nrow(tfr),
             spearman_all_vs_high = cor(tfr$rho_all, tfr$rho_high, use = "complete.obs"),
             frac_same_sign = mean(tfr$same_sign, na.rm = TRUE),
             frac_same_sign_strong = mean(tfr$same_sign[abs(tfr$rho_all) >= 0.20], na.rm = TRUE)),
  data.frame(layer = "PROGENy", n_features = nrow(pwr),
             spearman_all_vs_high = cor(pwr$rho_all, pwr$rho_high, use = "complete.obs"),
             frac_same_sign = mean(pwr$same_sign, na.rm = TRUE),
             frac_same_sign_strong = mean(pwr$same_sign[abs(pwr$rho_all) >= 0.20], na.rm = TRUE)))
print(as.data.frame(sumr), row.names = FALSE, digits = 3)

p6_sec("5. Patient-level direction consistency")
pc <- rbind(
  tf |> transmute(layer = "TF (CollecTRI)", n_patients_concordant, pooled_rho,
                  recurrent_direction),
  pw |> transmute(layer = paste0("Pathway (", layer, ")"), n_patients_concordant,
                  pooled_rho, recurrent_direction)) |>
  group_by(layer) |>
  summarise(n_pairs = n(),
            n_strong_pooled = sum(abs(pooled_rho) >= 0.20, na.rm = TRUE),
            n_strong_and_concordant_ge3 = sum(recurrent_direction, na.rm = TRUE),
            frac_strong_that_are_concordant = sum(recurrent_direction, na.rm = TRUE) /
              max(1, sum(abs(pooled_rho) >= 0.20, na.rm = TRUE)), .groups = "drop")
print(as.data.frame(pc), row.names = FALSE, digits = 3)
p6_msg("A strong pooled association that is NOT concordant across >= 3 patients is")
p6_msg("exactly what n = 4 with sample_id = patient = dataset produces by accident.")

p6_sec("6. PHASE6_ROBUSTNESS.tsv")
rob <- rbind(
  cs |> transmute(test = "clone-size threshold",
                  setting = paste0(">= ", min_clone_size, " cells"),
                  metric = "median eta^2", value = median_eta2, n = n_clones_used,
                  note = paste0(n_clones_not_evaluable, " clones NOT EVALUABLE at this threshold")),
  dt |> transmute(test = "dominant-program margin",
                  setting = paste0("margin >= ", margin_threshold),
                  metric = "median effective programs per clone",
                  value = median_effective_n, n = n_clones,
                  note = sprintf("%.0f%% of clones span >= 2 programs; median %.0f%% of cells Mixed",
                                 100 * frac_clones_multi, 100 * median_frac_mixed)),
  sumr |> transmute(test = "high-confidence malignant only", setting = layer,
                    metric = "spearman(rho_all, rho_high)",
                    value = spearman_all_vs_high, n = n_features,
                    note = sprintf("%.0f%% of strong associations keep their sign",
                                   100 * frac_same_sign_strong)),
  pc |> transmute(test = "patient-level direction consistency", setting = layer,
                  metric = "fraction of strong pooled associations concordant in >= 3 patients",
                  value = frac_strong_that_are_concordant, n = n_pairs,
                  note = sprintf("%d strong pooled, %d also concordant",
                                 n_strong_pooled, n_strong_and_concordant_ge3)))
if (!is.null(conf))
  rob <- rbind(rob, data.frame(test = "high-confidence malignant only",
    setting = "clone-program eta^2", metric = "spearman(eta2_all, eta2_high)",
    value = suppressWarnings(cor(conf$eta2_all, conf$eta2_high,
                                 use = "complete.obs", method = "spearman")),
    n = nrow(conf), note = "per patient x program"))
p6_tsv(rob, file.path(P6_TAB, "PHASE6_ROBUSTNESS.tsv"))
print(as.data.frame(rob), row.names = FALSE, digits = 3)

facts$clone_size_sensitivity <- cs
facts$dominant_margin_sensitivity <- dt
facts$high_confidence <- sumr
facts$patient_consistency <- pc
facts$eta2_high_vs_all <- conf
facts$loo_note <- "leave-one-patient-out is deliberately NOT applied to the clone analyses: they are inherently within-patient, so withholding a patient deletes the analysis instead of testing it"
p6_json(facts, file.path(P6_VAL, "m49_robustness_facts.json"))
p6_sec("M49 complete")
