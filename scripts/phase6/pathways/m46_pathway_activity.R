#!/usr/bin/env Rscript
# =============================================================================
# Phase 6 - M46 - pathway activity across malignant cells.
#
# Two independent, transparent layers:
#   PROGENy (footprint-based, via decoupleR)  - responsive pathway activity
#   MSigDB Hallmark (gene-set scoring)        - programme-level gene sets
#
# They are reported side by side. They are NOT merged into an invented
# composite score and NOT collapsed into a single numeric rank - the brief
# forbids it and the two layers answer different questions.
#
# As in M45, patient-level direction and concordance carry the evidence, not a
# pooled cell-level p-value.
# =============================================================================
suppressPackageStartupMessages({ library(decoupleR); library(Matrix) })
source("scripts/phase6/utils/phase6_common.R")
set.seed(42)
facts <- list(milestone = "M46", generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
NET_CACHE <- "external/networks"; dir.create(NET_CACHE, showWarnings = FALSE, recursive = TRUE)

p6_sec("1. PROGENy network")
prg_f <- file.path(NET_CACHE, "progeny_human_top500.tsv")
if (file.exists(prg_f)) {
  prg <- read.delim(prg_f, check.names = FALSE)
  p6_msg("PROGENy loaded from cache")
} else {
  prg <- as.data.frame(decoupleR::get_progeny(organism = "human", top = 500))
  write.table(prg, prg_f, sep = "\t", quote = FALSE, row.names = FALSE)
  p6_msg("PROGENy downloaded and cached: %s", prg_f)
}
p6_msg("PROGENy: %d edges across %d pathways", nrow(prg), length(unique(prg$source)))
p6_msg("pathways: %s", paste(sort(unique(prg$source)), collapse = ", "))

p6_sec("2. Expression input")
lg <- readRDS(file.path(P6_CLONE, "phase6_malignant_lognorm.rds"))
keep <- Matrix::rowSums(lg > 0) >= 0.01 * ncol(lg)
lg <- lg[keep, , drop = FALSE]
p6_msg("genes: %d | malignant cells: %d", nrow(lg), ncol(lg))

p6_sec("3. PROGENy activity per cell (run_mlm)")
pa <- decoupleR::run_mlm(mat = as.matrix(lg), network = prg, .source = "source",
                         .target = "target", .mor = "weight", minsize = 5)
P <- pa |> filter(statistic == "mlm") |> select(source, condition, score) |>
  pivot_wider(names_from = source, values_from = score) |> as.data.frame()
rownames(P) <- P$condition; P$condition <- NULL
P <- as.matrix(P)[colnames(lg), , drop = FALSE]
p6_msg("PROGENy activity matrix: %d cells x %d pathways", nrow(P), ncol(P))
saveRDS(P, file.path(P6_PATH, "m46_progeny_activity_matrix.rds"))

p6_sec("4. Hallmark gene-set scores per cell")
read_gmt <- function(f) {
  ln <- strsplit(readLines(f), "\t")
  setNames(lapply(ln, function(x) unique(x[-c(1, 2)])), vapply(ln, `[`, "", 1))
}
hm <- read_gmt("external/genesets/h.all.v2024.1.Hs.symbols.gmt")
hm <- lapply(hm, function(g) intersect(g, rownames(lg)))
hm <- hm[lengths(hm) >= 10]
p6_msg("Hallmark sets with >= 10 detected genes: %d", length(hm))
# mean z of the set's genes, so a score is comparable across sets
Z <- as.matrix(lg)
Z <- (Z - rowMeans(Z)) / (matrixStats::rowSds(Z) + 1e-9)
H <- vapply(hm, function(g) colMeans(Z[g, , drop = FALSE]), numeric(ncol(Z)))
rownames(H) <- colnames(lg)
p6_msg("Hallmark score matrix: %d cells x %d sets", nrow(H), ncol(H))
saveRDS(H, file.path(P6_PATH, "m46_hallmark_score_matrix.rds"))
rm(Z); invisible(gc())

p6_sec("5. Pathway activity versus each malignant program, per patient")
md <- readRDS(file.path(P6_CLONE, "phase6_cell_metadata.rds"))
mm <- md[colnames(lg), ]
PROGS <- grep("^program_P[0-9]+_score$", colnames(mm), value = TRUE)
PID <- sub("^program_(P[0-9]+)_score$", "\\1", PROGS)
lab <- read.delim("results/phase5/tables/final/MALIGNANT_PROGRAMS.tsv")
LABEL <- setNames(lab$program_label, lab$program_id)

assoc_layer <- function(M, layer) {
  do.call(rbind, lapply(seq_along(PROGS), function(i) {
    y <- mm[[PROGS[i]]]
    per <- vapply(SAMPLES, function(s) {
      k <- mm$sample_id == s
      if (sum(k) < 50) return(setNames(rep(NA_real_, ncol(M)), colnames(M)))
      suppressWarnings(cor(M[k, , drop = FALSE], y[k], method = "spearman")[, 1])
    }, numeric(ncol(M)))
    pooled <- suppressWarnings(cor(M, y, method = "spearman")[, 1])
    ok <- !is.na(per)
    npos <- rowSums(per > 0 & ok); nneg <- rowSums(per < 0 & ok)
    out <- data.frame(program = PID[i], program_label = unname(LABEL[PID[i]]),
      layer = layer, pathway = colnames(M), pooled_rho = pooled,
      as.data.frame(per), n_patients_evaluable = rowSums(ok),
      n_patients_concordant = pmax(npos, nneg),
      direction = ifelse(pooled > 0, "positive", "negative"), row.names = NULL)
    names(out)[names(out) %in% SAMPLES] <- paste0("rho_", SAMPLES)
    out
  }))
}
pw <- rbind(assoc_layer(P, "PROGENy"), assoc_layer(H, "Hallmark"))
pw$recurrent_direction <- pw$n_patients_concordant >= 3 & abs(pw$pooled_rho) >= 0.20
pw <- pw[order(pw$program, pw$layer, -abs(pw$pooled_rho)), ]
p6_tsv(pw, file.path(P6_TAB, "PROGRAM_PATHWAY_ACTIVITY.tsv"))

for (p in PID) {
  d <- pw[pw$program == p & pw$layer == "PROGENy", ]
  d <- d[order(-abs(d$pooled_rho)), ][1:5, ]
  p6_msg("  %-4s %-28s PROGENy: %s", p, LABEL[p],
         paste(sprintf("%s(%.2f%s)", d$pathway, d$pooled_rho,
               ifelse(d$recurrent_direction, "*", "")), collapse = " "))
  d2 <- pw[pw$program == p & pw$layer == "Hallmark", ]
  d2 <- d2[order(-d2$pooled_rho), ][1:4, ]
  p6_msg("       %-28s Hallmark: %s", "",
         paste(sprintf("%s(%.2f%s)", sub("^HALLMARK_", "", d2$pathway),
               d2$pooled_rho, ifelse(d2$recurrent_direction, "*", "")),
               collapse = " "))
}

p6_sec("6. Three independent evidence layers, kept separate")
tf <- read.delim(file.path(P6_TAB, "PROGRAM_TF_ACTIVITY.tsv"))
ev <- rbind(
  tf |> filter(recurrent_direction) |>
    transmute(program, program_label, layer = "TF activity (CollecTRI)",
              feature = tf, pooled_rho, n_patients_concordant, direction),
  pw |> filter(recurrent_direction) |>
    transmute(program, program_label, layer = paste0("Pathway (", layer, ")"),
              feature = pathway, pooled_rho, n_patients_concordant, direction))
ev <- ev[order(ev$program, ev$layer, -abs(ev$pooled_rho)), ]
ev$note <- "layers are reported side by side; no composite score is computed and no single numeric rank is imposed"
p6_tsv(ev, file.path(P6_TAB, "PROGRAM_REGULATORY_EVIDENCE.tsv"))
p6_msg("cross-patient-concordant evidence rows: %d (%d TF, %d pathway)",
       nrow(ev), sum(grepl("^TF", ev$layer)), sum(grepl("^Pathway", ev$layer)))

facts$progeny <- list(pathways = ncol(P), method = "decoupleR run_mlm, PROGENy top 500")
facts$hallmark <- list(sets = ncol(H), method = "mean z-score of set genes, MSigDB v2024.1.Hs")
facts$n_recurrent_pathway_pairs <- sum(pw$recurrent_direction, na.rm = TRUE)
facts$integration_rule <- "program, TF activity and pathway activity are three independent layers; they are never merged into a weighted composite"
p6_json(facts, file.path(P6_VAL, "m46_pathway_facts.json"))
p6_sec("M46 complete")
