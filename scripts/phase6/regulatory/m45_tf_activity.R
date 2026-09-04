#!/usr/bin/env Rscript
# =============================================================================
# Phase 6 - M45 - transcription-factor / regulon activity.
#
# decoupleR + CollecTRI, both ALREADY present in the frozen R_env, so no
# package is installed for this step. pySCENIC is deliberately NOT introduced:
# it would destabilise the frozen stack for no gain on the question asked here.
#
# Input is gene expression (the frozen RNA log-normalised layer) - never
# Harmony, never UMAP, never the SCEVAN CNA matrix.
#
# A TF-program relationship is NOT called recurrent because a pooled cell-level
# association is strong. It is reported per patient, with the number of
# patients showing a concordant direction and the dominant-patient fraction.
# =============================================================================
suppressPackageStartupMessages({ library(decoupleR); library(Matrix) })
source("scripts/phase6/utils/phase6_common.R")
set.seed(42)
facts <- list(milestone = "M45", generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))

NET_CACHE <- "external/networks"
dir.create(NET_CACHE, showWarnings = FALSE, recursive = TRUE)
MIN_TARGETS <- 5L

p6_sec("1. Regulatory network: CollecTRI")
# decoupleR::get_collectri() and OmnipathR::collectri() both fail in the frozen
# stack: OmnipathR 3.14.0's unnest_evidences() errors on the CollecTRI static
# table ("if (.keep) . else select(., -!!evs_col)"). The fix is NOT to upgrade
# OmnipathR - that would move a frozen dependency for one convenience wrapper.
# The same data is fetched from OmniPath's documented REST endpoint and cached
# with its provenance, so the network is identical and fully inspectable.
ctri_raw <- file.path(NET_CACHE, "collectri_omnipath_raw.tsv")
ctri_f   <- file.path(NET_CACHE, "collectri_human.tsv")
CTRI_URL <- paste0("https://omnipathdb.org/interactions?datasets=collectri",
                   "&genesymbols=yes&organisms=9606&fields=sources,references,curation_effort")
if (!file.exists(ctri_f)) {
  if (!file.exists(ctri_raw)) {
    p6_msg("fetching CollecTRI from the OmniPath REST endpoint")
    ok <- utils::download.file(CTRI_URL, ctri_raw, quiet = TRUE)
    stopifnot(ok == 0, file.exists(ctri_raw))
  }
  raw <- read.delim(ctri_raw, check.names = FALSE, quote = "")
  p6_msg("raw CollecTRI interactions: %d", nrow(raw))
  # mor follows the CollecTRI convention: -1 only for purely repressive edges.
  net <- data.frame(
    source = raw$source_genesymbol, target = raw$target_genesymbol,
    mor = ifelse(raw$is_stimulation %in% c("True", TRUE), 1,
          ifelse(raw$is_inhibition %in% c("True", TRUE), -1, 1)),
    stringsAsFactors = FALSE)
  net <- net[net$source != "" & net$target != "" &
               !grepl("COMPLEX", net$source) & !grepl("COMPLEX", net$target), ]
  net <- net[!duplicated(paste(net$source, net$target)), ]
  write.table(net, ctri_f, sep = "\t", quote = FALSE, row.names = FALSE)
  p6_msg("CollecTRI built and cached: %s", ctri_f)
} else {
  net <- read.delim(ctri_f, check.names = FALSE)
  p6_msg("CollecTRI loaded from cache: %s", ctri_f)
}
p6_msg("CollecTRI: %d edges, %d TFs, %d targets (%d activating, %d repressing)",
       nrow(net), length(unique(net$source)), length(unique(net$target)),
       sum(net$mor > 0), sum(net$mor < 0))
facts$network <- list(
  resource = "CollecTRI, fetched from the OmniPath REST endpoint",
  url = CTRI_URL, cache = ctri_f, raw_cache = ctri_raw,
  workaround = "decoupleR::get_collectri() and OmnipathR::collectri() both error in OmnipathR 3.14.0 (unnest_evidences); OmnipathR was NOT upgraded - a frozen dependency is not moved for a convenience wrapper - and the identical data was taken from the documented REST endpoint instead",
  decoupleR = as.character(packageVersion("decoupleR")),
  OmnipathR = as.character(packageVersion("OmnipathR")),
  edges = nrow(net), tfs = length(unique(net$source)),
  min_targets = MIN_TARGETS, method = "run_ulm (univariate linear model)")

p6_sec("2. Expression input")
lg <- readRDS(file.path(P6_CLONE, "phase6_malignant_lognorm.rds"))
p6_msg("frozen RNA log-normalised data: %d genes x %d malignant cells",
       nrow(lg), ncol(lg))
keep <- Matrix::rowSums(lg > 0) >= 0.01 * ncol(lg)
lg <- lg[keep, , drop = FALSE]
p6_msg("genes detected in >= 1%% of malignant cells: %d", nrow(lg))
net <- net[net$target %in% rownames(lg), ]
tf_n <- table(net$source)
net <- net[net$source %in% names(tf_n)[tf_n >= MIN_TARGETS], ]
p6_msg("network after gene filtering: %d edges, %d TFs with >= %d targets",
       nrow(net), length(unique(net$source)), MIN_TARGETS)

p6_sec("3. TF activity per cell (decoupleR run_ulm)")
t0 <- Sys.time()
acts <- decoupleR::run_ulm(mat = as.matrix(lg), network = net,
                           .source = "source", .target = "target",
                           .mor = "mor", minsize = MIN_TARGETS)
p6_msg("run_ulm finished in %.1f min; %d rows",
       as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(acts))
A <- acts |> filter(statistic == "ulm") |>
  select(source, condition, score) |>
  pivot_wider(names_from = source, values_from = score) |>
  as.data.frame()
rownames(A) <- A$condition; A$condition <- NULL
A <- as.matrix(A)[colnames(lg), , drop = FALSE]
p6_msg("TF activity matrix: %d cells x %d TFs", nrow(A), ncol(A))
saveRDS(A, file.path(P6_REG, "m45_tf_activity_matrix.rds"))

p6_sec("4. TF activity versus each malignant program, per patient")
md <- readRDS(file.path(P6_CLONE, "phase6_cell_metadata.rds"))
mm <- md[rownames(A), ]
PROGS <- grep("^program_P[0-9]+_score$", colnames(mm), value = TRUE)
PID <- sub("^program_(P[0-9]+)_score$", "\\1", PROGS)
lab <- read.delim("results/phase5/tables/final/MALIGNANT_PROGRAMS.tsv")
LABEL <- setNames(lab$program_label, lab$program_id)

# per-patient Spearman correlation between TF activity and program usage.
# Correlation is an effect size; patient concordance is the evidence.
res <- list()
for (i in seq_along(PROGS)) {
  y_all <- mm[[PROGS[i]]]
  per <- lapply(SAMPLES, function(s) {
    k <- mm$sample_id == s
    if (sum(k) < 50) return(setNames(rep(NA_real_, ncol(A)), colnames(A)))
    suppressWarnings(cor(A[k, , drop = FALSE], y_all[k], method = "spearman")[, 1])
  })
  names(per) <- SAMPLES
  P <- do.call(cbind, per)
  pooled <- suppressWarnings(cor(A, y_all, method = "spearman")[, 1])
  sgn <- sign(P); ok <- !is.na(P)
  n_conc_pos <- rowSums(sgn > 0 & ok); n_conc_neg <- rowSums(sgn < 0 & ok)
  res[[i]] <- data.frame(
    program = PID[i], program_label = unname(LABEL[PID[i]]),
    tf = colnames(A), pooled_rho = pooled,
    as.data.frame(P), n_patients_evaluable = rowSums(ok),
    n_patients_positive = n_conc_pos, n_patients_negative = n_conc_neg,
    n_patients_concordant = pmax(n_conc_pos, n_conc_neg),
    direction = ifelse(pooled > 0, "positive", "negative"),
    min_abs_rho_across_patients = apply(abs(P), 1, function(v)
      if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE)),
    row.names = NULL)
}
tfp <- do.call(rbind, res)
names(tfp)[names(tfp) %in% SAMPLES] <- paste0("rho_", SAMPLES)
# recurrent = concordant direction in >= 3 evaluable patients AND a pooled
# effect that is not trivial
tfp$recurrent_direction <- tfp$n_patients_concordant >= 3 & abs(tfp$pooled_rho) >= 0.20
p6_tsv(tfp[order(tfp$program, -abs(tfp$pooled_rho)), ],
       file.path(P6_TAB, "PROGRAM_TF_ACTIVITY.tsv"))

p6_sec("5. Top regulators per program")
for (p in PID) {
  d <- tfp[tfp$program == p, ]
  up <- d[order(-d$pooled_rho), ][1:8, ]
  dn <- d[order(d$pooled_rho), ][1:8, ]
  p6_msg("  %-4s %s", p, LABEL[p])
  p6_msg("     positive: %s", paste(sprintf("%s(%.2f%s)", up$tf, up$pooled_rho,
         ifelse(up$recurrent_direction, "*", "")), collapse = " "))
  p6_msg("     negative: %s", paste(sprintf("%s(%.2f%s)", dn$tf, dn$pooled_rho,
         ifelse(dn$recurrent_direction, "*", "")), collapse = " "))
}
p6_msg("* = concordant direction in >= 3 of the evaluable patients AND |pooled rho| >= 0.20")

facts$n_tfs <- ncol(A); facts$n_cells <- nrow(A)
facts$recurrence_rule <- "a TF-program relationship is called recurrent only when its direction is concordant in >= 3 evaluable patients AND |pooled rho| >= 0.20; a strong pooled association alone is never sufficient"
facts$n_recurrent_pairs <- sum(tfp$recurrent_direction, na.rm = TRUE)
p6_json(facts, file.path(P6_VAL, "m45_tf_activity_facts.json"))
p6_sec("M45 complete")
