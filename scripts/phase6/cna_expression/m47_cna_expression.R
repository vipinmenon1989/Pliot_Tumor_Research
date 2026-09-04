#!/usr/bin/env Rscript
# =============================================================================
# Phase 6 - M47 - transcriptional consequences of broad inferred CNA events.
#
# THE LIMITATION THAT DEFINES THIS MILESTONE: SCEVAN infers copy number FROM
# EXPRESSION. An association between an inferred broad event and the expression
# of genes on the affected chromosome is therefore an INTERNAL TRANSCRIPTIONAL
# CONSISTENCY analysis - it is NOT independent validation of the copy-number
# call, and it is never written as such.
#
# Claims are restricted to chromosome / arm / large-segment level. No
# single-gene deletion or amplification is asserted anywhere. Where a locus is
# named it is described only as "a broad segment whose interval contains X".
# =============================================================================
source("scripts/phase6/utils/phase6_common.R")
set.seed(42)
facts <- list(milestone = "M47", generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
MIN_SEG_MB <- 10; MIN_CELLS <- 20L; MIN_GENES_SEG <- 30L

p6_sec("1. Frozen Phase 4 broad clonal segments")
seg <- read.delim("results/phase4/tables/final/SCEVAN_CNV_SUMMARY.tsv",
                  check.names = FALSE)
seg <- seg |> filter(level == "clonal", broad %in% c(TRUE, "TRUE"),
                     event != "neutral", width_mb >= MIN_SEG_MB,
                     Chr %in% 1:22, sample_id %in% CLONE_RELIABLE)
p6_msg("broad clonal non-neutral segments (>= %d Mb) in the reliable patients: %d",
       MIN_SEG_MB, nrow(seg))
print(seg |> count(sample_id, event) |> as.data.frame(), row.names = FALSE)
p6_msg("MPNST_3 excluded from every clone/CNA analysis: %s", CLONE_EXCLUSION_REASON)

p6_sec("2. Per-subclone segments, so an event can be tied to specific clones")
sub <- read.delim("results/phase4/tables/final/SCEVAN_CNV_SUMMARY.tsv",
                  check.names = FALSE) |>
  filter(grepl("^subclone", level), broad %in% c(TRUE, "TRUE"),
         event != "neutral", width_mb >= MIN_SEG_MB, Chr %in% 1:22,
         sample_id %in% CLONE_RELIABLE) |>
  mutate(subclone_n = as.integer(sub("^subclone", "", level)),
         clone = sprintf("%s_clone%d", sample_id, subclone_n))
p6_msg("broad subclonal segments: %d across %d clone labels", nrow(sub),
       length(unique(sub$clone)))

p6_sec("3. Gene -> chromosome map, from the stored SCEVAN annotation")
ann <- NULL
for (s in CLONE_RELIABLE) {
  f <- sprintf("results/phase4/scevan/by_sample/%s/primary/output/%s_primary_count_mtx_annot.RData", s, s)
  if (!file.exists(f)) next
  e <- new.env(); load(f, envir = e)
  a <- get("count_mtx_annot", e)
  ann <- rbind(ann, data.frame(gene = a$gene_name, chr = as.character(a$seqnames),
                               start = a$start, end = a$end))
}
ann <- ann[!duplicated(ann$gene) & ann$chr %in% as.character(1:22), ]
p6_msg("gene coordinates available for %d genes on chr1-22", nrow(ann))

p6_sec("4. Expression input")
lg <- readRDS(file.path(P6_CLONE, "phase6_malignant_lognorm.rds"))
md <- readRDS(file.path(P6_CLONE, "phase6_cell_metadata.rds"))
mm <- md[colnames(lg), ]
keep <- Matrix::rowSums(lg > 0) >= 0.05 * ncol(lg)
lg <- lg[keep, , drop = FALSE]
ann <- ann[ann$gene %in% rownames(lg), ]
p6_msg("genes retained (detected in >= 5%% of malignant cells) with coordinates: %d",
       nrow(ann))

p6_sec("5. Segment-level expression, per patient, clone-carrier vs non-carrier")
# For each broad SUBCLONAL segment, compare the mean expression of the genes in
# that interval between the clones that carry the event and the clones in the
# SAME patient that do not. Segment-level, never single-gene.
res <- list()
for (s in CLONE_RELIABLE) {
  cells <- mm$cell_id[mm$sample_id == s & !is.na(mm$tumor_clone_phase4)]
  if (!length(cells)) next
  cl <- mm$tumor_clone_phase4[match(cells, mm$cell_id)]
  big <- names(which(table(cl) >= MIN_CELLS))
  ss <- sub[sub$sample_id == s, ]
  # collapse identical intervals across subclones into one event
  ss$key <- sprintf("chr%s:%.1f-%.1fMb_%s", ss$Chr, ss$Pos / 1e6, ss$End / 1e6, ss$event)
  for (k in unique(ss$key)) {
    e <- ss[ss$key == k, ]
    carriers <- intersect(unique(e$clone), big)
    noncar <- setdiff(big, unique(e$clone))
    g <- ann$gene[ann$chr == as.character(e$Chr[1]) &
                    ann$start >= e$Pos[1] & ann$end <= e$End[1]]
    if (length(g) < MIN_GENES_SEG || !length(carriers) || !length(noncar)) {
      res[[length(res) + 1L]] <- data.frame(
        sample_id = s, segment = k, chr = e$Chr[1], event = e$event[1],
        width_mb = e$width_mb[1], n_genes = length(g),
        n_carrier_clones = length(carriers), n_noncarrier_clones = length(noncar),
        n_carrier_cells = sum(cl %in% carriers), n_noncarrier_cells = sum(cl %in% noncar),
        mean_carrier = NA_real_, mean_noncarrier = NA_real_, delta = NA_real_,
        cohens_d = NA_real_, direction_matches_event = NA,
        evaluable = FALSE,
        reason = if (length(g) < MIN_GENES_SEG) sprintf("< %d genes in the interval", MIN_GENES_SEG)
                 else if (!length(carriers)) "no carrier clone above the cell minimum"
                 else "no non-carrier clone in this patient to compare against")
      next
    }
    sc <- colMeans(as.matrix(lg[g, cells, drop = FALSE]))
    a <- sc[cl %in% carriers]; b <- sc[cl %in% noncar]
    sp <- sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) /
                 (length(a) + length(b) - 2))
    d <- (mean(a) - mean(b)) / (sp + 1e-12)
    res[[length(res) + 1L]] <- data.frame(
      sample_id = s, segment = k, chr = e$Chr[1], event = e$event[1],
      width_mb = e$width_mb[1], n_genes = length(g),
      n_carrier_clones = length(carriers), n_noncarrier_clones = length(noncar),
      n_carrier_cells = length(a), n_noncarrier_cells = length(b),
      mean_carrier = mean(a), mean_noncarrier = mean(b),
      delta = mean(a) - mean(b), cohens_d = d,
      direction_matches_event = (e$event[1] == "gain" & d > 0) |
        (e$event[1] == "loss" & d < 0),
      evaluable = TRUE, reason = "evaluable")
  }
}
cx <- do.call(rbind, res)
p6_tsv(cx, file.path(P6_TAB, "BROAD_CNA_EXPRESSION_EFFECTS.tsv"))
ce <- cx[cx$evaluable, ]
p6_msg("evaluable broad events: %d of %d", nrow(ce), nrow(cx))
if (nrow(ce)) {
  p6_msg("direction matches the inferred event in %d of %d (%.0f%%)",
         sum(ce$direction_matches_event), nrow(ce),
         100 * mean(ce$direction_matches_event))
  print(as.data.frame(ce |> arrange(sample_id, desc(abs(cohens_d))) |>
    select(sample_id, segment, event, n_genes, n_carrier_cells,
           n_noncarrier_cells, cohens_d, direction_matches_event) |> head(20)),
    row.names = FALSE, digits = 3)
}
p6_msg("READ THIS AS INTERNAL CONSISTENCY, NOT VALIDATION: SCEVAN inferred these")
p6_msg("events from expression in the first place, so agreement is expected and")
p6_msg("is not independent evidence for the copy-number call.")

p6_sec("6. Loci of interest, described only as broad segments")
loci <- data.frame(gene = c("NF1", "NF2"), chr = c("17", "22"),
                   pos = c(31094927, 29603556))
lo <- do.call(rbind, lapply(seq_len(nrow(loci)), function(i) {
  hit <- seg[seg$Chr == as.integer(loci$chr[i]) & seg$Pos <= loci$pos[i] &
               seg$End >= loci$pos[i], ]
  if (!nrow(hit)) return(NULL)
  data.frame(locus = loci$gene[i], sample_id = hit$sample_id, chr = hit$Chr,
             segment_mb = sprintf("%.1f-%.1f", hit$Pos / 1e6, hit$End / 1e6),
             width_mb = hit$width_mb, inferred_event = hit$event, CN = hit$CN,
             permitted_wording = sprintf(
               "a broad inferred %s segment of %.1f Mb on chr%s whose interval contains the %s locus",
               hit$event, hit$width_mb, hit$Chr, loci$gene[i]),
             prohibited_wording = sprintf("%s deletion / %s-deleted cells",
                                          loci$gene[i], loci$gene[i]))
}))
if (!is.null(lo)) { print(as.data.frame(lo[, 1:7]), row.names = FALSE)
  p6_tsv(lo, file.path(P6_CNA, "M47_BROAD_SEGMENTS_CONTAINING_LOCI.tsv")) }

p6_sec("7. Broad CNA event <-> program association, within each patient")
PROGS <- grep("^program_P[0-9]+_score$", colnames(mm), value = TRUE)
PID <- sub("^program_(P[0-9]+)_score$", "\\1", PROGS)
lab <- read.delim("results/phase5/tables/final/MALIGNANT_PROGRAMS.tsv")
LABEL <- setNames(lab$program_label, lab$program_id)
pa <- list()
for (s in CLONE_RELIABLE) {
  cells <- mm$cell_id[mm$sample_id == s & !is.na(mm$tumor_clone_phase4)]
  cl <- mm$tumor_clone_phase4[match(cells, mm$cell_id)]
  big <- names(which(table(cl) >= MIN_CELLS))
  ss <- sub[sub$sample_id == s, ]
  ss$key <- sprintf("chr%s:%.1f-%.1fMb_%s", ss$Chr, ss$Pos / 1e6, ss$End / 1e6, ss$event)
  for (k in unique(ss$key)) {
    e <- ss[ss$key == k, ]
    carriers <- intersect(unique(e$clone), big); noncar <- setdiff(big, unique(e$clone))
    if (!length(carriers) || !length(noncar)) next
    ic <- cl %in% carriers; inc <- cl %in% noncar
    for (i in seq_along(PROGS)) {
      y <- mm[[PROGS[i]]][match(cells, mm$cell_id)]
      a <- y[ic]; b <- y[inc]
      sp <- sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) /
                   (length(a) + length(b) - 2))
      pa[[length(pa) + 1L]] <- data.frame(
        sample_id = s, segment = k, chr = e$Chr[1], event = e$event[1],
        program = PID[i], program_label = unname(LABEL[PID[i]]),
        n_carrier_cells = length(a), n_noncarrier_cells = length(b),
        median_carrier = median(a), median_noncarrier = median(b),
        cohens_d = (mean(a) - mean(b)) / (sp + 1e-12))
    }
  }
}
if (length(pa)) {
  pa <- do.call(rbind, pa)
  pa$strong <- abs(pa$cohens_d) >= 0.5
  pa <- pa[order(-abs(pa$cohens_d)), ]
  p6_tsv(pa, file.path(P6_TAB, "CNA_PROGRAM_ASSOCIATION.tsv"))
  p6_msg("clone-level CNA <-> program comparisons: %d (|d| >= 0.5 in %d)",
         nrow(pa), sum(pa$strong))
  print(as.data.frame(head(pa, 15)), row.names = FALSE, digits = 3)
} else p6_msg("no evaluable CNA <-> program comparison")

facts$segments <- list(broad_clonal_reliable = nrow(seg),
                       broad_subclonal = nrow(sub),
                       evaluable_events = if (exists("ce")) nrow(ce) else 0L,
                       direction_match_rate = if (exists("ce") && nrow(ce))
                         mean(ce$direction_matches_event) else NA_real_)
facts$non_independence_warning <- paste(
  "SCEVAN infers copy number FROM EXPRESSION. Coordinated expression on an",
  "affected chromosome is therefore an internal transcriptional consistency",
  "check, NOT orthogonal validation of the copy-number call.")
facts$cna_claim_level <- "chromosome / arm / large-segment only; no single-gene deletion or amplification is asserted"
p6_json(facts, file.path(P6_VAL, "m47_cna_expression_facts.json"))
p6_sec("M47 complete")
