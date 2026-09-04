#!/usr/bin/env Rscript
# =============================================================================
# Phase 5 - M38 - program annotation, pathways, recurrence and the relationship
# to the Phase 4 discrete tumour states.
#
# Reads the cNMF consensus solution chosen by the pre-declared M37d rule.
# Nothing here re-runs factorization and nothing re-selects K.
#
# Labels are CAUTIOUS and evidence-printed: each program's label is derived
# from its own top genes and its own top gene-set hit, and the evidence is
# written next to the label so a reader can disagree with it.
# =============================================================================
suppressPackageStartupMessages({ library(fgsea); library(ggplot2) })
source("scripts/phase5/utils/phase5_common.R")
source("scripts/phase5/utils/phase5_plot_utils.R")
set.seed(42)

args  <- commandArgs(trailingOnly = TRUE)
KSTAR <- as.integer(args[1])
DT    <- if (length(args) > 1) args[2] else "0_1"
stopifnot(!is.na(KSTAR))
p5_msg("selected K = %d, local-density-threshold tag = %s", KSTAR, DT)

CNMF <- "results/phase5/programs/cnmf"
GS   <- "external/genesets"
facts <- list(milestone = "M38", K = KSTAR, dt = DT,
              generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))

cnmf_path <- function(run, what, k = KSTAR, dt = DT)
  file.path(CNMF, run, sprintf("%s.%s.k_%d.dt_%s.%s", run, what, k, dt,
                               if (what == "usages") "consensus.txt" else "txt"))
read_tsv0 <- function(p) {
  x <- read.delim(p, row.names = 1, check.names = FALSE)
  as.matrix(x)
}

p5_sec("1. Load the chosen consensus solution")
U <- read_tsv0(cnmf_path("primary", "usages"))            # cells x K
Z <- read_tsv0(cnmf_path("primary", "gene_spectra_score")) # K x genes
TPM <- read_tsv0(cnmf_path("primary", "gene_spectra_tpm"))
colnames(U) <- paste0("P", seq_len(ncol(U)))
rownames(Z) <- colnames(U); rownames(TPM) <- colnames(U)
p5_msg("usages %d cells x %d programs | spectra %d programs x %d genes",
       nrow(U), ncol(U), nrow(Z), ncol(Z))

md <- readRDS(file.path(P5_PROG, "phase5_cell_metadata.rds"))
mmd <- md[rownames(U), ]
stopifnot(nrow(mmd) == nrow(U), all(mmd$malignancy_refined == "Malignant"))
REL <- U / rowSums(U)                                     # relative usage
PROGS <- colnames(U)

p5_sec("2. Per-cell program scores")
cs <- data.frame(cell_id = rownames(U), sample_id = mmd$sample_id,
                 malignancy_confidence = mmd$malignancy_confidence,
                 tumor_state_phase4 = mmd$tumor_state_phase4,
                 tumor_clone_phase4 = mmd$tumor_clone_phase4,
                 hist_annotation = mmd$annotation_ccc_phase3,
                 as.data.frame(REL), check.names = FALSE)
cs$dominant_program <- PROGS[max.col(REL, ties.method = "first")]
cs$dominant_program_score <- apply(REL, 1, max)
p5_tsv(cs, file.path(P5_TAB, "MALIGNANT_PROGRAM_CELL_SCORES.tsv"))

p5_sec("3. Top genes per program")
TOPN <- 50L
tg <- do.call(rbind, lapply(PROGS, function(p) {
  o <- order(Z[p, ], decreasing = TRUE)[seq_len(TOPN)]
  data.frame(program = p, rank = seq_len(TOPN), gene = colnames(Z)[o],
             spectra_score = Z[p, o], tpm = TPM[p, colnames(Z)[o]])
}))
p5_tsv(tg, file.path(P5_TAB, "MALIGNANT_PROGRAM_TOP_GENES.tsv"))
for (p in PROGS)
  p5_msg("  %-4s %s", p, paste(tg$gene[tg$program == p][1:15], collapse = " "))

p5_sec("4. Pathway enrichment (fgsea on the full ranked spectra score)")
read_gmt <- function(f) {
  ln <- strsplit(readLines(f), "\t")
  setNames(lapply(ln, function(x) unique(x[-c(1, 2)])), vapply(ln, `[`, "", 1))
}
gsets <- list(
  Hallmark = read_gmt(file.path(GS, "h.all.v2024.1.Hs.symbols.gmt")),
  Reactome = read_gmt(file.path(GS, "c2.cp.reactome.v2024.1.Hs.symbols.gmt")))
for (n in names(gsets)) p5_msg("  %s: %d sets", n, length(gsets[[n]]))

pw <- list()
for (coll in names(gsets)) {
  sets <- lapply(gsets[[coll]], function(g) intersect(g, colnames(Z)))
  sets <- sets[lengths(sets) >= 10 & lengths(sets) <= 500]
  for (p in PROGS) {
    st <- sort(Z[p, ], decreasing = TRUE)
    r <- suppressWarnings(fgsea(sets, st, minSize = 10, maxSize = 500,
                                nPermSimple = 10000))
    r <- r[order(r$padj, -abs(r$NES)), ]
    r <- as.data.frame(r[r$padj < 0.05 & r$NES > 0, ])
    if (!nrow(r)) next
    r$leadingEdge <- vapply(r$leadingEdge, function(x)
      paste(utils::head(x, 15), collapse = ","), character(1))
    pw[[length(pw) + 1L]] <- data.frame(program = p, collection = coll,
      pathway = r$pathway, NES = r$NES, pval = r$pval, padj = r$padj,
      size = r$size, leading_edge = r$leadingEdge)
  }
}
pw <- do.call(rbind, pw)
pw <- pw[order(pw$program, pw$padj, -pw$NES), ]
p5_tsv(pw, file.path(P5_TAB, "MALIGNANT_PROGRAM_PATHWAYS.tsv"))
for (p in PROGS) {
  s <- pw[pw$program == p, ]
  p5_msg("  %-4s %d enriched sets | top: %s", p, nrow(s),
         if (nrow(s)) paste(utils::head(s$pathway, 3), collapse = " | ") else "none")
}

p5_sec("5. Cautious program labels")
# Declared marker families, carried over verbatim from Phase 4 M31 so Phase 4
# and Phase 5 speak the same biological vocabulary. A label is assigned only
# when a family is clearly the best match; otherwise the program keeps a
# pathway-derived or explicitly Unassigned label.
FAM <- list(
  Schwann_like     = c("S100B","PLP1","MPZ","PMP22","SOX10","MBP","CNP","GFRA3","CDH19",
                       "NGFR","L1CAM","MIA","SEMA3B","PLEKHB1","CRYAB"),
  NeuralCrest_like = c("SOX9","TWIST1","PAX3","ZIC1","FOXD3","ETS1","NES","ERBB3","ABCB5",
                       "SOX4","ID1","TFAP2A","MSX1"),
  Mesenchymal_ECM  = c("COL1A1","COL1A2","COL3A1","COL5A1","COL6A3","FN1","POSTN","SPARC",
                       "THBS2","LOX","TNC","FBN1","VCAN","TAGLN","ACTA2"),
  Cycling          = c("MKI67","TOP2A","UBE2C","NEK2","KIF20A","DLGAP5","HJURP","ASPM",
                       "CCNB1","CDK1","PLK1","BIRC5","TYMS","RRM2"),
  Stress           = c("HSPA1A","HSPA1B","HSPB1","DNAJB1","HSPH1","JUN","JUNB","FOS",
                       "FOSB","EGR1","ATF3","DDIT3","HSPA6"),
  Interferon       = c("ISG15","IFI6","IFI27","IFIT1","IFIT3","MX1","MX2","OAS1","OASL",
                       "STAT1","IRF7","B2M","BST2"),
  Hypoxia_Angio    = c("VEGFA","ANGPT2","ANGPTL4","ADM","NDRG1","SLC2A1","HIF1A","EGLN3",
                       "PGK1","LDHA","CA9","BNIP3"),
  AntigenPresent   = c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","B2M","TAP1","TAP2","PSMB9",
                       "CD74","HLA-DRA","HLA-DRB1"))
TOP200 <- lapply(PROGS, function(p) colnames(Z)[order(Z[p, ], decreasing = TRUE)[1:200]])
names(TOP200) <- PROGS
famhit <- do.call(rbind, lapply(PROGS, function(p) {
  ov <- vapply(FAM, function(g) length(intersect(g, TOP200[[p]])), integer(1))
  fr <- ov / lengths(FAM)
  data.frame(program = p, family = names(FAM), n_overlap = ov, frac_family = fr,
             row.names = NULL)
}))
p5_tsv(famhit, file.path(P5_VAL, "M38_PROGRAM_MARKER_FAMILY_OVERLAP.tsv"))

# A factorization of malignant cells can still capture technical structure:
# ribosomal-protein programs, pseudogene/lncRNA programs, and ambient myeloid
# RNA. Those are detected explicitly and named for what they are, because a
# ribosomal factor labelled "Mesenchymal_ECM" on a 0.20 family overlap would
# mislead every downstream reader.
is_ribo   <- function(g) grepl("^RP[LS][0-9]|^RPLP[0-9]|^RPSA$", g)
is_pseudo <- function(g) grepl("^RP[0-9]+-|^AC[0-9]{6}|^AL[0-9]{6}|^CTD-|^CTA-|^CTC-|^LINC[0-9]|^MIR[0-9]|-AS[0-9]?$|^RP[LS][0-9]+P[0-9]+$|^AP[0-9]{6}", g)
MYELOID <- c("C1QA","C1QB","C1QC","CD14","TYROBP","AIF1","CSF1R","LYZ","FCGR3A",
             "FCER1G","MS4A6A","MS4A4A","ITGAM","PTPRC","CD68","CD163","FOLR2",
             "MNDA","F13A1","RNASE1","CCL3")
tech <- do.call(rbind, lapply(PROGS, function(p) {
  g50 <- tg$gene[tg$program == p][1:50]
  data.frame(program = p,
             frac_ribosomal = mean(is_ribo(g50)),
             frac_pseudogene_lncRNA = mean(is_pseudo(g50)),
             n_myeloid_markers = sum(g50 %in% MYELOID))
}))
tech$likely_technical <- tech$frac_ribosomal >= 0.30 |
  tech$frac_pseudogene_lncRNA >= 0.30 | tech$n_myeloid_markers >= 4
tech$technical_class <- with(tech, ifelse(
  frac_ribosomal >= 0.30, "Translation_ribosomal",
  ifelse(frac_pseudogene_lncRNA >= 0.30, "Pseudogene_lncRNA",
  ifelse(n_myeloid_markers >= 4, "Myeloid_ambient_like", NA_character_))))
print(tech, row.names = FALSE, digits = 3)
p5_tsv(tech, file.path(P5_VAL, "M38_PROGRAM_TECHNICAL_CONTENT.tsv"))

lab <- do.call(rbind, lapply(PROGS, function(p) {
  f <- famhit[famhit$program == p, ]
  f <- f[order(-f$frac_family), ]
  best <- f$family[1]; bf <- f$frac_family[1]; second <- f$frac_family[2]
  tp <- pw$pathway[pw$program == p][1]
  tc <- tech[tech$program == p, ]
  cautious <- if (isTRUE(tc$likely_technical)) tc$technical_class
              else if (bf >= 0.20 && bf >= 2 * second) best
              else if (bf >= 0.20) paste0(best, "_mixed")
              else if (!is.na(tp)) paste0("Unassigned_", sub("^(HALLMARK|REACTOME)_", "", tp))
              else "Unassigned"
  basis <- if (isTRUE(tc$likely_technical))
    sprintf("technical-dominated: %.0f%% of the top 50 genes are ribosomal, %.0f%% pseudogene/lncRNA, %d myeloid markers; top set %s. The marker-family match (%s at %.2f) is not treated as the label.",
            100 * tc$frac_ribosomal, 100 * tc$frac_pseudogene_lncRNA,
            tc$n_myeloid_markers, ifelse(is.na(tp), "none", tp), best, bf)
    else sprintf("top-200 genes overlap %s at %.2f of the family (next %.2f); top set %s",
                 best, bf, second, ifelse(is.na(tp), "none", tp))
  data.frame(program = p, program_label = cautious,
             best_family = best, best_family_frac = bf,
             second_family_frac = second, top_pathway = tp,
             likely_technical = tc$likely_technical,
             technical_class = tc$technical_class,
             frac_ribosomal_top50 = tc$frac_ribosomal,
             frac_pseudogene_top50 = tc$frac_pseudogene_lncRNA,
             n_myeloid_markers_top50 = tc$n_myeloid_markers,
             label_basis = basis)
}))
print(lab[, c("program","program_label","best_family","best_family_frac","top_pathway")],
      row.names = FALSE)

p5_sec("6. Patient distribution and recurrence (pre-declared criteria)")
ACT <- 0.20; MINCELL <- 10L; MINFRAC <- 0.05; MINPAT <- 3L
p5_msg("declared: activity >= %.2f relative usage; a patient carries a program with >= %d active cells AND >= %.0f%% of its malignant cells; recurrent needs >= %d patients",
       ACT, MINCELL, 100 * MINFRAC, MINPAT)
npat <- table(mmd$sample_id)
pd <- do.call(rbind, lapply(PROGS, function(p) {
  act <- REL[, p] >= ACT
  do.call(rbind, lapply(SAMPLES, function(s) {
    m <- mmd$sample_id == s
    data.frame(program = p, sample_id = s, n_malignant = sum(m),
               n_active = sum(act & m), frac_of_patient = sum(act & m) / sum(m),
               median_rel_usage = median(REL[m, p]),
               mean_rel_usage = mean(REL[m, p]),
               n_dominant = sum(cs$dominant_program == p & m),
               carries = sum(act & m) >= MINCELL &&
                 (sum(act & m) / sum(m)) >= MINFRAC)
  }))
}))
p5_tsv(pd, file.path(P5_TAB, "MALIGNANT_PROGRAM_PATIENT_DISTRIBUTION.tsv"))

rec <- do.call(rbind, lapply(PROGS, function(p) {
  s <- pd[pd$program == p, ]
  act <- REL[, p] >= ACT
  na <- sum(act)
  tp <- if (na) sort(table(mmd$sample_id[act]), decreasing = TRUE) else integer(0)
  nc <- sum(s$carries)
  data.frame(program = p, n_active_cells = na,
             carrying_patients = paste(s$sample_id[s$carries], collapse = ","),
             n_carrying = nc,
             dominant_patient = if (na) names(tp)[1] else NA_character_,
             dominant_patient_fraction = if (na) unname(tp[1] / na) else NA_real_,
             n_dominant_cells = sum(cs$dominant_program == p),
             recurrence_status = if (nc >= MINPAT) "recurrent"
                                 else if (nc == 2) "shared-limited"
                                 else if (nc == 1) "patient-private" else "uncertain")
}))
print(rec, row.names = FALSE)
facts$recurrence_criteria <- list(activity_threshold = ACT, min_cells = MINCELL,
                                  min_frac_of_patient = MINFRAC,
                                  min_patients = MINPAT,
                                  declared = "before program labels were assigned")

p5_sec("7. Balanced-run support")
# Does each primary program survive when the largest patient is cut down to the
# size of the smallest sufficiently sized one?
bal_ok <- file.exists(cnmf_path("balanced", "gene_spectra_score"))
if (bal_ok) {
  Zb <- read_tsv0(cnmf_path("balanced", "gene_spectra_score"))
  rownames(Zb) <- paste0("B", seq_len(nrow(Zb)))
  g <- intersect(colnames(Z), colnames(Zb))
  cosm <- function(a, b) sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
  bal <- do.call(rbind, lapply(PROGS, function(p) {
    sims <- vapply(rownames(Zb), function(q) cosm(Z[p, g], Zb[q, g]), numeric(1))
    jac <- vapply(rownames(Zb), function(q) {
      a <- colnames(Z)[order(Z[p, ], decreasing = TRUE)[1:100]]
      b <- colnames(Zb)[order(Zb[q, ], decreasing = TRUE)[1:100]]
      length(intersect(a, b)) / length(union(a, b))
    }, numeric(1))
    i <- which.max(sims)
    data.frame(program = p, best_balanced_match = names(sims)[i],
               cosine = unname(sims[i]), top100_jaccard = unname(jac[i]),
               balanced_support = unname(sims[i]) >= 0.60 | unname(jac[i]) >= 0.25)
  }))
  print(bal, row.names = FALSE)
} else {
  p5_msg("  balanced consensus not available - recorded as NA")
  bal <- data.frame(program = PROGS, best_balanced_match = NA_character_,
                    cosine = NA_real_, top100_jaccard = NA_real_,
                    balanced_support = NA)
}
p5_tsv(bal, file.path(P5_TAB, "PROGRAM_PATIENT_BALANCE_SENSITIVITY.tsv"))

p5_sec("8. Programs versus the Phase 4 discrete states")
# The specific question: are the five patient-private Mesenchymal_ECM-* states
# manifestations of one shared continuous program, several, or none?
st <- do.call(rbind, lapply(PROGS, function(p) {
  do.call(rbind, lapply(sort(unique(mmd$tumor_state_phase4)), function(s) {
    m <- mmd$tumor_state_phase4 == s
    data.frame(program = p, tumor_state_phase4 = s, n_cells = sum(m),
               median_rel_usage = median(REL[m, p]),
               mean_rel_usage = mean(REL[m, p]),
               frac_active = mean(REL[m, p] >= ACT),
               frac_dominant = mean(cs$dominant_program[m] == p),
               n_patients_in_state = length(unique(mmd$sample_id[m])))
  }))
}))
p5_tsv(st, file.path(P5_TAB, "PROGRAM_VS_PHASE4_STATE.tsv"))
ecm_states <- grep("^Mesenchymal_ECM", unique(mmd$tumor_state_phase4), value = TRUE)
p5_msg("Mesenchymal_ECM states: %s", paste(ecm_states, collapse = ", "))
for (s in ecm_states) {
  d <- st[st$tumor_state_phase4 == s, ]
  d <- d[order(-d$median_rel_usage), ]
  p5_msg("  %-18s top programs: %s", s,
         paste(sprintf("%s=%.3f", d$program[1:3], d$median_rel_usage[1:3]),
               collapse = "  "))
}
top_by_state <- st |> group_by(tumor_state_phase4) |>
  slice_max(median_rel_usage, n = 1, with_ties = FALSE) |> ungroup()
p5_msg("dominant program per Phase 4 state:")
print(as.data.frame(top_by_state[, c("tumor_state_phase4","program",
                                     "median_rel_usage","n_cells")]), row.names = FALSE)
ecm_top <- top_by_state$program[top_by_state$tumor_state_phase4 %in% ecm_states]
facts$ecm_state_programs <- list(
  states = ecm_states, dominant_program_per_state = ecm_top,
  n_distinct_programs = length(unique(ecm_top)),
  interpretation = if (length(unique(ecm_top)) == 1)
    "the five patient-private Mesenchymal_ECM states share a single dominant continuous program"
    else sprintf("the Mesenchymal_ECM states map onto %d distinct programs",
                 length(unique(ecm_top))))
p5_msg("=> %s", facts$ecm_state_programs$interpretation)

p5_sec("9. MALIGNANT_PROGRAMS.tsv")
progs <- lab |>
  left_join(rec, by = "program") |>
  left_join(bal |> select(program, balanced_cosine = cosine,
                          balanced_top100_jaccard = top100_jaccard,
                          balanced_support), by = "program") |>
  mutate(top_genes = vapply(program, function(p)
           paste(tg$gene[tg$program == p][1:20], collapse = ","), character(1)),
         dominant_pathway = top_pathway,
         patients_represented = vapply(program, function(p)
           sum(pd$carries[pd$program == p]), integer(1)),
         n_cells_dominant = n_dominant_cells,
         interpretation = sprintf(
           "%s; %s across patients (%d carrying); dominant patient %s at %.2f of active cells; balanced-run cosine %s",
           label_basis, recurrence_status, n_carrying, dominant_patient,
           dominant_patient_fraction,
           ifelse(is.na(balanced_cosine), "NA", sprintf("%.2f", balanced_cosine)))) |>
  select(program_id = program, program_label, top_genes, dominant_pathway,
         likely_technical, technical_class, frac_ribosomal_top50,
         frac_pseudogene_top50, n_myeloid_markers_top50,
         best_family, best_family_frac, n_active_cells, n_cells_dominant,
         patients_represented, carrying_patients, dominant_patient,
         dominant_patient_fraction, recurrence_status, balanced_cosine,
         balanced_top100_jaccard, balanced_support, label_basis, interpretation)
p5_tsv(progs, file.path(P5_TAB, "MALIGNANT_PROGRAMS.tsv"))
print(as.data.frame(progs[, c("program_id","program_label","recurrence_status",
                              "patients_represented","dominant_patient_fraction",
                              "balanced_cosine")]), row.names = FALSE, digits = 3)

facts$programs <- progs
facts$n_programs <- nrow(progs)
facts$n_recurrent <- sum(progs$recurrence_status == "recurrent")
facts$n_shared_limited <- sum(progs$recurrence_status == "shared-limited")
facts$n_patient_private <- sum(progs$recurrence_status == "patient-private")
facts$n_uncertain <- sum(progs$recurrence_status == "uncertain")
p5_msg("recurrent %d | shared-limited %d | patient-private %d | uncertain %d",
       facts$n_recurrent, facts$n_shared_limited, facts$n_patient_private,
       facts$n_uncertain)
facts$n_likely_technical <- sum(progs$likely_technical)
facts$technical_programs <- progs$program_id[progs$likely_technical]
p5_msg("programs flagged as technical-dominated: %d (%s)", facts$n_likely_technical,
       paste(facts$technical_programs, collapse = ", "))
p5_msg("biologically interpretable programs: %d", nrow(progs) - facts$n_likely_technical)
p5_json(facts, file.path(P5_VAL, "m38_program_annotation_facts.json"))
p5_sec("M38 complete")
