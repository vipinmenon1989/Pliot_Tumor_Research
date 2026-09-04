#!/usr/bin/env Rscript
# =============================================================================
# Phase 5 - M40c - a further DECLARED sensitivity gene universe.
#
# Why this run exists: with the primary gene universe, 3 of the 8 selected
# programs are technical-dominated - two ribosomal/pseudogene factors and one
# carrying 18 myeloid markers in its top 50 genes (ambient RNA in a
# malignant-only factorization). The only program carried by more than one
# patient is one of the ribosomal ones.
#
# The question this run answers, and it is a real scientific question rather
# than a cosmetic clean-up: once the technical gene CLASSES are removed, do
# biologically interpretable programs become shared across patients, or is the
# patient-private result unchanged?
#
# Removed: ribosomal-protein genes, ribosomal pseudogenes, lncRNA/clone-name
# genes, and the canonical myeloid markers that mark ambient contamination.
# ECM, HLA, interferon, Schwann, neural-crest, angiogenesis and cell-cycle
# genes are all RETAINED - they are the signal.
# =============================================================================
source("scripts/phase5/utils/phase5_common.R")
set.seed(42)
OUT <- file.path(P5_PROG, "cnmf_input")
facts <- list(milestone = "M40c", generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))

md <- readRDS(file.path(P5_PROG, "phase5_cell_metadata.rds"))
cnts <- readRDS(file.path(P5_PROG, "phase5_working_counts.rds"))
genes <- readLines(file.path(OUT, "primary", "genes.txt"))
mal <- md$cell_id[md$malignancy_refined == "Malignant"]

is_ribo   <- grepl("^RP[LS][0-9]|^RPLP[0-9]|^RPSA$", genes)
is_pseudo <- grepl("^RP[0-9]+-|^AC[0-9]{6}|^AL[0-9]{6}|^CTD-|^CTA-|^CTC-|^LINC[0-9]|^MIR[0-9]|-AS[0-9]?$|^RP[LS][0-9]+P[0-9]+$|^AP[0-9]{6}", genes)
MYELOID <- c("C1QA","C1QB","C1QC","CD14","TYROBP","AIF1","CSF1R","LYZ","FCGR3A",
             "FCER1G","MS4A6A","MS4A4A","ITGAM","PTPRC","CD68","CD163","FOLR2",
             "MNDA","F13A1","RNASE1","CCL3","CCL4","SPI1","CYBB","LST1","HLA-DRA")
is_mye <- genes %in% MYELOID
drop <- is_ribo | is_pseudo | is_mye
p5_msg("primary gene universe            %d", length(genes))
p5_msg("  ribosomal proteins removed     %d", sum(is_ribo))
p5_msg("  pseudogene / lncRNA removed    %d", sum(is_pseudo))
p5_msg("  canonical myeloid removed      %d", sum(is_mye))
p5_msg("  nortech gene universe          %d", sum(!drop))
keep_genes <- genes[!drop]
for (fam in list(ECM = c("COL1A1","COL1A2","COL3A1","FN1","POSTN","SPARC"),
                 Schwann = c("S100B","PLP1","SOX10","MPZ","CDH19"),
                 Cycling = c("MKI67","TOP2A","CDK1","CCNB1"),
                 IFN = c("ISG15","IFIT1","MX1","STAT1"),
                 HLA = c("HLA-A","HLA-B","HLA-E","B2M")))
  p5_msg("  retained: %s", paste(intersect(fam, keep_genes), collapse = " "))

d <- file.path(OUT, "nortech"); dir.create(d, showWarnings = FALSE, recursive = TRUE)
sub <- cnts[keep_genes, mal, drop = FALSE]
sub <- sub[Matrix::rowSums(sub) > 0, , drop = FALSE]
Matrix::writeMM(as(sub, "dgCMatrix"), file.path(d, "counts.mtx"))
writeLines(rownames(sub), file.path(d, "genes.txt"))
writeLines(colnames(sub), file.path(d, "cells.txt"))
p5_msg("  [nortech] %d genes x %d cells", nrow(sub), ncol(sub))
facts$nortech <- list(genes_in = length(genes), ribosomal_removed = sum(is_ribo),
                      pseudogene_removed = sum(is_pseudo), myeloid_removed = sum(is_mye),
                      genes_out = nrow(sub), cells = ncol(sub),
                      rationale = "tests whether the patient-private program result survives removal of the technical gene classes")
p5_json(facts, file.path(P5_VAL, "m40c_nortech_facts.json"))
p5_sec("M40c complete")
