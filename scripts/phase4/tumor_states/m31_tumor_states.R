#!/usr/bin/env Rscript
# =============================================================================
# Phase 4 · M31 — malignant-only tumour-state analysis
#
# Operates on malignancy_refined == "Malignant" only (§30). The broad Phase 2
# `MPNST-Tumor` label is NOT reused as the subset definition.
#
# REDUCTIONS (§31): the Phase 2 `pca` / `postint_harmony` / `postint_umap_harmony`
# reductions are FROZEN and are never overwritten. This script builds NEW,
# separately named malignant-only reductions and reports how they compare with
# simply reusing the Phase 2 embedding, so the choice is visible:
#     malignant_pca · malignant_harmony · malignant_umap
#
# CLONE vs STATE (§33): tumor_clone_phase4 is CNV architecture (from SCEVAN);
# tumor_state_phase4 is transcriptional phenotype. Their relationship is
# tabulated; one-to-one correspondence is never assumed.
# =============================================================================
options(stringsAsFactors = FALSE); options(future.globals.maxSize = 64 * 1024^3)
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix); library(harmony)
  library(ggplot2); library(dplyr); library(tidyr); library(jsonlite)
})
source("scripts/phase4/utils/phase4_plot_utils.R")
source("scripts/R/provenance_utils.R")
set.seed(42)

SC  <- "results/phase4/scevan/by_sample"
MAL <- "results/phase4/malignancy"
TS  <- "results/phase4/tumor_states"
TAB <- "results/phase4/tables"
FIG <- "results/phase4/figures"
SAMPLES <- c("MPNST_1","MPNST_2","MPNST_3","MPNST_4")
for (d in c(TS, TAB, FIG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
log_ <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
sec  <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")
CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))

# --- 1. malignant subset ---------------------------------------------------
sec("1. Malignant subset from malignancy_refined")
calls <- read.delim(file.path(MAL, "PHASE4_MALIGNANCY_CALLS.tsv"), stringsAsFactors = FALSE)
rownames(calls) <- calls$cell_id
mal <- calls[calls$malignancy_refined == "Malignant", , drop = FALSE]
log_("Malignant cells: ", nrow(mal), " of ", nrow(calls),
     sprintf(" (%.2f%%)", 100*nrow(mal)/nrow(calls)))
log_("by confidence:"); print(table(mal$malignancy_confidence))
log_("by patient:");    print(table(mal$sample_id))
log_("by Phase 2 annotation of origin:")
print(sort(table(mal$annotation_ccc_phase3), decreasing = TRUE))
if (nrow(mal) < 200) stop("STOP: fewer than 200 malignant cells; tumour-state analysis is not interpretable.")

# --- 2. expression: RNA counts -> LogNormalize -----------------------------
sec("2. Expression basis")
# Built from the M28-extracted RNA counts, the same basis Phase 3 and the Phase 4
# CCC work use. The 6 GB Phase 2 object is not reloaded and SCT (four models) is
# not used.
mats <- lapply(SAMPLES, function(s) {
  m <- readRDS(file.path(SC, s, "counts_raw.rds"))
  keep <- intersect(colnames(m), mal$cell_id)
  m[, keep, drop = FALSE] })
genes <- sort(Reduce(union, lapply(mats, rownames)))
pad <- function(m) { miss <- setdiff(genes, rownames(m))
  if (length(miss)) m <- rbind(m, Matrix(0, nrow = length(miss), ncol = ncol(m),
      sparse = TRUE, dimnames = list(miss, colnames(m))))
  m[genes, , drop = FALSE] }
cm <- do.call(cbind, lapply(mats, pad)); rm(mats); invisible(gc(FALSE))
log_("malignant counts: ", nrow(cm), " genes x ", ncol(cm), " cells")

obj <- CreateSeuratObject(counts = cm, project = "MPNST_malignant_phase4",
                          min.cells = 3, min.features = 0)
log_("after min.cells=3: ", nrow(obj), " genes x ", ncol(obj), " cells")
mm <- mal[colnames(obj), , drop = FALSE]
for (f in c("sample_id","annotation_ccc_phase3","annotation_ccc_refined","malignancy_refined",
            "malignancy_confidence","malignancy_rule","scevan_call","scevan_clone",
            "population_class","cnv_burden","primary_cluster")) {
  if (f %in% colnames(mm)) obj[[f]] <- mm[[f]]
}
obj$tumor_clone_phase4 <- mm$scevan_clone            # CNV architecture (§33)
obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                     scale.factor = 1e4, verbose = FALSE)

# --- 3. NEW malignant-only reductions -------------------------------------
sec("3. New malignant-only reductions (Phase 2 reductions untouched)")
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
NPCS <- min(50L, ncol(obj) - 1L)
obj <- RunPCA(obj, npcs = NPCS, reduction.name = "malignant_pca",
              reduction.key = "malPC_", verbose = FALSE)
sdev <- Stdev(obj, reduction = "malignant_pca")
varexp <- sdev^2 / sum(sdev^2)
cum <- cumsum(varexp)
NDIM <- min(30L, length(sdev) - 1L)
log_("malignant_pca: ", length(sdev), " PCs | dims 1:", NDIM, " capture ",
     sprintf("%.1f%%", 100*cum[NDIM]), " of variance",
     " (30 dims matches the Phase 2 clustering convention)")

# sample_id = patient = dataset, so malignant cells from four patients must be
# integrated or every cluster will be patient-private. Harmony with the SAME
# batch variable Phase 2 used, written to a NEW reduction name.
set.seed(42)
obj <- RunHarmony(obj, group.by.vars = "sample_id", reduction.use = "malignant_pca",
                  dims.use = 1:NDIM, reduction.save = "malignant_harmony", verbose = TRUE)
obj <- RunUMAP(obj, reduction = "malignant_harmony", dims = 1:NDIM,
               reduction.name = "malignant_umap", reduction.key = "malUMAP_",
               seed.use = 42, verbose = FALSE)
log_("reductions on the malignant object: ", paste(Reductions(obj), collapse = ", "))
log_("Phase 2 'pca'/'postint_harmony'/'postint_umap_harmony' are NOT present here and were never touched.")

# Does the malignant-only embedding actually differ from reusing Phase 2's? (§31)
h2 <- readRDS(file.path(MAL, "phase2_harmony_embedding.rds"))[colnames(obj), 1:NDIM]
hm <- Embeddings(obj, "malignant_harmony")[, 1:NDIM]
knn_overlap <- function(A, B, k = 20) {
  ia <- FNN::get.knn(A, k = k)$nn.index; ib <- FNN::get.knn(B, k = k)$nn.index
  mean(vapply(seq_len(nrow(ia)), function(i) length(intersect(ia[i,], ib[i,]))/k, numeric(1))) }
ov <- if (requireNamespace("FNN", quietly = TRUE)) knn_overlap(h2, hm) else NA_real_
log_(sprintf("kNN(k=20) overlap between Phase 2 Harmony and malignant-only Harmony: %.4f", ov))
log_("A low overlap is the justification for a new reduction: the Phase 2 embedding was")
log_("optimised to separate cell TYPES, not states within the malignant compartment.")

# --- 4. resolution sweep ---------------------------------------------------
sec("4. Clustering resolution sweep (no predetermined number of states, §32)")
obj <- FindNeighbors(obj, reduction = "malignant_harmony", dims = 1:NDIM,
                     k.param = 20, verbose = FALSE)
RES <- c(0.2, 0.3, 0.4, 0.5, 0.7, 1.0)
sw <- list()
for (r in RES) {
  obj <- FindClusters(obj, resolution = r, algorithm = 1, random.seed = 42, verbose = FALSE)
  cl <- obj@meta.data[[paste0("RNA_snn_res.", r)]]
  n <- nlevels(droplevels(factor(cl)))
  tiny <- sum(table(cl) < 30)
  # A state present in only one patient is exploratory, not recurrent biology (§36).
  multi <- sum(vapply(split(as.character(obj$sample_id), cl),
                      function(x) length(unique(x)) >= 2, logical(1)))
  sw[[length(sw)+1]] <- data.frame(resolution = r, n_clusters = n,
    n_clusters_lt30_cells = tiny, n_clusters_multipatient = multi,
    frac_multipatient = multi/n, min_size = min(table(cl)), max_size = max(table(cl)))
  log_(sprintf("  res %.1f -> %d clusters | %d with <30 cells | %d/%d multi-patient",
               r, n, tiny, multi, n))
}
sw <- bind_rows(sw)
write.table(sw, file.path(TAB, "TUMOR_STATE_RESOLUTION_SWEEP.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
# Selection rule, printed before it is applied: prefer the largest resolution with
# no cluster under 30 cells and every cluster multi-patient; if none qualifies,
# take the largest fraction of multi-patient clusters, then the smaller resolution.
ok <- sw[sw$n_clusters_lt30_cells == 0 & sw$frac_multipatient == 1, ]
CHOSEN <- if (nrow(ok)) max(ok$resolution) else {
  b <- sw[sw$frac_multipatient == max(sw$frac_multipatient), ]; min(b$resolution) }
log_("SELECTED resolution: ", CHOSEN, " (rule printed above, applied after the sweep)")
obj$tumor_cluster_phase4 <- factor(obj@meta.data[[paste0("RNA_snn_res.", CHOSEN)]])
Idents(obj) <- obj$tumor_cluster_phase4
log_("cluster sizes:"); print(table(obj$tumor_cluster_phase4))

# --- 5. programme scores (§34) --------------------------------------------
sec("5. Candidate biological programmes (candidates, not labels to force)")
PROG <- list(
  Schwann_like       = c("S100B","PLP1","MPZ","PMP22","SOX10","MBP","CNP","GFRA3","CDH19",
                         "NGFR","L1CAM","MIA","SEMA3B","PLEKHB1","CRYAB"),
  NeuralCrest_like   = c("SOX9","TWIST1","PAX3","ZIC1","FOXD3","ETS1","NES","ERBB3","ABCB5",
                         "SOX4","ID1","TFAP2A","MSX1"),
  Mesenchymal_ECM    = c("COL1A1","COL1A2","COL3A1","COL5A1","COL6A3","FN1","POSTN","SPARC",
                         "THBS2","LOX","TNC","FBN1","VCAN","TAGLN","ACTA2"),
  Cycling            = c("MKI67","TOP2A","UBE2C","NEK2","KIF20A","DLGAP5","HJURP","ASPM",
                         "CCNB1","CDK1","PLK1","BIRC5","TYMS","RRM2"),
  Stress             = c("HSPA1A","HSPA1B","HSPB1","DNAJB1","HSPH1","JUN","JUNB","FOS",
                         "FOSB","EGR1","ATF3","DDIT3","HSPA6"),
  Interferon         = c("ISG15","IFI6","IFI27","IFIT1","IFIT3","MX1","MX2","OAS1","OASL",
                         "STAT1","IRF7","B2M","BST2"),
  Angiogenic         = c("VEGFA","ANGPT2","ANGPTL4","ADM","NDRG1","SLC2A1","HIF1A","EGLN3",
                         "PGK1","LDHA","CA9","BNIP3"),
  Immune_interacting = c("APP","CD99","ANXA1","THBS1","HLA-A","HLA-B","HLA-C","HLA-E","HLA-F",
                         "B2M","CD47","LGALS1","LGALS3","MIF","CSF1","IL15"),
  Notch_perivascular = c("JAG1","JAG2","DLL4","DLL1","NOTCH1","NOTCH2","NOTCH3","HEY1","HEY2",
                         "HES1","RGS5","PDGFRB"))
present <- lapply(PROG, function(g) intersect(g, rownames(obj)))
for (n in names(present))
  log_(sprintf("  %-20s %2d/%2d genes present", n, length(present[[n]]), length(PROG[[n]])))
nbin <- 24; ctrl <- min(100L, floor(nrow(obj)/(nbin*3)))
log_("AddModuleScore nbin=", nbin, " ctrl=", ctrl)
obj <- AddModuleScore(obj, features = present, name = "PROG_", nbin = nbin,
                      ctrl = ctrl, seed = 42, verbose = FALSE)
pn <- paste0("PROG_", seq_along(present))
for (i in seq_along(present))
  obj[[paste0("prog_", names(present)[i])]] <- obj@meta.data[[pn[i]]]
prog_cols <- paste0("prog_", names(present))

pm <- obj@meta.data |> select(cluster = tumor_cluster_phase4, all_of(prog_cols)) |>
  group_by(cluster) |> summarise(across(everything(), median), .groups = "drop")
log_("median programme score per cluster:")
print(as.data.frame(pm), row.names = FALSE, digits = 3)
write.table(pm, file.path(TAB, "TUMOR_STATE_PROGRAM_SCORES.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- 6. markers ------------------------------------------------------------
sec("6. Cluster markers (characterisation for annotation, NOT condition DE)")
mk <- FindAllMarkers(obj, assay = "RNA", layer = "data", only.pos = TRUE,
                     min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
mk <- mk |> filter(p_val_adj < 0.05) |> arrange(cluster, desc(avg_log2FC))
write.table(mk, file.path(TAB, "TUMOR_STATE_MARKERS.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
log_("markers (p_adj<0.05): ", nrow(mk), " across ", length(unique(mk$cluster)), " clusters")
if (!nrow(mk))
  stop("STOP: no significant cluster markers. Tumour states would not be interpretable.")
top <- mk |> group_by(cluster) |> slice_head(n = 12) |> ungroup()
for (c in levels(obj$tumor_cluster_phase4))
  log_("  C", c, ": ", paste(top$gene[top$cluster == c], collapse = ", "))

# --- 7. state assignment ---------------------------------------------------
sec("7. tumor_state_phase4 assignment")
# Rule, printed before application: a cluster takes the name of its highest-scoring
# programme ONLY if that programme's median score is positive and exceeds the next
# best by >= 0.02; otherwise the state is "Uncertain" (§34).
zs <- as.data.frame(pm[, prog_cols]); rownames(zs) <- pm$cluster
assign_state <- function(v) {
  o <- order(v, decreasing = TRUE)
  if (v[o[1]] <= 0) return(c("Uncertain", "top programme score is not positive"))
  if ((v[o[1]] - v[o[2]]) < 0.02)
    return(c("Uncertain", sprintf("top two programmes within 0.02 (%s %.3f vs %s %.3f)",
             names(v)[o[1]], v[o[1]], names(v)[o[2]], v[o[2]])))
  c(sub("^prog_", "", names(v)[o[1]]),
    sprintf("%s %.3f, next %s %.3f", names(v)[o[1]], v[o[1]], names(v)[o[2]], v[o[2]]))
}
st <- t(apply(zs, 1, function(v) assign_state(setNames(as.numeric(v), colnames(zs)))))
state_map <- data.frame(cluster = rownames(zs), state = st[,1], evidence = st[,2],
                        stringsAsFactors = FALSE)
# Disambiguate repeated programme winners so each cluster keeps its own identity.
dup <- table(state_map$state)
for (s in names(dup)[dup > 1 & names(dup) != "Uncertain"]) {
  i <- which(state_map$state == s)
  state_map$state[i] <- paste0(s, "-", seq_along(i))
}
print(state_map, row.names = FALSE)
obj$tumor_state_phase4 <- state_map$state[match(as.character(obj$tumor_cluster_phase4),
                                                state_map$cluster)]
log_("tumor_state_phase4:"); print(table(obj$tumor_state_phase4))
write.table(state_map, file.path(TAB, "TUMOR_STATE_DEFINITIONS.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- 8. patient distribution (§36) ----------------------------------------
sec("8. Patient distribution per state")
pd <- obj@meta.data |> group_by(tumor_state_phase4) |>
  summarise(cluster = paste(sort(unique(as.character(tumor_cluster_phase4))), collapse = ","),
            n_cells = n(),
            fraction_of_malignant = n()/ncol(obj),
            patients_represented = length(unique(sample_id)),
            patient_breakdown = paste(sprintf("%s=%d", names(table(sample_id)),
                                              as.integer(table(sample_id))), collapse = "; "),
            dominant_patient = names(sort(table(sample_id), decreasing = TRUE))[1],
            dominant_patient_fraction = max(table(sample_id))/n(),
            scevan_clones = paste(sort(unique(na.omit(tumor_clone_phase4))), collapse = "; "),
            n_scevan_clones = length(unique(na.omit(tumor_clone_phase4))),
            median_cnv_burden = median(cnv_burden, na.rm = TRUE),
            phase2_origin = paste(sprintf("%s=%d",
              names(sort(table(annotation_ccc_phase3), decreasing = TRUE)),
              sort(table(annotation_ccc_phase3), decreasing = TRUE)), collapse = "; "),
            .groups = "drop") |>
  # §36: a state present in only one patient must be called patient-specific.
  # Mere PRESENCE in >=3 patients is not recurrence - a state with 1,587 cells from
  # one patient and 1, 3 and 7 from the others is patient-private with stragglers.
  # Recurrence therefore requires that no single patient dominates (<=80% of the
  # state's cells) AND that at least 3 patients each contribute a non-trivial share
  # (>=5% of the state, and at least 10 cells).
  mutate(
    n_patients_nontrivial = vapply(seq_len(n()), function(i) {
      br <- strsplit(patient_breakdown[i], "; ")[[1]]
      cn <- suppressWarnings(as.numeric(sub(".*=", "", br)))
      sum(cn >= 10 & cn >= 0.05 * n_cells[i], na.rm = TRUE) }, integer(1)),
    recurrence = case_when(
      patients_represented == 1 ~
        "PATIENT-SPECIFIC / exploratory - not recurrent MPNST biology",
      dominant_patient_fraction > 0.80 ~
        "PATIENT-DOMINATED (>80% one patient) / exploratory - not recurrent MPNST biology",
      n_patients_nontrivial >= 3 ~ "recurrent: >=3 patients each contribute >=5%",
      n_patients_nontrivial == 2 ~ "shared between 2 patients",
      TRUE ~ "PATIENT-DOMINATED / exploratory - not recurrent MPNST biology")) |>
  arrange(desc(n_cells))
print(as.data.frame(pd), row.names = FALSE)
write.table(pd, file.path(TAB, "TUMOR_STATE_PATIENT_DISTRIBUTION.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
n_recur <- sum(grepl("^recurrent", pd$recurrence))
n_shared <- sum(grepl("^shared", pd$recurrence))
n_private <- sum(grepl("PATIENT-SPECIFIC|PATIENT-DOMINATED", pd$recurrence))
log_("STATE RECURRENCE SUMMARY (§36)")
log_("  recurrent (>=3 patients, each >=5%) : ", n_recur, " of ", nrow(pd))
log_("  shared between 2 patients           : ", n_shared)
log_("  patient-specific or patient-dominated: ", n_private)
log_(sprintf("  fraction of malignant cells in patient-private states: %.3f",
             sum(pd$fraction_of_malignant[grepl("PATIENT-", pd$recurrence)])))
if (n_recur == 0)
  log_("  *** NO malignant transcriptional state is recurrent across patients. This is a")
  log_("      NEGATIVE RESULT and must be reported as one: with n = 4 and sample_id = patient")
  log_("      = dataset, residual patient structure dominates the malignant compartment. ***")

# --- 9. clone vs state (§33) ----------------------------------------------
sec("9. Clone versus transcriptional state")
cv <- obj@meta.data |> filter(!is.na(tumor_clone_phase4)) |>
  count(sample_id, tumor_clone_phase4, tumor_state_phase4)
if (nrow(cv)) {
  ct <- xtabs(n ~ tumor_clone_phase4 + tumor_state_phase4, data = cv)
  print(ct)
  # Cramer's V as a descriptive measure of association only - NOT a claim of
  # one-to-one correspondence, and not a hypothesis test.
  chi <- suppressWarnings(chisq.test(ct))
  V <- sqrt(as.numeric(chi$statistic) / (sum(ct) * (min(dim(ct)) - 1)))
  log_(sprintf("descriptive association (Cramer's V) between clone and state: %.3f", V))
  log_("This is descriptive. Clone and state are distinct concepts (§33) and a moderate")
  log_("association does not make them equivalent.")
  write.table(as.data.frame(ct), file.path(TAB, "CLONE_VS_TUMOR_STATE.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  clone_state_V <- V
} else { log_("no clone labels among malignant cells"); clone_state_V <- NA_real_ }

# --- 10. assignments + object ---------------------------------------------
sec("10. Outputs")
asg <- obj@meta.data |> mutate(cell_id = rownames(obj@meta.data)) |>
  select(cell_id, sample_id, annotation_ccc_phase3, annotation_ccc_refined,
         malignancy_confidence, malignancy_rule, primary_cluster,
         tumor_cluster_phase4, tumor_state_phase4, tumor_clone_phase4,
         cnv_burden, all_of(prog_cols))
write.table(asg, file.path(TAB, "TUMOR_STATE_ASSIGNMENTS.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(obj, file.path(TS, "malignant_only_object.rds"))
log_("wrote malignant-only object (NEW reductions only; Phase 2 reductions untouched)")

# --- 11. figures (§35) -----------------------------------------------------
sec("11. M31 figures")
em <- as.data.frame(Embeddings(obj, "malignant_umap")); colnames(em) <- c("U1","U2")
em <- cbind(em, obj@meta.data[rownames(em), ])
nstate <- length(unique(em$tumor_state_phase4))
spal <- setNames(colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(nstate),
                 sort(unique(em$tumor_state_phase4)))
base_cap <- paste(
  "Malignant cells only, defined by malignancy_refined - NOT by the broad Phase 2 MPNST-Tumor label.",
  "Coordinates are a NEW malignant-only reduction (malignant_pca -> malignant_harmony -> malignant_umap).",
  "The Phase 2 pca / postint_harmony / postint_umap_harmony reductions are frozen and were never modified.",
  sep = "\n")

g1 <- ggplot(em, aes(U1, U2, colour = tumor_state_phase4)) +
  geom_point(size = 0.6, alpha = 0.8, stroke = 0) +
  scale_colour_manual(values = spal, name = "tumor_state_phase4") +
  guides(colour = guide_legend(override.aes = list(size = 2.8, alpha = 1))) +
  coord_equal() + p4_theme() + theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(title = "MPNST malignant transcriptional states",
       subtitle = sprintf("%s malignant cells, resolution %.1f, %d states",
                          format(ncol(obj), big.mark = ","), CHOSEN, nstate),
       x = "malignant UMAP 1", y = "malignant UMAP 2", caption = base_cap)
save_fig(g1, file.path(FIG, "31_01_tumor_state_umap.pdf"), 8.6, 7.4)

g2 <- ggplot(em, aes(U1, U2, colour = sample_id)) +
  geom_point(size = 0.55, alpha = 0.8, stroke = 0) +
  scale_colour_manual(values = P4_SAMPLE, name = "patient") +
  guides(colour = guide_legend(override.aes = list(size = 2.8, alpha = 1))) +
  facet_wrap(~ tumor_state_phase4) + coord_equal() + p4_theme() +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(title = "Patient contribution to each malignant state",
       subtitle = "A state drawn from one patient is exploratory, not recurrent MPNST biology",
       x = "malignant UMAP 1", y = "malignant UMAP 2", caption = base_cap)
save_fig(g2, file.path(FIG, "31_02_tumor_state_patient_contribution.pdf"), 11.5, 9.2)

mc <- em |> filter(!is.na(tumor_clone_phase4))
if (nrow(mc)) {
  cpal <- setNames(colorRampPalette(RColorBrewer::brewer.pal(8,"Dark2"))(
            length(unique(mc$tumor_clone_phase4))), sort(unique(mc$tumor_clone_phase4)))
  g3 <- ggplot(em, aes(U1, U2)) + geom_point(colour = "grey88", size = 0.4, stroke = 0) +
    geom_point(data = mc, aes(colour = tumor_clone_phase4), size = 0.55, alpha = 0.85, stroke = 0) +
    scale_colour_manual(values = cpal, name = "SCEVAN clone") +
    guides(colour = guide_legend(override.aes = list(size = 2.8, alpha = 1), ncol = 1)) +
    coord_equal() + p4_theme() + theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    labs(title = "SCEVAN clones on the malignant-state UMAP",
         subtitle = "Clone = CNV architecture. State = transcriptional phenotype. These are different things (§33).",
         x = "malignant UMAP 1", y = "malignant UMAP 2", caption = base_cap)
  save_fig(g3, file.path(FIG, "31_03_scevan_clone_umap.pdf"), 9.2, 7.4)

  cs <- mc |> count(tumor_clone_phase4, tumor_state_phase4) |>
    group_by(tumor_clone_phase4) |> mutate(frac = n/sum(n)) |> ungroup()
  g4 <- ggplot(cs, aes(x = tumor_state_phase4, y = tumor_clone_phase4, fill = frac)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = n), size = 2.6,
              colour = ifelse(cs$frac > 0.55, "white", "grey15")) +
    scale_fill_gradient(low = "#F7FBFF", high = "#08519C",
                        labels = percent_format(accuracy = 1),
                        name = "fraction\nof clone") +
    p4_theme() + theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
    labs(title = "CNV clone versus transcriptional state",
         subtitle = sprintf("Cell counts printed. Descriptive association (Cramer's V) = %.3f — an association, not an equivalence.",
                            clone_state_V),
         x = "tumor_state_phase4", y = "tumor_clone_phase4")
  save_fig(g4, file.path(FIG, "31_04_clone_vs_tumor_state.pdf"), 10.5, 7.2)
}

tg <- top |> group_by(cluster) |> slice_head(n = 6) |> ungroup()
tg$label <- state_map$state[match(as.character(tg$cluster), state_map$cluster)]
g5 <- DotPlot(obj, features = unique(tg$gene), group.by = "tumor_state_phase4",
              assay = "RNA", cluster.idents = FALSE) +
  scale_colour_gradient2(low = "#2166AC", mid = "grey92", high = "#B2182B", midpoint = 0) +
  p4_theme() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6.4)) +
  labs(title = "Top markers per malignant state",
       subtitle = "Cluster characterisation for annotation — NOT condition-level differential expression",
       x = NULL, y = NULL, caption = base_cap)
save_fig(g5, file.path(FIG, "31_05_tumor_state_marker_dotplot.pdf"), 17.5, 5.8)

pl <- pm |> pivot_longer(-cluster, names_to = "programme", values_to = "score") |>
  mutate(programme = sub("^prog_", "", programme),
         state = state_map$state[match(as.character(cluster), state_map$cluster)])
g6 <- ggplot(pl, aes(x = programme, y = state, fill = score)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", score)), size = 2.5,
            colour = ifelse(pl$score > max(pl$score)*0.6, "white", "grey15")) +
  scale_fill_gradient2(low = "#2166AC", mid = "grey95", high = "#B2182B", midpoint = 0,
                       name = "median\nmodule score") +
  p4_theme() + theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(title = "Candidate biological programmes per malignant state",
       subtitle = paste("Programmes are CANDIDATE concepts tested against the data, not labels imposed on it.",
                        "\nA state whose top two programmes are within 0.02 is called Uncertain."),
       x = NULL, y = NULL, caption = base_cap)
save_fig(g6, file.path(FIG, "31_06_tumor_state_program_heatmap.pdf"), 11.5, 6.4)

topm <- mk |> group_by(cluster) |> slice_head(n = 8) |> ungroup()
oh <- subset(obj, features = unique(topm$gene))
oh <- ScaleData(oh, features = unique(topm$gene), verbose = FALSE)
g7 <- DoHeatmap(oh, features = unique(topm$gene), group.by = "tumor_state_phase4",
                size = 3, angle = 30) +
  scale_fill_gradient2(low = "#2166AC", mid = "grey95", high = "#B2182B", midpoint = 0,
                       name = "scaled\nexpression") +
  theme(axis.text.y = element_text(size = 5.4)) +
  labs(title = "Top marker heatmap by malignant state", caption = base_cap)
save_fig(g7, file.path(FIG, "31_07_tumor_state_marker_heatmap.pdf"), 13.5, 11.5)

pcomp <- pd |> select(tumor_state_phase4, patient_breakdown) |>
  separate_rows(patient_breakdown, sep = "; ") |>
  separate(patient_breakdown, into = c("sample_id","n"), sep = "=", convert = TRUE) |>
  group_by(tumor_state_phase4) |> mutate(frac = n/sum(n)) |> ungroup()
g8 <- ggplot(pcomp, aes(y = tumor_state_phase4, x = frac, fill = sample_id)) +
  geom_col(width = 0.78) +
  scale_fill_manual(values = P4_SAMPLE, name = "patient") +
  scale_x_continuous(labels = percent_format(accuracy = 1), expand = expansion(0)) +
  p4_theme() +
  labs(title = "Patient composition of each malignant state",
       subtitle = "Read alongside TUMOR_STATE_PATIENT_DISTRIBUTION.tsv, which flags single-patient states explicitly",
       x = "fraction of state", y = NULL, caption = base_cap)
save_fig(g8, file.path(FIG, "31_08_tumor_state_patient_composition.pdf"), 10.5, 5.4)

# --- 12. provenance -------------------------------------------------------
facts <- list(
  malignant_cells = ncol(obj), total_cells = nrow(calls),
  malignant_fraction = ncol(obj)/nrow(calls),
  by_confidence = as.list(table(mm$malignancy_confidence)),
  by_patient = as.list(table(obj$sample_id)),
  phase2_origin = as.list(sort(table(mm$annotation_ccc_phase3), decreasing = TRUE)),
  reductions_created = c("malignant_pca","malignant_harmony","malignant_umap"),
  phase2_reductions_modified = FALSE,
  dims_used = NDIM, npcs_computed = NPCS, variance_explained_at_ndim = cum[NDIM],
  knn_overlap_phase2_vs_malignant_harmony = ov,
  harmony_batch_variable = "sample_id",
  resolution_sweep = sw, resolution_selected = CHOSEN,
  n_states = nstate, states = as.list(table(obj$tumor_state_phase4)),
  state_recurrence_summary = list(recurrent = n_recur, shared_2_patients = n_shared,
    patient_private_or_dominated = n_private,
    fraction_malignant_in_patient_private_states =
      sum(pd$fraction_of_malignant[grepl("PATIENT-", pd$recurrence)])),
  state_definitions = state_map,
  clone_state_cramers_v = clone_state_V,
  programme_panels_present = present,
  markers_significant = nrow(mk),
  expression_basis = "RNA counts (M28 extraction) -> LogNormalize 1e4; SCT not used (four models)",
  slurm = list(job_id = Sys.getenv("SLURM_JOB_ID"), cpus = CORES),
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
write_json(facts, file.path(TS, "m31_tumor_state_facts.json"),
           auto_unbox = TRUE, pretty = TRUE, digits = 8, null = "null")
sec("M31 COMPLETE")
