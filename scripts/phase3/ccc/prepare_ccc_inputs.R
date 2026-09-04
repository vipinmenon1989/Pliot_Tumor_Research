# scripts/phase3/ccc/prepare_ccc_inputs.R
#
# Phase 3 / Milestone M19 - sample-aware CCC input preparation.
#
# Builds standardised, sample-aware inputs for every CCC framework from the frozen Phase 2
# final object, and establishes the documented minimum-cell policy that decides which
# sender/receiver pairs are EVALUABLE in which sample.
#
# EXPRESSION BASIS DECISION (documented, not silent):
#   CCC uses the RNA assay with JoinLayers() + LogNormalize, NOT the SCT assay.
#   Reason: the SCT assay carries FOUR SCTransform models (Phase 1 M8 ran SCTransform
#   layer-wise), so SCT values are not normalised on a common footing across samples.
#   Every CCC framework expects a single uniformly log-normalised matrix. LogNormalize on
#   joined RNA counts gives exactly that. Harmony coordinates are NEVER used as expression.
#
# Absence of evidence is not evidence of absence: a population below the minimum in a given
# sample is recorded as NOT EVALUABLE for that sample rather than as "no signalling".

options(stringsAsFactors=FALSE); options(future.globals.maxSize=+Inf)
suppressPackageStartupMessages({library(Seurat); library(SeuratObject); library(Matrix)
  library(ggplot2); library(patchwork); library(RColorBrewer); library(jsonlite); library(digest)})
source("scripts/R/logging_utils.R"); source("scripts/R/provenance_utils.R"); setup_strict_logging()
STAGE <- "phase3_m19_ccc_inputs"

args <- commandArgs(trailingOnly=TRUE)
input_rds  <- "results/phase2/phase2_final_object.rds"
out_dir    <- "results/phase3/ccc"
tab_dir    <- "results/phase3/tables"
fig_dir    <- "results/phase3/figures/M19"
cpdb_dir   <- "results/phase3/ccc/cellphonedb_inputs"
min_cells  <- 10L      # per population per sample
min_pops   <- 3L       # populations needed for a sample to be analysable
seed <- 42L; expected_cells <- 19716L
i <- 1
while (i <= length(args)) { a <- args[i]
  if (a=="--input"){input_rds<-args[i+1];i<-i+2} else if (a=="--out-dir"){out_dir<-args[i+1];i<-i+2}
  else if (a=="--tab-dir"){tab_dir<-args[i+1];i<-i+2} else if (a=="--fig-dir"){fig_dir<-args[i+1];i<-i+2}
  else if (a=="--cpdb-dir"){cpdb_dir<-args[i+1];i<-i+2}
  else if (a=="--min-cells"){min_cells<-as.integer(args[i+1]);i<-i+2}
  else if (a=="--random-seed"){seed<-as.integer(args[i+1]);i<-i+2}
  else stop(sprintf("Unknown argument: %s", a)) }

t0 <- Sys.time(); warns <- character(0)
wrn <- function(m){warns<<-c(warns,m); log_warn(m,stage=STAGE)}
fail <- function(m){log_error(m,stage=STAGE); stop(m,call.=FALSE)}
for (d in c(out_dir,tab_dir,fig_dir,cpdb_dir)) dir.create(d,recursive=TRUE,showWarnings=FALSE)
set.seed(seed)

log_info("============= M19 SAMPLE-AWARE CCC INPUT PREPARATION =============", stage=STAGE)
in_md5 <- digest(input_rds, file=TRUE, algo="md5")
log_info(sprintf("Input %s (md5 %s)", input_rds, in_md5), stage=STAGE)
if (!identical(in_md5, "153d5f6acc70f9c05aa48cabc4f4ac2d"))
  wrn(sprintf("Phase 2 final object md5 is %s, expected 153d5f6acc70f9c05aa48cabc4f4ac2d.", in_md5))
obj <- readRDS(input_rds)
log_info(sprintf("Loaded %d features x %d cells | assays %s | default %s",
  nrow(obj), ncol(obj), paste(Assays(obj),collapse=","), DefaultAssay(obj)), stage=STAGE)
if (ncol(obj) != expected_cells) fail("Cell count mismatch with Phase 2.")
for (cc in c("annotation_ccc","annotation_ccc_compartment","annotation_ccc_ccc_ready",
             "sample_id","postint_celltype_level2_refined"))
  if (!(cc %in% colnames(obj@meta.data))) fail(sprintf("Required column missing: %s", cc))

md <- obj@meta.data
ccc <- as.character(md$annotation_ccc); grp <- as.character(md$sample_id)
ready <- md$annotation_ccc_ccc_ready
lv <- sort(unique(grp))
log_info(sprintf("Samples (= datasets = patients): %s", paste(lv, collapse=", ")), stage=STAGE)
log_info(sprintf("CCC identities: %d | flagged not CCC-ready: %s", length(unique(ccc)),
  paste(sort(unique(ccc[!ready])), collapse=", ")), stage=STAGE)

# ---------------------------------------------------------------------------
# 1. Expression basis: RNA, joined, LogNormalized
# ---------------------------------------------------------------------------
log_info("Preparing the expression basis: RNA assay -> JoinLayers -> LogNormalize ...", stage=STAGE)
DefaultAssay(obj) <- "RNA"
lay_before <- SeuratObject::Layers(obj[["RNA"]])
log_info(sprintf("  RNA layers before join: %s", paste(lay_before, collapse=";")), stage=STAGE)
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
lay_after <- SeuratObject::Layers(obj[["RNA"]])
log_info(sprintf("  RNA layers after join : %s", paste(lay_after, collapse=";")), stage=STAGE)
if (!("counts" %in% lay_after)) fail("JoinLayers did not produce a single 'counts' layer.")
obj <- NormalizeData(obj, assay="RNA", normalization.method="LogNormalize",
                     scale.factor=1e4, verbose=FALSE)
log_info("  LogNormalize complete (scale.factor = 1e4).", stage=STAGE)
gene_space <- rownames(obj[["RNA"]])
log_info(sprintf("  RNA gene space: %d features", length(gene_space)), stage=STAGE)

# ---------------------------------------------------------------------------
# 2. Minimum-cell policy and the sample x population count table
# ---------------------------------------------------------------------------
log_info("=== MINIMUM-CELL POLICY ===", stage=STAGE)
log_info(sprintf("  (a) Populations flagged not CCC-ready in Phase 2 are EXCLUDED outright: %s",
  paste(sort(unique(ccc[!ready])), collapse=", ")), stage=STAGE)
log_info(sprintf("  (b) A population is EVALUABLE within a sample only if it has >= %d cells in that sample.", min_cells), stage=STAGE)
log_info("      Basis: %d matches CellChat's `min.cells` default and is the conventional floor for", stage=STAGE)
log_info("      per-group mean-expression estimates in LR inference. It is stated, not silently chosen.", stage=STAGE)
log_info(sprintf("  (c) A sample is ANALYSABLE only if >= %d populations are evaluable in it.", min_pops), stage=STAGE)
log_info("  (d) A population below the minimum is recorded NOT EVALUABLE, never as 'no signalling'.", stage=STAGE)

keep_pop <- sort(unique(ccc[ready]))
ctab <- table(factor(ccc, levels=sort(unique(ccc))), factor(grp, levels=lv))
counts_df <- do.call(rbind, lapply(rownames(ctab), function(p) do.call(rbind, lapply(lv, function(s) {
  n <- as.integer(ctab[p, s])
  data.frame(annotation_ccc=p, sample_id=s, patient=s, dataset=s, n_cells=n,
    phase2_ccc_ready=p %in% keep_pop,
    evaluable_in_sample=(p %in% keep_pop) && n >= min_cells,
    reason=if (!(p %in% keep_pop)) "excluded: flagged not CCC-ready in Phase 2"
           else if (n < min_cells) sprintf("not evaluable: %d < %d cells in this sample", n, min_cells)
           else "evaluable",
    stringsAsFactors=FALSE) }))))
write.table(counts_df, file.path(tab_dir,"CCC_SAMPLE_CELL_COUNTS.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

ev <- counts_df[counts_df$evaluable_in_sample, ]
per_sample_pops <- table(ev$sample_id)
analysable <- names(per_sample_pops)[per_sample_pops >= min_pops]
log_info(sprintf("Evaluable populations per sample: %s",
  paste(sprintf("%s=%d", names(per_sample_pops), as.integer(per_sample_pops)), collapse=", ")), stage=STAGE)
log_info(sprintf("Analysable samples: %s", paste(analysable, collapse=", ")), stage=STAGE)
if (!length(analysable)) fail("No sample is analysable under the minimum-cell policy.")

pop_summary <- do.call(rbind, lapply(keep_pop, function(p) {
  r <- counts_df[counts_df$annotation_ccc==p, ]
  data.frame(annotation_ccc=p, total_cells=sum(r$n_cells),
    samples_evaluable=sum(r$evaluable_in_sample), samples_total=length(lv),
    min_cells_in_a_sample=min(r$n_cells), max_cells_in_a_sample=max(r$n_cells),
    median_cells_per_sample=as.numeric(median(r$n_cells)),
    evaluable_samples=paste(r$sample_id[r$evaluable_in_sample], collapse=";"),
    not_evaluable_samples=paste(r$sample_id[!r$evaluable_in_sample], collapse=";"),
    stringsAsFactors=FALSE) }))
pop_summary <- pop_summary[order(-pop_summary$total_cells), ]
write.table(pop_summary, file.path(tab_dir,"CCC_POPULATION_EVALUABILITY.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
log_info("Population evaluability:", stage=STAGE)
for (r in seq_len(nrow(pop_summary)))
  log_info(sprintf("  %-30s total=%5d  evaluable in %d/%d samples  (min %d per sample)",
    pop_summary$annotation_ccc[r], pop_summary$total_cells[r], pop_summary$samples_evaluable[r],
    pop_summary$samples_total[r], pop_summary$min_cells_in_a_sample[r]), stage=STAGE)

# directed pair evaluability, tumour-centric first
pairs <- expand.grid(sender=keep_pop, receiver=keep_pop, stringsAsFactors=FALSE)
pairs <- pairs[pairs$sender != pairs$receiver, ]
pair_rows <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
  s <- pairs$sender[i]; r <- pairs$receiver[i]
  okv <- vapply(lv, function(x)
    isTRUE(counts_df$evaluable_in_sample[counts_df$annotation_ccc==s & counts_df$sample_id==x]) &&
    isTRUE(counts_df$evaluable_in_sample[counts_df$annotation_ccc==r & counts_df$sample_id==x]), logical(1))
  prio <- if (s=="MPNST-Tumor" && r!="MPNST-Tumor") "P1_tumor_to_TME"
          else if (r=="MPNST-Tumor" && s!="MPNST-Tumor") "P1_TME_to_tumor" else "P2_TME_to_TME"
  data.frame(sender=s, receiver=r, priority=prio,
    samples_evaluable=sum(okv), samples_total=length(lv),
    evaluable_samples=paste(lv[okv], collapse=";"), stringsAsFactors=FALSE) }))
pair_rows <- pair_rows[order(pair_rows$priority, -pair_rows$samples_evaluable), ]
write.table(pair_rows, file.path(tab_dir,"CCC_PAIR_EVALUABILITY.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
log_info(sprintf("Directed pairs: %d total | tumour->TME %d | TME->tumour %d | TME<->TME %d",
  nrow(pair_rows), sum(pair_rows$priority=="P1_tumor_to_TME"),
  sum(pair_rows$priority=="P1_TME_to_tumor"), sum(pair_rows$priority=="P2_TME_to_TME")), stage=STAGE)
log_info(sprintf("Pairs evaluable in all 4 samples: %d | in >=3: %d | in >=2: %d",
  sum(pair_rows$samples_evaluable==4), sum(pair_rows$samples_evaluable>=3),
  sum(pair_rows$samples_evaluable>=2)), stage=STAGE)

# ---------------------------------------------------------------------------
# 3. Gene-symbol / expression QC
# ---------------------------------------------------------------------------
log_info("Gene-symbol and expression QC ...", stage=STAGE)
dup_genes <- sum(duplicated(gene_space))
if (dup_genes) wrn(sprintf("%d duplicated gene symbols in the RNA assay.", dup_genes))
probe <- c("CSF1","CSF1R","TGFB1","TGFBR1","TGFBR2","VEGFA","KDR","FLT1","CXCL12","CXCR4",
           "IL34","SPP1","CD44","ITGB1","THBS1","PDGFB","PDGFRB","JAG1","NOTCH1","NOTCH3",
           "TNF","IFNG","IL1B","IL10","CCL2","CCR2","CXCL9","CXCL10","CXCR3","PDCD1","CD274",
           "LGALS9","HAVCR2","TIGIT","NECTIN2","MIF","CD74","APP","CD47","SIRPA","GAS6","AXL")
found <- intersect(probe, gene_space)
log_info(sprintf("  Ligand/receptor probe genes present: %d/%d (%s)", length(found), length(probe),
  paste(setdiff(probe, gene_space), collapse=",")), stage=STAGE)
if (length(found) < 0.7*length(probe)) wrn("Fewer than 70% of probe LR genes found; check gene symbols.")

expr_rows <- list()
dat <- GetAssayData(obj, assay="RNA", layer="data")
for (p in keep_pop) for (s in lv) {
  sel <- ccc==p & grp==s
  if (sum(sel) < min_cells) next
  sub <- dat[, sel, drop=FALSE]
  detected <- Matrix::rowSums(sub > 0)
  expr_rows[[paste(p,s)]] <- data.frame(annotation_ccc=p, sample_id=s, n_cells=sum(sel),
    genes_detected_in_ge10pct=sum(detected >= 0.10*sum(sel)),
    genes_detected_at_all=sum(detected > 0),
    median_genes_per_cell=as.numeric(median(Matrix::colSums(sub > 0))),
    stringsAsFactors=FALSE)
}
expr_qc <- do.call(rbind, expr_rows)
write.table(expr_qc, file.path(tab_dir,"CCC_INPUT_EXPRESSION_QC.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# ---------------------------------------------------------------------------
# 4. Build the CCC-ready object and per-sample exports
# ---------------------------------------------------------------------------
log_info("Building the CCC-ready object (CCC-ready populations only) ...", stage=STAGE)
keep_cells <- colnames(obj)[ccc %in% keep_pop]
ccc_obj <- subset(obj, cells=keep_cells)
ccc_obj$ccc_label <- factor(as.character(ccc_obj$annotation_ccc), levels=keep_pop)
ccc_obj$ccc_sample <- as.character(ccc_obj$sample_id)
DefaultAssay(ccc_obj) <- "RNA"
# CCC frameworks only need RNA; dropping SCT keeps the exported object small and prevents
# any framework from silently picking up multi-model SCT values.
if ("SCT" %in% Assays(ccc_obj)) ccc_obj[["SCT"]] <- NULL
for (r in setdiff(Reductions(ccc_obj), c("postint_harmony","postint_umap_harmony"))) ccc_obj[[r]] <- NULL
log_info(sprintf("CCC object: %d features x %d cells | %d populations | assays %s | reductions %s",
  nrow(ccc_obj), ncol(ccc_obj), nlevels(ccc_obj$ccc_label),
  paste(Assays(ccc_obj),collapse=","), paste(Reductions(ccc_obj),collapse=",")), stage=STAGE)
ccc_rds <- file.path(out_dir,"ccc_input_object.rds")
saveRDS(ccc_obj, ccc_rds)
ccc_md5 <- digest(ccc_rds, file=TRUE, algo="md5")
log_info(sprintf("Saved %s (%.2f GB, md5 %s)", ccc_rds, file.info(ccc_rds)$size/1024^3, ccc_md5), stage=STAGE)

log_info("Writing per-sample CellPhoneDB inputs (sparse MTX + meta) ...", stage=STAGE)
for (s in analysable) {
  sd <- file.path(cpdb_dir, s); dir.create(sd, recursive=TRUE, showWarnings=FALSE)
  cells_s <- colnames(ccc_obj)[ccc_obj$ccc_sample == s]
  lab_s <- as.character(ccc_obj$ccc_label[match(cells_s, colnames(ccc_obj))])
  okp <- names(which(table(lab_s) >= min_cells))
  cells_s <- cells_s[lab_s %in% okp]; lab_s <- lab_s[lab_s %in% okp]
  m <- GetAssayData(ccc_obj, assay="RNA", layer="data")[, cells_s, drop=FALSE]
  Matrix::writeMM(as(m, "dgCMatrix"), file.path(sd,"matrix.mtx"))
  writeLines(colnames(m), file.path(sd,"barcodes.tsv"))
  writeLines(rownames(m), file.path(sd,"features.tsv"))
  write.table(data.frame(Cell=cells_s, cell_type=lab_s), file.path(sd,"meta.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
  log_info(sprintf("  %s: %d cells x %d genes, %d populations -> %s", s, ncol(m), nrow(m), length(okp), sd), stage=STAGE)
}

# ---------------------------------------------------------------------------
# 5. FIGURES
# ---------------------------------------------------------------------------
log_info("Generating M19 QC figures ...", stage=STAGE)
frows <- list()
sf <- function(p,name,w,h,method,params){ for (ext in c("pdf","png")) {
  fp <- file.path(fig_dir, paste0(name,".",ext))
  ggsave(fp,p,width=w,height=h,units="in",dpi=200,device=ext,limitsize=FALSE)
  frows[[length(frows)+1]] <<- data.frame(figure_path=fp, phase="phase3", milestone="M19",
    analysis="ccc_input_qc", sender=NA_character_, receiver=NA_character_, method="input_prep",
    input=input_rds, input_checksum=in_md5, script="scripts/phase3/ccc/prepare_ccc_inputs.R",
    parameters=params, stringsAsFactors=FALSE) }
  log_info(sprintf("  wrote %s", name), stage=STAGE) }
th <- theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(), plot.title=element_text(face="bold"),
      plot.subtitle=element_text(size=9,colour="grey30"))
cd <- counts_df; cd$annotation_ccc <- factor(cd$annotation_ccc,
  levels=rev(c(keep_pop, setdiff(unique(cd$annotation_ccc), keep_pop))))
sf(ggplot(cd, aes(sample_id, annotation_ccc, fill=log10(n_cells+1)))+geom_tile(colour="white")+
   geom_text(aes(label=n_cells, colour=evaluable_in_sample), size=3.1, fontface="bold")+
   scale_fill_viridis_c(option="mako", name="log10(cells+1)")+
   scale_colour_manual(values=c(`TRUE`="white",`FALSE`="#FF6666"), name=sprintf("evaluable (>=%d cells)", min_cells))+
   labs(title="Cells per sample x CCC population",
        subtitle=sprintf("Red counts are NOT EVALUABLE. Populations below the line were excluded outright as not CCC-ready in Phase 2. sample = dataset = patient."),
        x=NULL,y=NULL)+th,
   "M19_01_sample_x_population_counts", 9, 8, "heatmap", sprintf("min_cells=%d", min_cells))
sf(ggplot(pop_summary, aes(reorder(annotation_ccc, total_cells), total_cells, fill=factor(samples_evaluable)))+
   geom_col()+coord_flip()+scale_fill_brewer(palette="Greens", name="samples\nevaluable")+
   geom_text(aes(label=total_cells), hjust=-0.1, size=3)+expand_limits(y=max(pop_summary$total_cells)*1.2)+
   labs(title="CCC-ready population sizes and sample coverage",
        subtitle="Only populations passing the Phase 2 readiness audit are shown", x=NULL,y="cells")+th,
   "M19_02_population_evaluability", 9, 6.5, "barplot", sprintf("min_cells=%d", min_cells))
pr <- pair_rows[pair_rows$priority != "P2_TME_to_TME", ]
sf(ggplot(pr, aes(sender, receiver, fill=factor(samples_evaluable)))+geom_tile(colour="white")+
   scale_fill_brewer(palette="Blues", name="samples\nevaluable")+
   facet_wrap(~priority, scales="free")+
   theme_bw(base_size=11)+theme(axis.text.x=element_text(angle=45,hjust=1), panel.grid=element_blank())+
   labs(title="Tumour-centric directed pair evaluability",
        subtitle="How many of the four samples support each directed sender -> receiver pair", x="sender",y="receiver"),
   "M19_03_tumor_pair_evaluability", 13, 7, "heatmap", sprintf("min_cells=%d", min_cells))
sf(ggplot(expr_qc, aes(n_cells, genes_detected_in_ge10pct, colour=sample_id, label=annotation_ccc))+
   geom_point(size=2.4)+scale_x_log10()+scale_colour_brewer(palette="Set1")+
   labs(title="Detected gene breadth versus population size",
        subtitle="Genes detected in >=10% of cells, per sample x population. Small populations detect fewer genes, which limits LR sensitivity.",
        x="cells in population (log10)", y="genes detected in >=10% of cells")+th,
   "M19_04_expression_breadth", 9, 6.5, "scatter", "detection>=10%")
write.table(do.call(rbind, frows), file.path(out_dir,"figure_index_m19.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

rec <- list(milestone="M19", phase="phase3", timestamp=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),
  script="scripts/phase3/ccc/prepare_ccc_inputs.R",
  input=list(path=input_rds, md5=in_md5, cells=expected_cells),
  expression_basis=list(assay="RNA", layers_joined=TRUE, layers_before=lay_before, layers_after=lay_after,
    normalization="LogNormalize", scale_factor=1e4,
    rationale="The SCT assay carries four SCTransform models so SCT values are not on a common footing across samples. CCC frameworks expect one uniformly log-normalised matrix. Harmony coordinates are never used as expression."),
  minimum_cell_policy=list(min_cells_per_population_per_sample=min_cells,
    min_populations_per_sample=min_pops,
    basis=sprintf("%d matches CellChat's min.cells default and the conventional floor for per-group mean-expression estimates in LR inference", min_cells),
    excluded_populations=sort(unique(ccc[!ready])),
    absence_rule="A population below the minimum is recorded NOT EVALUABLE for that sample, never as absence of signalling."),
  populations=keep_pop, n_populations=length(keep_pop),
  samples=lv, analysable_samples=analysable,
  pairs=list(total=nrow(pair_rows), tumor_to_TME=sum(pair_rows$priority=="P1_tumor_to_TME"),
    TME_to_tumor=sum(pair_rows$priority=="P1_TME_to_tumor"), TME_to_TME=sum(pair_rows$priority=="P2_TME_to_TME"),
    evaluable_in_all_4=sum(pair_rows$samples_evaluable==4)),
  outputs=list(ccc_object=ccc_rds, ccc_object_md5=ccc_md5,
    cellphonedb_inputs=cpdb_dir,
    tables=file.path(tab_dir, c("CCC_SAMPLE_CELL_COUNTS.tsv","CCC_POPULATION_EVALUABILITY.tsv",
      "CCC_PAIR_EVALUABILITY.tsv","CCC_INPUT_EXPRESSION_QC.tsv"))),
  gene_space=list(n_features=length(gene_space), duplicated_symbols=dup_genes,
    lr_probe_found=length(found), lr_probe_total=length(probe),
    lr_probe_missing=setdiff(probe, gene_space)),
  warnings=warns, slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM"),
  git_commit=tryCatch(trimws(system("git rev-parse HEAD",intern=TRUE)), error=function(e) NA_character_),
  elapsed_seconds=as.numeric(difftime(Sys.time(),t0,units="secs")), session_info=capture.output(sessionInfo()))
write_json(rec, file.path(out_dir,"m19_ccc_input_record.json"), auto_unbox=TRUE, pretty=TRUE, null="null", digits=NA)
prov <- record_provenance("phase3_m19_ccc_inputs", inputs=list(phase2_final=input_rds),
  outputs=list(ccc_object=ccc_rds, counts=file.path(tab_dir,"CCC_SAMPLE_CELL_COUNTS.tsv"),
    record=file.path(out_dir,"m19_ccc_input_record.json")),
  parameters=list(milestone="M19", min_cells=min_cells, seed=seed,
    slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM")), dataset="combined")
save_provenance_json(prov, file.path(out_dir,"prov_m19_ccc_inputs.json"))
log_info(sprintf("M19 COMPLETE in %.1f min. Warnings: %d", as.numeric(difftime(Sys.time(),t0,units="mins")), length(warns)), stage=STAGE)
log_system_usage(stage=STAGE)
