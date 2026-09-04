# scripts/phase3/ccc/run_cellchat.R
#
# Phase 3 / Milestone M20 - CellChat execution, sample-aware.
#
# CellChat is used as an INDEPENDENT pathway- and network-oriented framework, not as another
# LR list generator. Its distinctive outputs are pathway-level communication, outgoing and
# incoming signalling roles, and information flow.
#
# CellChat communication probability is NOT numerically comparable to a LIANA rank or a
# CellPhoneDB p-value. It contributes only a boolean support flag to the M21 concordance.
#
# Expression: RNA assay, joined, LogNormalize (prepared in M19).

options(stringsAsFactors=FALSE); options(future.globals.maxSize=+Inf)
suppressPackageStartupMessages({library(Seurat); library(SeuratObject); library(Matrix)
  library(CellChat); library(dplyr); library(jsonlite); library(digest); library(future)})
source("scripts/R/logging_utils.R"); source("scripts/R/provenance_utils.R"); setup_strict_logging()
STAGE <- "phase3_m20_cellchat"

args <- commandArgs(trailingOnly=TRUE)
input_rds <- "results/phase3/ccc/ccc_input_object.rds"
out_dir   <- "results/phase3/ccc/by_method"
min_cells <- 10L; seed <- 42L; nworkers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","4"))
only_samples <- NULL   # when set, process just these samples (one SLURM job per sample)
tag <- ""
i <- 1
while (i <= length(args)) { a <- args[i]
  if (a=="--input"){input_rds<-args[i+1];i<-i+2} else if (a=="--out-dir"){out_dir<-args[i+1];i<-i+2}
  else if (a=="--min-cells"){min_cells<-as.integer(args[i+1]);i<-i+2}
  else if (a=="--random-seed"){seed<-as.integer(args[i+1]);i<-i+2}
  else if (a=="--samples"){only_samples<-unlist(strsplit(args[i+1],","));i<-i+2}
  else if (a=="--tag"){tag<-args[i+1];i<-i+2}
  else stop(sprintf("Unknown argument: %s", a)) }

t0 <- Sys.time(); warns <- character(0)
wrn <- function(m){warns<<-c(warns,m); log_warn(m,stage=STAGE)}
fail <- function(m){log_error(m,stage=STAGE); stop(m,call.=FALSE)}
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
set.seed(seed)
future::plan("multisession", workers=max(1, nworkers-1))

log_info("================ M20 CellChat (sample-aware) ================", stage=STAGE)
log_info(sprintf("CellChat %s | workers %d | seed %d", packageVersion("CellChat"), nworkers-1, seed), stage=STAGE)
db <- CellChatDB.human
log_info(sprintf("CellChatDB.human: %d interactions, %d pathways, categories: %s",
  nrow(db$interaction), length(unique(db$interaction$pathway_name)),
  paste(names(table(db$interaction$annotation)), collapse=", ")), stage=STAGE)

in_md5 <- digest(input_rds, file=TRUE, algo="md5")
obj <- readRDS(input_rds); DefaultAssay(obj) <- "RNA"
log_info(sprintf("Input %s (md5 %s): %d x %d", input_rds, in_md5, nrow(obj), ncol(obj)), stage=STAGE)
lv <- sort(unique(as.character(obj$ccc_sample)))
if (!is.null(only_samples)) {
  lv <- intersect(lv, only_samples)
  if (!length(lv)) fail("None of the requested samples are present.")
  log_info(sprintf("Restricted to samples: %s (parallel per-sample execution)", paste(lv, collapse=", ")), stage=STAGE)
}
sfx <- if (nzchar(tag)) paste0("_", tag) else ""

res_all <- list(); path_all <- list(); centr_all <- list(); cc_objs <- list()
for (s in lv) {
  log_info(sprintf("--- sample %s ---", s), stage=STAGE)
  cells <- colnames(obj)[as.character(obj$ccc_sample)==s]
  sub <- subset(obj, cells=cells)
  keep <- names(which(table(as.character(sub$ccc_label)) >= min_cells))
  sub <- subset(sub, cells=colnames(sub)[as.character(sub$ccc_label) %in% keep])
  lab <- factor(as.character(sub$ccc_label), levels=sort(keep))
  dat <- GetAssayData(sub, assay="RNA", layer="data")
  meta <- data.frame(ccc_label=lab, row.names=colnames(sub))
  log_info(sprintf("  %d cells, %d populations", ncol(sub), length(keep)), stage=STAGE)

  cc <- tryCatch({
    x <- createCellChat(object=dat, meta=meta, group.by="ccc_label")
    x@DB <- db
    x <- subsetData(x)
    x <- identifyOverExpressedGenes(x)
    x <- identifyOverExpressedInteractions(x)
    # population.size=TRUE corrects for unequal population abundance, which matters here
    # because populations range from 106 to 5,064 cells.
    x <- computeCommunProb(x, raw.use=TRUE, population.size=TRUE, seed.use=seed)
    x <- filterCommunication(x, min.cells=min_cells)
    x <- computeCommunProbPathway(x)
    x <- aggregateNet(x)
    x <- netAnalysis_computeCentrality(x, slot.name="netP")
    x }, error=function(e) { wrn(sprintf("CellChat failed for %s: %s", s, conditionMessage(e))); NULL })
  if (is.null(cc)) next

  lrdf <- tryCatch(subsetCommunication(cc), error=function(e) NULL)
  if (!is.null(lrdf) && nrow(lrdf)) { lrdf$sample_id <- s; res_all[[s]] <- lrdf }
  pdf_ <- tryCatch(subsetCommunication(cc, slot.name="netP"), error=function(e) NULL)
  if (!is.null(pdf_) && nrow(pdf_)) { pdf_$sample_id <- s; path_all[[s]] <- pdf_ }
  # outgoing / incoming signalling roles per population
  cen <- tryCatch({
    m <- cc@netP$centr
    do.call(rbind, lapply(names(m), function(pw) data.frame(
      sample_id=s, pathway=pw, population=names(m[[pw]]$outdeg),
      outdeg=as.numeric(m[[pw]]$outdeg), indeg=as.numeric(m[[pw]]$indeg),
      stringsAsFactors=FALSE))) }, error=function(e) NULL)
  if (!is.null(cen)) centr_all[[s]] <- cen
  cc_objs[[s]] <- cc
  log_info(sprintf("  %d LR interactions, %d pathway-level rows",
    if (is.null(lrdf)) 0 else nrow(lrdf), if (is.null(pdf_)) 0 else nrow(pdf_)), stage=STAGE)
}
if (!length(res_all)) fail("CellChat produced no results for any sample.")
saveRDS(cc_objs, file.path(out_dir, sprintf("cellchat_objects%s.rds", sfx)))

cc_df <- dplyr::bind_rows(res_all)
colnames(cc_df)[colnames(cc_df)=="source"] <- "source"
cc_df$cellchat_supported <- TRUE   # subsetCommunication returns only significant links (pval < 0.05)
write.table(cc_df, file.path(out_dir, sprintf("all_CellChat_interactions%s.tsv", sfx)), sep="\t", row.names=FALSE, quote=FALSE)
log_info(sprintf("Pooled CellChat LR rows: %d across %d samples", nrow(cc_df), length(res_all)), stage=STAGE)
log_info(sprintf("Columns: %s", paste(colnames(cc_df), collapse=", ")), stage=STAGE)

if (length(path_all)) {
  pw_df <- dplyr::bind_rows(path_all)
  write.table(pw_df, file.path(out_dir, sprintf("CellChat_pathway_communication%s.tsv", sfx)), sep="\t", row.names=FALSE, quote=FALSE)
  log_info(sprintf("Pathway-level rows: %d | distinct pathways: %d", nrow(pw_df), length(unique(pw_df$pathway_name))), stage=STAGE)
}
if (length(centr_all)) {
  cen_df <- dplyr::bind_rows(centr_all)
  write.table(cen_df, file.path(out_dir, sprintf("CellChat_signaling_roles%s.tsv", sfx)), sep="\t", row.names=FALSE, quote=FALSE)
  tum <- cen_df[cen_df$population=="MPNST-Tumor", ]
  log_info(sprintf("MPNST-Tumor signalling role rows: %d | mean outdeg %.4f | mean indeg %.4f",
    nrow(tum), mean(tum$outdeg), mean(tum$indeg)), stage=STAGE)
}
per_sample <- as.data.frame(table(sample_id=cc_df$sample_id))
log_info(sprintf("Interactions per sample: %s", paste(sprintf("%s=%d", per_sample$sample_id, per_sample$Freq), collapse=", ")), stage=STAGE)
log_info(sprintf("Tumour-involving: %d (tumour->TME %d, TME->tumour %d)",
  sum(cc_df$source=="MPNST-Tumor" | cc_df$target=="MPNST-Tumor"),
  sum(cc_df$source=="MPNST-Tumor"), sum(cc_df$target=="MPNST-Tumor")), stage=STAGE)

rec <- list(milestone="M20", phase="phase3", method="CellChat",
  timestamp=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),
  script="scripts/phase3/ccc/run_cellchat.R",
  input=list(path=input_rds, md5=in_md5),
  parameters=list(cellchat_version=as.character(packageVersion("CellChat")),
    database="CellChatDB.human", db_interactions=nrow(db$interaction),
    db_pathways=length(unique(db$interaction$pathway_name)),
    assay="RNA", layer="data", normalization="LogNormalize (M19)",
    raw.use=TRUE, population.size=TRUE, min.cells=min_cells, seed=seed,
    population_size_rationale="Populations range from 106 to 5,064 cells, so population.size=TRUE is used to correct for unequal abundance."),
  comparability_caveat="CellChat communication probability is not numerically comparable to a LIANA rank or a CellPhoneDB p-value. It contributes only a boolean support flag to the concordance model.",
  results=list(samples_run=names(res_all), n_lr_rows=nrow(cc_df),
    per_sample=as.list(setNames(per_sample$Freq, per_sample$sample_id)),
    tumour_involving=sum(cc_df$source=="MPNST-Tumor" | cc_df$target=="MPNST-Tumor")),
  warnings=warns, slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM"),
  git_commit=tryCatch(trimws(system("git rev-parse HEAD",intern=TRUE)), error=function(e) NA_character_),
  elapsed_seconds=as.numeric(difftime(Sys.time(),t0,units="secs")), session_info=capture.output(sessionInfo()))
write_json(rec, file.path(out_dir, sprintf("m20_cellchat_record%s.json", sfx)), auto_unbox=TRUE, pretty=TRUE, null="null", digits=NA)
prov <- record_provenance("phase3_m20_cellchat", inputs=list(ccc_input=input_rds),
  outputs=list(all=file.path(out_dir, sprintf("all_CellChat_interactions%s.tsv", sfx)),
    record=file.path(out_dir, sprintf("m20_cellchat_record%s.json", sfx))),
  parameters=list(milestone="M20", method="CellChat",
    cellchat_version=as.character(packageVersion("CellChat")), seed=seed,
    slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM")), dataset="combined")
save_provenance_json(prov, file.path(out_dir, sprintf("prov_m20_cellchat%s.json", sfx)))
log_info(sprintf("CellChat COMPLETE in %.1f min. Warnings: %d", as.numeric(difftime(Sys.time(),t0,units="mins")), length(warns)), stage=STAGE)
log_system_usage(stage=STAGE)
