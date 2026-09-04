# scripts/phase3/ccc/run_liana.R
#
# Phase 3 / Milestone M20 - LIANA execution, sample-aware.
#
# LIANA is the primary multi-resource / multi-score CONSENSUS LR framework. It is run
# independently on each sample and then pooled, so recurrence across patients can be
# assessed rather than assumed.
#
# NOTE ON INDEPENDENCE: LIANA's method set includes a CellPhoneDB-style score. LIANA and the
# standalone CellPhoneDB run are therefore NOT fully independent. This is recorded here and
# carried into the M21 concordance report rather than glossed over.
#
# Expression: RNA assay, joined, LogNormalize (prepared in M19). Harmony coordinates are
# never used as expression.

options(stringsAsFactors=FALSE); options(future.globals.maxSize=+Inf)
suppressPackageStartupMessages({library(Seurat); library(SeuratObject); library(Matrix)
  library(SingleCellExperiment); library(liana); library(dplyr); library(tibble)
  library(jsonlite); library(digest)})
source("scripts/R/logging_utils.R"); source("scripts/R/provenance_utils.R"); setup_strict_logging()
STAGE <- "phase3_m20_liana"

args <- commandArgs(trailingOnly=TRUE)
input_rds <- "results/phase3/ccc/ccc_input_object.rds"
out_dir   <- "results/phase3/ccc/by_method"
min_cells <- 10L; seed <- 42L
resource  <- "Consensus"
i <- 1
while (i <= length(args)) { a <- args[i]
  if (a=="--input"){input_rds<-args[i+1];i<-i+2} else if (a=="--out-dir"){out_dir<-args[i+1];i<-i+2}
  else if (a=="--min-cells"){min_cells<-as.integer(args[i+1]);i<-i+2}
  else if (a=="--resource"){resource<-args[i+1];i<-i+2}
  else if (a=="--random-seed"){seed<-as.integer(args[i+1]);i<-i+2}
  else stop(sprintf("Unknown argument: %s", a)) }

t0 <- Sys.time(); warns <- character(0)
wrn <- function(m){warns<<-c(warns,m); log_warn(m,stage=STAGE)}
fail <- function(m){log_error(m,stage=STAGE); stop(m,call.=FALSE)}
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
set.seed(seed)

log_info("================== M20 LIANA (sample-aware) ==================", stage=STAGE)
log_info(sprintf("liana %s | OmnipathR %s | resource '%s' | seed %d",
  packageVersion("liana"), packageVersion("OmnipathR"), resource, seed), stage=STAGE)
in_md5 <- digest(input_rds, file=TRUE, algo="md5")
obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"
log_info(sprintf("Input %s (md5 %s): %d features x %d cells", input_rds, in_md5, nrow(obj), ncol(obj)), stage=STAGE)
if (!all(c("ccc_label","ccc_sample") %in% colnames(obj@meta.data))) fail("M19 columns ccc_label/ccc_sample missing.")
lv <- sort(unique(as.character(obj$ccc_sample)))
pops <- levels(obj$ccc_label)
log_info(sprintf("Samples: %s | populations: %d", paste(lv, collapse=", "), length(pops)), stage=STAGE)

liana_methods <- c("natmi","connectome","logfc","sca","cellphonedb")
log_info(sprintf("LIANA methods: %s", paste(liana_methods, collapse=", ")), stage=STAGE)
wrn("LIANA's method set includes a CellPhoneDB-style score, so LIANA and the standalone CellPhoneDB run are not fully independent. Recorded for the M21 concordance report.")

res_all <- list()
for (s in lv) {
  log_info(sprintf("--- sample %s ---", s), stage=STAGE)
  cells <- colnames(obj)[as.character(obj$ccc_sample) == s]
  sub <- subset(obj, cells=cells)
  keep <- names(which(table(as.character(sub$ccc_label)) >= min_cells))
  sub <- subset(sub, cells=colnames(sub)[as.character(sub$ccc_label) %in% keep])
  sub$ccc_label <- factor(as.character(sub$ccc_label), levels=sort(keep))
  Idents(sub) <- sub$ccc_label
  log_info(sprintf("  %d cells, %d populations (>= %d cells): %s",
    ncol(sub), length(keep), min_cells, paste(sort(keep), collapse=", ")), stage=STAGE)
  # Route through SingleCellExperiment. liana 0.1.14's Seurat path calls
  # GetAssayData(slot=), which is DEFUNCT in SeuratObject 5.3.0. Converting to SCE uses
  # liana's SCE code path instead and avoids downgrading SeuratObject, which would
  # destabilise the frozen Phase 2 software stack.
  sce <- SingleCellExperiment(
    assays = list(counts   = GetAssayData(sub, assay="RNA", layer="counts"),
                  logcounts = GetAssayData(sub, assay="RNA", layer="data")),
    colData = DataFrame(ccc_label = factor(as.character(sub$ccc_label), levels=sort(keep)))
  )
  set.seed(seed)
  tw <- Sys.time()
  lr <- tryCatch(
    liana_wrap(sce, idents_col="ccc_label",
               method=liana_methods, resource=resource, verbose=FALSE),
    error=function(e) { wrn(sprintf("liana_wrap failed for %s: %s", s, conditionMessage(e))); NULL })
  if (is.null(lr)) next
  agg <- tryCatch(liana_aggregate(lr, verbose=FALSE),
                  error=function(e) { wrn(sprintf("liana_aggregate failed for %s: %s", s, conditionMessage(e))); NULL })
  if (is.null(agg)) next
  agg <- as.data.frame(agg); agg$sample_id <- s
  # liana_aggregate reports complexes; expose canonical single-symbol columns too
  if (!("ligand" %in% colnames(agg)) && "ligand.complex" %in% colnames(agg)) agg$ligand <- agg$ligand.complex
  if (!("receptor" %in% colnames(agg)) && "receptor.complex" %in% colnames(agg)) agg$receptor <- agg$receptor.complex
  log_info(sprintf("  %d aggregated interactions in %.1f s", nrow(agg),
    as.numeric(difftime(Sys.time(), tw, units="secs"))), stage=STAGE)
  res_all[[s]] <- agg
  saveRDS(lr, file.path(out_dir, sprintf("liana_raw_%s.rds", s)))
}
if (!length(res_all)) fail("LIANA produced no results for any sample.")
liana_df <- dplyr::bind_rows(res_all)
log_info(sprintf("Pooled LIANA rows across %d samples: %d", length(res_all), nrow(liana_df)), stage=STAGE)
log_info(sprintf("Columns: %s", paste(colnames(liana_df), collapse=", ")), stage=STAGE)

write.table(liana_df, file.path(out_dir,"all_LIANA_interactions.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# LIANA support definition: aggregate_rank is an empirical p-value-like consensus statistic
# (lower = more consistently top-ranked across methods). The 0.05 cut is LIANA's own
# documented convention and is stated rather than tuned.
rank_col <- if ("aggregate_rank" %in% colnames(liana_df)) "aggregate_rank" else
            if ("mean_rank" %in% colnames(liana_df)) "mean_rank" else NA_character_
if (is.na(rank_col)) fail("No LIANA aggregate rank column found.")
LIANA_CUT <- 0.05
liana_df$liana_supported <- liana_df[[rank_col]] <= LIANA_CUT
log_info(sprintf("LIANA support: %s <= %.3f -> %d / %d rows supported", rank_col, LIANA_CUT,
  sum(liana_df$liana_supported), nrow(liana_df)), stage=STAGE)

sig <- liana_df[liana_df$liana_supported, ]
write.table(sig, file.path(out_dir,"LIANA_supported_interactions.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
per_sample <- as.data.frame(table(sample_id=sig$sample_id))
log_info(sprintf("Supported interactions per sample: %s",
  paste(sprintf("%s=%d", per_sample$sample_id, per_sample$Freq), collapse=", ")), stage=STAGE)
tum <- sig[sig$source=="MPNST-Tumor" | sig$target=="MPNST-Tumor", ]
log_info(sprintf("Tumour-involving supported interactions: %d (tumour->TME %d, TME->tumour %d)",
  nrow(tum), sum(sig$source=="MPNST-Tumor"), sum(sig$target=="MPNST-Tumor")), stage=STAGE)

rec <- list(milestone="M20", phase="phase3", method="LIANA",
  timestamp=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),
  script="scripts/phase3/ccc/run_liana.R",
  input=list(path=input_rds, md5=in_md5, cells=ncol(obj)),
  parameters=list(liana_version=as.character(packageVersion("liana")),
    omnipathr_version=as.character(packageVersion("OmnipathR")),
    resource=resource, methods=liana_methods, assay="RNA", layer="data",
    normalization="LogNormalize (prepared in M19)", organism="human",
    input_route="SingleCellExperiment (liana 0.1.14's Seurat path uses the defunct GetAssayData(slot=); SeuratObject was NOT downgraded)",
    min_cells_per_population_per_sample=min_cells, seed=seed,
    aggregation="liana_aggregate (robust rank aggregation)",
    support_rule=sprintf("%s <= %.3f", rank_col, LIANA_CUT),
    rank_column=rank_col),
  independence_caveat="LIANA's method set includes a CellPhoneDB-style score, so LIANA and the standalone CellPhoneDB run are not fully independent.",
  results=list(samples_run=names(res_all), n_rows_pooled=nrow(liana_df),
    n_supported=sum(liana_df$liana_supported),
    supported_per_sample=as.list(setNames(per_sample$Freq, per_sample$sample_id)),
    tumour_involving_supported=nrow(tum)),
  warnings=warns, slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM"),
  git_commit=tryCatch(trimws(system("git rev-parse HEAD",intern=TRUE)), error=function(e) NA_character_),
  elapsed_seconds=as.numeric(difftime(Sys.time(),t0,units="secs")), session_info=capture.output(sessionInfo()))
write_json(rec, file.path(out_dir,"m20_liana_record.json"), auto_unbox=TRUE, pretty=TRUE, null="null", digits=NA)
prov <- record_provenance("phase3_m20_liana", inputs=list(ccc_input=input_rds),
  outputs=list(all=file.path(out_dir,"all_LIANA_interactions.tsv"),
    supported=file.path(out_dir,"LIANA_supported_interactions.tsv"),
    record=file.path(out_dir,"m20_liana_record.json")),
  parameters=list(milestone="M20", method="LIANA", resource=resource,
    liana_version=as.character(packageVersion("liana")), seed=seed,
    slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM")), dataset="combined")
save_provenance_json(prov, file.path(out_dir,"prov_m20_liana.json"))
log_info(sprintf("LIANA COMPLETE in %.1f min. Warnings: %d", as.numeric(difftime(Sys.time(),t0,units="mins")), length(warns)), stage=STAGE)
log_system_usage(stage=STAGE)
