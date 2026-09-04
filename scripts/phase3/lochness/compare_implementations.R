# scripts/phase3/lochness/compare_implementations.R
#
# Phase 3 / M24 - LochNESS implementation comparison (spec SS36).
# Runs the official MMCA formulation (R, FNN) and the perturb-seq formulation (Python,
# sklearn) on the SAME deterministic subset and reports agreement.

options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(Seurat); library(SeuratObject); library(FNN)
  library(dplyr); library(ggplot2); library(patchwork); library(jsonlite); library(digest)})
source("scripts/R/logging_utils.R"); source("scripts/R/provenance_utils.R"); setup_strict_logging()
STAGE <- "phase3_m24_impl_comparison"
args <- commandArgs(trailingOnly=TRUE)
obj_rds <- "results/phase3/ccc/ccc_input_object.rds"
ctx_p   <- "results/phase3/lochness/lochness_context_labels.tsv"
out_dir <- "results/phase3/lochness"; fig_dir <- "results/phase3/figures/M24"
lineage <- "Macrophage"; n_sub <- 2000L; seed <- 42L
i <- 1
while (i <= length(args)) { a <- args[i]
  if (a=="--object"){obj_rds<-args[i+1];i<-i+2} else if (a=="--context"){ctx_p<-args[i+1];i<-i+2}
  else if (a=="--out-dir"){out_dir<-args[i+1];i<-i+2} else if (a=="--fig-dir"){fig_dir<-args[i+1];i<-i+2}
  else if (a=="--lineage"){lineage<-args[i+1];i<-i+2} else if (a=="--n-sub"){n_sub<-as.integer(args[i+1]);i<-i+2}
  else stop(sprintf("Unknown argument: %s", a)) }
t0 <- Sys.time(); warns <- character(0)
wrn <- function(m){warns<<-c(warns,m); log_warn(m,stage=STAGE)}
fail <- function(m){log_error(m,stage=STAGE); stop(m,call.=FALSE)}
for (d in c(out_dir,fig_dir)) dir.create(d,recursive=TRUE,showWarnings=FALSE)
set.seed(seed)

log_info("======= M24 LochNESS IMPLEMENTATION COMPARISON =======", stage=STAGE)
obj <- readRDS(obj_rds)
ctx <- read.delim(ctx_p, sep="\t")
ctx_map <- setNames(as.character(ctx$context), as.character(ctx$sample_id))
cells <- colnames(obj)[as.character(obj$ccc_label)==lineage]
if (!length(cells)) fail(sprintf("Lineage %s absent.", lineage))
sub <- subset(obj, cells=cells)
# deterministic stratified subsample (pattern adapted from perturbseq distance._sample_cell_indices)
smp_all <- as.character(sub$ccc_sample)
if (ncol(sub) > n_sub) {
  per <- max(1L, floor(n_sub/length(unique(smp_all))))
  set.seed(seed)
  keep <- unlist(lapply(unique(smp_all), function(s) { cc <- colnames(sub)[smp_all==s]
    if (length(cc) > per) sample(cc, per) else cc }))
  sub <- subset(sub, cells=keep)
}
sub <- Seurat::L2Dim(sub, reduction="postint_harmony", new.dr="hL2", new.key="HL2_")
emb <- Embeddings(sub, "hL2"); smp <- as.character(sub$ccc_sample); ctxv <- ctx_map[smp]
N <- nrow(emb); kadj <- as.integer(round(0.5*sqrt(N)))
log_info(sprintf("Deterministic subset: lineage %s, %d cells (seed %d), k = round(0.5*sqrt(N)) = %d",
  lineage, N, seed, kadj), stage=STAGE)
log_info(sprintf("  per sample: %s", paste(sprintf("%s=%d", names(table(smp)), table(smp)), collapse=", ")), stage=STAGE)

# ---- R, official MMCA formulation ----
loch_R <- function(emb, samples, ctx_vec, k, target="high") {
  out <- rep(NA_real_, nrow(emb))
  for (j in unique(samples)) {
    qi <- which(samples==j); ri <- which(samples!=j)
    if (!length(qi) || length(ri) <= k) next
    knn <- FNN::get.knnx(emb[ri,,drop=FALSE], emb[qi,,drop=FALSE], k=k, algorithm="kd_tree")
    nt <- rowSums(matrix(ctx_vec[ri][knn$nn.index], nrow=length(qi)) == target)
    gf <- mean(ctx_vec[ri]==target); if (gf<=0 || gf>=1) next
    out[qi] <- (nt/k)/gf - 1 }
  out }
sc_R <- loch_R(emb, smp, ctxv, kadj)
log_info(sprintf("R (MMCA): %d scored | mean %+.5f | sd %.5f", sum(is.finite(sc_R)),
  mean(sc_R,na.rm=TRUE), sd(sc_R,na.rm=TRUE)), stage=STAGE)

# ---- export for Python ----
ep <- file.path(out_dir,"impl_cmp_embedding.tsv"); mp <- file.path(out_dir,"impl_cmp_meta.tsv")
write.table(data.frame(cell=rownames(emb), emb, check.names=FALSE), ep, sep="\t", row.names=FALSE, quote=FALSE)
write.table(data.frame(cell=rownames(emb), sample_id=smp, context=ctxv), mp, sep="\t", row.names=FALSE, quote=FALSE)

py <- "scripts/phase3/lochness/lochness_python_port.py"
res <- list()
for (mode in c("mmca","perturbseq")) {
  op <- file.path(out_dir, sprintf("impl_cmp_python_%s.tsv", mode))
  cmd <- sprintf("python %s --embedding %s --meta %s --out %s --k %d --mode %s", py, ep, mp, op, kadj, mode)
  log_info(sprintf("Running Python (%s): %s", mode, cmd), stage=STAGE)
  o <- suppressWarnings(system(cmd, intern=TRUE))
  log_info(sprintf("  python: %s", paste(tail(o,1), collapse=" ")), stage=STAGE)
  if (!file.exists(op)) { wrn(sprintf("Python %s mode produced no output", mode)); next }
  d <- read.delim(op, sep="\t")
  res[[mode]] <- d$lochness_python[match(rownames(emb), d$cell)]
}
if (!length(res)) fail("Python port produced no output in either mode.")

cmp_rows <- list()
for (mode in names(res)) {
  v <- res[[mode]]; ok <- is.finite(sc_R) & is.finite(v)
  pear <- suppressWarnings(cor(sc_R[ok], v[ok]))
  spear <- suppressWarnings(cor(sc_R[ok], v[ok], method="spearman"))
  mad_ <- max(abs(sc_R[ok]-v[ok])); rmse <- sqrt(mean((sc_R[ok]-v[ok])^2))
  cmp_rows[[mode]] <- data.frame(python_mode=mode, n_compared=sum(ok),
    r_mean=mean(sc_R[ok]), py_mean=mean(v[ok]), r_sd=sd(sc_R[ok]), py_sd=sd(v[ok]),
    pearson=pear, spearman=spear, max_abs_diff=mad_, rmse=rmse,
    same_sample_exclusion_python=(mode=="mmca"),
    denominator_python=ifelse(mode=="mmca","requested k","actual neighbour count"),
    global_fraction_scope_python=ifelse(mode=="mmca","reference set","all cells"),
    stringsAsFactors=FALSE)
  log_info(sprintf("  %-11s vs R(MMCA): n=%d | Pearson %.4f | Spearman %.4f | max|diff| %.3e | RMSE %.3e",
    mode, sum(ok), pear, spear, mad_, rmse), stage=STAGE)
}
cmp <- do.call(rbind, cmp_rows)
write.table(cmp, file.path(out_dir,"LOCHNESS_IMPLEMENTATION_COMPARISON.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

df <- data.frame(r=sc_R, mmca=res[["mmca"]] %||% NA_real_, perturbseq=res[["perturbseq"]] %||% NA_real_)
p1 <- ggplot(df, aes(r, mmca))+geom_point(size=0.5, alpha=0.5, colour="#1A9850")+
  geom_abline(slope=1, intercept=0, linetype="dashed")+
  labs(title="R (official MMCA) vs Python, MMCA-equivalent settings",
       subtitle=sprintf("Same k, same L2 embedding, same-sample exclusion, reference-set global fraction. Pearson %.4f",
         cmp$pearson[cmp$python_mode=="mmca"]),
       x="lochNESS (R, FNN)", y="lochNESS (Python, sklearn)")+theme_bw(base_size=11)
p2 <- ggplot(df, aes(r, perturbseq))+geom_point(size=0.5, alpha=0.5, colour="#C1452B")+
  geom_abline(slope=1, intercept=0, linetype="dashed")+
  labs(title="R (official MMCA) vs Python, perturb-seq formulation",
       subtitle=sprintf("No same-sample exclusion, actual-count denominator, whole-dataset global fraction. Pearson %.4f",
         cmp$pearson[cmp$python_mode=="perturbseq"]),
       x="lochNESS (R, MMCA)", y="lochNESS (Python, perturb-seq formulation)")+theme_bw(base_size=11)
pp <- (p1|p2)+plot_annotation(title="LochNESS implementation comparison",
  subtitle="Left: the two implementations agree when given mathematically equivalent inputs. Right: the perturb-seq formulation differs because it omits same-sample exclusion, which matters in MPNST where sample = patient = batch.")
for (ext in c("pdf","png")) ggsave(file.path(fig_dir, paste0("M24_06_implementation_comparison.",ext)),
  pp, width=13, height=6, units="in", dpi=200, device=ext)
log_info("  wrote M24_06_implementation_comparison", stage=STAGE)

rec <- list(milestone="M24", phase="phase3", component="implementation_comparison",
  timestamp=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),
  r_script="scripts/phase3/lochness/compare_implementations.R",
  python_script=py,
  subset=list(lineage=lineage, n_cells=N, k=kadj, seed=seed,
    sampling="deterministic stratified by sample"),
  comparison=lapply(seq_len(nrow(cmp)), function(r) as.list(cmp[r,])),
  conclusion=sprintf("With mathematically equivalent inputs the R and Python implementations agree (Pearson %.4f). The perturb-seq formulation as written diverges (Pearson %.4f) because it omits same-sample exclusion, uses the actual neighbour count as denominator, and takes the global fraction over all cells rather than the reference set. The R implementation is therefore used for Phase 3.",
    cmp$pearson[cmp$python_mode=="mmca"], cmp$pearson[cmp$python_mode=="perturbseq"]),
  warnings=warns, slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM"),
  git_commit=tryCatch(trimws(system("git rev-parse HEAD",intern=TRUE)), error=function(e) NA_character_),
  elapsed_seconds=as.numeric(difftime(Sys.time(),t0,units="secs")))
write_json(rec, file.path(out_dir,"m24_implementation_comparison_record.json"), auto_unbox=TRUE, pretty=TRUE, null="null", digits=NA)
log_info(sprintf("Implementation comparison COMPLETE in %.1f min.", as.numeric(difftime(Sys.time(),t0,units="mins"))), stage=STAGE)
