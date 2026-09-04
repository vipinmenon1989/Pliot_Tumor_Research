# scripts/phase3/lochness/lochness_mpnst.R
#
# Phase 3 / Milestone M24 - LochNESS receiver-state analysis for MPNST.
#
# Implements the OFFICIAL MMCA formulation (shendurelab/MMCA @ af629c49,
# Section_5_step_1_run_lochness_calculation.R + demo_lochness.R):
#     lochNESS_i = ( n_target_neighbours_i / k ) / global_fraction_target_in_reference - 1
#   * k = round(0.5 * sqrt(N))                       [MMCA rule, not the Python constant]
#   * L2-normalised PCA-space kNN via Seurat::L2Dim  [MMCA choice]
#   * exact FNN::get.knnx(algorithm="kd_tree")       [MMCA choice]
#   * query and reference sets DISJOINT by sample    [MMCA same-embryo exclusion]
#   * global fraction computed over the REFERENCE set only
#
# Core formula cross-checked against, but NOT copied from:
#   perturbseq-pipeline/src/perturbseq_pipeline/lochness.py (read-only reference)
# Corrections applied relative to that Python port: same-sample exclusion ADDED,
# k switched to round(0.5*sqrt(N)), L2 normalisation added.
#
# LochNESS IS NOT A LIGAND-RECEPTOR METHOD. No output of this script may be read as
# evidence of cell-cell communication. See reports/phase3/LOCHNESS_MPNST_DESIGN.md.

options(stringsAsFactors=FALSE); options(future.globals.maxSize=+Inf)
suppressPackageStartupMessages({library(Seurat); library(SeuratObject); library(Matrix)
  library(FNN); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(viridis); library(RColorBrewer); library(jsonlite); library(digest)})
source("scripts/R/logging_utils.R"); source("scripts/R/provenance_utils.R"); setup_strict_logging()
STAGE <- "phase3_m24_lochness"

args <- commandArgs(trailingOnly=TRUE)
input_rds   <- "results/phase3/ccc/ccc_input_object.rds"
context_tsv <- "results/phase3/lochness/lochness_context_labels.tsv"
out_dir     <- "results/phase3/lochness"
fig_dir     <- "results/phase3/figures/M24"
lineages    <- c("Macrophage","Fibroblast","Endothelial","CD8-T")
n_perm_B    <- 100L
min_per_sample <- 10L
seed <- 42L
i <- 1
while (i <= length(args)) { a <- args[i]
  if (a=="--input"){input_rds<-args[i+1];i<-i+2} else if (a=="--context"){context_tsv<-args[i+1];i<-i+2}
  else if (a=="--out-dir"){out_dir<-args[i+1];i<-i+2} else if (a=="--fig-dir"){fig_dir<-args[i+1];i<-i+2}
  else if (a=="--lineages"){lineages<-unlist(strsplit(args[i+1],","));i<-i+2}
  else if (a=="--n-perm"){n_perm_B<-as.integer(args[i+1]);i<-i+2}
  else if (a=="--random-seed"){seed<-as.integer(args[i+1]);i<-i+2}
  else stop(sprintf("Unknown argument: %s", a)) }

t0 <- Sys.time(); warns <- character(0)
wrn <- function(m){warns<<-c(warns,m); log_warn(m,stage=STAGE)}
fail <- function(m){log_error(m,stage=STAGE); stop(m,call.=FALSE)}
for (d in c(out_dir,fig_dir)) dir.create(d,recursive=TRUE,showWarnings=FALSE)
set.seed(seed)
derive_seed <- function(base, id) as.integer(abs(base + sum(utf8ToInt(as.character(id)))) %% 2147483647L)

log_info("================ M24 LochNESS (MPNST adaptation) ================", stage=STAGE)
log_info("Formulation: MMCA official (k = round(0.5*sqrt(N)); L2 PCA; exact kNN; same-sample exclusion).", stage=STAGE)
log_info("LochNESS is NOT a ligand-receptor method. Outputs describe receiver-state neighbourhood composition only.", stage=STAGE)

in_md5 <- digest(input_rds, file=TRUE, algo="md5")
obj <- readRDS(input_rds)
log_info(sprintf("Input %s (md5 %s): %d cells", input_rds, in_md5, ncol(obj)), stage=STAGE)
if (!("postint_harmony" %in% Reductions(obj))) fail("postint_harmony reduction absent.")

ctx <- read.delim(context_tsv, sep="\t")
if (!all(c("sample_id","context") %in% colnames(ctx))) fail("Context table needs sample_id and context columns.")
log_info(sprintf("Context labels (sample-level, derived from TUMOUR cells): %s",
  paste(sprintf("%s=%s", ctx$sample_id, ctx$context), collapse=", ")), stage=STAGE)
if (length(unique(ctx$context)) != 2) fail("Context label must have exactly two levels.")
HIGH <- "high"
ctx_map <- setNames(as.character(ctx$context), as.character(ctx$sample_id))

# ---------------------------------------------------------------------------
# Core LochNESS, following MMCA exactly
# ---------------------------------------------------------------------------
lochness_core <- function(emb, samples, ctx_vec, kadj) {
  # emb: cells x dims (already L2-normalised). samples: per-cell sample id.
  # ctx_vec: per-cell context label. Query = one sample; reference = all other samples.
  out <- rep(NA_real_, nrow(emb)); names(out) <- rownames(emb)
  for (j in unique(samples)) {
    qi <- which(samples == j); ri <- which(samples != j)
    if (!length(qi) || length(ri) <= kadj) next
    knn <- FNN::get.knnx(emb[ri, , drop=FALSE], emb[qi, , drop=FALSE], k=kadj, algorithm="kd_tree")
    nb_ctx <- matrix(ctx_vec[ri][knn$nn.index], nrow=length(qi))
    n_target <- rowSums(nb_ctx == HIGH)
    global_frac <- mean(ctx_vec[ri] == HIGH)     # over the REFERENCE set only, as in MMCA
    if (global_frac <= 0 || global_frac >= 1) next
    out[qi] <- (n_target / kadj) / global_frac - 1
  }
  out
}

all_scores <- list(); summaries <- list(); perm_summ <- list(); refvar <- list(); skipped <- list()

for (LIN in lineages) {
  log_info(sprintf("=== lineage %s ===", LIN), stage=STAGE)
  cells <- colnames(obj)[as.character(obj$ccc_label) == LIN]
  if (!length(cells)) { skipped[[LIN]] <- "lineage absent"; wrn(sprintf("%s absent; skipped", LIN)); next }
  sub <- subset(obj, cells=cells)
  smp <- as.character(sub$ccc_sample)
  per <- table(smp)
  if (any(per < min_per_sample) || length(per) < 2) {
    skipped[[LIN]] <- sprintf("insufficient per-sample cells: %s", paste(sprintf("%s=%d", names(per), per), collapse=";"))
    wrn(sprintf("%s skipped: %s", LIN, skipped[[LIN]])); next }
  N <- ncol(sub); kadj <- as.integer(round(0.5*sqrt(N)))
  min_ref <- min(N - as.integer(per))
  if (kadj >= min_ref) {
    skipped[[LIN]] <- sprintf("k=%d >= smallest reference set %d", kadj, min_ref)
    wrn(sprintf("%s skipped: %s", LIN, skipped[[LIN]])); next }
  unstable <- kadj < 10L
  log_info(sprintf("  N=%d  k=round(0.5*sqrt(N))=%d  per-sample: %s%s", N, kadj,
    paste(sprintf("%s=%d", names(per), per), collapse=", "),
    if (unstable) "  [PRE-FLAGGED UNSTABLE: k < 10]" else ""), stage=STAGE)
  if (unstable) wrn(sprintf("%s: k=%d < 10, estimates are unstable and are reported with that caveat.", LIN, kadj))

  # primary representation: L2-normalised Harmony embedding restricted to the lineage
  sub <- Seurat::L2Dim(sub, reduction="postint_harmony", new.dr="hL2", new.key="HL2_")
  emb <- Embeddings(sub, "hL2")
  ctx_vec <- ctx_map[smp]
  set.seed(seed)
  sc <- lochness_core(emb, smp, ctx_vec, kadj)

  # secondary representation: within-lineage PCA (sensitivity to the embedding choice)
  sub2 <- NormalizeData(sub, assay="RNA", verbose=FALSE)
  sub2 <- FindVariableFeatures(sub2, nfeatures=2000, verbose=FALSE)
  sub2 <- ScaleData(sub2, verbose=FALSE)
  npc <- min(30L, ncol(sub2)-1L)
  sub2 <- RunPCA(sub2, npcs=npc, verbose=FALSE, seed.use=seed)
  sub2 <- Seurat::L2Dim(sub2, reduction="pca", new.dr="pL2", new.key="PL2_")
  set.seed(seed)
  sc_pca <- lochness_core(Embeddings(sub2, "pL2"), smp, ctx_vec, kadj)
  rho_rep <- suppressWarnings(cor(sc, sc_pca, use="complete.obs", method="spearman"))
  log_info(sprintf("  representation sensitivity (Harmony vs within-lineage PCA): Spearman rho = %.3f", rho_rep), stage=STAGE)

  uxy <- Embeddings(sub, "postint_umap_harmony")
  df <- data.frame(cell=colnames(sub), lineage=LIN, sample_id=smp, context=ctx_vec,
    lochness=as.numeric(sc), lochness_within_lineage_pca=as.numeric(sc_pca),
    UMAP_1=uxy[,1], UMAP_2=uxy[,2], k=kadj, n_lineage_cells=N,
    unstable_k=unstable, stringsAsFactors=FALSE)
  all_scores[[LIN]] <- df

  obs_stat <- mean(sc, na.rm=TRUE)
  obs_by_ctx <- tapply(sc, ctx_vec, mean, na.rm=TRUE)
  summaries[[LIN]] <- data.frame(lineage=LIN, n_cells=N, k=kadj, unstable_k=unstable,
    mean_lochness=obs_stat, sd_lochness=sd(sc, na.rm=TRUE),
    mean_in_high_context=unname(obs_by_ctx[HIGH]),
    mean_in_low_context=unname(obs_by_ctx[setdiff(names(obs_by_ctx), HIGH)][1]),
    frac_positive=mean(sc > 0, na.rm=TRUE),
    representation_rho_harmony_vs_pca=rho_rep, stringsAsFactors=FALSE)

  # ---- Null A: sample-level EXACT permutation of the context label ----
  smps <- sort(unique(smp)); n_high <- sum(unique(ctx_map[smps]) == HIGH)
  n_high <- sum(ctx_map[smps] == HIGH)
  combos <- combn(smps, n_high, simplify=FALSE)
  nullA <- vapply(combos, function(hi) {
    cv <- ifelse(smp %in% hi, HIGH, "low")
    set.seed(derive_seed(seed, paste(hi, collapse="_")))
    mean(lochness_core(emb, smp, cv, kadj), na.rm=TRUE) }, numeric(1))
  labA <- vapply(combos, function(h) paste(h, collapse="+"), character(1))
  obs_idx <- which(vapply(combos, function(h) setequal(h, smps[ctx_map[smps]==HIGH]), logical(1)))
  rankA <- sum(abs(nullA) >= abs(obs_stat))
  pA <- rankA / length(nullA)
  log_info(sprintf("  Null A (exact, %d labelings): observed mean %.4f | rank %d/%d | descriptive p = %.3f (floor 1/%d = %.3f)",
    length(nullA), obs_stat, rankA, length(nullA), pA, length(nullA), 1/length(nullA)), stage=STAGE)

  # ---- Null B: cell-level context shuffle preserving per-sample sizes ----
  nullB <- numeric(n_perm_B)
  for (b in seq_len(n_perm_B)) {
    set.seed(derive_seed(seed + 1000L*b, LIN))
    perm_map <- setNames(sample(ctx_map[smps]), smps)
    nullB[b] <- mean(lochness_core(emb, smp, perm_map[smp], kadj), na.rm=TRUE)
  }
  pB <- mean(abs(nullB) >= abs(obs_stat))
  log_info(sprintf("  Null B (%d shuffles): mean %.4f sd %.4f | descriptive p = %.3f",
    n_perm_B, mean(nullB), sd(nullB), pB), stage=STAGE)
  perm_summ[[LIN]] <- data.frame(lineage=LIN, observed_mean_lochness=obs_stat,
    nullA_n_labelings=length(nullA), nullA_rank=rankA, nullA_p_descriptive=pA,
    nullA_p_floor=1/length(nullA), nullA_labelings=paste(labA, collapse=";"),
    nullA_values=paste(round(nullA,4), collapse=";"),
    nullB_n=n_perm_B, nullB_mean=mean(nullB), nullB_sd=sd(nullB), nullB_p_descriptive=pB,
    interpretation="Descriptive ranks only. With four samples no sample-level permutation test can reach conventional significance.",
    stringsAsFactors=FALSE)

  # ---- per-reference-sample variant (MMCA 'alternative implementation') ----
  rv <- list()
  for (hold in smps) {
    keep <- smp != hold
    if (sum(keep) <= kadj+1) next
    k2 <- as.integer(round(0.5*sqrt(sum(keep))))
    if (k2 < 3) next
    cv <- ctx_map[smp[keep]]
    if (length(unique(cv)) < 2) next
    set.seed(derive_seed(seed, hold))
    s2 <- lochness_core(emb[keep,,drop=FALSE], smp[keep], cv, k2)
    v <- rep(NA_real_, length(smp)); v[keep] <- s2; rv[[hold]] <- v
  }
  if (length(rv) >= 2) {
    M <- do.call(cbind, rv)
    cm <- suppressWarnings(cor(M, use="pairwise.complete.obs", method="spearman"))
    mean_off <- mean(cm[upper.tri(cm)], na.rm=TRUE)
    log_info(sprintf("  leave-one-sample-out: %d variants | mean pairwise Spearman rho = %.3f",
      ncol(M), mean_off), stage=STAGE)
    refvar[[LIN]] <- data.frame(lineage=LIN, n_variants=ncol(M),
      mean_pairwise_spearman=mean_off,
      min_pairwise_spearman=min(cm[upper.tri(cm)], na.rm=TRUE),
      held_out_samples=paste(colnames(M), collapse=";"), stringsAsFactors=FALSE)
  }
}

if (!length(all_scores)) fail("No lineage was analysable for LochNESS.")
scores <- dplyr::bind_rows(all_scores)
write.table(scores, file.path(out_dir,"lochness_scores.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
summ <- dplyr::bind_rows(summaries)
write.table(summ, file.path(out_dir,"lochness_summary.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
ps <- dplyr::bind_rows(perm_summ)
write.table(ps, file.path(out_dir,"lochness_permutation_summary.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
if (length(refvar)) write.table(dplyr::bind_rows(refvar),
  file.path(out_dir,"lochness_leave_one_sample_out.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
if (length(skipped)) write.table(data.frame(lineage=names(skipped), reason=unlist(skipped)),
  file.path(out_dir,"lochness_skipped_lineages.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
log_info("LochNESS summary:", stage=STAGE)
for (r in seq_len(nrow(summ)))
  log_info(sprintf("  %-14s n=%5d k=%3d  mean=%+.4f  high=%+.4f low=%+.4f  frac>0=%.3f  rho(rep)=%.3f%s",
    summ$lineage[r], summ$n_cells[r], summ$k[r], summ$mean_lochness[r],
    summ$mean_in_high_context[r], summ$mean_in_low_context[r], summ$frac_positive[r],
    summ$representation_rho_harmony_vs_pca[r], if (summ$unstable_k[r]) "  [UNSTABLE k]" else ""), stage=STAGE)

# ---------------------------------------------------------------------------
# FIGURES
# ---------------------------------------------------------------------------
th <- theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),
      plot.title=element_text(face="bold"), plot.subtitle=element_text(size=9,colour="grey30"))
NOTE <- "LochNESS describes neighbourhood composition in receiver-state space. It is NOT evidence of cell-cell communication."
frows <- list()
sf <- function(p,name,w,h,analysis,params){ for (ext in c("pdf","png")) {
  fp <- file.path(fig_dir, paste0(name,".",ext))
  ggsave(fp,p,width=w,height=h,units="in",dpi=200,device=ext,limitsize=FALSE)
  frows[[length(frows)+1]] <<- data.frame(figure_path=fp, phase="phase3", milestone="M24",
    analysis=analysis, sender="MPNST-Tumor(context)", receiver="receiver lineages",
    method="LochNESS (MMCA formulation)", input=input_rds, input_checksum=in_md5,
    script="scripts/phase3/lochness/lochness_mpnst.R", parameters=params, stringsAsFactors=FALSE) }
  log_info(sprintf("  wrote %s", name), stage=STAGE) }

lim <- max(abs(quantile(scores$lochness, c(0.01,0.99), na.rm=TRUE)))
sf(ggplot(scores, aes(UMAP_1, UMAP_2, colour=lochness))+geom_point(size=0.5, alpha=0.85, stroke=0)+
   scale_colour_gradient2(low="#2166AC", mid="white", high="#B2182B", midpoint=0,
     limits=c(-lim,lim), oob=scales::squish, name="lochNESS")+
   facet_wrap(~lineage, scales="free")+
   labs(title="LochNESS on the post-Harmony UMAP, per receiver lineage",
        subtitle=paste("Positive = neighbourhood over-represented for high-context samples.", NOTE),
        x="UMAP 1", y="UMAP 2")+th,
   "M24_01_lochness_umap", 12, 9, "lochness_umap", "k=round(0.5*sqrt(N));L2 Harmony;same-sample excluded")
sf(ggplot(scores, aes(lochness, fill=context))+geom_density(alpha=0.5)+
   geom_vline(xintercept=0, linetype="dashed", colour="grey30")+
   facet_wrap(~lineage, scales="free")+scale_fill_brewer(palette="Set1", name="sample context")+
   labs(title="LochNESS distribution by tumour-derived signalling context",
        subtitle=paste("Context is a SAMPLE-level label derived from tumour cells; the score is computed in receiver space.", NOTE),
        x="lochNESS", y="density")+th,
   "M24_02_lochness_distribution_by_context", 11, 8, "lochness_distribution", "context = sample-level high/low")
sf(ggplot(scores, aes(sample_id, lochness, fill=context))+
   geom_violin(scale="width", alpha=0.7)+geom_hline(yintercept=0, linetype="dashed")+
   facet_wrap(~lineage, scales="free_y")+scale_fill_brewer(palette="Set1", name="context")+
   labs(title="LochNESS by sample and lineage (cross-sample consistency)",
        subtitle=paste("A score driven by one patient shows up here as a single deviant violin.", NOTE),
        x="sample (= patient)", y="lochNESS")+th+
   theme(axis.text.x=element_text(angle=45,hjust=1)),
   "M24_03_lochness_cross_sample_consistency", 11, 8, "consistency", "per sample x lineage")
if (nrow(ps)) {
  nd <- do.call(rbind, lapply(seq_len(nrow(ps)), function(r) {
    v <- as.numeric(unlist(strsplit(ps$nullA_values[r], ";")))
    data.frame(lineage=ps$lineage[r], value=v, observed=ps$observed_mean_lochness[r], stringsAsFactors=FALSE) }))
  sf(ggplot(nd, aes(value))+geom_dotplot(binwidth=diff(range(nd$value, na.rm=TRUE))/40, fill="grey70")+
     geom_vline(aes(xintercept=observed), colour="#B2182B", linewidth=1)+
     facet_wrap(~lineage, scales="free")+
     labs(title="Null A: exact sample-level permutation distribution",
          subtitle="Grey dots = all possible context labelings of the four samples; red line = observed. With 4 samples the p-value floor is 1/6 = 0.17, so this is a descriptive rank, not a significance test.",
          x="mean lochNESS under the labeling", y=NULL)+th+
     theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()),
     "M24_04_lochness_nullA_exact_permutation", 11, 7, "permutation", "exact sample-level null")
}
sf(ggplot(scores, aes(lochness, lochness_within_lineage_pca))+
   geom_point(size=0.4, alpha=0.4, colour="#3B6EA5")+geom_abline(slope=1, intercept=0, linetype="dashed")+
   facet_wrap(~lineage, scales="free")+
   labs(title="Representation sensitivity: Harmony-space versus within-lineage PCA LochNESS",
        subtitle=paste("Agreement means the score does not hinge on the embedding choice.", NOTE),
        x="lochNESS (L2 Harmony space, primary)", y="lochNESS (within-lineage PCA, sensitivity)")+th,
   "M24_05_lochness_representation_sensitivity", 11, 8, "sensitivity", "Harmony vs within-lineage PCA")
write.table(do.call(rbind, frows), file.path(out_dir,"figure_index_m24.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

rec <- list(milestone="M24", phase="phase3", timestamp=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),
  script="scripts/phase3/lochness/lochness_mpnst.R",
  formulation="Official MMCA (shendurelab/MMCA @ af629c498f421d6f0bfb325e9fd6b2dab0e22837)",
  parameters=list(k_rule="round(0.5*sqrt(N))", representation="Seurat::L2Dim on postint_harmony (30 dims), restricted to lineage",
    secondary_representation="within-lineage PCA, L2-normalised (sensitivity)",
    knn="FNN::get.knnx(algorithm='kd_tree'), exact",
    same_sample_exclusion=TRUE,
    global_fraction_scope="reference set only (cells eligible to be neighbours)",
    lineages_requested=lineages, min_cells_per_sample=min_per_sample,
    nullA="exact sample-level permutation of the context label (all labelings preserving group sizes)",
    nullB=sprintf("cell-level context shuffle preserving per-sample sizes, %d iterations", n_perm_B),
    seed=seed),
  context_labels=as.list(setNames(as.character(ctx$context), as.character(ctx$sample_id))),
  scientific_guard="LochNESS is NOT a ligand-receptor method. No output may be read as evidence of cell-cell communication.",
  hard_limitation="With four samples the exact sample-level null has at most 6 labelings, so the minimum attainable one-sided p-value is about 0.17. All permutation results are descriptive ranks, never significance tests.",
  results=list(lineages_analysed=summ$lineage, skipped=skipped,
    summary=lapply(seq_len(nrow(summ)), function(r) as.list(summ[r,]))),
  warnings=warns, slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM"),
  git_commit=tryCatch(trimws(system("git rev-parse HEAD",intern=TRUE)), error=function(e) NA_character_),
  elapsed_seconds=as.numeric(difftime(Sys.time(),t0,units="secs")), session_info=capture.output(sessionInfo()))
write_json(rec, file.path(out_dir,"m24_lochness_record.json"), auto_unbox=TRUE, pretty=TRUE, null="null", digits=NA)
prov <- record_provenance("phase3_m24_lochness", inputs=list(ccc_input=input_rds, context=context_tsv),
  outputs=list(scores=file.path(out_dir,"lochness_scores.tsv"),
    summary=file.path(out_dir,"lochness_summary.tsv"),
    permutation=file.path(out_dir,"lochness_permutation_summary.tsv"),
    record=file.path(out_dir,"m24_lochness_record.json")),
  parameters=list(milestone="M24", k_rule="round(0.5*sqrt(N))", same_sample_exclusion=TRUE, seed=seed,
    slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM")), dataset="combined")
save_provenance_json(prov, file.path(out_dir,"prov_m24_lochness.json"))
log_info(sprintf("M24 COMPLETE in %.1f min. Warnings: %d", as.numeric(difftime(Sys.time(),t0,units="mins")), length(warns)), stage=STAGE)
log_system_usage(stage=STAGE)
