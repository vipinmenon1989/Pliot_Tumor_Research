# scripts/phase3/receiver_response/run_nichenet.R
#
# Phase 3 / Milestone M23 - NicheNet receiver-response analysis.
#
# NicheNet is NOT used here as a fourth ligand-receptor list generator. It answers a
# different question: does a tumour-derived candidate ligand explain the RECEIVER's
# transcriptional programme?
#
#     MPNST-Tumor  ->  candidate ligand  ->  receiver receptor  ->  receiver target genes
#
# Gene set of interest = the receiver population's own marker genes (receiver vs all other
# CCC populations, RNA LogNormalize, Wilcoxon). This is cluster characterisation, exactly as
# in Phase 2 M14 - NOT condition-level differential expression.
# Background = genes expressed in the receiver.
#
# NicheNet ligand activity (AUPR) is not comparable to a LIANA rank, a CellChat probability
# or a CellPhoneDB p-value. It enters the evidence matrix as ORTHOGONAL receiver-response
# support, never as a fourth LR vote.

options(stringsAsFactors=FALSE); options(future.globals.maxSize=+Inf)
suppressPackageStartupMessages({library(Seurat); library(SeuratObject); library(Matrix)
  library(nichenetr); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(viridis); library(RColorBrewer); library(jsonlite); library(digest)})
source("scripts/R/logging_utils.R"); source("scripts/R/provenance_utils.R"); setup_strict_logging()
STAGE <- "phase3_m23_nichenet"

args <- commandArgs(trailingOnly=TRUE)
obj_rds  <- "results/phase3/ccc/ccc_input_object.rds"
master_p <- "results/phase3/ccc/prioritized/MPNST_CCC_MASTER_TABLE.tsv"
nn_dir   <- "external/nichenet"
out_dir  <- "results/phase3/receiver_response"
fig_dir  <- "results/phase3/figures/M23"
receivers <- c("Macrophage","CD8-T","NK","Fibroblast","Endothelial","Monocyte","Dendritic","CD4-T")
sender <- "MPNST-Tumor"; top_n_ligands <- 20L; seed <- 42L
i <- 1
while (i <= length(args)) { a <- args[i]
  if (a=="--object"){obj_rds<-args[i+1];i<-i+2} else if (a=="--master"){master_p<-args[i+1];i<-i+2}
  else if (a=="--nn-dir"){nn_dir<-args[i+1];i<-i+2} else if (a=="--out-dir"){out_dir<-args[i+1];i<-i+2}
  else if (a=="--fig-dir"){fig_dir<-args[i+1];i<-i+2}
  else if (a=="--receivers"){receivers<-unlist(strsplit(args[i+1],","));i<-i+2}
  else if (a=="--random-seed"){seed<-as.integer(args[i+1]);i<-i+2}
  else stop(sprintf("Unknown argument: %s", a)) }
t0 <- Sys.time(); warns <- character(0)
wrn <- function(m){warns<<-c(warns,m); log_warn(m,stage=STAGE)}
fail <- function(m){log_error(m,stage=STAGE); stop(m,call.=FALSE)}
for (d in c(out_dir,fig_dir,nn_dir)) dir.create(d,recursive=TRUE,showWarnings=FALSE)
set.seed(seed)

log_info("============ M23 NicheNet RECEIVER-RESPONSE ============", stage=STAGE)
log_info(sprintf("nichenetr %s | sender %s", packageVersion("nichenetr"), sender), stage=STAGE)

# ---- prior model (Zenodo) ----
urls <- list(
  ligand_target_matrix = "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds",
  lr_network           = "https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds",
  weighted_networks    = "https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds")
mods <- list()
for (nm in names(urls)) {
  fp <- file.path(nn_dir, paste0(nm, ".rds"))
  if (!file.exists(fp)) {
    log_info(sprintf("Downloading NicheNet %s ...", nm), stage=STAGE)
    # R's default 60 s timeout is far too short for a 262 MB file; raise it and fall back to curl.
    old_to <- getOption("timeout"); options(timeout=7200)
    ok <- tryCatch({ download.file(urls[[nm]], fp, mode="wb", quiet=TRUE); TRUE },
                   error=function(e) { wrn(sprintf("download.file failed for %s: %s; trying curl", nm, conditionMessage(e)))
                     rc <- suppressWarnings(system(sprintf("curl -sS -L --retry 3 --max-time 3600 -o %s '%s'", fp, urls[[nm]])))
                     file.exists(fp) && file.info(fp)$size > 1e5 })
    options(timeout=old_to)
    if (!ok || !file.exists(fp)) fail(sprintf("NicheNet prior model %s unavailable; M23 cannot proceed.", nm))
  }
  mods[[nm]] <- readRDS(fp)
  log_info(sprintf("  %-22s %s (md5 %s)", nm,
    paste(dim(mods[[nm]]) %||% length(mods[[nm]]), collapse="x"), digest(fp, file=TRUE, algo="md5")), stage=STAGE)
}
ltm <- mods$ligand_target_matrix; lrn <- mods$lr_network; wn <- mods$weighted_networks
if (is.data.frame(lrn)) lrn <- lrn %>% distinct(from, to)
log_info(sprintf("Prior model: ligand-target %d x %d | lr_network %d rows",
  nrow(ltm), ncol(ltm), nrow(lrn)), stage=STAGE)

obj <- readRDS(obj_rds); DefaultAssay(obj) <- "RNA"
in_md5 <- digest(obj_rds, file=TRUE, algo="md5")
obj$ccc_label <- factor(as.character(obj$ccc_label))
Idents(obj) <- obj$ccc_label
dat <- GetAssayData(obj, assay="RNA", layer="data")
lab <- as.character(obj$ccc_label)
log_info(sprintf("Object: %d genes x %d cells", nrow(dat), ncol(dat)), stage=STAGE)

expressed_in <- function(pop, frac=0.10) {
  m <- dat[, lab==pop, drop=FALSE]; rownames(dat)[Matrix::rowSums(m>0)/ncol(m) >= frac] }
sender_expressed <- expressed_in(sender)
log_info(sprintf("Genes expressed in >=10%% of %s cells: %d", sender, length(sender_expressed)), stage=STAGE)

# prioritised tumour ligands, if the master table exists
prio_lig <- character(0)
if (file.exists(master_p)) {
  mt <- read.delim(master_p, sep="\t")
  prio_lig <- unique(unlist(strsplit(mt$ligand[mt$direction=="tumor_to_TME"], "_")))
  log_info(sprintf("Prioritised tumour-derived ligands from the master table: %d", length(prio_lig)), stage=STAGE)
} else wrn("Master table not found; ligand activity will use all expressed tumour ligands.")

act_all <- list(); tgt_all <- list(); summ <- list()
for (RC in receivers) {
  if (!(RC %in% lab)) { wrn(sprintf("%s absent; skipped", RC)); next }
  n_rc <- sum(lab==RC)
  log_info(sprintf("=== receiver %s (%d cells) ===", RC, n_rc), stage=STAGE)
  rc_expressed <- expressed_in(RC)
  background <- intersect(rc_expressed, rownames(ltm))
  if (length(background) < 200) { wrn(sprintf("%s: only %d background genes in the prior model; skipped", RC, length(background))); next }

  # gene set of interest = receiver marker genes (cluster characterisation, not condition DE)
  mk <- tryCatch(FindMarkers(obj, ident.1=RC, assay="RNA", slot="data", test.use="wilcox",
                             min.pct=0.10, logfc.threshold=0.25, only.pos=TRUE, verbose=FALSE),
                 error=function(e) { wrn(sprintf("FindMarkers failed for %s: %s", RC, conditionMessage(e))); NULL })
  if (is.null(mk) || !nrow(mk)) { wrn(sprintf("%s: no markers; skipped", RC)); next }
  mk <- mk[mk$p_val_adj < 0.05, , drop=FALSE]
  geneset <- intersect(rownames(mk), background)
  log_info(sprintf("  gene set of interest: %d receiver marker genes in the prior model (background %d)",
    length(geneset), length(background)), stage=STAGE)
  if (length(geneset) < 20) { wrn(sprintf("%s: gene set too small (%d); skipped", RC, length(geneset))); next }

  # potential ligands: expressed by the tumour AND with a receptor expressed in the receiver
  lr_sub <- lrn %>% filter(from %in% sender_expressed, to %in% rc_expressed)
  potential <- intersect(unique(lr_sub$from), colnames(ltm))
  if (length(prio_lig)) {
    pr <- intersect(potential, prio_lig)
    log_info(sprintf("  potential tumour ligands: %d (of which %d are in the prioritised set)",
      length(potential), length(pr)), stage=STAGE)
  } else log_info(sprintf("  potential tumour ligands: %d", length(potential)), stage=STAGE)
  if (length(potential) < 5) { wrn(sprintf("%s: only %d potential ligands; skipped", RC, length(potential))); next }

  act <- tryCatch(predict_ligand_activities(geneset=geneset, background_expressed_genes=background,
                    ligand_target_matrix=ltm, potential_ligands=potential),
                  error=function(e) { wrn(sprintf("predict_ligand_activities failed for %s: %s", RC, conditionMessage(e))); NULL })
  if (is.null(act)) next
  score_col <- if ("aupr_corrected" %in% names(act)) "aupr_corrected" else if ("aupr" %in% names(act)) "aupr" else "pearson"
  act <- act %>% arrange(desc(.data[[score_col]])) %>% mutate(rank=row_number(), receiver=RC,
    n_receiver_cells=n_rc, score_column=score_col,
    in_prioritised_set=test_ligand %in% prio_lig)
  act_all[[RC]] <- act
  top <- act %>% head(top_n_ligands)
  log_info(sprintf("  top ligands by %s: %s", score_col,
    paste(sprintf("%s(%.3f%s)", top$test_ligand, top[[score_col]],
      ifelse(top$in_prioritised_set,"*","")), collapse=", ")), stage=STAGE)
  log_info("  (* = also prioritised by the LR concordance analysis)", stage=STAGE)

  lt <- tryCatch({
    do.call(rbind, lapply(top$test_ligand, function(lg) {
      w <- get_weighted_ligand_target_links(ligand=lg, geneset=geneset,
             ligand_target_matrix=ltm, n=200)
      if (is.null(w) || !nrow(w)) return(NULL)
      w$receiver <- RC; as.data.frame(w) })) },
    error=function(e) { wrn(sprintf("ligand-target links failed for %s: %s", RC, conditionMessage(e))); NULL })
  if (!is.null(lt) && nrow(lt)) tgt_all[[RC]] <- lt

  summ[[RC]] <- data.frame(receiver=RC, n_cells=n_rc, n_background=length(background),
    n_geneset=length(geneset), n_potential_ligands=length(potential),
    score_column=score_col, best_ligand=act$test_ligand[1], best_score=act[[score_col]][1],
    n_top_also_prioritised=sum(top$in_prioritised_set),
    top10_ligands=paste(head(act$test_ligand,10), collapse=";"), stringsAsFactors=FALSE)
}
if (!length(act_all)) fail("NicheNet produced no ligand activities for any receiver.")

acts <- dplyr::bind_rows(act_all)
write.table(acts, file.path(out_dir,"NicheNet_ligand_activity.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
if (length(tgt_all)) {
  tg <- dplyr::bind_rows(tgt_all)
  write.table(tg, file.path(out_dir,"receiver_target_programs.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
  log_info(sprintf("Ligand-target links: %d rows across %d receivers", nrow(tg), length(tgt_all)), stage=STAGE)
}
sm <- dplyr::bind_rows(summ)
write.table(sm, file.path(out_dir,"NicheNet_receiver_summary.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
log_info("NicheNet summary:", stage=STAGE)
for (r in seq_len(nrow(sm)))
  log_info(sprintf("  %-14s cells=%5d geneset=%4d ligands=%4d best=%-10s %s=%.4f  top20 overlapping prioritised: %d",
    sm$receiver[r], sm$n_cells[r], sm$n_geneset[r], sm$n_potential_ligands[r],
    sm$best_ligand[r], sm$score_column[r], sm$best_score[r], sm$n_top_also_prioritised[r]), stage=STAGE)

# receiver-response support flag for the evidence matrix: a tumour ligand is supported for a
# receiver if it ranks in that receiver's top N ligand activities.
rr <- acts %>% group_by(receiver) %>% slice_head(n=top_n_ligands) %>% ungroup() %>%
  transmute(receiver, ligand=test_ligand, nichenet_rank=rank,
            nichenet_score=.data[[acts$score_column[1]]], receiver_response_support=TRUE)
write.table(rr, file.path(out_dir,"NicheNet_receiver_response_support.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
log_info(sprintf("Receiver-response support flags: %d ligand x receiver pairs (top %d per receiver)",
  nrow(rr), top_n_ligands), stage=STAGE)

# ---- figures ----
th <- theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),
      plot.title=element_text(face="bold"), plot.subtitle=element_text(size=9,colour="grey30"))
CAV <- "NicheNet models regulatory potential from a prior network; it does not prove physical signalling."
frows <- list()
sf <- function(p,name,w,h,analysis,receiver,params){ for (ext in c("pdf","png")) {
  fp <- file.path(fig_dir, paste0(name,".",ext))
  ggsave(fp,p,width=w,height=h,units="in",dpi=200,device=ext,limitsize=FALSE)
  frows[[length(frows)+1]] <<- data.frame(figure_path=fp, phase="phase3", milestone="M23",
    analysis=analysis, sender=sender, receiver=receiver, method="NicheNet",
    input=obj_rds, input_checksum=in_md5,
    script="scripts/phase3/receiver_response/run_nichenet.R", parameters=params, stringsAsFactors=FALSE) }
  log_info(sprintf("  wrote %s", name), stage=STAGE) }
scol <- acts$score_column[1]
topa <- acts %>% group_by(receiver) %>% slice_head(n=15) %>% ungroup()
sf(ggplot(topa, aes(.data[[scol]], reorder(test_ligand, .data[[scol]]), fill=in_prioritised_set))+
   geom_col()+facet_wrap(~receiver, scales="free", ncol=3)+
   scale_fill_manual(values=c(`TRUE`="#B2182B",`FALSE`="grey65"),
     labels=c("NicheNet only","also LR-prioritised"), name=NULL)+
   labs(title=sprintf("Tumour-derived ligand activity per receiver (%s)", scol),
        subtitle=paste("Top 15 ligands predicted to explain each receiver's own marker programme.", CAV),
        x=scol, y=NULL)+th+theme(axis.text.y=element_text(size=6.5)),
   "M23_01_receiver_ligand_activity", 14, 10, "ligand_activity", "all", sprintf("top15;score=%s", scol))
if (length(tgt_all)) {
  tg <- dplyr::bind_rows(tgt_all)
  for (RC in names(tgt_all)[1:min(4,length(tgt_all))]) {
    d <- tg %>% filter(receiver==RC)
    tl <- d %>% group_by(ligand) %>% summarise(s=sum(weight), .groups="drop") %>% arrange(desc(s)) %>% head(12) %>% pull(ligand)
    tt <- d %>% filter(ligand %in% tl) %>% group_by(target) %>% summarise(s=sum(weight), .groups="drop") %>%
      arrange(desc(s)) %>% head(40) %>% pull(target)
    dd <- d %>% filter(ligand %in% tl, target %in% tt)
    if (!nrow(dd)) next
    sf(ggplot(dd, aes(target, ligand, fill=weight))+geom_tile(colour="white")+
       scale_fill_viridis_c(option="rocket", direction=-1, name="regulatory\npotential")+
       labs(title=sprintf("Ligand -> target regulatory potential: %s -> %s", sender, RC),
            subtitle=paste("Targets are that receiver's own marker genes.", CAV), x=NULL, y=NULL)+th+
       theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=6.5)),
       sprintf("M23_02_ligand_target_%s", gsub("[^A-Za-z0-9]","_",RC)),
       max(11,0.28*length(tt)), 6.5, "ligand_target", RC, "top12 ligands x top40 targets")
  }
}
sf(ggplot(sm, aes(reorder(receiver, n_top_also_prioritised), n_top_also_prioritised))+
   geom_col(fill="#3B6EA5")+geom_text(aes(label=n_top_also_prioritised), hjust=-0.2, size=3.4)+
   coord_flip(clip="off")+expand_limits(y=max(sm$n_top_also_prioritised)*1.25+1)+
   labs(title="Convergence of NicheNet and the LR concordance analysis",
        subtitle=sprintf("Of each receiver's top %d NicheNet ligands, how many were also prioritised by the LR frameworks. This is orthogonal agreement, not a shared vote.", top_n_ligands),
        x=NULL, y="overlapping ligands")+th,
   "M23_03_nichenet_lr_convergence", 9, 5.5, "convergence", "all", sprintf("top%d", top_n_ligands))
write.table(do.call(rbind, frows), file.path(out_dir,"figure_index_m23.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

rec <- list(milestone="M23", phase="phase3", timestamp=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),
  script="scripts/phase3/receiver_response/run_nichenet.R",
  nichenetr_version=as.character(packageVersion("nichenetr")),
  prior_model=lapply(names(urls), function(n) list(name=n, url=urls[[n]],
    file=file.path(nn_dir, paste0(n,".rds")), md5=digest(file.path(nn_dir, paste0(n,".rds")), file=TRUE, algo="md5"))),
  design=list(sender=sender, receivers=receivers,
    gene_set_of_interest="receiver population marker genes (receiver vs all other CCC populations; RNA LogNormalize; Wilcoxon; adj p < 0.05; only.pos) - cluster characterisation, NOT condition-level DE",
    background="genes detected in >=10% of receiver cells and present in the prior ligand-target matrix",
    potential_ligands="expressed in >=10% of MPNST-Tumor cells AND with a receptor expressed in the receiver",
    score=if (exists("scol")) scol else NA_character_, top_n_ligands=top_n_ligands, seed=seed),
  role="Orthogonal receiver-response evidence. NicheNet is NOT counted as a fourth LR framework and its AUPR is never averaged with LR scores.",
  results=lapply(seq_len(nrow(sm)), function(r) as.list(sm[r,])),
  warnings=warns, slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM"),
  git_commit=tryCatch(trimws(system("git rev-parse HEAD",intern=TRUE)), error=function(e) NA_character_),
  elapsed_seconds=as.numeric(difftime(Sys.time(),t0,units="secs")), session_info=capture.output(sessionInfo()))
write_json(rec, file.path(out_dir,"m23_nichenet_record.json"), auto_unbox=TRUE, pretty=TRUE, null="null", digits=NA)
prov <- record_provenance("phase3_m23_nichenet", inputs=list(object=obj_rds, master=master_p),
  outputs=list(activity=file.path(out_dir,"NicheNet_ligand_activity.tsv"),
    targets=file.path(out_dir,"receiver_target_programs.tsv"),
    record=file.path(out_dir,"m23_nichenet_record.json")),
  parameters=list(milestone="M23", nichenetr=as.character(packageVersion("nichenetr")), seed=seed,
    slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM")), dataset="combined")
save_provenance_json(prov, file.path(out_dir,"prov_m23_nichenet.json"))
log_info(sprintf("M23 COMPLETE in %.1f min. Warnings: %d", as.numeric(difftime(Sys.time(),t0,units="mins")), length(warns)), stage=STAGE)
log_system_usage(stage=STAGE)
