# scripts/phase3/concordance/integrate_evidence.R
#
# Phase 3 / Milestone M25 - CCC + receiver-response + LochNESS integration.
#
# Produces an EVIDENCE MATRIX, not a composite score. Phase 3 SS25/SS38 are explicit:
# these evidence streams are not commensurable, so they are displayed side by side and
# counted, never weighted-averaged into one number.
#
#   LR concordance | sample recurrence | expression | NicheNet receiver response | LochNESS
#
# LochNESS contributes ONLY receiver-state context. It is never converted into an LR score.

options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(RColorBrewer); library(viridis); library(igraph); library(jsonlite); library(digest)})
source("scripts/R/logging_utils.R"); source("scripts/R/provenance_utils.R"); setup_strict_logging()
STAGE <- "phase3_m25_integration"
args <- commandArgs(trailingOnly=TRUE)
master_p <- "results/phase3/ccc/prioritized/MPNST_CCC_MASTER_TABLE.tsv"
nn_p     <- "results/phase3/receiver_response/NicheNet_receiver_response_support.tsv"
loch_p   <- "results/phase3/lochness/lochness_summary.tsv"
lperm_p  <- "results/phase3/lochness/lochness_permutation_summary.tsv"
out_dir  <- "results/phase3/ccc/prioritized"; tab_dir <- "results/phase3/tables"
fig_dir  <- "results/phase3/figures/M25"
i <- 1
while (i <= length(args)) { a <- args[i]
  if (a=="--master"){master_p<-args[i+1];i<-i+2} else if (a=="--nichenet"){nn_p<-args[i+1];i<-i+2}
  else if (a=="--lochness"){loch_p<-args[i+1];i<-i+2} else if (a=="--loch-perm"){lperm_p<-args[i+1];i<-i+2}
  else if (a=="--out-dir"){out_dir<-args[i+1];i<-i+2} else if (a=="--tab-dir"){tab_dir<-args[i+1];i<-i+2}
  else if (a=="--fig-dir"){fig_dir<-args[i+1];i<-i+2} else stop(sprintf("Unknown argument: %s", a)) }
t0 <- Sys.time(); warns <- character(0)
wrn <- function(m){warns<<-c(warns,m); log_warn(m,stage=STAGE)}
for (d in c(out_dir,tab_dir,fig_dir)) dir.create(d,recursive=TRUE,showWarnings=FALSE)

log_info("=========== M25 EVIDENCE INTEGRATION ===========", stage=STAGE)
log_info("Evidence streams are displayed side by side and COUNTED. No weighted composite score is computed:", stage=STAGE)
log_info("  a CellChat probability, a CellPhoneDB p-value, a LIANA rank, a NicheNet AUPR and a LochNESS", stage=STAGE)
log_info("  enrichment are different quantities on different scales and are not commensurable.", stage=STAGE)

m <- read.delim(master_p, sep="\t", check.names=FALSE)
log_info(sprintf("Master table: %d rows", nrow(m)), stage=STAGE)

# ---- NicheNet receiver-response support (ligand x receiver) ----
if (file.exists(nn_p)) {
  nn <- read.delim(nn_p, sep="\t")
  m$ligand_first <- vapply(strsplit(m$ligand,"_"), function(z) z[1], character(1))
  m$NicheNet_ligand_support <- mapply(function(l, r) {
    any(nn$receiver==r & nn$ligand==l) }, m$ligand_first, m$receiver)
  m$receiver_response_support <- m$NicheNet_ligand_support
  log_info(sprintf("NicheNet receiver-response support: %d / %d master rows (ligand in that receiver's top NicheNet ligands)",
    sum(m$NicheNet_ligand_support), nrow(m)), stage=STAGE)
} else { wrn("NicheNet support table absent; receiver-response evidence unavailable.")
  m$NicheNet_ligand_support <- NA; m$receiver_response_support <- NA }

# ---- LochNESS receiver-state support (receiver level, NOT an LR score) ----
if (file.exists(loch_p)) {
  ls_ <- read.delim(loch_p, sep="\t")
  lp <- if (file.exists(lperm_p)) read.delim(lperm_p, sep="\t") else NULL
  ls_$lochness_structure <- abs(ls_$mean_lochness) >= 0.05 & !ls_$unstable_k
  if (!is.null(lp)) ls_ <- ls_ %>% left_join(lp %>% select(lineage, nullA_p_descriptive, nullB_p_descriptive), by="lineage")
  m$lochness_support <- ls_$lochness_structure[match(m$receiver, ls_$lineage)]
  m$lochness_mean_receiver <- ls_$mean_lochness[match(m$receiver, ls_$lineage)]
  m$lochness_k_unstable <- ls_$unstable_k[match(m$receiver, ls_$lineage)]
  log_info(sprintf("LochNESS receiver-state evidence available for %d receivers: %s",
    nrow(ls_), paste(sprintf("%s(mean %+.3f%s)", ls_$lineage, ls_$mean_lochness,
      ifelse(ls_$unstable_k," UNSTABLE","")), collapse=", ")), stage=STAGE)
  log_info("  NOTE: this is receiver-LINEAGE-level context, attached to every interaction with that receiver.", stage=STAGE)
  log_info("  It is NOT interaction-specific and is NEVER treated as ligand-receptor evidence.", stage=STAGE)
} else { wrn("LochNESS summary absent; receiver-state evidence unavailable.")
  m$lochness_support <- NA; m$lochness_mean_receiver <- NA; m$lochness_k_unstable <- NA }

# ---- evidence count: LR frameworks + recurrence + expression + receiver response ----
m$n_evidence_streams <- rowSums(cbind(
  LR2plus       = m$n_LR_methods_supported >= 2,
  recurrent3    = m$samples_supported >= 3,
  expression    = ifelse(is.na(m$expression_support), FALSE, m$expression_support),
  receiver_resp = ifelse(is.na(m$receiver_response_support), FALSE, m$receiver_response_support)), na.rm=TRUE)
m$evidence_profile <- with(m, paste0(
  ifelse(n_LR_methods_supported>=2,"LR",""), ifelse(samples_supported>=3,"+REC",""),
  ifelse(!is.na(expression_support)&expression_support,"+EXPR",""),
  ifelse(!is.na(receiver_response_support)&receiver_response_support,"+RR","")))
m$evidence_profile[m$evidence_profile==""] <- "none"
m <- m %>% arrange(desc(n_evidence_streams), overall_priority, desc(samples_supported))
write.table(m, file.path(out_dir,"MPNST_CCC_MASTER_TABLE.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(m, file.path(tab_dir,"MPNST_CCC_MASTER_TABLE.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
es <- as.data.frame(table(n_evidence_streams=m$n_evidence_streams))
log_info("Evidence-stream counts (max 4: >=2 LR frameworks, >=3/4 samples, expression, receiver response):", stage=STAGE)
for (r in seq_len(nrow(es))) log_info(sprintf("  %s streams : %6d interactions", es$n_evidence_streams[r], es$Freq[r]), stage=STAGE)

tc <- m %>% filter(direction != "TME_to_TME")
top <- tc %>% filter(n_evidence_streams >= 3) %>% arrange(desc(n_evidence_streams), desc(samples_supported))
log_info(sprintf("Tumour-centric interactions with >=3 evidence streams: %d", nrow(top)), stage=STAGE)
if (nrow(top)) for (r in seq_len(min(25,nrow(top))))
  log_info(sprintf("  [%d streams|%s] %s -> %s : %s -> %s (%d frameworks, %d/4 samples%s)",
    top$n_evidence_streams[r], top$evidence_profile[r], top$sender[r], top$receiver[r],
    top$ligand[r], top$receptor[r], top$n_LR_methods_supported[r], top$samples_supported[r],
    ifelse(isTRUE(top$receiver_response_support[r]), ", NicheNet RR", "")), stage=STAGE)
write.table(top, file.path(tab_dir,"CCC_INTEGRATED_EVIDENCE_TOP.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# ---- figures ----
th <- theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),
      plot.title=element_text(face="bold"), plot.subtitle=element_text(size=9,colour="grey30"))
CAV <- "Inferred communication potential; not evidence of physical adjacency or direct signalling."
frows <- list()
sf <- function(p,name,w,h,analysis,params){ for (ext in c("pdf","png")) {
  fp <- file.path(fig_dir, paste0(name,".",ext))
  ggsave(fp,p,width=w,height=h,units="in",dpi=200,device=ext,limitsize=FALSE)
  frows[[length(frows)+1]] <<- data.frame(figure_path=fp, phase="phase3", milestone="M25",
    analysis=analysis, sender="MPNST-Tumor", receiver="TME", method="integration",
    input=master_p, input_checksum=digest(master_p,file=TRUE,algo="md5"),
    script="scripts/phase3/concordance/integrate_evidence.R", parameters=params, stringsAsFactors=FALSE) }
  log_info(sprintf("  wrote %s", name), stage=STAGE) }

# evidence matrix heatmap for the top tumour-centric interactions
tm <- tc %>% arrange(desc(n_evidence_streams), desc(samples_supported)) %>% head(45) %>%
  mutate(lab=paste0(sender," -> ",receiver," | ",ligand,"->",receptor))
if (nrow(tm)) {
  em <- tm %>% transmute(lab,
    `LIANA`=LIANA_support, `CellChat`=CellChat_support, `CellPhoneDB`=CellPhoneDB_support,
    `>=3/4 samples`=samples_supported>=3, `expression`=ifelse(is.na(expression_support),FALSE,expression_support),
    `NicheNet RR`=ifelse(is.na(receiver_response_support),FALSE,receiver_response_support),
    `LochNESS state`=ifelse(is.na(lochness_support),FALSE,lochness_support)) %>%
    pivot_longer(-lab, names_to="evidence", values_to="present")
  em$evidence <- factor(em$evidence, levels=c("LIANA","CellChat","CellPhoneDB",">=3/4 samples","expression","NicheNet RR","LochNESS state"))
  sf(ggplot(em, aes(evidence, factor(lab, levels=rev(tm$lab)), fill=present))+geom_tile(colour="white")+
     scale_fill_manual(values=c(`TRUE`="#1A9850",`FALSE`="grey90"), labels=c("absent","present"), name=NULL)+
     labs(title="Integrated evidence matrix for the top tumour-centric interactions",
          subtitle=paste("Streams are shown side by side and counted, never averaged. 'LochNESS state' is receiver-lineage context, not LR evidence.", CAV),
          x=NULL, y=NULL)+th+
     theme(axis.text.y=element_text(size=6.4), axis.text.x=element_text(angle=35,hjust=1)),
     "M25_01_integrated_evidence_matrix", 10.5, 11, "evidence_matrix", "top45 tumour-centric")
}
sf(ggplot(tc, aes(factor(n_evidence_streams), fill=evidence_profile))+geom_bar()+
   facet_wrap(~direction)+
   labs(title="Evidence-stream depth for tumour-centric interactions",
        subtitle="Max 4 streams: >=2 LR frameworks, recurrence in >=3/4 samples, expression support, NicheNet receiver response.",
        x="evidence streams supporting", y="interactions", fill="profile")+th,
   "M25_02_evidence_depth", 12, 6, "evidence_depth", "counts by direction")

# integrated tumour-centric network, edges annotated by evidence depth
# exclude self-loops: geom_curve cannot draw identical endpoints
net <- tc %>% filter(sender != receiver, n_evidence_streams >= 2) %>% count(sender, receiver, name="n_interactions")
if (nrow(net)) {
  g <- graph_from_data_frame(net, directed=TRUE)
  set.seed(42); L <- layout_in_circle(g)
  vd <- data.frame(name=V(g)$name, x=L[,1], y=L[,2], stringsAsFactors=FALSE)
  IM <- c("Macrophage","Monocyte","Dendritic","Plasmacytoid-DC","CD8-T","CD4-T","NK","T-cell-other","B-cell","Plasma-cell")
  vd$compartment <- ifelse(vd$name=="MPNST-Tumor","Tumour", ifelse(vd$name %in% IM,"Immune",
                    ifelse(vd$name=="Endothelial","Endothelial","Stromal")))
  ed <- as.data.frame(as_edgelist(g)); names(ed) <- c("from","to"); ed$n <- E(g)$n_interactions
  ed <- ed %>% left_join(vd %>% select(name,x,y), by=c("from"="name")) %>% rename(x0=x,y0=y) %>%
    left_join(vd %>% select(name,x,y), by=c("to"="name")) %>% rename(x1=x,y1=y)
  sf(ggplot()+
     geom_curve(data=ed, aes(x0,y0,xend=x1,yend=y1,linewidth=n,colour=n), curvature=0.2, alpha=0.8,
                arrow=grid::arrow(length=unit(0.18,"cm"), type="closed"))+
     scale_linewidth_continuous(range=c(0.4,3), name="interactions\n(>=2 streams)")+
     scale_colour_viridis_c(option="rocket", direction=-1, name="interactions\n(>=2 streams)")+
     geom_point(data=vd, aes(x,y,fill=compartment), size=8, shape=21, colour="grey20")+
     scale_fill_manual(values=c(Tumour="#B2182B",Immune="#2166AC",Stromal="#D9A441",Endothelial="#1B7837"))+
     ggrepel::geom_text_repel(data=vd, aes(x,y,label=name), size=3.4, fontface="bold", box.padding=0.7)+
     coord_equal()+theme_void(base_size=12)+
     labs(title="Integrated MPNST tumour-microenvironment interaction map",
          subtitle=paste("Edges retained only where >=2 independent evidence streams support at least one interaction for that directed pair.", CAV)),
     "M25_03_integrated_tme_network", 10.5, 9.5, "network", ">=2 evidence streams")
}
write.table(do.call(rbind, frows), file.path(out_dir,"figure_index_m25.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

rec <- list(milestone="M25", phase="phase3", timestamp=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),
  script="scripts/phase3/concordance/integrate_evidence.R",
  design="Evidence matrix, not a composite score. Streams are counted, never weighted-averaged, because they are not commensurable.",
  evidence_streams=c(">=2 LR frameworks","recurrence in >=3/4 samples","ligand+receptor expression support","NicheNet receiver-response"),
  lochness_role="Receiver-LINEAGE-level context attached to interactions with that receiver. Not interaction-specific and never treated as ligand-receptor evidence.",
  results=list(master_rows=nrow(m), evidence_counts=as.list(setNames(es$Freq, es$n_evidence_streams)),
    tumour_centric=nrow(tc), tumour_centric_ge3_streams=nrow(top)),
  warnings=warns, slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM"),
  git_commit=tryCatch(trimws(system("git rev-parse HEAD",intern=TRUE)), error=function(e) NA_character_),
  elapsed_seconds=as.numeric(difftime(Sys.time(),t0,units="secs")), session_info=capture.output(sessionInfo()))
write_json(rec, file.path(out_dir,"m25_integration_record.json"), auto_unbox=TRUE, pretty=TRUE, null="null", digits=NA)
prov <- record_provenance("phase3_m25_integration",
  inputs=list(master=master_p, nichenet=nn_p, lochness=loch_p),
  outputs=list(master=file.path(out_dir,"MPNST_CCC_MASTER_TABLE.tsv"),
    top=file.path(tab_dir,"CCC_INTEGRATED_EVIDENCE_TOP.tsv"),
    record=file.path(out_dir,"m25_integration_record.json")),
  parameters=list(milestone="M25", slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM")), dataset="combined")
save_provenance_json(prov, file.path(out_dir,"prov_m25_integration.json"))
log_info(sprintf("M25 COMPLETE in %.1f min. Warnings: %d", as.numeric(difftime(Sys.time(),t0,units="mins")), length(warns)), stage=STAGE)
