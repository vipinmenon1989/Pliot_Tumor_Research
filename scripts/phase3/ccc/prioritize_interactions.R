# scripts/phase3/ccc/prioritize_interactions.R
#
# Phase 3 / Milestone M22 - biological prioritisation of MPNST tumour <-> TME interactions.
#
# Priority is an EVIDENCE MATRIX, not a single algorithmic score. Where an ordering is
# needed it is lexicographic on (concordance class, sample recurrence, expression support),
# and the rule is printed. No weighted composite is invented (Phase 3 SS24, SS69).
#
# Also emits the sample-level tumour-derived signalling CONTEXT labels that M24 LochNESS
# consumes. Those labels are computed from MPNST-Tumor cells only, never from receiver cells.

options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(Seurat); library(SeuratObject); library(Matrix)
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork); library(RColorBrewer)
  library(jsonlite); library(digest)})
source("scripts/R/logging_utils.R"); source("scripts/R/provenance_utils.R"); setup_strict_logging()
STAGE <- "phase3_m22_prioritization"

args <- commandArgs(trailingOnly=TRUE)
obj_rds <- "results/phase3/ccc/ccc_input_object.rds"
conc_p  <- "results/phase3/ccc/concordance/CCC_CONCORDANCE.tsv"
out_dir <- "results/phase3/ccc/prioritized"
tab_dir <- "results/phase3/tables"
fig_dir <- "results/phase3/figures/M22"
loch_dir<- "results/phase3/lochness"
i <- 1
while (i <= length(args)) { a <- args[i]
  if (a=="--object"){obj_rds<-args[i+1];i<-i+2} else if (a=="--concordance"){conc_p<-args[i+1];i<-i+2}
  else if (a=="--out-dir"){out_dir<-args[i+1];i<-i+2} else if (a=="--tab-dir"){tab_dir<-args[i+1];i<-i+2}
  else if (a=="--fig-dir"){fig_dir<-args[i+1];i<-i+2} else if (a=="--loch-dir"){loch_dir<-args[i+1];i<-i+2}
  else stop(sprintf("Unknown argument: %s", a)) }
t0 <- Sys.time(); warns <- character(0)
wrn <- function(m){warns<<-c(warns,m); log_warn(m,stage=STAGE)}
fail <- function(m){log_error(m,stage=STAGE); stop(m,call.=FALSE)}
for (d in c(out_dir,tab_dir,fig_dir,loch_dir)) dir.create(d,recursive=TRUE,showWarnings=FALSE)

log_info("=========== M22 BIOLOGICAL PRIORITISATION ===========", stage=STAGE)
conc <- read.delim(conc_p, sep="\t", check.names=FALSE)
sup <- conc %>% filter(n_LR_methods_supported >= 1)
log_info(sprintf("Concordance keys with support: %d", nrow(sup)), stage=STAGE)

obj <- readRDS(obj_rds); DefaultAssay(obj) <- "RNA"
in_md5 <- digest(obj_rds, file=TRUE, algo="md5")
dat <- GetAssayData(obj, assay="RNA", layer="data")
lab <- as.character(obj$ccc_label); smp <- as.character(obj$ccc_sample)
pops <- sort(unique(lab)); samples <- sort(unique(smp))
genes <- rownames(dat)
log_info(sprintf("Expression matrix: %d genes x %d cells | %d populations | %d samples",
  nrow(dat), ncol(dat), length(pops), length(samples)), stage=STAGE)

# per-population mean expression and detection rate (whole cohort and per sample)
log_info("Computing per-population expression support ...", stage=STAGE)
mean_expr <- sapply(pops, function(p) Matrix::rowMeans(dat[, lab==p, drop=FALSE]))
det_rate  <- sapply(pops, function(p) { m <- dat[, lab==p, drop=FALSE]; Matrix::rowSums(m>0)/ncol(m) })
rownames(mean_expr) <- genes; rownames(det_rate) <- genes

# complexes: take the MINIMUM across subunits (a complex needs all subunits present)
cplx_stat <- function(cplx, mat, pop) {
  parts <- unlist(strsplit(cplx, "_")); parts <- intersect(parts, rownames(mat))
  if (!length(parts)) return(NA_real_)
  min(mat[parts, pop])
}
log_info("Annotating sender ligand and receiver receptor expression (complex = min across subunits) ...", stage=STAGE)
sup$sender_ligand_mean_expr   <- mapply(function(l,s) cplx_stat(l, mean_expr, s), sup$ligand_n, sup$sender)
sup$sender_ligand_detect_frac <- mapply(function(l,s) cplx_stat(l, det_rate,  s), sup$ligand_n, sup$sender)
sup$receiver_receptor_mean_expr   <- mapply(function(r,s) cplx_stat(r, mean_expr, s), sup$receptor_n, sup$receiver)
sup$receiver_receptor_detect_frac <- mapply(function(r,s) cplx_stat(r, det_rate,  s), sup$receptor_n, sup$receiver)
sup$expression_support <- with(sup,
  !is.na(sender_ligand_detect_frac) & !is.na(receiver_receptor_detect_frac) &
  sender_ligand_detect_frac >= 0.10 & receiver_receptor_detect_frac >= 0.10)
log_info(sprintf("Expression support (ligand and receptor each detected in >=10%% of the relevant population): %d / %d",
  sum(sup$expression_support, na.rm=TRUE), nrow(sup)), stage=STAGE)
log_info(sprintf("  genes not found in the RNA assay for %d ligand and %d receptor entries",
  sum(is.na(sup$sender_ligand_mean_expr)), sum(is.na(sup$receiver_receptor_mean_expr))), stage=STAGE)

# ---- priority tiers: transparent, lexicographic, no invented composite ----
cls_rank <- c("High concordance"=3L, "Moderate concordance"=2L,
              "Single-method (only one framework could test it)"=1L,
              "Discordant/ambiguous"=1L, "Not supported"=0L)
sup$concordance_rank <- unname(cls_rank[sup$concordance_class]); sup$concordance_rank[is.na(sup$concordance_rank)] <- 0L
sup$overall_priority <- with(sup, case_when(
  concordance_rank >= 3 & samples_supported == 4 & expression_support ~ "P1_highest",
  concordance_rank >= 3 & samples_supported >= 3 & expression_support ~ "P2_high",
  concordance_rank >= 2 & samples_supported >= 3 & expression_support ~ "P3_moderate",
  concordance_rank >= 2 & samples_supported >= 2 & expression_support ~ "P4_moderate_low",
  samples_supported >= 2 & expression_support ~ "P5_low",
  TRUE ~ "P6_exploratory"))
log_info("Priority rule (printed, not hidden):", stage=STAGE)
log_info("  P1_highest      : all testable frameworks agree AND supported in 4/4 samples AND expression support", stage=STAGE)
log_info("  P2_high         : all testable frameworks agree AND >=3/4 samples AND expression support", stage=STAGE)
log_info("  P3_moderate     : >=2 frameworks AND >=3/4 samples AND expression support", stage=STAGE)
log_info("  P4_moderate_low : >=2 frameworks AND >=2/4 samples AND expression support", stage=STAGE)
log_info("  P5_low          : >=2/4 samples AND expression support", stage=STAGE)
log_info("  P6_exploratory  : everything else retained but flagged", stage=STAGE)
pt <- as.data.frame(table(sup$overall_priority))
for (r in seq_len(nrow(pt))) log_info(sprintf("  %-16s %6d", pt$Var1[r], pt$Freq[r]), stage=STAGE)

master <- sup %>% transmute(
  sender, receiver, ligand=ligand_n, receptor=receptor_n,
  pathway=if ("cellchat_pathway" %in% names(sup)) cellchat_pathway else NA_character_,
  classification=if ("cpdb_classification" %in% names(sup)) cpdb_classification else NA_character_,
  direction,
  LIANA_support=liana_sup, CellChat_support=cellchat_sup, CellPhoneDB_support=cpdb_sup,
  LIANA_testable=liana_testable, CellChat_testable=cellchat_testable, CellPhoneDB_testable=cpdb_testable,
  n_LR_methods_supported, n_methods_testable, concordance_class,
  NicheNet_ligand_support=NA, receiver_response_support=NA, lochness_support=NA,
  sender_expression=sender_ligand_mean_expr, sender_detect_frac=sender_ligand_detect_frac,
  receiver_expression=receiver_receptor_mean_expr, receiver_detect_frac=receiver_receptor_detect_frac,
  expression_support,
  samples_evaluable, samples_supported, sample_recurrence_fraction, supported_samples,
  patients_evaluable, patients_supported, patient_recurrence_fraction,
  liana_best_rank, cpdb_min_pvalue=if ("cpdb_min_pvalue" %in% names(sup)) cpdb_min_pvalue else NA_real_,
  literature_support=NA, overall_priority,
  notes="scRNA-seq LR analysis infers communication POTENTIAL; it does not establish physical adjacency or direct signalling."
) %>% arrange(overall_priority, desc(n_LR_methods_supported), desc(samples_supported), liana_best_rank)
write.table(master, file.path(out_dir,"MPNST_CCC_MASTER_TABLE.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(master, file.path(tab_dir,"MPNST_CCC_MASTER_TABLE.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
log_info(sprintf("Master table written: %d rows", nrow(master)), stage=STAGE)

for (dd in c("tumor_to_TME","TME_to_tumor")) {
  m <- master %>% filter(direction==dd)
  top <- m %>% filter(overall_priority %in% c("P1_highest","P2_high","P3_moderate"))
  log_info(sprintf("%s: %d supported | %d in P1-P3", dd, nrow(m), nrow(top)), stage=STAGE)
  if (nrow(top)) { tt <- top %>% arrange(overall_priority, desc(samples_supported)) %>% head(20)
    for (r in seq_len(nrow(tt))) log_info(sprintf("    [%s] %s -> %s : %s -> %s  (%d frameworks, %d/4 samples)",
      tt$overall_priority[r], tt$sender[r], tt$receiver[r], tt$ligand[r], tt$receptor[r],
      tt$n_LR_methods_supported[r], tt$samples_supported[r]), stage=STAGE) }
}

# ---- LochNESS context labels: sample-level, from TUMOUR cells only ----
log_info("Deriving the sample-level tumour-derived signalling context for M24 LochNESS ...", stage=STAGE)
t2i <- master %>% filter(direction=="tumor_to_TME",
  overall_priority %in% c("P1_highest","P2_high","P3_moderate"),
  receiver %in% c("Macrophage","CD8-T","Fibroblast","Endothelial"))
if (!nrow(t2i)) { t2i <- master %>% filter(direction=="tumor_to_TME") %>%
    arrange(desc(n_LR_methods_supported), desc(samples_supported)) %>% head(50)
  wrn("No P1-P3 tumour->receiver interactions; falling back to the top 50 tumour->TME rows for the context label.") }
lig_rank <- t2i %>% group_by(ligand) %>%
  summarise(n_receivers=dplyr::n_distinct(receiver), n_frameworks=max(n_LR_methods_supported),
            n_samples=max(samples_supported), .groups="drop") %>%
  arrange(desc(n_frameworks), desc(n_samples), desc(n_receivers))
ctx_ligand <- NA_character_
for (g in lig_rank$ligand) { parts <- intersect(unlist(strsplit(g,"_")), genes)
  if (length(parts)) { ctx_ligand <- g; break } }
if (is.na(ctx_ligand)) fail("No prioritised tumour ligand is present in the RNA assay.")
log_info(sprintf("Context ligand selected: %s (top of the prioritised tumour->receiver ligand ranking)", ctx_ligand), stage=STAGE)
parts <- intersect(unlist(strsplit(ctx_ligand,"_")), genes)
tum_cells <- lab == "MPNST-Tumor"
ctx_tbl <- do.call(rbind, lapply(samples, function(s) {
  sel <- tum_cells & smp == s
  v <- if (length(parts)==1) dat[parts, sel, drop=TRUE] else Matrix::colMeans(dat[parts, sel, drop=FALSE])
  data.frame(sample_id=s, n_tumor_cells=sum(sel), ligand=ctx_ligand,
    mean_ligand_expr=mean(v), detect_frac=mean(v>0), stringsAsFactors=FALSE) }))
ctx_tbl <- ctx_tbl %>% arrange(desc(mean_ligand_expr))
n_high <- floor(nrow(ctx_tbl)/2)
ctx_tbl$context <- c(rep("high", n_high), rep("low", nrow(ctx_tbl)-n_high))
write.table(ctx_tbl, file.path(loch_dir,"lochness_context_labels.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
log_info("Context labels (from MPNST-Tumor cells only; receiver cells never contribute):", stage=STAGE)
for (r in seq_len(nrow(ctx_tbl)))
  log_info(sprintf("  %-9s mean %s expr in tumour = %.4f (detected in %.1f%% of %d tumour cells) -> %s",
    ctx_tbl$sample_id[r], ctx_ligand, ctx_tbl$mean_ligand_expr[r], 100*ctx_tbl$detect_frac[r],
    ctx_tbl$n_tumor_cells[r], ctx_tbl$context[r]), stage=STAGE)

# ---- figures ----
th <- theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),
      plot.title=element_text(face="bold"), plot.subtitle=element_text(size=9,colour="grey30"))
CAV <- "Inferred communication potential; not evidence of physical adjacency or direct signalling."
frows <- list()
sf <- function(p,name,w,h,analysis,params){ for (ext in c("pdf","png")) {
  fp <- file.path(fig_dir, paste0(name,".",ext))
  ggsave(fp,p,width=w,height=h,units="in",dpi=200,device=ext,limitsize=FALSE)
  frows[[length(frows)+1]] <<- data.frame(figure_path=fp, phase="phase3", milestone="M22",
    analysis=analysis, sender="MPNST-Tumor", receiver="TME", method="prioritisation",
    input=conc_p, input_checksum=digest(conc_p,file=TRUE,algo="md5"),
    script="scripts/phase3/ccc/prioritize_interactions.R", parameters=params, stringsAsFactors=FALSE) }
  log_info(sprintf("  wrote %s", name), stage=STAGE) }

sf(ggplot(master %>% filter(direction!="TME_to_TME"), aes(overall_priority, fill=concordance_class))+
   geom_bar()+facet_wrap(~direction)+scale_fill_brewer(palette="Set2", name="concordance")+
   labs(title="Prioritisation tiers for tumour-centric interactions",
        subtitle=paste("Tiers are lexicographic on concordance, sample recurrence and expression support. No weighted composite score was invented.", CAV),
        x=NULL, y="interactions")+th+theme(axis.text.x=element_text(angle=45,hjust=1)),
   "M22_01_priority_tiers", 11, 6, "priority", "lexicographic tiers")
topm <- master %>% filter(direction!="TME_to_TME", overall_priority %in% c("P1_highest","P2_high","P3_moderate")) %>%
  mutate(lab=paste0(sender," -> ",receiver," | ",ligand,"->",receptor)) %>%
  arrange(overall_priority, desc(samples_supported)) %>% head(40)
if (nrow(topm)) sf(ggplot(topm, aes(samples_supported, factor(lab, levels=rev(lab)),
    colour=overall_priority, size=n_LR_methods_supported))+geom_point()+
   scale_colour_brewer(palette="Dark2", name="priority")+
   scale_size_continuous(range=c(2,5.5), name="frameworks", breaks=1:3)+
   scale_x_continuous(breaks=1:4, limits=c(0.5,4.5))+
   labs(title="Top prioritised tumour-centric interactions",
        subtitle=paste("Ordered by tier then sample recurrence.", CAV), x="samples supporting (of 4)", y=NULL)+th+
   theme(axis.text.y=element_text(size=7)),
   "M22_02_top_prioritised_interactions", 11, 10, "priority", "P1-P3, top 40")
sf(ggplot(ctx_tbl, aes(reorder(sample_id, -mean_ligand_expr), mean_ligand_expr, fill=context))+
   geom_col()+geom_text(aes(label=sprintf("%.3f", mean_ligand_expr)), vjust=-0.4, size=3.3)+
   scale_fill_brewer(palette="Set1", name="context")+
   labs(title=sprintf("Sample-level tumour-derived signalling context: %s", ctx_ligand),
        subtitle="Mean expression in MPNST-Tumor cells only. Receiver cells never contribute to this label, which is what keeps the M24 LochNESS design non-circular.",
        x="sample (= patient)", y=sprintf("mean %s expression in tumour cells", ctx_ligand))+th,
   "M22_03_lochness_context_label", 8, 5.5, "context", sprintf("ligand=%s", ctx_ligand))
write.table(do.call(rbind, frows), file.path(out_dir,"figure_index_m22.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

rec <- list(milestone="M22", phase="phase3", timestamp=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),
  script="scripts/phase3/ccc/prioritize_interactions.R",
  inputs=list(object=obj_rds, object_md5=in_md5, concordance=conc_p),
  priority_rule=list(type="lexicographic tiers, no weighted composite",
    P1_highest="all testable frameworks agree AND 4/4 samples AND expression support",
    P2_high="all testable frameworks agree AND >=3/4 samples AND expression support",
    P3_moderate=">=2 frameworks AND >=3/4 samples AND expression support",
    P4_moderate_low=">=2 frameworks AND >=2/4 samples AND expression support",
    P5_low=">=2/4 samples AND expression support", P6_exploratory="retained but flagged"),
  expression_support_rule="ligand and receptor each detected in >=10% of the relevant population; complexes scored as the MINIMUM across subunits",
  results=list(supported_keys=nrow(sup), priority_counts=as.list(setNames(pt$Freq, pt$Var1)),
    expression_supported=sum(sup$expression_support, na.rm=TRUE),
    tumor_to_TME=sum(master$direction=="tumor_to_TME"), TME_to_tumor=sum(master$direction=="TME_to_tumor")),
  lochness_context=list(ligand=ctx_ligand,
    derivation="mean expression of the prioritised tumour-derived ligand in MPNST-Tumor cells only, per sample; top half labelled high",
    non_circularity="Label uses tumour cells only; the LochNESS score is computed in receiver cell space. Sender and receiver cell sets are disjoint.",
    labels=lapply(seq_len(nrow(ctx_tbl)), function(r) as.list(ctx_tbl[r,]))),
  warnings=warns, slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM"),
  git_commit=tryCatch(trimws(system("git rev-parse HEAD",intern=TRUE)), error=function(e) NA_character_),
  elapsed_seconds=as.numeric(difftime(Sys.time(),t0,units="secs")), session_info=capture.output(sessionInfo()))
write_json(rec, file.path(out_dir,"m22_prioritization_record.json"), auto_unbox=TRUE, pretty=TRUE, null="null", digits=NA)
prov <- record_provenance("phase3_m22_prioritization", inputs=list(object=obj_rds, concordance=conc_p),
  outputs=list(master=file.path(out_dir,"MPNST_CCC_MASTER_TABLE.tsv"),
    context=file.path(loch_dir,"lochness_context_labels.tsv"),
    record=file.path(out_dir,"m22_prioritization_record.json")),
  parameters=list(milestone="M22", context_ligand=ctx_ligand,
    slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM")), dataset="combined")
save_provenance_json(prov, file.path(out_dir,"prov_m22_prioritization.json"))
log_info(sprintf("M22 COMPLETE in %.1f min. Warnings: %d", as.numeric(difftime(Sys.time(),t0,units="mins")), length(warns)), stage=STAGE)
