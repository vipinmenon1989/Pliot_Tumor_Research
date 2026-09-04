# scripts/phase3/concordance/robustness.R
#
# Phase 3 / Milestone M26 - robustness and replication.
#
# Not a benchmarking study. The question is whether the Phase 3 conclusions survive dropping
# any single patient, and where the methods genuinely disagree. Negative results are kept.

options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(RColorBrewer); library(jsonlite); library(digest)})
source("scripts/R/logging_utils.R"); source("scripts/R/provenance_utils.R"); setup_strict_logging()
STAGE <- "phase3_m26_robustness"
args <- commandArgs(trailingOnly=TRUE)
master_p <- "results/phase3/ccc/prioritized/MPNST_CCC_MASTER_TABLE.tsv"
loo_loch <- "results/phase3/lochness/lochness_leave_one_sample_out.tsv"
out_dir  <- "results/phase3/tables"; fig_dir <- "results/phase3/figures/M26"
i <- 1
while (i <= length(args)) { a <- args[i]
  if (a=="--master"){master_p<-args[i+1];i<-i+2} else if (a=="--out-dir"){out_dir<-args[i+1];i<-i+2}
  else if (a=="--fig-dir"){fig_dir<-args[i+1];i<-i+2} else if (a=="--loo-lochness"){loo_loch<-args[i+1];i<-i+2}
  else stop(sprintf("Unknown argument: %s", a)) }
t0 <- Sys.time(); warns <- character(0)
wrn <- function(m){warns<<-c(warns,m); log_warn(m,stage=STAGE)}
for (d in c(out_dir,fig_dir)) dir.create(d,recursive=TRUE,showWarnings=FALSE)

log_info("=========== M26 ROBUSTNESS AND REPLICATION ===========", stage=STAGE)
m <- read.delim(master_p, sep="\t", check.names=FALSE)
samples <- c("MPNST_1","MPNST_2","MPNST_3","MPNST_4")
sup_list <- strsplit(ifelse(is.na(m$supported_samples), "", m$supported_samples), ";")

# ---- leave-one-sample-out: does the interaction still meet a >=2-sample bar? ----
log_info("Leave-one-sample-out on CCC support ...", stage=STAGE)
loo <- do.call(rbind, lapply(samples, function(s) {
  remaining <- vapply(sup_list, function(z) sum(setdiff(z, s) %in% setdiff(samples,s)), integer(1))
  data.frame(dropped_sample=s,
    n_interactions_total=nrow(m),
    n_supported_ge2_full=sum(m$samples_supported>=2),
    n_supported_ge2_after_drop=sum(remaining>=2),
    n_lost=sum(m$samples_supported>=2 & remaining<2),
    frac_retained=round(sum(remaining>=2)/max(1,sum(m$samples_supported>=2)),4),
    stringsAsFactors=FALSE) }))
write.table(loo, file.path(out_dir,"CCC_LEAVE_ONE_SAMPLE_OUT.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
for (r in seq_len(nrow(loo)))
  log_info(sprintf("  drop %-9s : %d/%d interactions with >=2-sample support retained (%.1f%%), %d lost",
    loo$dropped_sample[r], loo$n_supported_ge2_after_drop[r], loo$n_supported_ge2_full[r],
    100*loo$frac_retained[r], loo$n_lost[r]), stage=STAGE)

# ---- single-sample dominance ----
dom <- m %>% mutate(single_sample_only = samples_supported==1,
  driving_sample = vapply(sup_list, function(z) if (length(z)==1) z[1] else NA_character_, character(1)))
by_s <- dom %>% filter(single_sample_only) %>% count(driving_sample, name="n_single_sample_interactions")
log_info(sprintf("Interactions supported in exactly ONE sample: %d (%.1f%% of all supported)",
  sum(dom$single_sample_only), 100*mean(dom$single_sample_only)), stage=STAGE)
if (nrow(by_s)) for (r in seq_len(nrow(by_s)))
  log_info(sprintf("  driven only by %-9s : %d", by_s$driving_sample[r], by_s$n_single_sample_interactions[r]), stage=STAGE)
write.table(by_s, file.path(out_dir,"CCC_SINGLE_SAMPLE_DOMINANCE.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
tc_dom <- dom %>% filter(direction!="TME_to_TME", single_sample_only) %>%
  select(sender, receiver, ligand, receptor, driving_sample, n_LR_methods_supported, concordance_class, overall_priority)
write.table(tc_dom, file.path(out_dir,"CCC_TUMOUR_SINGLE_SAMPLE_INTERACTIONS.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
log_info(sprintf("  of which tumour-centric: %d (reported, not discarded - these are the least generalisable findings)", nrow(tc_dom)), stage=STAGE)

# ---- method disagreement, preserved ----
dis <- m %>% filter(concordance_class=="Discordant/ambiguous")
log_info(sprintf("Discordant/ambiguous interactions (testable by >=2 frameworks, supported by only 1): %d", nrow(dis)), stage=STAGE)
dsum <- dis %>% mutate(only_method=ifelse(LIANA_support,"LIANA only",
  ifelse(CellChat_support,"CellChat only", ifelse(CellPhoneDB_support,"CellPhoneDB only","none")))) %>%
  count(only_method, name="n") %>% arrange(desc(n))
for (r in seq_len(nrow(dsum))) log_info(sprintf("  %-20s %6d", dsum$only_method[r], dsum$n[r]), stage=STAGE)
write.table(dis, file.path(out_dir,"CCC_METHOD_DISAGREEMENT.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(dsum, file.path(out_dir,"CCC_METHOD_DISAGREEMENT_SUMMARY.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# ---- rare population sensitivity ----
rare <- c("NK","Monocyte","B-cell","T-cell-other","Plasmacytoid-DC")
rr <- m %>% filter(sender %in% rare | receiver %in% rare) %>%
  group_by(population=ifelse(sender %in% rare, sender, receiver)) %>%
  summarise(n_interactions=n(), n_ge2_frameworks=sum(n_LR_methods_supported>=2),
            n_ge2_samples=sum(samples_supported>=2), median_samples=median(samples_supported), .groups="drop")
write.table(rr, file.path(out_dir,"CCC_RARE_POPULATION_SENSITIVITY.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
log_info("Rare-population sensitivity (smallest CCC-ready populations):", stage=STAGE)
for (r in seq_len(nrow(rr)))
  log_info(sprintf("  %-16s interactions=%5d  >=2 frameworks=%4d  >=2 samples=%4d  median samples=%.1f",
    rr$population[r], rr$n_interactions[r], rr$n_ge2_frameworks[r], rr$n_ge2_samples[r], rr$median_samples[r]), stage=STAGE)

if (file.exists(loo_loch)) {
  ll <- read.delim(loo_loch, sep="\t")
  log_info("LochNESS leave-one-sample-out (from M24):", stage=STAGE)
  for (r in seq_len(nrow(ll)))
    log_info(sprintf("  %-14s %d variants, mean pairwise Spearman rho = %.3f (min %.3f)",
      ll$lineage[r], ll$n_variants[r], ll$mean_pairwise_spearman[r], ll$min_pairwise_spearman[r]), stage=STAGE)
} else wrn("LochNESS leave-one-sample-out table absent.")

# ---- figures ----
th <- theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),
      plot.title=element_text(face="bold"), plot.subtitle=element_text(size=9,colour="grey30"))
frows <- list()
sf <- function(p,name,w,h,analysis,params){ for (ext in c("pdf","png")) {
  fp <- file.path(fig_dir, paste0(name,".",ext))
  ggsave(fp,p,width=w,height=h,units="in",dpi=200,device=ext,limitsize=FALSE)
  frows[[length(frows)+1]] <<- data.frame(figure_path=fp, phase="phase3", milestone="M26",
    analysis=analysis, sender="all", receiver="all", method="robustness", input=master_p,
    input_checksum=digest(master_p,file=TRUE,algo="md5"),
    script="scripts/phase3/concordance/robustness.R", parameters=params, stringsAsFactors=FALSE) }
  log_info(sprintf("  wrote %s", name), stage=STAGE) }
sf(ggplot(loo, aes(dropped_sample, 100*frac_retained))+geom_col(fill="#3B6EA5")+
   geom_text(aes(label=sprintf("%.1f%%",100*frac_retained)), vjust=-0.4, size=3.4)+
   ylim(0,110)+labs(title="Leave-one-sample-out robustness of CCC support",
     subtitle="Percentage of interactions with >=2-sample support that retain >=2-sample support after dropping each patient.",
     x="patient dropped", y="% retained")+th,
   "M26_01_leave_one_sample_out", 8.5, 5.5, "loo", ">=2-sample bar")
if (nrow(by_s)) sf(ggplot(by_s, aes(reorder(driving_sample,-n_single_sample_interactions), n_single_sample_interactions))+
   geom_col(fill="#C1452B")+geom_text(aes(label=n_single_sample_interactions), vjust=-0.4, size=3.4)+
   labs(title="Interactions supported by only one patient",
        subtitle="These are the least generalisable findings. They are reported, not discarded.",
        x="the single supporting patient", y="interactions")+th,
   "M26_02_single_sample_dominance", 8.5, 5.5, "dominance", "single-sample support")
if (nrow(dsum)) sf(ggplot(dsum, aes(reorder(only_method,-n), n))+geom_col(fill="#F46D43")+
   geom_text(aes(label=n), vjust=-0.4, size=3.4)+
   labs(title="Method disagreement, preserved",
        subtitle="Interactions testable by >=2 frameworks but supported by only one. Disagreement is a result, not noise to be averaged away.",
        x=NULL, y="interactions")+th,
   "M26_03_method_disagreement", 8.5, 5.5, "disagreement", "discordant class")
sf(ggplot(m %>% filter(direction!="TME_to_TME"), aes(factor(samples_supported), fill=concordance_class))+
   geom_bar()+facet_wrap(~direction)+scale_fill_brewer(palette="Set2", name="concordance")+
   labs(title="Sample recurrence of tumour-centric interactions",
        subtitle="How many of the four patients support each interaction, split by concordance class.",
        x="patients supporting (of 4)", y="interactions")+th,
   "M26_04_sample_recurrence_distribution", 11, 6, "recurrence", "by direction and class")
write.table(do.call(rbind, frows), file.path(out_dir,"figure_index_m26.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

rec <- list(milestone="M26", phase="phase3", timestamp=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),
  script="scripts/phase3/concordance/robustness.R",
  scope="Robustness of Phase 3 conclusions to dropping any single patient, plus preservation of method disagreement. NOT a benchmarking study.",
  results=list(
    leave_one_sample_out=lapply(seq_len(nrow(loo)), function(r) as.list(loo[r,])),
    single_sample_supported=sum(dom$single_sample_only),
    single_sample_fraction=round(mean(dom$single_sample_only),4),
    tumour_centric_single_sample=nrow(tc_dom),
    discordant=nrow(dis),
    discordant_breakdown=as.list(setNames(dsum$n, dsum$only_method)),
    rare_population_sensitivity=lapply(seq_len(nrow(rr)), function(r) as.list(rr[r,]))),
  warnings=warns, slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM"),
  git_commit=tryCatch(trimws(system("git rev-parse HEAD",intern=TRUE)), error=function(e) NA_character_),
  elapsed_seconds=as.numeric(difftime(Sys.time(),t0,units="secs")), session_info=capture.output(sessionInfo()))
write_json(rec, file.path(out_dir,"m26_robustness_record.json"), auto_unbox=TRUE, pretty=TRUE, null="null", digits=NA)
prov <- record_provenance("phase3_m26_robustness", inputs=list(master=master_p),
  outputs=list(loo=file.path(out_dir,"CCC_LEAVE_ONE_SAMPLE_OUT.tsv"),
    disagreement=file.path(out_dir,"CCC_METHOD_DISAGREEMENT.tsv"),
    record=file.path(out_dir,"m26_robustness_record.json")),
  parameters=list(milestone="M26", slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM")), dataset="combined")
save_provenance_json(prov, file.path(out_dir,"prov_m26_robustness.json"))
log_info(sprintf("M26 COMPLETE in %.1f min. Warnings: %d", as.numeric(difftime(Sys.time(),t0,units="mins")), length(warns)), stage=STAGE)
