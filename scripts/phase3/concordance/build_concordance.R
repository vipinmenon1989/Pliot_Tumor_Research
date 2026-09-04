# scripts/phase3/concordance/build_concordance.R
#
# Phase 3 / Milestone M21 - CCC concordance across independent frameworks.
#
# DESIGN PRINCIPLES
#  1. Raw scores are NEVER averaged across tools. A CellChat probability, a CellPhoneDB
#     p-value and a LIANA rank are different quantities on different scales. Each method
#     contributes a boolean SUPPORT flag under its own documented rule.
#  2. Concordance is scored against what each method could actually TEST. The three
#     frameworks use different LR resources (OmniPath/Consensus, CellChatDB, CellPhoneDB v5),
#     so "1/3 support" is meaningless unless we know whether the other two ever tested the
#     interaction. Every row therefore carries `testable_*` alongside `supported_*`.
#  3. Method DISAGREEMENT is preserved as an explicit `Discordant` class, never averaged away.
#  4. NicheNet is NOT counted here as a fourth LR framework - it is orthogonal
#     receiver-response evidence and is joined in M23/M25.
#
# KNOWN PARTIAL NON-INDEPENDENCE: LIANA's method set includes a CellPhoneDB-style score, so
# LIANA and the standalone CellPhoneDB run are not fully independent. Reported, not hidden.

options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(RColorBrewer); library(jsonlite); library(digest)})
source("scripts/R/logging_utils.R"); source("scripts/R/provenance_utils.R"); setup_strict_logging()
STAGE <- "phase3_m21_concordance"

args <- commandArgs(trailingOnly=TRUE)
by_dir  <- "results/phase3/ccc/by_method"
out_dir <- "results/phase3/ccc/concordance"
tab_dir <- "results/phase3/tables"
fig_dir <- "results/phase3/figures/M21"
i <- 1
while (i <= length(args)) { a <- args[i]
  if (a=="--by-dir"){by_dir<-args[i+1];i<-i+2} else if (a=="--out-dir"){out_dir<-args[i+1];i<-i+2}
  else if (a=="--tab-dir"){tab_dir<-args[i+1];i<-i+2} else if (a=="--fig-dir"){fig_dir<-args[i+1];i<-i+2}
  else stop(sprintf("Unknown argument: %s", a)) }
t0 <- Sys.time(); warns <- character(0)
wrn <- function(m){warns<<-c(warns,m); log_warn(m,stage=STAGE)}
fail <- function(m){log_error(m,stage=STAGE); stop(m,call.=FALSE)}
for (d in c(out_dir,tab_dir,fig_dir)) dir.create(d,recursive=TRUE,showWarnings=FALSE)

log_info("=============== M21 CCC CONCORDANCE ===============", stage=STAGE)

# ---- gene / complex symbol normalisation (adapted in spirit from perturbseq gene_sets.py) ----
norm_complex <- function(x) {
  x <- toupper(trimws(as.character(x)))
  vapply(x, function(v) {
    if (is.na(v) || !nzchar(v)) return(NA_character_)
    parts <- unlist(strsplit(v, "[_+&:/]"))
    parts <- trimws(parts); parts <- parts[nzchar(parts)]
    parts <- sub("^COMPLEX$", "", parts); parts <- parts[nzchar(parts)]
    if (!length(parts)) return(NA_character_)
    paste(sort(unique(parts)), collapse="_")
  }, character(1), USE.NAMES=FALSE)
}
mk_key <- function(s,r,l,rc) paste(s, r, norm_complex(l), norm_complex(rc), sep="|")

# =========================== LIANA ===========================
lp <- file.path(by_dir,"all_LIANA_interactions.tsv.gz")
if (!file.exists(lp)) lp <- file.path(by_dir,"all_LIANA_interactions.tsv")
if (!file.exists(lp)) fail("LIANA output not found.")
log_info(sprintf("Reading LIANA: %s", lp), stage=STAGE)
li <- read.delim(lp, sep="\t", check.names=FALSE)
rank_col <- if ("aggregate_rank" %in% names(li)) "aggregate_rank" else "mean_rank"
LIANA_CUT <- 0.05
li <- li %>% mutate(sender=source, receiver=target,
  ligand_n=norm_complex(ligand.complex), receptor_n=norm_complex(receptor.complex),
  key=paste(sender,receiver,ligand_n,receptor_n,sep="|"),
  liana_sup=.data[[rank_col]] <= LIANA_CUT)
log_info(sprintf("  LIANA rows %d | distinct keys %d | supported rows %d (%s <= %.2f)",
  nrow(li), dplyr::n_distinct(li$key), sum(li$liana_sup), rank_col, LIANA_CUT), stage=STAGE)
liana_by_key <- li %>% group_by(key, sender, receiver, ligand_n, receptor_n) %>%
  summarise(liana_testable=TRUE,
            liana_samples_supported=sum(liana_sup),
            liana_samples_tested=dplyr::n_distinct(sample_id),
            liana_best_rank=min(.data[[rank_col]], na.rm=TRUE),
            liana_supported_samples=paste(sort(unique(sample_id[liana_sup])), collapse=";"),
            .groups="drop") %>%
  mutate(liana_sup=liana_samples_supported > 0)

# =========================== CellChat ===========================
# CellChat was run as one SLURM job per sample (parallelism, identical results), so the
# per-sample tables are globbed and bound here.
cc_files <- sort(Sys.glob(file.path(by_dir,"all_CellChat_interactions*.tsv")))
cp <- paste(basename(cc_files), collapse=";")
cc_by_key <- NULL; cellchat_ok <- FALSE
if (length(cc_files)) {
  log_info(sprintf("Reading CellChat from %d file(s): %s", length(cc_files), cp), stage=STAGE)
  cc <- dplyr::bind_rows(lapply(cc_files, function(f) read.delim(f, sep="\t", check.names=FALSE)))
  log_info(sprintf("  CellChat samples present: %s", paste(sort(unique(cc$sample_id)), collapse=", ")), stage=STAGE)
  lig_col <- if ("ligand" %in% names(cc)) "ligand" else NA
  rec_col <- if ("receptor" %in% names(cc)) "receptor" else NA
  if (is.na(lig_col) || is.na(rec_col)) { wrn("CellChat table lacks ligand/receptor columns; CellChat excluded.")
  } else {
    cc <- cc %>% mutate(sender=source, receiver=target,
      ligand_n=norm_complex(.data[[lig_col]]), receptor_n=norm_complex(.data[[rec_col]]),
      key=paste(sender,receiver,ligand_n,receptor_n,sep="|"))
    # subsetCommunication() returns only significant links, so presence = support.
    cc_by_key <- cc %>% group_by(key) %>%
      summarise(cellchat_testable=TRUE, cellchat_samples_supported=dplyr::n_distinct(sample_id),
                cellchat_supported_samples=paste(sort(unique(sample_id)), collapse=";"),
                cellchat_pathway=paste(sort(unique(na.omit(pathway_name))), collapse=";"),
                cellchat_max_prob=if ("prob" %in% names(cc)) max(prob, na.rm=TRUE) else NA_real_,
                .groups="drop") %>% mutate(cellchat_sup=TRUE)
    cellchat_ok <- TRUE
    log_info(sprintf("  CellChat rows %d | distinct keys %d | pathways %d",
      nrow(cc), nrow(cc_by_key), dplyr::n_distinct(cc$pathway_name)), stage=STAGE)
  }
} else wrn("CellChat output not found; the LR-framework denominator drops from 3 to 2.")

# =========================== CellPhoneDB ===========================
pp <- file.path(by_dir,"all_CellPhoneDB_interactions.tsv.gz")
if (!file.exists(pp)) pp <- file.path(by_dir,"all_CellPhoneDB_interactions.tsv")
cpdb_by_key <- NULL; cpdb_ok <- FALSE
if (file.exists(pp)) {
  log_info(sprintf("Reading CellPhoneDB: %s", pp), stage=STAGE)
  cd <- read.delim(pp, sep="\t", check.names=FALSE)
  # gene_a/gene_b are populated for simple pairs; fall back to splitting interacting_pair
  ip <- strsplit(as.character(cd$interacting_pair), "_")
  ga <- ifelse(!is.na(cd$gene_a) & nzchar(as.character(cd$gene_a)), as.character(cd$gene_a),
               vapply(ip, function(z) if (length(z)) z[1] else NA_character_, character(1)))
  gb <- ifelse(!is.na(cd$gene_b) & nzchar(as.character(cd$gene_b)), as.character(cd$gene_b),
               vapply(ip, function(z) if (length(z)>1) paste(z[-1], collapse="_") else NA_character_, character(1)))
  # pandas writes booleans as the strings "True"/"False"; coerce to logical explicitly.
  as_logical_py <- function(x) {
    if (is.logical(x)) return(x)
    tolower(trimws(as.character(x))) %in% c("true","t","1")
  }
  cd$cellphonedb_supported <- as_logical_py(cd$cellphonedb_supported)
  cd$pvalue <- suppressWarnings(as.numeric(cd$pvalue))
  cd$mean_expr <- suppressWarnings(as.numeric(cd$mean_expr))
  log_info(sprintf("  CellPhoneDB support column coerced to logical: %d TRUE of %d rows",
    sum(cd$cellphonedb_supported, na.rm=TRUE), nrow(cd)), stage=STAGE)
  cd <- cd %>% mutate(sender=source, receiver=target,
    ligand_n=norm_complex(ga), receptor_n=norm_complex(gb),
    key=paste(sender,receiver,ligand_n,receptor_n,sep="|"))
  cpdb_by_key <- cd %>% group_by(key) %>%
    summarise(cpdb_testable=TRUE,
              cpdb_samples_supported=sum(cellphonedb_supported, na.rm=TRUE),
              cpdb_samples_tested=dplyr::n_distinct(sample_id),
              cpdb_min_pvalue=suppressWarnings(min(pvalue, na.rm=TRUE)),
              cpdb_supported_samples=paste(sort(unique(sample_id[cellphonedb_supported])), collapse=";"),
              cpdb_classification=paste(sort(unique(na.omit(classification))), collapse=";"),
              .groups="drop") %>% mutate(cpdb_sup=cpdb_samples_supported > 0)
  cpdb_ok <- TRUE
  log_info(sprintf("  CellPhoneDB rows %d | distinct keys %d | keys with support %d",
    nrow(cd), nrow(cpdb_by_key), sum(cpdb_by_key$cpdb_sup)), stage=STAGE)
} else wrn("CellPhoneDB output not found; the LR-framework denominator drops.")

n_frameworks <- 1L + as.integer(cellchat_ok) + as.integer(cpdb_ok)
log_info(sprintf("LR frameworks contributing: %d (LIANA%s%s)", n_frameworks,
  if (cellchat_ok) ", CellChat" else "", if (cpdb_ok) ", CellPhoneDB" else ""), stage=STAGE)

# =========================== JOIN ===========================
conc <- liana_by_key %>% select(key, sender, receiver, ligand_n, receptor_n,
  liana_testable, liana_sup, liana_samples_supported, liana_best_rank, liana_supported_samples)
if (cellchat_ok) conc <- conc %>% full_join(cc_by_key, by="key")
if (cpdb_ok)     conc <- conc %>% full_join(cpdb_by_key, by="key")
# rebuild identity columns for keys contributed only by CellChat/CellPhoneDB
kp <- do.call(rbind, strsplit(conc$key, "|", fixed=TRUE))
conc$sender     <- ifelse(is.na(conc$sender), kp[,1], conc$sender)
conc$receiver   <- ifelse(is.na(conc$receiver), kp[,2], conc$receiver)
conc$ligand_n   <- ifelse(is.na(conc$ligand_n), kp[,3], conc$ligand_n)
conc$receptor_n <- ifelse(is.na(conc$receptor_n), kp[,4], conc$receptor_n)
z <- function(x) ifelse(is.na(x), FALSE, x)
n0 <- function(x) ifelse(is.na(x), 0L, as.integer(x))
conc <- conc %>% mutate(
  liana_testable=z(liana_testable), liana_sup=z(liana_sup),
  cellchat_testable=if (cellchat_ok) z(cellchat_testable) else FALSE,
  cellchat_sup=if (cellchat_ok) z(cellchat_sup) else FALSE,
  cpdb_testable=if (cpdb_ok) z(cpdb_testable) else FALSE,
  cpdb_sup=if (cpdb_ok) z(cpdb_sup) else FALSE,
  liana_samples_supported=n0(liana_samples_supported),
  cellchat_samples_supported=if (cellchat_ok) n0(cellchat_samples_supported) else 0L,
  cpdb_samples_supported=if (cpdb_ok) n0(cpdb_samples_supported) else 0L)
conc <- conc %>% mutate(
  n_methods_testable = as.integer(liana_testable) + as.integer(cellchat_testable) + as.integer(cpdb_testable),
  n_LR_methods_supported = as.integer(liana_sup) + as.integer(cellchat_sup) + as.integer(cpdb_sup))

# ---- concordance classes, defined transparently from what was testable ----
conc <- conc %>% mutate(concordance_class = case_when(
  n_LR_methods_supported >= 3 ~ "High concordance",
  n_LR_methods_supported == 2 & n_methods_testable == 2 ~ "High concordance",
  n_LR_methods_supported == 2 ~ "Moderate concordance",
  n_LR_methods_supported == 1 & n_methods_testable == 1 ~ "Single-method (only one framework could test it)",
  n_LR_methods_supported == 1 & n_methods_testable >= 2 ~ "Discordant/ambiguous",
  n_LR_methods_supported == 0 ~ "Not supported",
  TRUE ~ "Unclassified"))
log_info("Concordance class definitions:", stage=STAGE)
log_info("  High concordance      : supported by all frameworks that could test it (>=2)", stage=STAGE)
log_info("  Moderate concordance  : supported by 2 of 3 testable frameworks", stage=STAGE)
log_info("  Single-method         : supported by the ONLY framework whose resource contained it", stage=STAGE)
log_info("  Discordant/ambiguous  : testable by >=2 frameworks but supported by only 1 -> genuine disagreement", stage=STAGE)
log_info("  Not supported         : present in a resource but not passing any framework's rule", stage=STAGE)

# ---- sample / patient recurrence (union across supporting frameworks) ----
sup_samples <- function(a,b,c) {
  f <- function(x) if (is.na(x) || !nzchar(x)) character(0) else unlist(strsplit(x,";"))
  mapply(function(x,y,zz) paste(sort(unique(c(f(x),f(y),f(zz)))), collapse=";"),
         a, b, c, USE.NAMES=FALSE)
}
conc$supported_samples <- sup_samples(conc$liana_supported_samples,
  if (cellchat_ok) conc$cellchat_supported_samples else rep(NA_character_, nrow(conc)),
  if (cpdb_ok) conc$cpdb_supported_samples else rep(NA_character_, nrow(conc)))
conc$samples_supported <- vapply(conc$supported_samples,
  function(x) if (!nzchar(x)) 0L else length(unlist(strsplit(x,";"))), integer(1), USE.NAMES=FALSE)
conc$samples_evaluable <- 4L      # M19: every population and pair is evaluable in all 4 samples
conc$sample_recurrence_fraction <- round(conc$samples_supported / conc$samples_evaluable, 3)
conc$patients_supported <- conc$samples_supported          # one sample = one patient
conc$patients_evaluable <- 4L
conc$patient_recurrence_fraction <- conc$sample_recurrence_fraction
conc$direction <- ifelse(conc$sender=="MPNST-Tumor" & conc$receiver!="MPNST-Tumor", "tumor_to_TME",
                  ifelse(conc$receiver=="MPNST-Tumor" & conc$sender!="MPNST-Tumor", "TME_to_tumor", "TME_to_TME"))

conc <- conc %>% arrange(desc(n_LR_methods_supported), desc(samples_supported), liana_best_rank)
write.table(conc, file.path(out_dir,"CCC_CONCORDANCE.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(conc, file.path(tab_dir,"CCC_CONCORDANCE.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

cls <- as.data.frame(table(concordance_class=conc$concordance_class))
log_info("Concordance class counts:", stage=STAGE)
for (r in seq_len(nrow(cls))) log_info(sprintf("  %-50s %6d", cls$concordance_class[r], cls$Freq[r]), stage=STAGE)
sup <- conc %>% filter(n_LR_methods_supported >= 1)
log_info(sprintf("Total distinct interaction keys: %d | with >=1 framework support: %d", nrow(conc), nrow(sup)), stage=STAGE)
log_info(sprintf("Resource overlap: testable by 3 frameworks %d | by 2 %d | by 1 %d",
  sum(conc$n_methods_testable==3), sum(conc$n_methods_testable==2), sum(conc$n_methods_testable==1)), stage=STAGE)

# directional tables
for (dd in c("tumor_to_TME","TME_to_tumor","TME_to_TME")) {
  sub <- conc %>% filter(direction==dd, n_LR_methods_supported >= 1)
  fn <- switch(dd, tumor_to_TME="tumor_to_immune_interactions.tsv",
                   TME_to_tumor="immune_to_tumor_interactions.tsv",
                   TME_to_TME="tme_to_tme_interactions.tsv")
  write.table(sub, file.path(tab_dir, fn), sep="\t", row.names=FALSE, quote=FALSE)
  log_info(sprintf("  %s: %d supported keys -> %s", dd, nrow(sub), fn), stage=STAGE)
}
strom <- conc %>% filter(n_LR_methods_supported >= 1,
  (sender=="MPNST-Tumor" & receiver %in% c("Fibroblast","Pericyte-VSMC","Endothelial")) |
  (receiver=="MPNST-Tumor" & sender %in% c("Fibroblast","Pericyte-VSMC","Endothelial")))
write.table(strom, file.path(tab_dir,"tumor_to_stromal_interactions.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
log_info(sprintf("  tumour<->stromal/endothelial supported keys: %d", nrow(strom)), stage=STAGE)

# recurrence tables
rec_s <- sup %>% group_by(sender, receiver, ligand_n, receptor_n, concordance_class) %>%
  summarise(samples_supported=max(samples_supported), supported_samples=first(supported_samples),
            sample_recurrence_fraction=max(sample_recurrence_fraction), .groups="drop") %>%
  arrange(desc(samples_supported))
write.table(rec_s, file.path(tab_dir,"sample_recurrence.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(rec_s, file.path(tab_dir,"patient_recurrence.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
log_info(sprintf("Interactions supported in all 4 samples: %d | >=3: %d | >=2: %d | 1 only: %d",
  sum(sup$samples_supported==4), sum(sup$samples_supported>=3),
  sum(sup$samples_supported>=2), sum(sup$samples_supported==1)), stage=STAGE)

# method overlap counts (for the UpSet-equivalent figure)
ov <- sup %>% mutate(combo=paste0(ifelse(liana_sup,"LIANA",""),
                                  ifelse(cellchat_sup, ifelse(liana_sup,"+CellChat","CellChat"), ""),
                                  ifelse(cpdb_sup, ifelse(liana_sup|cellchat_sup,"+CellPhoneDB","CellPhoneDB"), ""))) %>%
  count(combo, name="n_interactions") %>% arrange(desc(n_interactions))
write.table(ov, file.path(out_dir,"method_overlap_counts.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
log_info("Method overlap:", stage=STAGE)
for (r in seq_len(nrow(ov))) log_info(sprintf("  %-32s %6d", ov$combo[r], ov$n_interactions[r]), stage=STAGE)

rec <- list(milestone="M21", phase="phase3", timestamp=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),
  script="scripts/phase3/concordance/build_concordance.R",
  frameworks=list(n_contributing=n_frameworks, liana=TRUE, cellchat=cellchat_ok, cellphonedb=cpdb_ok),
  support_rules=list(LIANA=sprintf("%s <= %.2f in >=1 sample", rank_col, LIANA_CUT),
    CellChat="present in subsetCommunication() output, which returns only significant links (p < 0.05)",
    CellPhoneDB="permutation p < 0.05 AND mean expression > 0 in >=1 sample"),
  no_score_averaging="Raw scores are never averaged across tools. Each framework contributes a boolean support flag only.",
  independence_caveat="LIANA's method set includes a CellPhoneDB-style score, so LIANA and the standalone CellPhoneDB run are not fully independent.",
  resource_awareness="Concordance is scored against what each framework could TEST. Every row carries testable_* alongside supported_*, so 'single-method' distinguishes genuine disagreement from resource non-overlap.",
  nichenet_handling="NicheNet is NOT counted as a fourth LR framework; it is orthogonal receiver-response evidence joined in M23/M25.",
  results=list(total_keys=nrow(conc), keys_with_support=nrow(sup),
    class_counts=as.list(setNames(cls$Freq, cls$concordance_class)),
    testable_by_3=sum(conc$n_methods_testable==3), testable_by_2=sum(conc$n_methods_testable==2),
    testable_by_1=sum(conc$n_methods_testable==1),
    supported_in_4_samples=sum(sup$samples_supported==4),
    supported_in_ge3_samples=sum(sup$samples_supported>=3),
    tumor_to_TME=sum(sup$direction=="tumor_to_TME"), TME_to_tumor=sum(sup$direction=="TME_to_tumor")),
  warnings=warns, slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM"),
  git_commit=tryCatch(trimws(system("git rev-parse HEAD",intern=TRUE)), error=function(e) NA_character_),
  elapsed_seconds=as.numeric(difftime(Sys.time(),t0,units="secs")), session_info=capture.output(sessionInfo()))
write_json(rec, file.path(out_dir,"m21_concordance_record.json"), auto_unbox=TRUE, pretty=TRUE, null="null", digits=NA)
prov <- record_provenance("phase3_m21_concordance",
  inputs=list(liana=lp, cellchat=if (cellchat_ok) paste(cc_files, collapse=";") else "NOT_AVAILABLE", cellphonedb=if (cpdb_ok) pp else "NOT_AVAILABLE"),
  outputs=list(concordance=file.path(out_dir,"CCC_CONCORDANCE.tsv"),
    overlap=file.path(out_dir,"method_overlap_counts.tsv"),
    record=file.path(out_dir,"m21_concordance_record.json")),
  parameters=list(milestone="M21", n_frameworks=n_frameworks,
    slurm_job_id=Sys.getenv("SLURM_JOB_ID","NOT_IN_SLURM")), dataset="combined")
save_provenance_json(prov, file.path(out_dir,"prov_m21_concordance.json"))
log_info(sprintf("M21 COMPLETE in %.1f min. Warnings: %d", as.numeric(difftime(Sys.time(),t0,units="mins")), length(warns)), stage=STAGE)
