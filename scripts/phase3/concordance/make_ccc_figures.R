# scripts/phase3/concordance/make_ccc_figures.R
#
# Phase 3 / Milestones M20-M21 - CCC overview and concordance figure suite.
# Reads only finished method outputs and the M21 concordance table; performs no inference.

options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(RColorBrewer); library(viridis); library(igraph)
  library(jsonlite); library(digest)})
source("scripts/R/logging_utils.R"); source("scripts/R/provenance_utils.R"); setup_strict_logging()
STAGE <- "phase3_ccc_figures"

args <- commandArgs(trailingOnly=TRUE)
by_dir  <- "results/phase3/ccc/by_method"
conc_p  <- "results/phase3/ccc/concordance/CCC_CONCORDANCE.tsv"
fig20   <- "results/phase3/figures/M20"
fig21   <- "results/phase3/figures/M21"
out_dir <- "results/phase3/ccc/concordance"
i <- 1
while (i <= length(args)) { a <- args[i]
  if (a=="--by-dir"){by_dir<-args[i+1];i<-i+2} else if (a=="--concordance"){conc_p<-args[i+1];i<-i+2}
  else if (a=="--fig20"){fig20<-args[i+1];i<-i+2} else if (a=="--fig21"){fig21<-args[i+1];i<-i+2}
  else if (a=="--out-dir"){out_dir<-args[i+1];i<-i+2} else stop(sprintf("Unknown argument: %s", a)) }
t0 <- Sys.time(); warns <- character(0)
wrn <- function(m){warns<<-c(warns,m); log_warn(m,stage=STAGE)}
for (d in c(fig20,fig21,out_dir)) dir.create(d,recursive=TRUE,showWarnings=FALSE)

log_info("=========== CCC FIGURE SUITE (M20 / M21) ===========", stage=STAGE)
conc <- read.delim(conc_p, sep="\t", check.names=FALSE)
in_md5 <- digest(conc_p, file=TRUE, algo="md5")
sup <- conc %>% filter(n_LR_methods_supported >= 1)
log_info(sprintf("Concordance table: %d keys | %d with support", nrow(conc), nrow(sup)), stage=STAGE)

TUM <- "MPNST-Tumor"
IMMUNE <- c("Macrophage","Monocyte","Dendritic","Plasmacytoid-DC","CD8-T","CD4-T","NK","T-cell-other","B-cell","Plasma-cell")
STROM <- c("Fibroblast","Pericyte-VSMC"); ENDO <- "Endothelial"
pop_order <- c(TUM, IMMUNE, STROM, ENDO)
compartment_of <- function(x) ifelse(x==TUM,"Tumour", ifelse(x %in% IMMUNE,"Immune",
                             ifelse(x %in% STROM,"Stromal", ifelse(x==ENDO,"Endothelial","Other"))))
th <- theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),
      plot.title=element_text(face="bold"), plot.subtitle=element_text(size=9,colour="grey30"))
CAVEAT <- "Inferred communication potential from expression; not evidence of physical adjacency or direct signalling."

frows <- list()
sf <- function(p,dir,name,w,h,analysis,sender,receiver,method,params){
  for (ext in c("pdf","png")) { fp <- file.path(dir, paste0(name,".",ext))
    ggsave(fp,p,width=w,height=h,units="in",dpi=200,device=ext,limitsize=FALSE)
    frows[[length(frows)+1]] <<- data.frame(figure_path=fp, phase="phase3",
      milestone=if (identical(dir,fig20)) "M20" else "M21", analysis=analysis,
      sender=sender, receiver=receiver, method=method, input=conc_p, input_checksum=in_md5,
      script="scripts/phase3/concordance/make_ccc_figures.R", parameters=params,
      stringsAsFactors=FALSE) }
  log_info(sprintf("  wrote %s/%s", basename(dir), name), stage=STAGE) }

ord <- function(x) factor(x, levels=intersect(pop_order, unique(x)))

# ---- Figure 1: global sender x receiver interaction heatmap ----
gm <- sup %>% count(sender, receiver, name="n_interactions")
sf(ggplot(gm, aes(ord(receiver), ord(sender), fill=n_interactions))+geom_tile(colour="white")+
   geom_text(aes(label=n_interactions), size=2.7, colour="grey15")+
   scale_fill_viridis_c(option="mako", direction=-1, name="supported\ninteractions")+
   labs(title="Global sender x receiver communication potential",
        subtitle=paste("Supported ligand-receptor interactions per directed pair (>=1 framework).", CAVEAT),
        x="receiver", y="sender")+th+
   theme(axis.text.x=element_text(angle=45,hjust=1)),
   fig20, "M20_01_global_sender_receiver_heatmap", 11, 9,
   "global_network","all","all","LIANA+CellChat+CellPhoneDB","supported >=1 framework")

# ---- Figure 2/3: MPNST-Tumor outgoing and incoming ----
outg <- sup %>% filter(sender==TUM, receiver!=TUM) %>% count(receiver, concordance_class, name="n")
sf(ggplot(outg, aes(reorder(receiver, n, sum), n, fill=concordance_class))+geom_col()+coord_flip()+
   scale_fill_brewer(palette="Set2", name="concordance")+
   labs(title="MPNST-Tumor OUTGOING signalling",
        subtitle=paste("Supported interactions sent from MPNST-Tumor to each TME population.", CAVEAT),
        x=NULL, y="supported interactions")+th,
   fig20, "M20_02_mpnst_outgoing_signaling", 10, 6.5,
   "outgoing", TUM, "TME", "all", "supported >=1 framework")
inc <- sup %>% filter(receiver==TUM, sender!=TUM) %>% count(sender, concordance_class, name="n")
sf(ggplot(inc, aes(reorder(sender, n, sum), n, fill=concordance_class))+geom_col()+coord_flip()+
   scale_fill_brewer(palette="Set2", name="concordance")+
   labs(title="MPNST-Tumor INCOMING signalling",
        subtitle=paste("Supported interactions received by MPNST-Tumor from each TME population.", CAVEAT),
        x=NULL, y="supported interactions")+th,
   fig20, "M20_03_mpnst_incoming_signaling", 10, 6.5,
   "incoming", "TME", TUM, "all", "supported >=1 framework")

# ---- Figure 4/5: tumour->immune and immune->tumour LR heatmaps (top pairs) ----
lr_heat <- function(df, title, xlab_pop, fname, dir) {
  if (!nrow(df)) { wrn(sprintf("No rows for %s", fname)); return(invisible(NULL)) }
  df$lr <- paste0(df$ligand_n, " -> ", df$receptor_n)
  top <- df %>% group_by(lr) %>%
    summarise(score=sum(n_LR_methods_supported) + sum(samples_supported), .groups="drop") %>%
    arrange(desc(score)) %>% head(40) %>% pull(lr)
  d <- df %>% filter(lr %in% top)
  p <- ggplot(d, aes(ord(if (xlab_pop=="receiver") receiver else sender), factor(lr, levels=rev(top)),
                     fill=n_LR_methods_supported))+
    geom_tile(colour="white")+
    geom_point(aes(size=samples_supported), colour="grey15", alpha=0.55)+
    scale_fill_viridis_c(option="rocket", direction=-1, name="frameworks\nsupporting", breaks=1:3)+
    scale_size_continuous(range=c(0.4,3.2), name="samples\nsupporting", breaks=1:4)+
    labs(title=title, subtitle=paste("Top 40 ligand->receptor pairs by combined framework and sample support.", CAVEAT),
         x=xlab_pop, y=NULL)+th+
    theme(axis.text.y=element_text(size=6.6), axis.text.x=element_text(angle=45,hjust=1))
  sf(p, dir, fname, 11, 11, "LR_heatmap",
     if (xlab_pop=="receiver") TUM else "TME", if (xlab_pop=="receiver") "TME" else TUM,
     "all", "top40 by framework+sample support")
}
lr_heat(sup %>% filter(sender==TUM, receiver!=TUM), "Tumour -> immune / TME ligand-receptor map",
        "receiver", "M20_04_tumor_to_immune_LR_heatmap", fig20)
lr_heat(sup %>% filter(receiver==TUM, sender!=TUM), "Immune / TME -> tumour ligand-receptor map",
        "sender", "M20_05_immune_to_tumor_LR_heatmap", fig20)

# ---- Figure 6: method concordance heatmap for the top tumour-centric interactions ----
tc <- sup %>% filter(sender==TUM | receiver==TUM) %>%
  mutate(lab=paste0(sender," -> ",receiver," | ",ligand_n,"->",receptor_n)) %>%
  arrange(desc(n_LR_methods_supported), desc(samples_supported)) %>% head(45)
if (nrow(tc)) {
  mh <- tc %>% select(lab, LIANA=liana_sup, CellChat=cellchat_sup, CellPhoneDB=cpdb_sup) %>%
    pivot_longer(-lab, names_to="method", values_to="supported")
  th2 <- tc %>% select(lab, LIANA=liana_testable, CellChat=cellchat_testable, CellPhoneDB=cpdb_testable) %>%
    pivot_longer(-lab, names_to="method", values_to="testable")
  mh <- mh %>% left_join(th2, by=c("lab","method")) %>%
    mutate(state=ifelse(supported,"supported", ifelse(testable,"testable, not supported","not in this resource")))
  sf(ggplot(mh, aes(method, factor(lab, levels=rev(tc$lab)), fill=state))+geom_tile(colour="white")+
     scale_fill_manual(values=c(`supported`="#1A9850", `testable, not supported`="#FEE08B",
                                `not in this resource`="grey88"), name=NULL)+
     labs(title="Method concordance for the top tumour-centric interactions",
          subtitle=paste("Grey = the framework's LR resource never contained this pair, so it could not test it.", CAVEAT),
          x=NULL, y=NULL)+th+theme(axis.text.y=element_text(size=6.4)),
     fig21, "M21_06_method_concordance", 10, 11, "concordance","tumour-centric","tumour-centric","all","top45")
}

# ---- Figure 7: method overlap (UpSet-equivalent bar chart; no extra package installed) ----
ovp <- file.path(out_dir,"method_overlap_counts.tsv")
if (file.exists(ovp)) {
  ov <- read.delim(ovp, sep="\t")
  ov$combo <- factor(ov$combo, levels=ov$combo[order(ov$n_interactions)])
  sf(ggplot(ov, aes(combo, n_interactions))+geom_col(fill="#3B6EA5")+
     geom_text(aes(label=n_interactions), hjust=-0.1, size=3.2)+coord_flip(clip="off")+
     expand_limits(y=max(ov$n_interactions)*1.2)+
     labs(title="Framework support overlap",
          subtitle="An UpSet-equivalent bar chart. No package was installed purely for decoration.",
          x=NULL, y="supported interactions")+th,
     fig21, "M21_07_method_overlap", 9.5, 6, "overlap","all","all","all","combination counts")
}

# ---- Figure 8: sample recurrence matrix for top tumour-centric interactions ----
tr <- sup %>% filter(sender==TUM | receiver==TUM) %>%
  arrange(desc(samples_supported), desc(n_LR_methods_supported)) %>% head(40) %>%
  mutate(lab=paste0(sender," -> ",receiver," | ",ligand_n,"->",receptor_n))
if (nrow(tr)) {
  samples <- c("MPNST_1","MPNST_2","MPNST_3","MPNST_4")
  rr <- do.call(rbind, lapply(seq_len(nrow(tr)), function(i) {
    ss <- if (is.na(tr$supported_samples[i]) || !nzchar(tr$supported_samples[i])) character(0)
          else unlist(strsplit(tr$supported_samples[i], ";"))
    data.frame(lab=tr$lab[i], sample=samples, supported=samples %in% ss, stringsAsFactors=FALSE) }))
  sf(ggplot(rr, aes(sample, factor(lab, levels=rev(tr$lab)), fill=supported))+geom_tile(colour="white")+
     scale_fill_manual(values=c(`TRUE`="#1A9850",`FALSE`="grey90"),
                       labels=c("not supported","supported"), name=NULL)+
     labs(title="Sample / patient recurrence of top tumour-centric interactions",
          subtitle=paste("Every population and pair was evaluable in all four samples (M19), so a blank cell is a real negative, not a missing test.", CAVEAT),
          x="sample (= patient = dataset)", y=NULL)+th+theme(axis.text.y=element_text(size=6.6)),
     fig21, "M21_08_sample_recurrence", 8.5, 10, "recurrence","tumour-centric","tumour-centric","all","top40")
}

# ---- Figure: concordance class distribution by direction ----
cd <- sup %>% count(direction, concordance_class, name="n")
sf(ggplot(cd, aes(direction, n, fill=concordance_class))+geom_col(position="fill")+
   scale_fill_brewer(palette="Set2", name="concordance")+
   scale_y_continuous(labels=scales::percent)+
   labs(title="Concordance class composition by direction",
        subtitle="Proportion of supported interactions in each concordance class", x=NULL, y=NULL)+th,
   fig21, "M21_09_concordance_by_direction", 9, 5.5, "concordance","all","all","all","fill by class")

# ---- Figure: tumour-centric network graph (igraph layout drawn with ggplot2) ----
# Self-loops (autocrine, sender == receiver) cannot be drawn with geom_curve and are
# summarised separately rather than dropped silently.
auto <- sup %>% filter(sender==TUM, receiver==TUM, n_LR_methods_supported >= 2)
if (nrow(auto)) log_info(sprintf("  autocrine MPNST-Tumor interactions with >=2-framework support: %d (excluded from the network layout, reported in the caption)", nrow(auto)), stage=STAGE)
net <- sup %>% filter(sender != receiver, sender==TUM | receiver==TUM, n_LR_methods_supported >= 2) %>%
  count(sender, receiver, name="weight")
if (nrow(net)) {
  g <- graph_from_data_frame(net, directed=TRUE)
  set.seed(42); L <- layout_in_circle(g)
  vd <- data.frame(name=V(g)$name, x=L[,1], y=L[,2], stringsAsFactors=FALSE)
  vd$compartment <- compartment_of(vd$name)
  ed <- as.data.frame(as_edgelist(g)); names(ed) <- c("from","to"); ed$weight <- E(g)$weight
  ed <- ed %>% left_join(vd %>% select(name,x,y), by=c("from"="name")) %>%
    rename(x0=x, y0=y) %>% left_join(vd %>% select(name,x,y), by=c("to"="name")) %>% rename(x1=x, y1=y)
  sf(ggplot()+
     geom_curve(data=ed, aes(x=x0,y=y0,xend=x1,yend=y1,linewidth=weight, colour=weight),
                curvature=0.18, alpha=0.75, arrow=grid::arrow(length=unit(0.16,"cm"), type="closed"))+
     scale_linewidth_continuous(range=c(0.3,2.6), name="interactions")+
     scale_colour_viridis_c(option="rocket", direction=-1, name="interactions")+
     geom_point(data=vd, aes(x,y,fill=compartment), size=7, shape=21, colour="grey20")+
     scale_fill_manual(values=c(Tumour="#B2182B", Immune="#2166AC", Stromal="#D9A441", Endothelial="#1B7837"))+
     ggrepel::geom_text_repel(data=vd, aes(x,y,label=name), size=3.3, fontface="bold", box.padding=0.6)+
     coord_equal()+theme_void(base_size=12)+
     labs(title="Tumour-centric communication network",
          subtitle=paste(sprintf("Edges: directed pairs with >=2-framework support; width and colour = number of such interactions. Autocrine MPNST-Tumor interactions (n = %d) are excluded from the layout.", nrow(auto)), CAVEAT)),
     fig21, "M21_10_tumour_centric_network", 10, 9, "network", TUM, "TME", "all", ">=2 frameworks")
}
write.table(do.call(rbind, frows), file.path(out_dir,"figure_index_ccc.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
log_info(sprintf("CCC figure suite complete: %d files in %.1f min. Warnings: %d",
  length(frows), as.numeric(difftime(Sys.time(),t0,units="mins")), length(warns)), stage=STAGE)
