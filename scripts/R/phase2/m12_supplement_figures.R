# scripts/R/phase2/m12_supplement_figures.R
# M12 supplement: the individually-named, full-page pre/post Harmony metadata panels
# requested in the Phase 2 completion authorization. No new analysis is performed;
# these re-render the M12 comparison at single-panel size for readability.
options(stringsAsFactors = FALSE); options(future.globals.maxSize = +Inf)
suppressPackageStartupMessages({library(Seurat); library(SeuratObject); library(ggplot2)
  library(RColorBrewer); library(jsonlite); library(digest)})
source("scripts/R/logging_utils.R"); source("scripts/R/provenance_utils.R"); setup_strict_logging()
STAGE <- "phase2_m12_supplement"
IN  <- "results/phase2/harmony/phase2_harmony_integrated.rds"
UMAP_POST <- "results/phase2/harmony/evaluation/umap_harmony_m12_embedding.rds"
FIG <- "results/phase2/figures/M12"
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
SEED <- 42

log_info("Loading M11 object (read-only) ...", stage = STAGE)
in_md5 <- digest(IN, file = TRUE, algo = "md5")
obj <- readRDS(IN)
stopifnot(ncol(obj) == 19716L)
pre_xy  <- Embeddings(obj, "umap_preintegration")
post_xy <- readRDS(UMAP_POST)
stopifnot(identical(rownames(pre_xy), rownames(post_xy)))
md <- obj@meta.data
rm(obj); invisible(gc(verbose = FALSE))

grp <- as.character(md$sample_id)
lv  <- sort(unique(grp))
pal <- setNames(RColorBrewer::brewer.pal(max(3, length(lv)), "Set1")[seq_along(lv)], lv)
set.seed(SEED); ord <- sample.int(length(grp))

th <- theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 11, colour = "grey30"),
        legend.position = "right", legend.key.size = unit(0.7, "cm"))

rows <- list()
panel <- function(xy, title, sub, fname, role) {
  d <- data.frame(x = xy[ord, 1], y = xy[ord, 2], sample_id = grp[ord])
  p <- ggplot(d, aes(x, y, colour = sample_id)) +
    geom_point(size = 0.30, alpha = 0.8, stroke = 0) +
    scale_colour_manual(values = pal) +
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    labs(title = title, subtitle = sub, x = "UMAP 1", y = "UMAP 2", colour = role) + th
  for (ext in c("pdf","png")) {
    fp <- file.path(FIG, paste0(fname, ".", ext))
    ggsave(fp, p, width = 9, height = 7.5, units = "in", dpi = 220, device = ext)
    rows[[length(rows)+1]] <<- data.frame(
      figure_path = fp, phase = "phase2", milestone = "M12", dataset = "combined",
      processing_stage = "pre_post_harmony_evaluation",
      analysis_method = "UMAP_scatter_single_panel",
      parameters = "dims=1:30;seed=42;n.neighbors=30;min.dist=0.3;metric=cosine;point=0.30",
      input_object_checksum = in_md5,
      generating_script = "scripts/R/phase2/m12_supplement_figures.R",
      stringsAsFactors = FALSE)
  }
  log_info(sprintf("  wrote %s(.pdf/.png)", fname), stage = STAGE)
}

# sample_id is simultaneously dataset, sample, patient and the Harmony grouping
# variable (verified cell-for-cell identical to orig.ident in M10). The requested
# per-role panels are therefore the same comparison rendered under each role name,
# and each figure says so on its face. No 'condition' panel: no such field exists.
note <- "sample_id is simultaneously the dataset, the sample, the patient and the Harmony grouping variable - all one 4-level field."
PRE  <- "PRE-Harmony (pca 1:30 -> umap_preintegration, frozen Phase 1 baseline)"
POST <- "POST-Harmony (postint_harmony 1:30 -> umap_harmony_m12)"

panel(pre_xy,  PRE,  note, "M12_01_pre_harmony_by_dataset",  "dataset")
panel(post_xy, POST, note, "M12_02_post_harmony_by_dataset", "dataset")
panel(pre_xy,  PRE,  note, "M12_03_pre_harmony_by_sample",   "sample")
panel(post_xy, POST, note, "M12_04_post_harmony_by_sample",  "sample")
panel(pre_xy,  PRE,  paste(note, "One sample = one patient; no separate patient field exists."),
      "M12_05_pre_harmony_by_patient",  "patient")
panel(post_xy, POST, paste(note, "One sample = one patient; no separate patient field exists."),
      "M12_06_post_harmony_by_patient", "patient")
panel(pre_xy,  PRE,  note, "M12_09_pre_harmony_by_harmony_variable",  "Harmony group")
panel(post_xy, POST, note, "M12_10_post_harmony_by_harmony_variable", "Harmony group")

writeLines(c(
  "M12_07_pre_harmony_by_condition and M12_08_post_harmony_by_condition were NOT generated.",
  "",
  "No biological-condition field exists anywhere in the Phase 1 handoff object. The M10 audit",
  "of all 93 metadata columns found no treatment arm, disease stage, anatomical site, primary/",
  "metastatic status, NF1 status, tumour grade, age or sex. Fabricating a condition variable to",
  "satisfy a filename would be scientifically invalid, so these two figures are deliberately absent.",
  "",
  "See reports/phase2/HARMONY_VARIABLE_DECISION.md section 2 and results/phase2/handoff/handoff_metadata_inventory.tsv."
), file.path(FIG, "M12_07_08_condition_figures_NOT_APPLICABLE.txt"))

write.table(do.call(rbind, rows), file.path(FIG, "figure_index_m12_supplement.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
prov <- record_provenance("phase2_m12_supplement",
  inputs = list(m11_rds = IN, umap_post = UMAP_POST),
  outputs = list(figure_index = file.path(FIG, "figure_index_m12_supplement.tsv")),
  parameters = list(milestone = "M12_supplement", seed = SEED,
                    slurm_job_id = Sys.getenv("SLURM_JOB_ID", "NOT_IN_SLURM")),
  dataset = "combined")
save_provenance_json(prov, file.path(FIG, "prov_m12_supplement.json"))
log_info("M12 supplement figures complete.", stage = STAGE)
