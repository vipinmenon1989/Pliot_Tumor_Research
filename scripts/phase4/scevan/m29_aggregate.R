#!/usr/bin/env Rscript
# =============================================================================
# Phase 4 · M29 aggregation — collect the four per-patient SCEVAN runs
#
# Produces:
#   SCEVAN_CELL_CLASSIFICATION_RAW.tsv   per-cell calls, both runs, all patients
#   SCEVAN_CLONES.tsv                    clone/subclone inventory per patient
#   SCEVAN_CNV_SUMMARY.tsv               broad CN segments per clone
#   SCEVAN_ARM_LEVEL_EVENTS.tsv          cytoband-level gains/losses (oncoHeat)
#   SCEVAN_RUN_ACCOUNTING.tsv            what each run did, including failures
#   figures: UMAP overlays + clone composition; SCEVAN's own CNA heatmaps indexed
#
# Wording rule (§50): CNA patterns shared BETWEEN patients are described as
# "shared" or "recurrent", never as "the same clone".
# =============================================================================
suppressPackageStartupMessages({
  library(jsonlite); library(ggplot2); library(dplyr); library(tidyr)
})
source("scripts/phase4/utils/phase4_plot_utils.R")

SC  <- "results/phase4/scevan/by_sample"
TAB <- "results/phase4/tables"
FIG <- "results/phase4/figures"
MAL <- "results/phase4/malignancy"
SAMPLES <- c("MPNST_1","MPNST_2","MPNST_3","MPNST_4")
dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
log_ <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
sec  <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

md <- read.delim(file.path(MAL, "phase4_cell_metadata.tsv.gz"), stringsAsFactors = FALSE)
rownames(md) <- md$cell_id
UX <- "postint_umap_harmony_1"; UY <- "postint_umap_harmony_2"

# --- 1. run accounting -----------------------------------------------------
sec("1. Run accounting")
acc <- list(); cls_all <- list()
for (s in SAMPLES) {
  jf <- file.path(SC, s, paste0("m29_", s, "_summary.json"))
  if (!file.exists(jf)) { log_("MISSING summary for ", s); next }
  j <- fromJSON(jf, simplifyVector = TRUE)
  for (tag in c("primary","sensitivity")) {
    r <- j[[tag]]
    acc[[length(acc)+1]] <- data.frame(
      sample_id = s, run = tag, status = r$status,
      elapsed_min = if (is.null(r$elapsed_min)) NA_real_ else r$elapsed_min,
      cells_in = if (is.null(r$cells_in)) NA_integer_ else r$cells_in,
      cells_returned = if (is.null(r$cells_returned)) NA_integer_ else r$cells_returned,
      n_tumor = if (is.null(r$class_table$tumor)) 0L else r$class_table$tumor,
      n_normal = if (is.null(r$class_table$normal)) 0L else r$class_table$normal,
      n_filtered = if (is.null(r$class_table$filtered)) 0L else r$class_table$filtered,
      confident_normal = if (is.null(r$confident_normal)) NA_integer_ else r$confident_normal,
      n_subclones = if (is.null(r$n_subclones)) NA_integer_ else r$n_subclones,
      norm_cell_n = if (is.null(r$norm_cell_n)) NA_integer_ else r$norm_cell_n,
      n_files = if (is.null(r$files_written)) 0L else length(r$files_written),
      stringsAsFactors = FALSE)
  }
  if (!is.null(j$concordance))
    log_(sprintf("%s: primary/sensitivity agreement %.4f over %d cells",
                 s, j$concordance$agreement, j$concordance$cells_common))
  f <- file.path(SC, s, "primary",
                 sprintf("scevan_classification_%s_primary.rds", s))
  if (file.exists(f)) cls_all[[s]] <- readRDS(f)
}
acc <- bind_rows(acc)
write.table(acc, file.path(TAB, "SCEVAN_RUN_ACCOUNTING.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(acc, row.names = FALSE)

cls <- bind_rows(cls_all)
log_("total cells classified (primary): ", nrow(cls))
write.table(cls, file.path(TAB, "SCEVAN_CELL_CLASSIFICATION_RAW.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- 2. clone inventory ----------------------------------------------------
sec("2. Clone / subclone inventory")
cls$clone <- ifelse(is.na(cls$scevan_subclone), NA_character_,
                    paste0(cls$sample_id, "_clone", cls$scevan_subclone))
mdx <- md[match(cls$cell_id, md$cell_id), ]
cls$annotation_ccc <- mdx$annotation_ccc
cls$primary_cluster <- mdx$postint_harmony_primary_cluster
clone_tab <- cls |> filter(!is.na(clone)) |>
  group_by(sample_id, clone, scevan_subclone) |>
  summarise(n_cells = n(),
            frac_of_sample_tumor = n() / sum(cls$sample_id == sample_id[1] &
                                             cls$scevan_call == "tumor"),
            dominant_annotation_ccc = names(sort(table(annotation_ccc), decreasing = TRUE))[1],
            n_annotations = length(unique(annotation_ccc)),
            annotation_breakdown = paste(sprintf("%s=%d",
              names(sort(table(annotation_ccc), decreasing = TRUE)),
              sort(table(annotation_ccc), decreasing = TRUE)), collapse = "; "),
            dominant_cluster = names(sort(table(primary_cluster), decreasing = TRUE))[1],
            .groups = "drop") |>
  arrange(sample_id, desc(n_cells))
write.table(clone_tab, file.path(TAB, "SCEVAN_CLONES.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(as.data.frame(clone_tab), row.names = FALSE)

# --- 3. broad CN segments per clone ---------------------------------------
sec("3. Broad copy-number segments (clonal and subclonal .seg files)")
segs <- list()
for (s in SAMPLES) {
  od <- file.path(SC, s, "primary", "output")
  if (!dir.exists(od)) next
  for (f in list.files(od, pattern = "_CN\\.seg$", full.names = TRUE)) {
    b <- basename(f)
    lvl <- if (grepl("_Clonal_CN\\.seg$", b)) "clonal"
           else sub(".*_subclone([0-9]+)_CN\\.seg$", "subclone\\1", b)
    d <- read.delim(f, stringsAsFactors = FALSE)
    d$sample_id <- s; d$level <- lvl
    d$width_mb <- (d$End - d$Pos) / 1e6
    segs[[length(segs)+1]] <- d
  }
}
if (length(segs)) {
  sg <- bind_rows(segs)
  # Only broad events are reported. Single-gene CNV calls from scRNA-derived
  # inference are not asserted (§51, §65B).
  sg$event <- ifelse(sg$CN > 2, "gain", ifelse(sg$CN < 2, "loss", "neutral"))
  sg$broad <- sg$width_mb >= 10
  write.table(sg, file.path(TAB, "SCEVAN_CNV_SUMMARY.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  log_("segments: ", nrow(sg), " | broad (>=10 Mb): ", sum(sg$broad))
  print(sg |> filter(broad, event != "neutral") |>
        count(sample_id, level, event) |> as.data.frame(), row.names = FALSE)
  # recurrence across patients, at chromosome-arm resolution
  rec <- sg |> filter(broad, event != "neutral", level == "clonal") |>
    count(Chr, event, sample_id) |> count(Chr, event, name = "n_patients") |>
    arrange(desc(n_patients), Chr)
  write.table(rec, file.path(TAB, "SCEVAN_RECURRENT_BROAD_EVENTS.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  log_("recurrent broad clonal events (chromosome x direction x n patients):")
  print(head(as.data.frame(rec), 25), row.names = FALSE)
} else log_("no .seg files found")

# --- 4. cytoband-level events from oncoHeat -------------------------------
sec("4. Cytoband-level events (oncoHeat)")
oh <- list()
for (s in SAMPLES) {
  f <- file.path(SC, s, "primary", "output", paste0(s, "_primary", "PlotOncoHeat.RData"))
  if (!file.exists(f)) { log_("no oncoHeat for ", s); next }
  e <- new.env(); load(f, envir = e)
  if (!"oncoHeat" %in% ls(e)) next
  o <- as.data.frame(get("oncoHeat", e))
  o$clone <- rownames(o)
  ol <- pivot_longer(o, cols = -clone, names_to = "cytoband", values_to = "code")
  ol$sample_id <- s
  ol$event <- c("-2" = "loss", "0" = "neutral", "2" = "gain")[as.character(ol$code)]
  oh[[s]] <- ol
}
if (length(oh)) {
  ohl <- bind_rows(oh) |> filter(!is.na(event))
  write.table(ohl, file.path(TAB, "SCEVAN_ARM_LEVEL_EVENTS.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  log_("oncoHeat rows: ", nrow(ohl))
  print(ohl |> filter(event != "neutral") |> count(sample_id, event) |> as.data.frame(),
        row.names = FALSE)
} else log_("no oncoHeat objects found")

# --- 5. figures ------------------------------------------------------------
sec("5. M29 figures")
md$scevan_call <- cls$scevan_call[match(md$cell_id, cls$cell_id)]
md$scevan_call[is.na(md$scevan_call) | md$scevan_call == "filtered"] <- "not-assessed"
md$scevan_call <- c(tumor = "malignant", normal = "non-malignant",
                    `not-assessed` = "not-assessed")[md$scevan_call]
md$scevan_call <- factor(md$scevan_call,
                         levels = c("malignant","non-malignant","not-assessed"))
md$clone <- cls$clone[match(md$cell_id, cls$cell_id)]

# 29-01 per-patient UMAP overlay of the SCEVAN call
f1 <- umap_layer(md, UX, UY, "scevan_call", P4_SCEVAN,
  "SCEVAN malignancy call per patient, on the frozen Phase 2 Harmony UMAP",
  paste("SCEVAN was run SEPARATELY for each patient, because sample_id = patient = dataset.",
        "\nThe UMAP is for visualisation only; classification comes from raw counts and inferred copy number."),
  P4_CAPTION_CNV, highlight = "malignant", legend_title = "SCEVAN call") +
  facet_wrap(~ sample_id, nrow = 1)
save_fig(f1, file.path(FIG, "29_01_scevan_call_umap_by_patient.pdf"), 16.5, 5.4)

# 29-02 clone overlay (clones are patient-scoped by construction)
mc <- md |> filter(!is.na(clone))
if (nrow(mc)) {
  ncl <- length(unique(mc$clone))
  pal <- setNames(colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(ncl),
                  sort(unique(mc$clone)))
  f2 <- ggplot(md, aes(x = .data[[UX]], y = .data[[UY]])) +
    geom_point(colour = "grey88", size = 0.28, stroke = 0) +
    geom_point(data = mc, aes(colour = clone), size = 0.42, alpha = 0.85, stroke = 0) +
    scale_colour_manual(values = pal, name = "SCEVAN clone") +
    guides(colour = guide_legend(override.aes = list(size = 2.6, alpha = 1), ncol = 1)) +
    facet_wrap(~ sample_id, nrow = 1) + coord_equal() + p4_theme() +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    labs(title = "SCEVAN subclones, per patient",
         subtitle = paste("Clone labels are PATIENT-SCOPED. A clone in one patient is not the same",
                          "clone as one in another;\nsimilar CNA architecture across patients is a",
                          "SHARED or RECURRENT pattern, not a shared clone."),
         x = "UMAP 1 (Phase 2 Harmony)", y = "UMAP 2 (Phase 2 Harmony)",
         caption = P4_CAPTION_CNV)
  save_fig(f2, file.path(FIG, "29_02_scevan_clones_umap_by_patient.pdf"), 17.5, 5.8)

  # 29-03 clone composition vs Phase 2 annotation
  d3 <- mc |> count(sample_id, clone, annotation_ccc) |>
    group_by(clone) |> mutate(frac = n/sum(n), tot = sum(n)) |> ungroup()
  f3 <- ggplot(d3, aes(y = clone, x = frac, fill = annotation_ccc)) +
    geom_col(width = 0.78) +
    scale_x_continuous(labels = percent_format(accuracy = 1), expand = expansion(0)) +
    scale_fill_brewer(palette = "Paired", name = "Phase 2 annotation_ccc") +
    facet_grid(sample_id ~ ., scales = "free_y", space = "free_y") +
    labs(title = "What each SCEVAN clone is made of, in Phase 2 terms",
         subtitle = paste("A clone drawn largely from Fibroblast is the signature of a",
                          "malignant Mes-NC-like population\nthat marker-based annotation read as stroma."),
         x = "fraction of clone", y = NULL, caption = P4_CAPTION_CNV) + p4_theme()
  save_fig(f3, file.path(FIG, "29_03_clone_composition_by_annotation.pdf"), 11.5, 7.6)
} else log_("no clones reported; clone figures skipped")

# 29-04 index SCEVAN's own native CNA heatmaps rather than redrawing them
sec("6. SCEVAN native figure index")
nat <- list()
for (s in SAMPLES) for (tag in c("primary","sensitivity")) {
  od <- file.path(SC, s, tag, "output")
  if (!dir.exists(od)) next
  for (p in list.files(od, pattern = "\\.png$", full.names = TRUE)) {
    b <- basename(p)
    kind <- if (grepl("onlytumorheatmap", b)) "CNA heatmap, tumour cells only"
      else if (grepl("heatmap_subclones", b)) "CNA heatmap with subclone annotation track"
      else if (grepl("heatmap", b))           "CNA heatmap, all cells, with malignant/normal track"
      else if (grepl("OncoHeat", b))          "cytoband-level onco-heatmap by subclone"
      else if (grepl("consensus", b))         "consensus clonal CN profile"
      else if (grepl("ClonalCNProfile", b))   "clonal CN profile along the genome"
      else if (grepl("umap_CNA", b))          "UMAP of CNA space, coloured by subclone"
      else if (grepl("umap_scRNA", b))        "UMAP of expression space, tumour vs normal"
      else "SCEVAN output"
    nat[[length(nat)+1]] <- data.frame(sample_id = s, run = tag, file = p,
                                       figure_type = kind, stringsAsFactors = FALSE)
  }
}
nat <- bind_rows(nat)
write.table(nat, file.path(TAB, "SCEVAN_NATIVE_FIGURE_INDEX.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
log_("SCEVAN native figures indexed: ", nrow(nat))
print(nat |> count(run, figure_type) |> as.data.frame(), row.names = FALSE)

sec("M29 AGGREGATION COMPLETE")
