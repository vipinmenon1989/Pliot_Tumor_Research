#!/usr/bin/env Rscript
# =============================================================================
# Phase 4 · M32 — targeted Phase 3 CCC sensitivity analysis
#
# The question is NOT "what are the interactions" (Phase 3 answered that). It is:
#   does the Phase 3 biological architecture survive a better tumour definition?
#
# Method: rerun ONLY the tumour-centric part of the CCC space with
# annotation_ccc_refined in place of the Phase 3 annotation_ccc, using the
# SAME framework versions, resources, thresholds, expression basis and
# sample-aware design, then classify each interaction's change.
#
# Change classes (§42) are transparent categories. NO weighted score is invented.
# =============================================================================
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(jsonlite); library(ggplot2)
})
source("scripts/phase4/utils/phase4_plot_utils.R")
set.seed(42)

P3_MASTER <- "results/phase3/ccc/prioritized/MPNST_CCC_MASTER_TABLE.tsv"
P3_CONC   <- "results/phase3/tables/final/CCC_CONCORDANCE.tsv"
P4_DIR    <- "results/phase4/ccc_refinement"
TAB       <- "results/phase4/tables"
FIG       <- "results/phase4/figures"
for (d in c(TAB, FIG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
log_ <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
sec  <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

# --- the prioritized Phase 3 axes M32 must adjudicate (§38) ---------------
FOCUS <- data.frame(
  ligand = c("APP","CD99","ANXA1","THBS1","HLA-F","HLA-F","HLA-E","BAG6","CD58",
             "VEGFA","VEGFA","VEGFA","JAG1","JAG1","JAG1","JAG2","DLL4","FGF2",
             "SLIT2","COL1A1","COL1A2","COL6A2","FN1"),
  receptor = c("CD74","PILRA","FPR1","CD36","LILRB1","LILRB2","KLRC1","NCR3","CD2",
               "KDR","FLT1","NRP1","NOTCH3","NOTCH4","NOTCH2","NOTCH2","NOTCH2","FGFR1",
               "ROBO1","ITGAV_ITGB8","ITGAV_ITGB8","ITGAV_ITGB8","ITGAV_ITGB8"),
  axis = c("tumour->myeloid antigen-presentation","tumour->myeloid inhibitory",
           "tumour->myeloid chemotaxis","tumour->myeloid scavenger",
           "tumour->myeloid/lymphoid inhibitory","tumour->myeloid inhibitory",
           "tumour->NK inhibitory","tumour->NK activating",
           "tumour->lymphoid adhesion/costimulation",
           "tumour->endothelial angiogenic","tumour->endothelial angiogenic",
           "tumour->endothelial angiogenic","tumour->pericyte Notch",
           "tumour->endothelial Notch","TME->tumour Notch","TME->tumour Notch",
           "TME->tumour Notch","tumour->stromal growth factor",
           "tumour->stromal guidance","ECM->tumour integrin","ECM->tumour integrin",
           "ECM->tumour integrin","ECM->tumour integrin"),
  stringsAsFactors = FALSE)

# --- 1. Phase 3 baseline ---------------------------------------------------
sec("1. Phase 3 baseline")
p3m <- read.delim(P3_MASTER, stringsAsFactors = FALSE)
log_("Phase 3 master table (prioritized): ", nrow(p3m), " rows x ", ncol(p3m), " cols")
p3 <- read.delim(P3_CONC, stringsAsFactors = FALSE)
log_("Phase 3 concordance table: ", nrow(p3), " rows x ", ncol(p3), " cols")

# --- 2. Phase 4 concordance output ----------------------------------------
sec("2. Phase 4 refined concordance")
# Produced by the SAME unmodified scripts/phase3/concordance/build_concordance.R,
# so the schema is identical and the comparison is column-for-column fair.
p4f <- file.path(P4_DIR, "concordance", "CCC_CONCORDANCE.tsv")
if (!file.exists(p4f)) p4f <- file.path(TAB, "CCC_CONCORDANCE_REFINED.tsv")
stopifnot(file.exists(p4f))
p4 <- read.delim(p4f, stringsAsFactors = FALSE)
log_("Phase 4 refined concordance: ", nrow(p4), " rows x ", ncol(p4), " cols")
log_("schemas identical: ", identical(sort(colnames(p3)), sort(colnames(p4))))

as_lgl <- function(v) { if (is.logical(v)) return(ifelse(is.na(v), FALSE, v))
  v <- tolower(trimws(as.character(v))); !is.na(v) & v %in% c("true","t","yes","1") }

# Column names are taken from the build_concordance.R schema, not guessed:
#   sender receiver ligand_n receptor_n
#   liana_sup cellchat_sup cpdb_sup   liana_testable cellchat_testable cpdb_testable
#   n_LR_methods_supported n_methods_testable concordance_class
#   samples_supported patients_supported
norm_tbl <- function(df, prefix) {
  need <- c("sender","receiver","ligand_n","receptor_n","liana_sup","cellchat_sup",
            "cpdb_sup","n_LR_methods_supported","n_methods_testable",
            "concordance_class","patients_supported","samples_supported")
  miss <- setdiff(need, colnames(df))
  if (length(miss)) stop("missing expected concordance columns: ", paste(miss, collapse = ", "))
  out <- data.frame(
    sender = as.character(df$sender), receiver = as.character(df$receiver),
    ligand = as.character(df$ligand_n), receptor = as.character(df$receptor_n),
    n_frameworks = as.integer(df$n_LR_methods_supported),
    n_testable   = as.integer(df$n_methods_testable),
    concordance  = as.character(df$concordance_class),
    patients_supported = suppressWarnings(as.numeric(df$patients_supported)),
    samples_supported  = suppressWarnings(as.numeric(df$samples_supported)),
    liana = as_lgl(df$liana_sup), cellchat = as_lgl(df$cellchat_sup),
    cpdb = as_lgl(df$cpdb_sup), stringsAsFactors = FALSE)
  out$support <- out$n_frameworks >= 1
  nm <- c("n_frameworks","n_testable","concordance","patients_supported",
          "samples_supported","liana","cellchat","cpdb","support")
  names(out)[match(nm, names(out))] <- paste0(prefix, "_", nm)
  out
}
t3 <- norm_tbl(p3, "phase3")
t4 <- norm_tbl(p4, "phase4")
log_("Phase 3 supported: ", sum(t3$phase3_support), " | Phase 4 supported: ", sum(t4$phase4_support))

# --- 3. label-lineage mapping ---------------------------------------------
sec("3. Mapping Phase 3 populations onto Phase 4 populations")
# Phase 4 collapses newly-confirmed malignant cells into MPNST-Tumor, so a Phase 3
# interaction whose sender was Fibroblast may reappear with sender MPNST-Tumor.
# That is a SENDER REASSIGNMENT, not a loss - which is why change classes include
# "Sender-reassigned" (§42).
canon <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
t3$key  <- paste(canon(t3$sender), canon(t3$receiver), canon(t3$ligand), canon(t3$receptor), sep = "|")
t4$key  <- paste(canon(t4$sender), canon(t4$receiver), canon(t4$ligand), canon(t4$receptor), sep = "|")
t3$lrk  <- paste(canon(t3$ligand), canon(t3$receptor), sep = "|")
t4$lrk  <- paste(canon(t4$ligand), canon(t4$receptor), sep = "|")

# canon() strips punctuation, so an upstream resource recording the same complex
# both as "CD8 RECEPTOR" and "CD8_RECEPTOR" collapses to one key. That is correct
# behaviour - they ARE the same interaction - but it leaves duplicate rows and a
# many-to-many join. Collapse them, keeping the strongest support, and report how
# many were affected rather than silencing the warning.
collapse_dups <- function(d, prefix) {
  dup <- sum(duplicated(d$key))
  if (!dup) return(d)
  log_(prefix, ": collapsing ", dup, " duplicate key(s) created by punctuation variants ",
       "in the upstream resource (e.g. 'CD8 RECEPTOR' vs 'CD8_RECEPTOR')")
  num <- grep("_(n_frameworks|n_testable|patients_supported|samples_supported)$",
              names(d), value = TRUE)
  lgl <- grep("_(liana|cellchat|cpdb|support)$", names(d), value = TRUE)
  spl <- split(seq_len(nrow(d)), d$key)
  keep <- vapply(spl, function(i) i[1], integer(1))
  out <- d[keep, , drop = FALSE]
  multi <- spl[lengths(spl) > 1]
  for (k in names(multi)) {
    i <- multi[[k]]; j <- which(out$key == k)
    for (cc in num) out[j, cc] <- suppressWarnings(max(d[i, cc], na.rm = TRUE))
    for (cc in lgl) out[j, cc] <- any(d[i, cc], na.rm = TRUE)
  }
  for (cc in num) out[[cc]][!is.finite(out[[cc]])] <- NA_real_
  out
}
t3 <- collapse_dups(t3, "Phase 3")
t4 <- collapse_dups(t4, "Phase 4")
stopifnot(!any(duplicated(t3$key)), !any(duplicated(t4$key)))
log_("distinct keys - Phase 3: ", nrow(t3), " | Phase 4: ", nrow(t4))

m <- full_join(t3 |> select(-sender,-receiver,-ligand,-receptor),
               t4 |> select(-sender,-receiver,-ligand,-receptor), by = c("key","lrk"))
ann <- bind_rows(t3 |> select(key, sender, receiver, ligand, receptor),
                 t4 |> select(key, sender, receiver, ligand, receptor)) |> distinct(key, .keep_all = TRUE)
m <- left_join(m, ann, by = "key")
log_("union of interaction keys: ", nrow(m))

# --- 4. change classification (§42) ---------------------------------------
sec("4. Change classification")
z <- function(v) ifelse(is.na(v), 0, v)
m$phase3_support <- ifelse(is.na(m$phase3_support), FALSE, m$phase3_support)
m$phase4_support <- ifelse(is.na(m$phase4_support), FALSE, m$phase4_support)

# Sender reassignment: the same ligand-receptor pair, same receiver, supported in
# both phases, but the sender identity moved into MPNST-Tumor.
reasg <- m |> filter(phase4_support) |>
  select(lrk, receiver, p4_sender = sender) |> distinct() |>
  inner_join(t3 |> filter(phase3_support) |> select(lrk, receiver, p3_sender = sender) |> distinct(),
             by = c("lrk","receiver")) |>
  filter(p4_sender != p3_sender, p4_sender == "MPNST-Tumor") |>
  mutate(rk = paste(lrk, canon(receiver), sep = "||")) |> pull(rk) |> unique()
m$rk <- paste(m$lrk, canon(m$receiver), sep = "||")

m$change_class <- with(m, case_when(
  !phase3_support & !phase4_support                    ~ "Not-supported-either",
   phase3_support & !phase4_support &  rk %in% reasg   ~ "Sender-reassigned",
  !phase3_support &  phase4_support                    ~ "Newly-supported",
   phase3_support & !phase4_support                    ~ "Lost",
   z(phase4_n_frameworks) >  z(phase3_n_frameworks)    ~ "Strengthened",
   z(phase4_n_frameworks) <  z(phase3_n_frameworks)    ~ "Weakened",
   z(phase4_patients_supported) > z(phase3_patients_supported) ~ "Strengthened",
   z(phase4_patients_supported) < z(phase3_patients_supported) ~ "Weakened",
   TRUE                                                ~ "Stable"))
# Interactions untestable in Phase 4 because a population stopped existing are
# Ambiguous, not Lost: absence of evidence is not evidence of absence.
p4_pops <- unique(c(t4$sender, t4$receiver))
m$change_class[m$change_class == "Lost" &
  (!(m$sender %in% p4_pops) | !(m$receiver %in% p4_pops))] <- "Ambiguous"
log_("change_class distribution:")
print(sort(table(m$change_class), decreasing = TRUE))

out <- m |> transmute(sender, receiver, ligand, receptor,
  phase3_support, phase4_support,
  phase3_patients_supported, phase4_patients_supported,
  phase3_concordance, phase4_concordance,
  phase3_n_frameworks, phase4_n_frameworks,
  phase3_n_testable, phase4_n_testable,
  change_class,
  notes = case_when(
    change_class == "Sender-reassigned" ~
      "same ligand-receptor and receiver; sender identity moved into MPNST-Tumor after malignancy refinement",
    change_class == "Ambiguous" ~
      "a participating population does not exist under annotation_ccc_refined; NOT EVALUABLE, not absent",
    change_class == "Newly-supported" ~
      "supported only after refinement; the refined tumour compartment includes cells Phase 3 labelled stroma",
    TRUE ~ "")) |>
  arrange(factor(change_class, levels = c("Sender-reassigned","Strengthened","Newly-supported",
                                          "Stable","Weakened","Lost","Ambiguous",
                                          "Not-supported-either")),
          desc(phase4_n_frameworks), desc(phase4_patients_supported))
f1 <- file.path(P4_DIR, "PHASE3_VS_PHASE4_CCC.tsv")
dir.create(dirname(f1), showWarnings = FALSE, recursive = TRUE)
write.table(out, f1, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(out, file.path(TAB, "PHASE3_VS_PHASE4_CCC.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
log_("wrote ", f1, " (", nrow(out), " rows)")

# --- 5. the named axes (§38) ----------------------------------------------
sec("5. Verdict on each prioritized Phase 3 axis")
FOCUS$lrk <- paste(canon(FOCUS$ligand), canon(FOCUS$receptor), sep = "|")
fa <- FOCUS |> left_join(
  m |> group_by(lrk) |> summarise(
    n_keys = n(),
    p3_supported_keys = sum(phase3_support),
    p4_supported_keys = sum(phase4_support),
    p3_max_frameworks = suppressWarnings(max(phase3_n_frameworks, na.rm = TRUE)),
    p4_max_frameworks = suppressWarnings(max(phase4_n_frameworks, na.rm = TRUE)),
    p3_max_patients = suppressWarnings(max(phase3_patients_supported, na.rm = TRUE)),
    p4_max_patients = suppressWarnings(max(phase4_patients_supported, na.rm = TRUE)),
    p3_senders = paste(sort(unique(sender[phase3_support])), collapse = "; "),
    p4_senders = paste(sort(unique(sender[phase4_support])), collapse = "; "),
    p3_receivers = paste(sort(unique(receiver[phase3_support])), collapse = "; "),
    p4_receivers = paste(sort(unique(receiver[phase4_support])), collapse = "; "),
    classes = paste(sprintf("%s=%d", names(sort(table(change_class), decreasing = TRUE)),
                            sort(table(change_class), decreasing = TRUE)), collapse = "; "),
    .groups = "drop"), by = "lrk")
fix <- function(v) ifelse(is.finite(v), v, NA_real_)
fa <- fa |> mutate(across(starts_with("p3_max"), fix), across(starts_with("p4_max"), fix),
                   p3_supported_keys = ifelse(is.na(p3_supported_keys), 0L, p3_supported_keys),
                   p4_supported_keys = ifelse(is.na(p4_supported_keys), 0L, p4_supported_keys))
fa$verdict <- with(fa, case_when(
  p3_supported_keys == 0 & p4_supported_keys == 0 ~ "not supported in either phase",
  p3_supported_keys == 0 & p4_supported_keys  > 0 ~ "newly supported",
  p4_supported_keys == 0                          ~ "lost",
  p4_supported_keys > p3_supported_keys |
    (!is.na(p4_max_frameworks) & !is.na(p3_max_frameworks) & p4_max_frameworks > p3_max_frameworks) ~ "strengthened",
  p4_supported_keys < p3_supported_keys |
    (!is.na(p4_max_frameworks) & !is.na(p3_max_frameworks) & p4_max_frameworks < p3_max_frameworks) ~ "weakened",
  TRUE ~ "stable"))
fa$sender_changed <- with(fa, !is.na(p3_senders) & !is.na(p4_senders) & p3_senders != p4_senders)
print(as.data.frame(fa |> select(ligand, receptor, axis, p3_supported_keys, p4_supported_keys,
                                 p3_max_frameworks, p4_max_frameworks, verdict, sender_changed)),
      row.names = FALSE)
write.table(fa, file.path(TAB, "PHASE3_AXIS_SENSITIVITY.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- 6. the fibroblast question (§43) -------------------------------------
sec("6. Were apparent fibroblast interactions actually tumour interactions?")
fib3 <- t3 |> filter(phase3_support, sender == "Fibroblast" | receiver == "Fibroblast")
log_("Phase 3 supported interactions involving Fibroblast: ", nrow(fib3))
if (nrow(fib3)) {
  fib3$rk <- paste(fib3$lrk, canon(fib3$receiver), sep = "||")
  fib3$still <- fib3$key %in% m$key[m$phase4_support]
  fib3$reassigned <- fib3$rk %in% reasg
  s <- fib3 |> summarise(
    n = n(),
    retained_as_fibroblast = sum(still),
    reassigned_to_tumour = sum(!still & reassigned),
    lost_or_ambiguous = sum(!still & !reassigned))
  print(as.data.frame(s), row.names = FALSE)
  log_(sprintf("=> %.1f%% of Phase 3 fibroblast interactions are reassigned to the refined tumour compartment",
               100*s$reassigned_to_tumour/max(1,s$n)))
  write.table(fib3, file.path(TAB, "FIBROBLAST_INTERACTION_REASSIGNMENT.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  fib_summary <- as.list(s)
} else fib_summary <- list(n = 0L)

# --- 7. figures ------------------------------------------------------------
sec("7. Figures")
cl_ord <- c("Sender-reassigned","Strengthened","Newly-supported","Stable",
            "Weakened","Lost","Ambiguous","Not-supported-either")
cpal <- setNames(c("#6A3D9A","#B2182B","#E7298A","#4D9221","#F0A202","#999999",
                   "#80B1D3","grey85"), cl_ord)
d <- m |> filter(change_class != "Not-supported-either") |>
  count(change_class) |> mutate(change_class = factor(change_class, levels = cl_ord))
g1 <- ggplot(d, aes(y = change_class, x = n, fill = change_class)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = format(n, big.mark = ",")), hjust = -0.12, size = 3) +
  scale_fill_manual(values = cpal, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
  scale_y_discrete(limits = rev(cl_ord)) +
  p4_theme() +
  labs(title = "How Phase 3 interactions change after malignancy refinement",
       subtitle = paste("Transparent categories only - no weighted composite score was invented (§42).",
                        "\n'Ambiguous' means a participating population no longer exists: NOT EVALUABLE, not absent."),
       x = "interaction keys", y = NULL,
       caption = paste("Same framework versions, thresholds, resources, expression basis and sample-aware design",
                       "as Phase 3.\nOnly the ccc_label definition differs, so this is a label-sensitivity analysis.",
                       sep = "\n"))
save_fig(g1, file.path(FIG, "32_01_phase3_vs_phase4_ccc.pdf"), 10.5, 5.4)

fd <- fa |> mutate(pair = paste0(ligand, " → ", receptor),
                   verdict = factor(verdict, levels = c("strengthened","newly supported","stable",
                                                        "weakened","lost","not supported in either phase")))
vpal <- c(strengthened = "#B2182B", `newly supported` = "#E7298A", stable = "#4D9221",
          weakened = "#F0A202", lost = "#999999", `not supported in either phase` = "grey85")
g2 <- ggplot(fd, aes(y = reorder(pair, as.integer(verdict)), x = verdict, fill = verdict)) +
  geom_tile(colour = "white", linewidth = 0.5, width = 0.9, height = 0.85) +
  scale_fill_manual(values = vpal, guide = "none") +
  facet_grid(axis ~ ., scales = "free_y", space = "free_y", switch = "y") +
  p4_theme() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 6.6),
        strip.placement = "outside") +
  labs(title = "Sensitivity of each prioritized Phase 3 axis to malignancy refinement",
       subtitle = "The axes named in the Phase 4 authorization (§38), grouped by biological axis",
       x = NULL, y = NULL)
save_fig(g2, file.path(FIG, "32_02_axis_sensitivity.pdf"), 11.5, 9.6)

facts <- list(
  phase3_master = list(path = P3_MASTER, rows = nrow(p3m)),
  phase3_concordance = list(path = P3_CONC, rows = nrow(p3)),
  phase4_concordance = list(path = p4f, rows = nrow(p4)),
  union_keys = nrow(m),
  change_class_counts = as.list(sort(table(m$change_class), decreasing = TRUE)),
  sender_reassignment_groups = length(reasg),
  axis_verdicts = fa |> select(ligand, receptor, axis, verdict, sender_changed,
                               p3_supported_keys, p4_supported_keys),
  fibroblast_question = fib_summary,
  method_reuse = paste("scripts/phase3/ccc/run_liana.R, run_cellchat.R and",
    "run_cellphonedb.py reused UNMODIFIED at Phase 3 versions (liana 0.1.14,",
    "CellChat 2.2.0.9001, CellPhoneDB 5.0.1); only ccc_label differs."),
  no_composite_score = TRUE,
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  slurm = list(job_id = Sys.getenv("SLURM_JOB_ID")))
write_json(facts, file.path(P4_DIR, "m32_sensitivity_facts.json"),
           auto_unbox = TRUE, pretty = TRUE, digits = 8, null = "null")
sec("M32 SENSITIVITY ANALYSIS COMPLETE")
