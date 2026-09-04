#!/usr/bin/env Rscript
# =============================================================================
# Phase 6 - M42 - reconstruction and feasibility.
#
# Verifies the Phase 5 object by md5 AND sha256 against its manifest, confirms
# every field Phase 6 needs, re-asserts the frozen Phase 4 counts, and exports
# the small durable Phase 6 working artefacts so the large object is opened
# exactly once.
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject) })
source("scripts/phase6/utils/phase6_common.R")
set.seed(42)
facts <- list(milestone = "M42", generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))

p6_sec("1. Verify the Phase 5 object against its manifest")
mf <- fromJSON("results/phase5/phase5_manifest.json", simplifyVector = FALSE)
exp_md5 <- mf$phase5_object$md5; exp_sha <- mf$phase5_object$sha256
md5 <- p6_md5(P5_OBJECT); sha <- p6_sha256(P5_OBJECT)
p6_msg("manifest md5    %s", exp_md5); p6_msg("observed md5    %s", md5)
p6_msg("manifest sha256 %s", exp_sha); p6_msg("observed sha256 %s", sha)
stopifnot(md5 == exp_md5, sha == exp_sha)
p6_msg("  [OK] Phase 5 object verified")
facts$phase5_object <- list(path = P5_OBJECT, md5 = md5, sha256 = sha,
                            size_bytes = file.info(P5_OBJECT)$size,
                            manifest_verified = TRUE)
KSTAR <- mf$selected_K
p6_msg("Phase 5 selected K = %s", KSTAR)

p6_sec("2. Load")
obj <- readRDS(P5_OBJECT)
md <- obj@meta.data; md$cell_id <- rownames(md)
p6_msg("cells %d | metadata columns %d | assays %s", ncol(obj), ncol(md) - 1L,
       paste(Assays(obj), collapse = ", "))
stopifnot(ncol(obj) == 19716L)

p6_sec("3. Frozen Phase 4 counts still hold inside the Phase 5 object")
ref <- table(md$malignancy_refined)
FROZ <- c(Malignant = 6434L, `Non-malignant` = 9078L, Ambiguous = 3766L,
          `Excluded-low-quality` = 438L)
for (k in names(FROZ))
  p6_msg("  [%s] refined_%-22s expected %-5d observed %d",
         if (ref[[k]] == FROZ[[k]]) "OK  " else "FAIL", k, FROZ[[k]], ref[[k]])
stopifnot(all(vapply(names(FROZ), function(k) ref[[k]] == FROZ[[k]], logical(1))))

p6_sec("4. Fields Phase 6 requires")
PROGS <- grep("^program_P[0-9]+_score$", colnames(md), value = TRUE)
REQ <- c("tumor_clone_phase4", "scevan_clone", "scevan_sample_reliable",
         "cnv_burden", "cnv_frac_gain", "cnv_frac_loss", "cnv_mean_abs",
         "sample_id", "malignancy_confidence", "malignancy_refined",
         "tumor_state_phase4", "dominant_malignant_program",
         "dominant_malignant_program_label", "dominant_program_score",
         "dominant_program_margin", "malignant_program_confidence",
         "program_score_source", "nCount_RNA", "nFeature_RNA", PROGS)
fld <- do.call(rbind, lapply(REQ, function(f) data.frame(
  field = f, present = f %in% colnames(md),
  n_non_na = if (f %in% colnames(md)) sum(!is.na(md[[f]])) else NA_integer_)))
print(fld, row.names = FALSE, max = 400)
p6_tsv(fld, file.path(P6_VAL, "M42_REQUIRED_FIELD_AUDIT.tsv"))
stopifnot(all(fld$present))
p6_msg("  [OK] all %d required fields present (%d program-score fields)",
       nrow(fld), length(PROGS))
facts$programs <- PROGS; facts$n_programs <- length(PROGS)

p6_sec("5. MPNST_3 clone reliability - re-asserted, not assumed")
rel <- md |> group_by(sample_id) |>
  summarise(scevan_sample_reliable = unique(as.character(scevan_sample_reliable))[1],
            n_clones = length(unique(na.omit(scevan_clone))),
            n_malignant = sum(malignancy_refined == "Malignant"),
            n_malignant_with_clone = sum(malignancy_refined == "Malignant" &
                                           !is.na(tumor_clone_phase4)),
            .groups = "drop")
print(as.data.frame(rel), row.names = FALSE)
p6_tsv(rel, file.path(P6_VAL, "M42_CLONE_RELIABILITY_BY_PATIENT.tsv"))
unrel <- rel$sample_id[rel$scevan_sample_reliable %in% c(FALSE, "FALSE")]
p6_msg("unreliable by the frozen Phase 4 flag: %s", paste(unrel, collapse = ", "))
stopifnot(identical(as.character(unrel), CLONE_EXCLUDED))
p6_msg("  [OK] the frozen flag independently reproduces the declared exclusion")
p6_msg("%s", CLONE_EXCLUSION_REASON)
facts$clone_reliability <- rel
facts$clone_reliable_patients <- CLONE_RELIABLE
facts$clone_excluded_patients <- CLONE_EXCLUDED
facts$clone_exclusion_reason <- CLONE_EXCLUSION_REASON

p6_sec("6. Clone x program evaluability")
mal <- md[md$malignancy_refined == "Malignant", ]
cl <- mal |> filter(!is.na(tumor_clone_phase4)) |>
  group_by(sample_id, tumor_clone_phase4) |>
  summarise(n_cells = n(), .groups = "drop") |>
  mutate(clone_reliable = sample_id %in% CLONE_RELIABLE,
         ge20 = n_cells >= 20, ge50 = n_cells >= 50, ge100 = n_cells >= 100)
print(as.data.frame(cl), row.names = FALSE)
p6_tsv(cl, file.path(P6_VAL, "M42_CLONE_SIZES.tsv"))
p6_msg("clones in the reliable patients: %d (>=20 cells: %d, >=50: %d, >=100: %d)",
       sum(cl$clone_reliable), sum(cl$clone_reliable & cl$ge20),
       sum(cl$clone_reliable & cl$ge50), sum(cl$clone_reliable & cl$ge100))
p6_msg("malignant cells WITHOUT a clone label: %d (Phase 4 assigned clones only to cells its SCEVAN run placed)",
       sum(is.na(mal$tumor_clone_phase4)))
facts$clone_sizes <- cl
facts$malignant_without_clone <- sum(is.na(mal$tumor_clone_phase4))

p6_sec("7. Export Phase 6 working artefacts")
dir.create(P6_CLONE, showWarnings = FALSE, recursive = TRUE)
saveRDS(md, file.path(P6_CLONE, "phase6_cell_metadata.rds"))
p6_msg("  [rds] phase6_cell_metadata.rds (%d x %d)", nrow(md), ncol(md))
rna <- JoinLayers(obj[["RNA"]])
malids <- mal$cell_id
lg <- LayerData(rna, layer = "data")[, malids, drop = FALSE]
saveRDS(lg, file.path(P6_CLONE, "phase6_malignant_lognorm.rds"))
p6_msg("  [rds] phase6_malignant_lognorm.rds (%d genes x %d malignant cells)",
       nrow(lg), ncol(lg))
rm(rna, obj); invisible(gc())

p6_sec("8. Phase 4 and Phase 5 objects untouched")
stopifnot(p6_md5(P4_OBJECT) == P4_MD5, p6_md5(P5_OBJECT) == md5)
p6_msg("Phase 4 md5 %s unchanged | Phase 5 md5 %s unchanged", P4_MD5, md5)
facts$upstream_unchanged <- TRUE
p6_json(facts, file.path(P6_VAL, "m42_feasibility_facts.json"))
p6_sec("M42 complete")
