#!/bin/bash
#SBATCH --job-name=p3_install3
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=logs/phase3/M18/install_retry2_%j.out
#SBATCH --error=logs/phase3/M18/install_retry2_%j.err
# Phase 3 dependency retry #2 — two precisely diagnosed blockers remain:
#   CellChat  : svglite 2.2.2 requires systemfonts >= 1.3.0; conda gave 1.2.3.
#   nichenetr : 'shadowtext' unavailable (and units/sf/ggiraph/gdtools need system
#               libraries: udunits2, GDAL/PROJ/GEOS).
# FIX: install those specific packages as prebuilt conda-forge R 4.4 binaries, with r-base
# PINNED to 4.4.3 so the interpreter cannot move under the frozen Phase 2 stack.
set -uo pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase3/M18
echo "=== retry2 start: $(date) | Job ${SLURM_JOB_ID:-NA} | $(hostname) ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
export MAKEFLAGS="-j8"

echo "--- baseline BEFORE ---"
Rscript -e 'for (p in c("Seurat","SeuratObject","harmony","Matrix","sctransform","liana","systemfonts")) cat(sprintf("%s %s\n", p, tryCatch(as.character(packageVersion(p)), error=function(e)"NA")))'

echo "=== conda-forge: systemfonts>=1.3.2 + nichenetr system-library deps (r-base pinned) ==="
conda install -y -c conda-forge \
  "r-base==4.4.3" "r-systemfonts>=1.3.2" r-shadowtext r-units r-sf r-ggiraph r-gdtools \
  r-fdrtool r-rocr r-catools r-hmisc 2>&1 | tail -30 || echo "CONDA INSTALL RETURNED NONZERO"

Rscript -e '
need <- c("systemfonts","svglite","shadowtext","units","sf","ggiraph","gdtools","fdrtool","ROCR","caTools","Hmisc")
for (p in need) cat(sprintf("  %-14s %s\n", p, tryCatch(as.character(packageVersion(p)), error=function(e) "STILL_MISSING")))'

echo "=== RETRY: CellChat ==="
Rscript -e '
options(repos=c(CRAN="https://cloud.r-project.org"), Ncpus=8)
if (!requireNamespace("CellChat", quietly=TRUE))
  tryCatch(remotes::install_github("jinworks/CellChat", upgrade="never", dependencies=TRUE),
           error=function(e) cat("CELLCHAT ERROR:", conditionMessage(e), "\n"))
cat("CellChat:", tryCatch(as.character(packageVersion("CellChat")), error=function(e) "FAILED"), "\n")' 2>&1 | tail -20

echo "=== RETRY: nichenetr (mlrMBO chain may be heavy) ==="
Rscript -e '
options(repos=c(CRAN="https://cloud.r-project.org"), Ncpus=8)
for (p in c("mlrMBO","ParamHelpers","smoof","mlr","lhs","emoa")) {
  if (!requireNamespace(p, quietly=TRUE)) try(install.packages(p), silent=TRUE)
  cat(sprintf("  %-14s %s\n", p, tryCatch(as.character(packageVersion(p)), error=function(e)"MISSING")))
}
if (!requireNamespace("nichenetr", quietly=TRUE))
  tryCatch(remotes::install_github("saeyslab/nichenetr", upgrade="never", dependencies=TRUE),
           error=function(e) cat("NICHENET ERROR:", conditionMessage(e), "\n"))
cat("nichenetr:", tryCatch(as.character(packageVersion("nichenetr")), error=function(e) "FAILED"), "\n")' 2>&1 | tail -25

echo "--- baseline AFTER (must be unchanged) ---"
Rscript -e 'for (p in c("Seurat","SeuratObject","harmony","Matrix","sctransform","liana")) cat(sprintf("%s %s\n", p, tryCatch(as.character(packageVersion(p)), error=function(e)"NA")))'
Rscript -e 'cat("R version:", as.character(getRversion()), "\n")'
echo "=== FINAL ROSTER ==="
Rscript -e '
pk <- c("liana","CellChat","nichenetr","OmnipathR","decoupleR","NMF","circlize","ComplexHeatmap",
        "igraph","ggraph","ggpubr","svglite","systemfonts","FNN","RANN","Seurat","SeuratObject",
        "harmony","Matrix","sctransform","presto","aricode")
for (p in pk) cat(sprintf("%-16s %s\n", p, tryCatch(as.character(packageVersion(p)), error=function(e) "NOT_INSTALLED")))'
echo "=== retry2 end: $(date) ==="
