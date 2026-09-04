#!/bin/bash
#SBATCH --job-name=p3_install2
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=logs/phase3/M18/install_retry_%j.out
#SBATCH --error=logs/phase3/M18/install_retry_%j.err
# Phase 3 dependency retry.
# ROOT CAUSE of the first attempt's CellChat and nichenetr failures: source compilation of
# `systemfonts` and `Deriv` failed, which cascaded to svglite, doBy, pbkrtest, car, rstatix
# and ggpubr. CellChat and nichenetr both hard-depend on that chain.
# FIX: install those specific packages as prebuilt conda-forge R 4.4 binaries instead of
# compiling from source, with --freeze-installed so no existing package (Seurat, Matrix,
# harmony, sctransform) can be altered. Then retry CellChat and nichenetr.
set -uo pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase3/M18
echo "=== retry start: $(date) | Job ${SLURM_JOB_ID:-NA} | $(hostname) ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
export MAKEFLAGS="-j8"

echo "--- baseline BEFORE (must not change) ---"
Rscript -e 'for (p in c("Seurat","SeuratObject","harmony","Matrix","sctransform","liana")) cat(sprintf("%s %s\n", p, tryCatch(as.character(packageVersion(p)), error=function(e)"NA")))'

echo "=== FIX: prebuilt conda-forge binaries for the failing compile chain ==="
conda install -y -c conda-forge --freeze-installed \
  r-systemfonts r-svglite r-deriv r-doby r-pbkrtest r-car r-rstatix r-ggpubr \
  r-ggraph r-tidygraph r-diagrammer r-caret r-randomforest r-e1071 2>&1 | tail -25 \
  || echo "CONDA FREEZE-INSTALL FAILED; retrying without --freeze-installed for this set only"

Rscript -e '
need <- c("systemfonts","svglite","Deriv","doBy","pbkrtest","car","rstatix","ggpubr","ggraph","tidygraph","DiagrammeR","caret","randomForest","e1071")
for (p in need) cat(sprintf("  %-16s %s\n", p, tryCatch(as.character(packageVersion(p)), error=function(e) "STILL_MISSING")))'

echo "=== RETRY: CellChat ==="
Rscript -e '
options(repos=c(CRAN="https://cloud.r-project.org"), Ncpus=8)
if (!requireNamespace("CellChat", quietly=TRUE))
  tryCatch(remotes::install_github("jinworks/CellChat", upgrade="never", dependencies=TRUE),
           error=function(e) cat("CELLCHAT RETRY ERROR:", conditionMessage(e), "\n"))
cat("CellChat:", tryCatch(as.character(packageVersion("CellChat")), error=function(e) "FAILED"), "\n")' 2>&1 | tail -25

echo "=== RETRY: nichenetr ==="
Rscript -e '
options(repos=c(CRAN="https://cloud.r-project.org"), Ncpus=8)
if (!requireNamespace("nichenetr", quietly=TRUE))
  tryCatch(remotes::install_github("saeyslab/nichenetr", upgrade="never", dependencies=TRUE),
           error=function(e) cat("NICHENET RETRY ERROR:", conditionMessage(e), "\n"))
cat("nichenetr:", tryCatch(as.character(packageVersion("nichenetr")), error=function(e) "FAILED"), "\n")' 2>&1 | tail -25

echo "--- baseline AFTER (must be unchanged) ---"
Rscript -e 'for (p in c("Seurat","SeuratObject","harmony","Matrix","sctransform","liana")) cat(sprintf("%s %s\n", p, tryCatch(as.character(packageVersion(p)), error=function(e)"NA")))'

echo "=== FINAL Phase 3 package roster ==="
Rscript -e '
pk <- c("liana","CellChat","nichenetr","OmnipathR","decoupleR","basilisk","NMF","circlize",
        "ComplexHeatmap","igraph","ggraph","tidygraph","ggpubr","svglite","FNN","RANN",
        "Seurat","SeuratObject","harmony","Matrix","sctransform","presto","aricode","reshape2")
for (p in pk) cat(sprintf("%-18s %s\n", p, tryCatch(as.character(packageVersion(p)), error=function(e) "NOT_INSTALLED")))'
echo "=== retry end: $(date) ==="
