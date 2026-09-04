#!/bin/bash
#SBATCH --job-name=p3_install
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=logs/phase3/M18/install_%j.out
#SBATCH --error=logs/phase3/M18/install_%j.err
# Phase 3 dependency installation. Staged and fault-tolerant: each stage is
# attempted independently and its success/failure recorded, so a single failing
# package cannot abort the whole install. Foundational packages (Seurat,
# SeuratObject, harmony, Matrix) are pinned and must NOT change.
set -uo pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase3/M18 reports/phase3/environment external
echo "=== Phase 3 install start: $(date) | Job ${SLURM_JOB_ID:-NA} | $(hostname) ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
export MAKEFLAGS="-j8"
R --version | head -1

echo "--- baseline versions (must not change) ---"
Rscript -e 'for (p in c("Seurat","SeuratObject","harmony","Matrix","sctransform")) cat(sprintf("%s %s\n", p, packageVersion(p)))'

# ---------- Stage 1: CRAN support packages ----------
echo "=== STAGE 1: CRAN support packages ==="
Rscript -e '
options(repos=c(CRAN="https://cloud.r-project.org"), Ncpus=8)
need <- c("NMF","ggalluvial","svglite","ggraph","ggnetwork","ggrepel","ggpubr",
          "future.apply","BiocManager","remotes","sna","network","statnet.common",
          "FNN","RANN","irlba","reshape2","cowplot","patchwork","circlize","psych","caret","randomForest","e1071","DiagrammeR")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly=TRUE)]
cat("missing:", paste(miss, collapse=", "), "\n")
if (length(miss)) install.packages(miss)
for (p in need) cat(sprintf("  %-20s %s\n", p, tryCatch(as.character(packageVersion(p)), error=function(e) "FAILED")))
' 2>&1 | tail -40

# ---------- Stage 2: Bioconductor packages ----------
echo "=== STAGE 2: Bioconductor packages ==="
Rscript -e '
options(repos=c(CRAN="https://cloud.r-project.org"), Ncpus=8)
library(BiocManager)
need <- c("OmnipathR","decoupleR","BiocNeighbors","ComplexHeatmap","basilisk","SingleCellExperiment","SummarizedExperiment","scater","scran","limma","edgeR","AnnotationDbi","org.Hs.eg.db")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly=TRUE)]
cat("missing:", paste(miss, collapse=", "), "\n")
if (length(miss)) BiocManager::install(miss, ask=FALSE, update=FALSE)
for (p in need) cat(sprintf("  %-24s %s\n", p, tryCatch(as.character(packageVersion(p)), error=function(e) "FAILED")))
' 2>&1 | tail -40

# ---------- Stage 3: LIANA ----------
echo "=== STAGE 3: LIANA (saezlab/liana) ==="
Rscript -e '
options(repos=c(CRAN="https://cloud.r-project.org"), Ncpus=8)
if (!requireNamespace("liana", quietly=TRUE))
  tryCatch(remotes::install_github("saezlab/liana", upgrade="never", dependencies=TRUE),
           error=function(e) cat("LIANA INSTALL ERROR:", conditionMessage(e), "\n"))
cat("liana:", tryCatch(as.character(packageVersion("liana")), error=function(e) "FAILED"), "\n")
' 2>&1 | tail -30

# ---------- Stage 4: CellChat ----------
echo "=== STAGE 4: CellChat (jinworks/CellChat) ==="
Rscript -e '
options(repos=c(CRAN="https://cloud.r-project.org"), Ncpus=8)
if (!requireNamespace("CellChat", quietly=TRUE))
  tryCatch(remotes::install_github("jinworks/CellChat", upgrade="never", dependencies=TRUE),
           error=function(e) cat("CELLCHAT INSTALL ERROR:", conditionMessage(e), "\n"))
cat("CellChat:", tryCatch(as.character(packageVersion("CellChat")), error=function(e) "FAILED"), "\n")
' 2>&1 | tail -30

# ---------- Stage 5: NicheNet ----------
echo "=== STAGE 5: nichenetr (saeyslab/nichenetr) ==="
Rscript -e '
options(repos=c(CRAN="https://cloud.r-project.org"), Ncpus=8)
if (!requireNamespace("nichenetr", quietly=TRUE))
  tryCatch(remotes::install_github("saeyslab/nichenetr", upgrade="never", dependencies=TRUE),
           error=function(e) cat("NICHENET INSTALL ERROR:", conditionMessage(e), "\n"))
cat("nichenetr:", tryCatch(as.character(packageVersion("nichenetr")), error=function(e) "FAILED"), "\n")
' 2>&1 | tail -30

# ---------- Stage 6: vendor MMCA reference (scripts repo, NOT an R package) ----------
echo "=== STAGE 6: vendor MMCA LochNESS reference ==="
if [ ! -d external/MMCA/.git ]; then
  rm -rf external/MMCA
  git clone --depth 50 https://github.com/shendurelab/MMCA.git external/MMCA 2>&1 | tail -5 || echo "MMCA CLONE FAILED"
fi
if [ -d external/MMCA/.git ]; then
  ( cd external/MMCA && echo "MMCA commit: $(git rev-parse HEAD)" && git log -1 --format='%H %ad %s' )
  find external/MMCA -iname "*lochness*" | head -20
fi

# ---------- Stage 7: CellPhoneDB in an isolated env ----------
echo "=== STAGE 7: CellPhoneDB (isolated env cpdb_env) ==="
# Isolated deliberately: cellphonedb pins numpy/pandas versions that would
# destabilise R_env (numpy 2.3.5, reticulate). R_env is left untouched.
if ! conda env list | grep -qE "^cpdb_env "; then
  conda create -y -n cpdb_env python=3.10 2>&1 | tail -3
fi
conda activate cpdb_env
pip install --quiet "cellphonedb" 2>&1 | tail -10 || echo "CELLPHONEDB PIP FAILED"
python -c "
import importlib
for m in ['cellphonedb','anndata','scanpy','numpy','pandas']:
    try:
        mod=importlib.import_module(m); print(f'  {m:16s} {getattr(mod,\"__version__\",\"?\")}')
    except Exception as e: print(f'  {m:16s} FAILED ({type(e).__name__})')
"
conda activate R_env

echo "--- baseline versions AFTER install (must be unchanged) ---"
Rscript -e 'for (p in c("Seurat","SeuratObject","harmony","Matrix","sctransform")) cat(sprintf("%s %s\n", p, packageVersion(p)))'
echo "=== Phase 3 install end: $(date) ==="
