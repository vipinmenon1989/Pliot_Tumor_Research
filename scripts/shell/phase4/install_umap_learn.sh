#!/bin/bash
#SBATCH --job-name=p4_umapenv
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/phase4/M28/install_umap_%j.out
#SBATCH --error=logs/phase4/M28/install_umap_%j.err
# Phase 4 · M28 — provide Python umap-learn for SCEVAN's subclone stage.
#
# WHY THIS IS NEEDED (root cause, SLURM 19896623):
#   SCEVAN:::subcloneAnalysisPipeline calls plotTSNE(), which calls
#   umap::umap(..., method = "umap-learn") — i.e. the PYTHON umap-learn
#   implementation via reticulate. plotTSNE() runs BEFORE the line that writes
#   subclone assignments into classDf, so its failure loses the clone labels
#   entirely. It is not an optional cosmetic step.
#
# WHY AN ISOLATED ENV:
#   Installing into R_env risks the conda solver moving R packages, and Phase 2/3
#   froze that stack (Seurat 5.4.0 / SeuratObject 5.3.0 / harmony 1.2.4 /
#   Matrix 1.7.4 / r-base 4.4.3). A separate Python-only env cannot touch them.
#   This follows the Phase 3 precedent of isolating CellPhoneDB in cpdb_env.
#   R_env is NOT modified by this job.
# NOTE: no "set -u". Conda's qt-main activate.d hook dereferences an unset
# QT_XCB_GL_INTEGRATION, so "set -u" aborts the script on conda activate/deactivate
# (root cause of SLURM 19896654). Conda's own guidance is to disable nounset.
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase4/M28 reports/phase4/environment
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh

echo "=== R_env BEFORE (must be identical after) ==="
conda activate R_env
Rscript -e 'for (p in c("Seurat","SeuratObject","harmony","Matrix","sctransform","SCEVAN","yaGST","umap","reticulate")) cat(sprintf("%-14s %s\n", p, tryCatch(as.character(packageVersion(p)), error=function(e)"NA")))'
R --version | head -1
conda deactivate

echo
echo "=== creating isolated p4_umap_env (python only) ==="
if conda env list | grep -qE "^p4_umap_env[[:space:]]"; then
  echo "p4_umap_env already exists; leaving it alone"
else
  conda create -y -n p4_umap_env -c conda-forge python=3.11 "umap-learn>=0.5" numpy scipy scikit-learn numba pynndescent
fi

conda activate p4_umap_env
python -c "import umap, numpy, sklearn, numba; print('umap-learn', umap.__version__); print('numpy', numpy.__version__); print('sklearn', sklearn.__version__); print('numba', numba.__version__)"
PY_BIN="$(which python)"
echo "PY_BIN=${PY_BIN}"
conda env export > reports/phase4/environment/p4_umap_env_PHASE4.yml
conda deactivate

echo
echo "=== R_env AFTER (verify nothing moved) ==="
conda activate R_env
Rscript -e 'for (p in c("Seurat","SeuratObject","harmony","Matrix","sctransform","SCEVAN","yaGST","umap","reticulate")) cat(sprintf("%-14s %s\n", p, tryCatch(as.character(packageVersion(p)), error=function(e)"NA")))'
R --version | head -1

echo
echo "=== reticulate can reach the isolated interpreter ==="
RETICULATE_PYTHON="${PY_BIN}" Rscript -e '
library(reticulate)
cat("python:", py_config()$python, "\n")
cat("umap module available:", py_module_available("umap"), "\n")
if (!py_module_available("umap")) quit(status = 1)
# exercise the exact call SCEVAN makes
library(umap); set.seed(1)
m <- matrix(rnorm(400 * 30), nrow = 400)
u <- umap(m, method = "umap-learn", n_components = 2, n_neighbors = 15,
          min_dist = 0.1, metric = "euclidean", seed = 1)
cat("umap-learn layout dims:", paste(dim(u$layout), collapse = " x "), "\n")
cat("SCEVAN plotTSNE dependency SATISFIED\n")'
rc=$?
echo "PY_BIN_FOR_M29=${PY_BIN}"
echo "=== install end $(date) | exit ${rc} ==="
exit $rc
