#!/bin/bash
#SBATCH --job-name=mpnst_smoke
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/slurm/smoke_job_%j.out
#SBATCH --error=logs/slurm/smoke_job_%j.err

set -e

echo "=== SLURM Job Start ==="
date
hostname

echo "=== Loading Conda Environment ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env

echo "=== Verifying Paths and Versions ==="
which R
R --version
which python
python --version
which snakemake
snakemake --version

echo "=== Running R Library Verification ==="
Rscript -e '
library(Seurat)
library(SeuratObject)
cat("Seurat and SeuratObject loaded successfully on compute node!\n")
'

echo "=== SLURM Job Finished Successfully ==="
date
