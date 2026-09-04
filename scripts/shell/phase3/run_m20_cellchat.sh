#!/bin/bash
#SBATCH --job-name=p3_m20_cchat
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=logs/phase3/M20/cellchat_%j.out
#SBATCH --error=logs/phase3/M20/cellchat_%j.err
set -euo pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase3/M20 results/phase3/ccc/by_method
echo "=== M20 CellChat start: $(date) | Job ${SLURM_JOB_ID:-NA} | $(hostname) | CPUs ${SLURM_CPUS_PER_TASK:-1} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
Rscript -e 'cat(sprintf("CellChat %s | Seurat %s | NMF %s\n", packageVersion("CellChat"), packageVersion("Seurat"), packageVersion("NMF")))'
Rscript scripts/phase3/ccc/run_cellchat.R \
  --input results/phase3/ccc/ccc_input_object.rds \
  --out-dir results/phase3/ccc/by_method \
  --min-cells 10 --random-seed 42
echo "=== End: $(date) ==="
