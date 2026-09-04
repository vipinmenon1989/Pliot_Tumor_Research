#!/bin/bash
#SBATCH --job-name=p3_m20_liana
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=logs/phase3/M20/liana_%j.out
#SBATCH --error=logs/phase3/M20/liana_%j.err
set -euo pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase3/M20 results/phase3/ccc/by_method
echo "=== M20 LIANA start: $(date) | Job ${SLURM_JOB_ID:-NA} | $(hostname) ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
Rscript -e 'cat(sprintf("liana %s | OmnipathR %s | decoupleR %s | Seurat %s\n", packageVersion("liana"), packageVersion("OmnipathR"), packageVersion("decoupleR"), packageVersion("Seurat")))'
Rscript scripts/phase3/ccc/run_liana.R \
  --input results/phase3/ccc/ccc_input_object.rds \
  --out-dir results/phase3/ccc/by_method \
  --min-cells 10 --resource Consensus --random-seed 42
echo "=== End: $(date) ==="
