#!/bin/bash
#SBATCH --job-name=p2_m12_suppl
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=logs/phase2/slurm/m12_supplement_%j.out
#SBATCH --error=logs/phase2/slurm/m12_supplement_%j.err
set -euo pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase2/slurm results/phase2/figures/M12
echo "=== M12 supplement start: $(date) | Job ${SLURM_JOB_ID:-NA} | $(hostname) ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
Rscript scripts/R/phase2/m12_supplement_figures.R
echo "=== End: $(date) ==="
