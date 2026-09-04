#!/bin/bash
#SBATCH --job-name=p2_m17_freeze
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --output=logs/phase2/slurm/m17_freeze_%j.out
#SBATCH --error=logs/phase2/slurm/m17_freeze_%j.err
set -euo pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase2/slurm results/phase2/tables/final
echo "=== M17 start: $(date) | Job ${SLURM_JOB_ID:-NA} | $(hostname) ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
R --version | head -1
Rscript scripts/R/phase2/validate_and_freeze.R \
  --input results/phase2/annotation/phase2_harmony_ccc.rds \
  --final-rds results/phase2/phase2_final_object.rds \
  --tables-dir results/phase2/tables/final \
  --manifest results/phase2/phase2_manifest.json \
  --expected-cells 19716
echo "=== End: $(date) ==="
