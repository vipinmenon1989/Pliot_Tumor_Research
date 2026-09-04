#!/bin/bash
#SBATCH --job-name=p2_m16_comp
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --output=logs/phase2/slurm/m16_composition_%j.out
#SBATCH --error=logs/phase2/slurm/m16_composition_%j.err
set -euo pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase2/slurm results/phase2/composition results/phase2/figures/M16
echo "=== M16 start: $(date) | Job ${SLURM_JOB_ID:-NA} | $(hostname) ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
Rscript scripts/R/phase2/refine_and_compose.R \
  --input results/phase2/annotation/phase2_harmony_annotated.rds \
  --out-rds results/phase2/composition/phase2_harmony_refined.rds \
  --out-dir results/phase2/composition --fig-dir results/phase2/figures/M16 \
  --random-seed 42 --expected-cells 19716
echo "=== End: $(date) ==="
