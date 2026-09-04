#!/bin/bash
#SBATCH --job-name=p2_m15_annot
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --output=logs/phase2/slurm/m15_annotation_%j.out
#SBATCH --error=logs/phase2/slurm/m15_annotation_%j.err
set -euo pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase2/slurm results/phase2/annotation results/phase2/figures/M15
echo "=== M15 start: $(date) | Job ${SLURM_JOB_ID:-NA} | $(hostname) ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
Rscript scripts/R/phase2/annotate_celltypes.R \
  --input results/phase2/clustering/phase2_harmony_clustered.rds \
  --map config/phase2/annotation_map_M15.tsv \
  --out-rds results/phase2/annotation/phase2_harmony_annotated.rds \
  --out-dir results/phase2/annotation --fig-dir results/phase2/figures/M15 \
  --random-seed 42 --expected-cells 19716
echo "=== End: $(date) ==="
