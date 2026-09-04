#!/bin/bash
#SBATCH --job-name=p2_m15a_ccc
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --output=logs/phase2/slurm/m15a_ccc_%j.out
#SBATCH --error=logs/phase2/slurm/m15a_ccc_%j.err
# Phase 2 / M15A - CCC-oriented annotation layer (approved amendment).
# Consumes the validated M16 object. Harmony, M12, M13 and M14 are NOT rerun.
set -euo pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase2/slurm results/phase2/annotation results/phase2/figures/CCC_annotation
echo "=== M15A start: $(date) | Job ${SLURM_JOB_ID:-NA} | $(hostname) | CPUs ${SLURM_CPUS_PER_TASK:-1} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
R --version | head -1
Rscript -e 'cat(sprintf("Seurat %s | SeuratObject %s | ggplot2 %s | reshape2 %s\n", packageVersion("Seurat"), packageVersion("SeuratObject"), packageVersion("ggplot2"), packageVersion("reshape2")))'
Rscript scripts/R/phase2/build_ccc_annotation.R \
  --input results/phase2/composition/phase2_harmony_refined.rds \
  --map config/phase2/ccc_annotation_map.tsv \
  --out-rds results/phase2/annotation/phase2_harmony_ccc.rds \
  --out-dir results/phase2/annotation \
  --comp-dir results/phase2/composition \
  --fig-dir results/phase2/figures/CCC_annotation \
  --random-seed 42 --expected-cells 19716
echo "=== End: $(date) ==="
