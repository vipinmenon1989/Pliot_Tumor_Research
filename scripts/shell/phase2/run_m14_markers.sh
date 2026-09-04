#!/bin/bash
#SBATCH --job-name=p2_m14_markers
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=12
#SBATCH --mem=250G
#SBATCH --time=12:00:00
#SBATCH --output=logs/phase2/slurm/m14_markers_%j.out
#SBATCH --error=logs/phase2/slurm/m14_markers_%j.err
# Phase 2 / M14 - cluster marker discovery. Memory headroom is for
# PrepSCTFindMarkers(), which recorrects counts across the four SCT models.
set -euo pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase2/slurm results/phase2/markers results/phase2/figures/M14
echo "=== M14 start: $(date) | Job ${SLURM_JOB_ID:-NA} | $(hostname) | CPUs ${SLURM_CPUS_PER_TASK:-1} | Mem 250G ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
R --version | head -1
Rscript -e 'cat(sprintf("Seurat %s | SeuratObject %s | presto %s | ggplot2 %s\n", packageVersion("Seurat"), packageVersion("SeuratObject"), packageVersion("presto"), packageVersion("ggplot2")))'
Rscript scripts/R/phase2/discover_markers_harmony.R \
  --input results/phase2/clustering/phase2_harmony_clustered.rds \
  --out-dir results/phase2/markers \
  --fig-dir results/phase2/figures/M14 \
  --assay SCT --min-pct 0.25 --logfc 0.25 --random-seed 42 --expected-cells 19716
echo "=== End: $(date) ==="
