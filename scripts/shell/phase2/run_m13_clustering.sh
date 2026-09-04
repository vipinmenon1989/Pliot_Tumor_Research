#!/bin/bash
#SBATCH --job-name=p2_m13_cluster
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=logs/phase2/slurm/m13_clustering_%j.out
#SBATCH --error=logs/phase2/slurm/m13_clustering_%j.err
# Phase 2 / M13 - post-Harmony neighbours, UMAP and Louvain resolution sweep 0.1-1.0.
# Resources from M12 accounting (8.14 GiB) plus the SNN graph, 10 silhouette
# computations over one 19,716 x 19,716 distance matrix, and a ~6 GB object save.
set -euo pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase2/slurm results/phase2/clustering results/phase2/figures/M13
echo "=== M13 start: $(date) | Job ${SLURM_JOB_ID:-NA} | $(hostname) | CPUs ${SLURM_CPUS_PER_TASK:-1} | Mem 128G ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
echo "which R: $(which R)"; R --version | head -1
Rscript -e 'cat(sprintf("Seurat %s | SeuratObject %s | harmony %s | aricode %s | cluster %s | ggplot2 %s\n", packageVersion("Seurat"), packageVersion("SeuratObject"), packageVersion("harmony"), packageVersion("aricode"), packageVersion("cluster"), packageVersion("ggplot2")))'
Rscript scripts/R/phase2/cluster_sweep_harmony.R \
  --input results/phase2/harmony/phase2_harmony_integrated.rds \
  --out-rds results/phase2/clustering/phase2_harmony_clustered.rds \
  --out-dir results/phase2/clustering \
  --fig-dir results/phase2/figures/M13 \
  --dims 30 --k-param 20 --random-seed 42 --expected-cells 19716
echo "=== End: $(date) ==="
