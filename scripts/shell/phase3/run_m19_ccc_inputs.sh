#!/bin/bash
#SBATCH --job-name=p3_m19_inputs
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=06:00:00
#SBATCH --output=logs/phase3/M19/m19_inputs_%j.out
#SBATCH --error=logs/phase3/M19/m19_inputs_%j.err
# Phase 3 / M19 - sample-aware CCC input preparation.
# Memory headroom is for JoinLayers + LogNormalize on the full 31,764 x 19,716 RNA assay
# (Phase 2 M14 peaked at 16.4 GiB on a comparable operation).
set -euo pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase3/M19 results/phase3/figures/M19 results/phase3/tables
echo "=== M19 start: $(date) | Job ${SLURM_JOB_ID:-NA} | $(hostname) | CPUs ${SLURM_CPUS_PER_TASK:-1} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
R --version | head -1
Rscript scripts/phase3/ccc/prepare_ccc_inputs.R \
  --input results/phase2/phase2_final_object.rds \
  --out-dir results/phase3/ccc \
  --tab-dir results/phase3/tables \
  --fig-dir results/phase3/figures/M19 \
  --cpdb-dir results/phase3/ccc/cellphonedb_inputs \
  --min-cells 10 --random-seed 42
echo "=== End: $(date) ==="
