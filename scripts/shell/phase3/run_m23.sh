#!/bin/bash
#SBATCH --job-name=p3_m23_nn
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=10:00:00
#SBATCH --output=logs/phase3/M23/m23_%j.out
#SBATCH --error=logs/phase3/M23/m23_%j.err
set -euo pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase3/M23 results/phase3/figures/M23
echo "=== M23 NicheNet start: $(date) | Job ${SLURM_JOB_ID:-NA} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
Rscript -e 'cat(sprintf("nichenetr %s | Seurat %s\n", packageVersion("nichenetr"), packageVersion("Seurat")))'
Rscript scripts/phase3/receiver_response/run_nichenet.R \
  --object results/phase3/ccc/ccc_input_object.rds \
  --master results/phase3/ccc/prioritized/MPNST_CCC_MASTER_TABLE.tsv \
  --nn-dir external/nichenet --out-dir results/phase3/receiver_response \
  --fig-dir results/phase3/figures/M23 \
  --receivers Macrophage,CD8-T,NK,Fibroblast,Endothelial,Monocyte,Dendritic,CD4-T \
  --random-seed 42
echo "=== End: $(date) ==="
