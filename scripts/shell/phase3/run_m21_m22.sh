#!/bin/bash
#SBATCH --job-name=p3_m21_m22
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=logs/phase3/M21/m21_m22_%j.out
#SBATCH --error=logs/phase3/M21/m21_m22_%j.err
set -euo pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase3/M21 logs/phase3/M22 results/phase3/figures/{M20,M21,M22}
echo "=== M21+M22 start: $(date) | Job ${SLURM_JOB_ID:-NA} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
echo "--- M21 concordance ---"
Rscript scripts/phase3/concordance/build_concordance.R \
  --by-dir results/phase3/ccc/by_method --out-dir results/phase3/ccc/concordance \
  --tab-dir results/phase3/tables --fig-dir results/phase3/figures/M21
echo "--- CCC figure suite (M20/M21) ---"
Rscript scripts/phase3/concordance/make_ccc_figures.R \
  --by-dir results/phase3/ccc/by_method \
  --concordance results/phase3/ccc/concordance/CCC_CONCORDANCE.tsv \
  --fig20 results/phase3/figures/M20 --fig21 results/phase3/figures/M21 \
  --out-dir results/phase3/ccc/concordance
echo "--- M22 prioritisation ---"
Rscript scripts/phase3/ccc/prioritize_interactions.R \
  --object results/phase3/ccc/ccc_input_object.rds \
  --concordance results/phase3/ccc/concordance/CCC_CONCORDANCE.tsv \
  --out-dir results/phase3/ccc/prioritized --tab-dir results/phase3/tables \
  --fig-dir results/phase3/figures/M22 --loch-dir results/phase3/lochness
echo "=== End: $(date) ==="
