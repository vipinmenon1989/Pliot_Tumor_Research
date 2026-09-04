#!/bin/bash
#SBATCH --job-name=p3_m25_m26
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=logs/phase3/M25/m25_m26_%j.out
#SBATCH --error=logs/phase3/M25/m25_m26_%j.err
set -euo pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase3/M25 logs/phase3/M26 results/phase3/figures/{M25,M26}
echo "=== M25+M26 start: $(date) | Job ${SLURM_JOB_ID:-NA} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
echo "--- M25 evidence integration ---"
Rscript scripts/phase3/concordance/integrate_evidence.R \
  --master results/phase3/ccc/prioritized/MPNST_CCC_MASTER_TABLE.tsv \
  --nichenet results/phase3/receiver_response/NicheNet_receiver_response_support.tsv \
  --lochness results/phase3/lochness/lochness_summary.tsv \
  --loch-perm results/phase3/lochness/lochness_permutation_summary.tsv \
  --out-dir results/phase3/ccc/prioritized --tab-dir results/phase3/tables \
  --fig-dir results/phase3/figures/M25
echo "--- M26 robustness ---"
Rscript scripts/phase3/concordance/robustness.R \
  --master results/phase3/ccc/prioritized/MPNST_CCC_MASTER_TABLE.tsv \
  --loo-lochness results/phase3/lochness/lochness_leave_one_sample_out.tsv \
  --out-dir results/phase3/tables --fig-dir results/phase3/figures/M26
echo "=== End: $(date) ==="
