#!/bin/bash
#SBATCH --job-name=p5_m38b
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/phase5/M38/m38b_%j.out
#SBATCH --error=logs/phase5/M38/m38b_%j.err
# Phase 5 - M38 re-run with explicit technical-content detection in the program
# labels. Nothing about the factorization, K selection or recurrence rule
# changes - only how a program is NAMED once its own top genes are inspected.
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
K=$(cat results/phase5/validation/selected_K.txt)
echo "=== M38 relabel start $(date) | Job ${SLURM_JOB_ID} | K=${K} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
Rscript scripts/phase5/programs/m38_annotate_programs.R "${K}" 0_1; rc=$?
echo "=== end $(date) | rc=${rc} ==="
exit $rc
