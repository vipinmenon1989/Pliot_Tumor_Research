#!/bin/bash
#SBATCH --job-name=p5_m38m39
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/phase5/M38/m38m39_%j.out
#SBATCH --error=logs/phase5/M38/m38m39_%j.err
# =============================================================================
# Phase 5 - M37d (K selection) + M38 (program annotation) + M39 (malignant ECM).
#
# 64 G: the largest object touched is the 31,764 x 7,462 working count matrix
# and its dense TPM projection (~1.8 GB); fgsea over Hallmark + Reactome for
# each program is the main CPU cost. Evidence-based, not the 450 G maximum.
# =============================================================================
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase5/M38
echo "=== M37d/M38/M39 start $(date) | Job ${SLURM_JOB_ID} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh

echo; echo "=== M37d: measure every K and apply the declared selection rule ==="
conda activate p5_cnmf_env
python scripts/phase5/programs/m37_k_selection.py; rc=$?
conda deactivate
[ $rc -ne 0 ] && { echo "M37d failed rc=$rc"; exit $rc; }

K=$(python3 -c "import json;print(json.load(open('results/phase5/validation/m37d_k_selection_decision.json'))['decision']['selected_K'])")
echo "SELECTED K = ${K}"
[ -z "$K" ] || [ "$K" = "None" ] && { echo "no K selected - STOP"; exit 1; }
echo "$K" > results/phase5/validation/selected_K.txt

conda activate R_env
echo; echo "=== M38: program annotation, pathways, Phase 4 state comparison ==="
Rscript scripts/phase5/programs/m38_annotate_programs.R "${K}" 0_1; rc1=$?
echo; echo "=== M39: malignant ECM-like cells vs true fibroblasts ==="
Rscript scripts/phase5/malignant_ecm/m39_malignant_ecm.R "${K}" 0_1; rc2=$?
echo "=== end $(date) | m38=${rc1} m39=${rc2} ==="
[ $rc1 -ne 0 ] && exit $rc1
exit $rc2
