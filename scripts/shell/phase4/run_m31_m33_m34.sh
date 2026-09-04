#!/bin/bash
#SBATCH --job-name=p4_m31_34
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/phase4/M31/m31_m33_m34_%j.out
#SBATCH --error=logs/phase4/M31/m31_m33_m34_%j.err
# Phase 4 · M31 (tumour states) -> M33 (state -> TME model) -> M34 (robustness).
# Chained in one job because M33 consumes M31's object and M34 consumes both plus
# the M32 sensitivity table.
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase4/M31
echo "=== M31/M33/M34 start $(date) | Job ${SLURM_JOB_ID} | ${SLURM_CPUS_PER_TASK} cpus ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
export RETICULATE_PYTHON=/local/projects-t3/lilab/vmenon/anaconda3/envs/p4_umap_env/bin/python
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo; echo "=== M31 — malignant-only tumour states ==="
Rscript scripts/phase4/tumor_states/m31_tumor_states.R;      rc1=$?
echo; echo "=== M33 — tumour-state -> TME evidence model ==="
Rscript scripts/phase4/tumor_states/m33_state_tme_model.R;   rc2=$?
echo; echo "=== M34 — robustness ==="
Rscript scripts/phase4/malignancy/m34_robustness.R;          rc3=$?
echo "=== end $(date) | m31=${rc1} m33=${rc2} m34=${rc3} ==="
[ $rc1 -ne 0 ] && exit $rc1
[ $rc2 -ne 0 ] && exit $rc2
exit $rc3
