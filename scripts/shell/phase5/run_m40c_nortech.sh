#!/bin/bash
#SBATCH --job-name=p5_m40c
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=logs/phase5/M40/m40c_%j.out
#SBATCH --error=logs/phase5/M40/m40c_%j.err
# Phase 5 - M40c - full K grid on the technical-gene-free universe, so this run
# selects its own K by the same declared rule rather than inheriting K = 8.
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase5/M40
echo "=== M40c start $(date) | Job ${SLURM_JOB_ID} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
Rscript scripts/phase5/validation/m40_prepare_nortech.R || exit 1
conda deactivate
conda activate p5_cnmf_env
OUTDIR=results/phase5/programs/cnmf; IN=results/phase5/programs/cnmf_input
KGRID="4 5 6 7 8 9 10 11 12 13 14 15"; W=8
python scripts/phase5/programs/m37_build_anndata.py "${IN}/nortech" \
    "${IN}/nortech/counts.h5ad" "${IN}/malignant_cell_covariates.tsv" || exit 1
cnmf prepare --output-dir "${OUTDIR}" --name nortech -c "${IN}/nortech/counts.h5ad" \
    -k ${KGRID} --n-iter 100 --seed 42 --numgenes 2000 --total-workers ${W} || exit 1
for w in $(seq 0 $((W-1))); do
  cnmf factorize --output-dir "${OUTDIR}" --name nortech --worker-index ${w} --total-workers ${W} &
done
wait
cnmf combine --output-dir "${OUTDIR}" --name nortech || exit 1
cnmf k_selection_plot --output-dir "${OUTDIR}" --name nortech || exit 1
for K in ${KGRID}; do
  cnmf consensus --output-dir "${OUTDIR}" --name nortech --components ${K} \
      --local-density-threshold 0.10 || exit 1
done
echo "=== M40c end $(date) ==="
