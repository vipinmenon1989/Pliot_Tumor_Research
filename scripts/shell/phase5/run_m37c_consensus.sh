#!/bin/bash
#SBATCH --job-name=p5_m37c
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/phase5/M37/m37c_%j.out
#SBATCH --error=logs/phase5/M37/m37c_%j.err
# =============================================================================
# Phase 5 - M37c - consensus spectra and usages for EVERY K in the grid.
#
# Consensus is run across the whole grid, not only at a pre-chosen K, so the
# K-selection table in M37d is built from measured properties of each solution
# (stability, error, duplication, uniqueness, coherence, patient
# representation) rather than from one solution inspected after the fact.
#
# local-density-threshold 0.10 is the cNMF-standard outlier filter and is
# declared here; dt = 2.00 (no filtering) is also produced for comparison.
# 32 G / 4 CPUs: consensus is a clustering of 100 x K spectra vectors.
# =============================================================================
set -o pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p logs/phase5/M37
echo "=== M37c start $(date) | Job ${SLURM_JOB_ID} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate p5_cnmf_env

OUTDIR=results/phase5/programs/cnmf
KGRID="4 5 6 7 8 9 10 11 12 13 14 15"
rc=0
for RUN in primary balanced nocc; do
  for K in ${KGRID}; do
    echo "--- ${RUN} K=${K} dt=0.10 ---"
    cnmf consensus --output-dir "${OUTDIR}" --name "${RUN}" --components ${K} \
        --local-density-threshold 0.10 --show-clustering || rc=$?
  done
  echo "--- ${RUN} K grid, dt=2.00 (unfiltered, for comparison) ---"
  for K in ${KGRID}; do
    cnmf consensus --output-dir "${OUTDIR}" --name "${RUN}" --components ${K} \
        --local-density-threshold 2.00 || rc=$?
  done
done
echo "=== M37c end $(date) | rc=${rc} ==="
exit $rc
