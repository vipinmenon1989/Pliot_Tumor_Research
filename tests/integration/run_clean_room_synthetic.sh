#!/bin/bash
#SBATCH --job-name=clean_room_synth
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=logs/slurm/clean_room_synth_%j.out
#SBATCH --error=logs/slurm/clean_room_synth_%j.err

set -e

echo "=== Clean-Room Synthetic Workflow Execution ==="
date
hostname

# 1. Clean previous clean-room runs
rm -rf test_clean_room
mkdir -p test_clean_room/reports/milestones
mkdir -p test_clean_room/reports/audits

# 2. Symlink required directories and files
ln -sfn ../scripts test_clean_room/scripts
ln -sfn ../tests test_clean_room/tests
ln -sfn ../config test_clean_room/config
ln -sfn ../workflow test_clean_room/workflow
ln -sfn ../PROJECT.md test_clean_room/PROJECT.md
ln -sfn ../README.md test_clean_room/README.md
ln -sfn ../PROGRESS.md test_clean_room/PROGRESS.md
ln -sfn ../CHANGELOG.md test_clean_room/CHANGELOG.md

# 3. Copy static report files needed for validation
cp reports/milestones/M0_REPORT.md test_clean_room/reports/milestones/
cp reports/milestones/M1_REPORT.md test_clean_room/reports/milestones/
cp reports/milestones/M2_REPORT.md test_clean_room/reports/milestones/
cp reports/milestones/M3_REPORT.md test_clean_room/reports/milestones/
cp reports/milestones/M9_REPORT.md test_clean_room/reports/milestones/
cp reports/audits/M0_M8_RECONCILIATION.md test_clean_room/reports/audits/
cp reports/audits/M0_M8_RECONCILIATION.tsv test_clean_room/reports/audits/
cp reports/PHASE1_HANDOFF.md test_clean_room/reports/PHASE1_HANDOFF.md

# 4. Load conda environment
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env

# 5. Run Snakemake workflow
cd test_clean_room
snakemake -s workflow/Snakefile \
  --configfile config/config.test.yaml \
  --cores 4 \
  synthetic_complete \
  --rerun-triggers mtime

echo "=== Clean-Room Synthetic Workflow Finished Successfully ==="
date
