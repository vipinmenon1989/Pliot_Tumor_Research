#!/bin/bash
#SBATCH --job-name=p3_nnfetch2
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=logs/phase3/M23/nnfetch2_%j.out
#SBATCH --error=logs/phase3/M23/nnfetch2_%j.err
# Retry the NicheNet ligand-target matrix with curl. ROOT CAUSE of the first failure:
# R download.file()'s default 60 s timeout against a 262 MB Zenodo file.
set -uo pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p external/nichenet
F=external/nichenet/ligand_target_matrix.rds
echo "=== retry start: $(date) | Job ${SLURM_JOB_ID:-NA} ==="
if [ ! -s "$F" ] || [ "$(stat -c%s "$F")" -lt 100000000 ]; then
  curl -sS -L --retry 3 --max-time 3600 \
    -o "$F" "https://zenodo.org/api/records/7074291/files/ligand_target_matrix_nsga2r_final.rds/content"
fi
ls -la external/nichenet/
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
Rscript -e '
for (f in list.files("external/nichenet", full.names=TRUE)) {
  x <- tryCatch(readRDS(f), error=function(e) NULL)
  cat(sprintf("%-40s %8.1f MB  %-14s md5 %s\n", basename(f), file.info(f)$size/1e6,
    if (is.null(x)) "UNREADABLE" else paste(dim(x) %||% length(x), collapse="x"), tools::md5sum(f)))
}'
echo "=== End: $(date) ==="
