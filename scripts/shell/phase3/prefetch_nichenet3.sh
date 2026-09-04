#!/bin/bash
#SBATCH --job-name=p3_nnfetch3
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=logs/phase3/M23/nnfetch3_%j.out
#SBATCH --error=logs/phase3/M23/nnfetch3_%j.err
# Third attempt at the NicheNet ligand-target matrix.
# ROOT CAUSE of attempt 2: the transfer was truncated at 112.7 MB of the expected 262.1 MB.
# FIX: resumable download (curl -C -), retry on transient errors, and an explicit size check
# against the Zenodo-reported 262116605 bytes before accepting the file.
set -uo pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
F=external/nichenet/ligand_target_matrix.rds
URL="https://zenodo.org/api/records/7074291/files/ligand_target_matrix_nsga2r_final.rds/content"
EXPECT=262116605
echo "=== attempt 3 start: $(date) | Job ${SLURM_JOB_ID:-NA} ==="
for try in 1 2 3 4; do
  SZ=$(stat -c%s "$F" 2>/dev/null || echo 0)
  if [ "$SZ" -eq "$EXPECT" ]; then echo "size OK ($SZ) after $((try-1)) retries"; break; fi
  echo "try $try: current size $SZ / $EXPECT — resuming"
  curl -sS -L -C - --retry 5 --retry-delay 5 --retry-all-errors \
       --speed-time 120 --speed-limit 10000 --max-time 3000 -o "$F" "$URL" || echo "  curl exit $?"
done
SZ=$(stat -c%s "$F" 2>/dev/null || echo 0)
echo "final size: $SZ (expected $EXPECT)"
[ "$SZ" -eq "$EXPECT" ] || echo "WARNING: size still mismatched"
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
Rscript -e '
for (f in list.files("external/nichenet", full.names=TRUE, pattern="\\.rds$")) {
  x <- tryCatch(readRDS(f), error=function(e) NULL)
  cat(sprintf("%-34s %8.1f MB  %-16s md5 %s\n", basename(f), file.info(f)$size/1e6,
    if (is.null(x)) "UNREADABLE" else paste(dim(x) %||% length(x), collapse="x"), tools::md5sum(f)))
}'
echo "=== End: $(date) ==="
