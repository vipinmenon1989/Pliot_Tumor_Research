#!/bin/bash
#SBATCH --job-name=p3_nnfetch
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/phase3/M23/nnfetch_%j.out
#SBATCH --error=logs/phase3/M23/nnfetch_%j.err
# Prefetch the NicheNet prior model from Zenodo so it is off the M23 critical path.
set -uo pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
mkdir -p external/nichenet logs/phase3/M23
echo "=== NicheNet prefetch start: $(date) | Job ${SLURM_JOB_ID:-NA} ==="
source /local/projects-t3/lilab/vmenon/anaconda3/etc/profile.d/conda.sh
conda activate R_env
Rscript -e '
d <- "external/nichenet"
urls <- list(
  ligand_target_matrix = "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds",
  lr_network           = "https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds",
  weighted_networks    = "https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds")
for (nm in names(urls)) {
  fp <- file.path(d, paste0(nm,".rds"))
  if (file.exists(fp)) { cat(sprintf("%-22s already present (%.1f MB)\n", nm, file.info(fp)$size/1e6)); next }
  cat(sprintf("downloading %s ...\n", nm))
  ok <- tryCatch({ download.file(urls[[nm]], fp, mode="wb", quiet=TRUE); TRUE },
                 error=function(e) { cat("  FAILED:", conditionMessage(e), "\n"); FALSE })
  if (ok) {
    x <- tryCatch(readRDS(fp), error=function(e) NULL)
    cat(sprintf("%-22s %.1f MB | %s | md5 %s\n", nm, file.info(fp)$size/1e6,
      if (is.null(x)) "UNREADABLE" else paste(dim(x) %||% length(x), collapse="x"),
      tools::md5sum(fp)))
  }
}'
echo "=== End: $(date) ==="
