#!/bin/bash
#SBATCH --job-name=p3_netprobe
#SBATCH --account=ihc
#SBATCH --partition=ihc
#SBATCH --nodelist=ihc-grid-1-1-1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/phase3/M18/netprobe_%j.out
#SBATCH --error=logs/phase3/M18/netprobe_%j.err
set -uo pipefail
cd /local/projects-t3/lilab/vmenon/Tumor_research/Pilot_MPNST
echo "host: $(hostname)  date: $(date)"
for url in https://cran.r-project.org https://api.github.com https://bioconductor.org https://pypi.org https://zenodo.org; do
  printf "%-34s " "$url"
  timeout 15 curl -sS -o /dev/null -w "HTTP %{http_code}\n" "$url" 2>&1 | tail -1 || echo "UNREACHABLE"
done
