# MPNST Phase 1 Conda Environments Specification

This directory contains specifications and logs of the Conda environment used to execute Phase 1 analysis.

## Files

1. **Portable Environment Specification**:
   - `R_env_portable.yaml`: A machine-agnostic representation of the environment, created by stripping the machine-specific `prefix:` field from the final export. This file is intended for researchers to recreate or clone the environment on other clusters or environments.

2. **Raw Provenance Snapshots**:
   - These are frozen snapshots of the exact Conda state before and after installing Snakemake. They preserve absolute local directory paths and installation paths as immutable execution history.
   - `R_env_baseline.yaml` / `R_env_baseline.explicit.txt`: State before installing Snakemake.
   - `R_env_final.yaml` / `R_env_final.explicit.txt`: State after installing Snakemake.
   - `R_env_info.txt` / `R_env_final_info.txt`: Output of `conda info` containing system architecture, UID, GID, platform, and packages cache directories.
   - `R_env_list.txt` / `R_env_final_list.txt`: Output of `conda list` showing active packages.
   - `R_env_final_sessionInfo.txt`: Session information from R showing active library namespaces and attachments.
