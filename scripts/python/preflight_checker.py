#!/usr/bin/env python3
"""
Preflight Checker for MPNST Phase 1 Single-Cell Analysis Workflow
Verifies directories, configurations, python and R package availability, and system dependencies.
"""

import os
import sys
import argparse
import subprocess
import yaml
from jsonschema import validate, ValidationError

def log_success(message):
    print(f"\033[92m[PASS]\033[0m {message}")

def log_failure(message, fix_hint=None):
    print(f"\033[91m[FAIL]\033[0m {message}")
    if fix_hint:
        print(f"       \033[93mHINT:\033[0m {fix_hint}")
    print()

def log_info(message):
    print(f"\033[94m[INFO]\033[0m {message}")

def run_command(cmd):
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
        return True, res.stdout.strip()
    except Exception as e:
        return False, str(e)

def main():
    parser = argparse.ArgumentParser(description="Preflight Checker for MPNST Phase 1 Workflow")
    parser.add_argument(
        "--mode",
        choices=["synthetic", "real"],
        required=True,
        help="Execution mode: 'synthetic' (test/CI execution on synthetic data) or 'real' (execution on real patient tumor data)."
    )
    parser.add_argument(
        "--config",
        help="Path to alternative configuration file. If not specified, defaults based on mode."
    )
    parser.add_argument(
        "--input-rds",
        help="CLI override path to the input RDS file. Overrides the config file setting."
    )
    parser.add_argument(
        "--require-slurm",
        action="store_true",
        help="If set, requires the SLURM scheduler (sbatch) to be available (fails otherwise)."
    )
    parser.add_argument(
        "--skip-r-packages",
        action="store_true",
        help="Skip the slow check for critical R packages (Seurat, sctransform, scDblFinder)."
    )
    
    args = parser.parse_args()
    mode = args.mode
    
    # 1. Establish default config path based on mode
    if args.config:
        config_path = args.config
    else:
        if mode == "synthetic":
            config_path = "config/config.test.yaml"
        else:
            config_path = "config/config.yaml"

    print(f"=== Starting Preflight Validation Checker (Mode: {mode.upper()}) ===")
    print(f"Config path: {config_path}")
    errors = 0

    # 2. Check directories
    required_dirs = ["workflow", "config", "scripts", "tests"]
    for d in required_dirs:
        if os.path.isdir(d):
            log_success(f"Directory '{d}/' exists.")
        else:
            log_failure(f"Directory '{d}/' is missing.", f"Create the directory '{d}/' or restore it from the repository.")
            errors += 1

    # 3. Check Python version and dependencies
    py_ok, py_ver = run_command([sys.executable, "--version"])
    if py_ok:
        log_success(f"Python is available: {py_ver}")
    else:
        log_failure("Python is not available or executable.", "Activate the correct Conda environment (e.g. conda activate R_env).")
        errors += 1

    # 4. Check R availability
    r_ok, r_ver = run_command(["Rscript", "--version"])
    if r_ok:
        log_success(f"Rscript is available: {r_ver}")
    else:
        log_failure("Rscript is not available in the current PATH.", "Ensure R is installed and Rscript is in your PATH, or activate R_env.")
        errors += 1

    # 5. Check Snakemake availability
    sm_ok, sm_ver = run_command(["snakemake", "--version"])
    if sm_ok:
        log_success(f"Snakemake is available: version {sm_ver}")
    else:
        log_failure("Snakemake is not available in the current PATH.", "Install snakemake (pip install snakemake or conda install snakemake).")
        errors += 1

    # 6. Check critical R packages (without loading the large RDS file)
    if args.skip_r_packages:
        log_info("Skipping critical R packages check (--skip-r-packages).")
    else:
        r_pkg_cmd = ["Rscript", "-e", "libs <- c('Seurat', 'sctransform', 'scDblFinder'); for (l in libs) { library(l, character.only=TRUE) }; cat('OK')"]
        r_pkg_ok, r_pkg_out = run_command(r_pkg_cmd)
        if r_pkg_ok and "OK" in r_pkg_out:
            log_success("Critical R packages are available: Seurat, sctransform, scDblFinder.")
        else:
            log_failure("Failed to load critical R packages in R.", f"Install packages Seurat, sctransform, or scDblFinder. Detailed error: {r_pkg_out}")
            errors += 1

    # 7. Read and validate configuration file against schema
    schema_path = "config/schemas/config.schema.yaml"
    resolved_rds_path = None
    if os.path.exists(config_path) and os.path.exists(schema_path):
        try:
            with open(config_path, "r") as cf:
                config_data = yaml.safe_load(cf)
            with open(schema_path, "r") as sf:
                schema_data = yaml.safe_load(sf)
            validate(instance=config_data, schema=schema_data)
            log_success(f"Configuration file {config_path} is valid against schemas.")
            
            # 8. Path resolution and precedence
            if args.input_rds:
                resolved_rds_path = args.input_rds
                log_info(f"Resolved input path (via CLI override): '{resolved_rds_path}'")
            else:
                resolved_rds_path = config_data.get("input_rds")
                log_info(f"Resolved input path (via config): '{resolved_rds_path}'")
            
            # Check if input rds exists (file existence check only)
            if resolved_rds_path:
                if os.path.exists(resolved_rds_path):
                    log_success(f"Input dataset exists: '{resolved_rds_path}'.")
                else:
                    if mode == "real":
                        log_failure(f"Input immutable dataset not found at: '{resolved_rds_path}'.", "Verify that the processed_mpnst.rds file is located at the specified path.")
                        errors += 1
                    else:
                        log_info(f"Synthetic/test input dataset not found at: '{resolved_rds_path}'. (This is expected as it will be generated dynamically by the Snakemake workflow).")
            else:
                log_failure(f"Configuration file {config_path} is missing 'input_rds' field.")
                errors += 1

        except ValidationError as ve:
            log_failure(f"Configuration file validation failed: {ve.message}", f"Correct {config_path} format.")
            errors += 1
        except Exception as e:
            log_failure(f"Error reading configuration/schema: {str(e)}")
            errors += 1
    else:
        log_failure(f"Missing configuration file at {config_path} or schema at {schema_path}.")
        errors += 1

    # 9. Check write permissions in outputs
    out_dirs = ["results", "reports", "logs"]
    for od in out_dirs:
        if os.path.exists(od):
            if os.access(od, os.W_OK):
                log_success(f"Write permissions verified for directory '{od}/'.")
            else:
                log_failure(f"No write permissions in directory '{od}/'.", f"Run chmod or change ownership to get write access to '{od}/'.")
                errors += 1
        else:
            # Check if parent is writable so Snakemake can create it
            if os.access(".", os.W_OK):
                log_success(f"Directory '{od}/' does not exist but parent directory is writable (Snakemake will create it).")
            else:
                log_failure(f"Directory '{od}/' does not exist and parent directory is not writable.", "Change write permissions of the project root.")
                errors += 1

    # 10. Check SLURM commands if available
    sbatch_ok, sbatch_ver = run_command(["sbatch", "--version"])
    if sbatch_ok:
        log_success(f"SLURM scheduler is available: {sbatch_ver}")
    else:
        if args.require_slurm:
            log_failure("SLURM scheduler (sbatch) is missing, and explicitly required (--require-slurm).", "Run in an environment/node with SLURM cluster access.")
            errors += 1
        else:
            log_info("SLURM is not available locally. (This is fine since --require-slurm was not specified).")

    print("=== Preflight Validation Summary ===")
    if errors == 0:
        print("\033[92mAll preflight checks PASSED.\033[0m Project is ready for execution.")
        sys.exit(0)
    else:
        print(f"\033[91mFailed {errors} preflight check(s).\033[0m Resolve the issues before proceeding.")
        sys.exit(1)

if __name__ == "__main__":
    main()
