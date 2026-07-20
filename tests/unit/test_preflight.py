#!/usr/bin/env python3
import subprocess
import os
import sys
import yaml
import tempfile

def run_preflight(args):
    # Always append --skip-r-packages to make tests extremely fast
    cmd = [sys.executable, "scripts/python/preflight_checker.py"] + args + ["--skip-r-packages"]
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return res.returncode, res.stdout, res.stderr

def print_failure_details(test_name, cmd_args, ret, stdout, stderr, extra=""):
    print(f"\n\033[91m========================================\033[0m")
    print(f"\033[91mFAIL: {test_name}\033[0m")
    print(f"Command executed: python scripts/python/preflight_checker.py {' '.join(cmd_args)}")
    print(f"Exit code: {ret}")
    if extra:
        print(f"Extra info: {extra}")
    print("--- STDOUT ---")
    print(stdout)
    print("--- STDERR ---")
    print(stderr)
    print(f"\033[91m========================================\033[0m\n")

def main():
    print("=== Running Comprehensive Preflight Checker Tests ===")
    
    # Check if sbatch is available on this machine
    sbatch_avail = subprocess.run(["which", "sbatch"], stdout=subprocess.PIPE, stderr=subprocess.PIPE).returncode == 0
    print(f"SLURM sbatch available on host: {sbatch_avail}")

    # Create a temporary directory for test config and mock files
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create a dummy config file
        temp_config_path = os.path.join(tmpdir, "config.yaml")
        with open("config/config.yaml", "r") as f:
            config_data = yaml.safe_load(f)
            
        # Create a dummy dataset file (non-RDS text file)
        dummy_rds_path = os.path.join(tmpdir, "dummy_mpnst.rds")
        with open(dummy_rds_path, "w") as f:
            f.write("This is a dummy text file to check that readRDS is NOT called during preflight.")
            
        # Write dummy config with real path configured
        config_data["input_rds"] = dummy_rds_path
        with open(temp_config_path, "w") as f:
            yaml.safe_dump(config_data, f)

        # TEST 1: Synthetic mode (No real RDS, No SLURM requirement) -> Expected: PASS
        ret, stdout, stderr = run_preflight(["--mode", "synthetic"])
        if ret == 0:
            print("[PASS] TEST 1: Synthetic mode passes.")
        else:
            print_failure_details("TEST 1: Synthetic mode should pass", ["--mode", "synthetic"], ret, stdout, stderr)
            sys.exit(1)
            
        # TEST 2: Real mode with configured RDS absent -> Expected: FAIL
        absent_config_path = os.path.join(tmpdir, "config_absent.yaml")
        config_data_absent = config_data.copy()
        config_data_absent["input_rds"] = os.path.join(tmpdir, "nonexistent_processed_mpnst.rds")
        with open(absent_config_path, "w") as f:
            yaml.safe_dump(config_data_absent, f)
            
        ret, stdout, stderr = run_preflight(["--mode", "real", "--config", absent_config_path])
        if ret != 0:
            print("[PASS] TEST 2: Real mode fails as expected when input RDS is absent.")
        else:
            print_failure_details("TEST 2: Real mode should fail when input RDS is absent", ["--mode", "real", "--config", absent_config_path], ret, stdout, stderr)
            sys.exit(1)
            
        # TEST 3: Real mode with temporary configured RDS exists and is readable, no SLURM requirement -> Expected: PASS
        # (This also checks TEST 5 implicitly: if it called readRDS, it would crash because dummy_rds_path is a text file!)
        ret, stdout, stderr = run_preflight(["--mode", "real", "--config", temp_config_path])
        if ret == 0:
            print("[PASS] TEST 3: Real mode passes when input RDS exists.")
        else:
            print_failure_details("TEST 3: Real mode should pass when input RDS exists", ["--mode", "real", "--config", temp_config_path], ret, stdout, stderr)
            sys.exit(1)
            
        # TEST 4: HPC/SLURM-required mode -> Expected behavior based on sbatch availability
        ret, stdout, stderr = run_preflight(["--mode", "real", "--config", temp_config_path, "--require-slurm"])
        if sbatch_avail:
            if ret == 0:
                print("[PASS] TEST 4: HPC/SLURM-required mode passes when sbatch is available.")
            else:
                print_failure_details("TEST 4: HPC/SLURM-required mode failed on host with sbatch", ["--mode", "real", "--config", temp_config_path, "--require-slurm"], ret, stdout, stderr)
                sys.exit(1)
        else:
            if ret != 0:
                print("[PASS] TEST 4: HPC/SLURM-required mode fails as expected when sbatch is unavailable.")
            else:
                print_failure_details("TEST 4: HPC/SLURM-required mode passed on host without sbatch!", ["--mode", "real", "--config", temp_config_path, "--require-slurm"], ret, stdout, stderr)
                sys.exit(1)
                
        # TEST 5: Verify that the preflight checker does not call readRDS / load RDS via R
        with open("scripts/python/preflight_checker.py", "r") as f:
            code = f.read()
        if "readRDS" not in code and "read_rds" not in code:
            print("[PASS] TEST 5: Preflight checker code does not contain load/deserialize commands (readRDS/read_rds).")
        else:
            print("[FAIL] TEST 5: Preflight checker code contains loading functions (readRDS/read_rds)!")
            sys.exit(1)

    print("=== All Preflight Tests Passed successfully! ===")

if __name__ == "__main__":
    main()
