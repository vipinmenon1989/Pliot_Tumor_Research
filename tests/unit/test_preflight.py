#!/usr/bin/env python3
import subprocess
import os
import sys
import yaml

def run_preflight(args):
    cmd = [sys.executable, "scripts/python/preflight_checker.py"] + args
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return res.returncode, res.stdout, res.stderr

def main():
    print("=== Running Preflight Checker Tests ===")
    
    # Check 1: Synthetic mode passes without real RDS requirement
    ret, stdout, stderr = run_preflight(["--mode", "synthetic"])
    print(f"Synthetic Mode Exit Code: {ret}")
    if ret == 0:
        print("[PASS] Synthetic mode passes successfully.")
    else:
        print(f"[FAIL] Synthetic mode failed. Stderr: {stderr}")
        sys.exit(1)
        
    # Check 2: Real mode fails when the configured real RDS is absent
    temp_config_path = "config/config.temp_test.yaml"
    with open("config/config.yaml", "r") as f:
        config_data = yaml.safe_load(f)
    
    config_data["input_rds"] = "nonexistent_file_path.rds"
    with open(temp_config_path, "w") as f:
        yaml.safe_dump(config_data, f)
        
    try:
        ret, stdout, stderr = run_preflight(["--mode", "real", "--config", temp_config_path])
        print(f"Real Mode (with absent RDS) Exit Code: {ret}")
        if ret != 0:
            print("[PASS] Real mode fails as expected when input RDS is absent.")
        else:
            print("[FAIL] Real mode passed when input RDS was absent!")
            sys.exit(1)
    finally:
        if os.path.exists(temp_config_path):
            os.remove(temp_config_path)
            
    # Check 3: Real mode passes when a valid configured input path exists
    # (Note: real mode also checks sbatch, which is available on this login node)
    ret, stdout, stderr = run_preflight(["--mode", "real"])
    print(f"Real Mode (with valid RDS) Exit Code: {ret}")
    if ret == 0:
        print("[PASS] Real mode passes when a valid configured input path exists.")
    else:
        print(f"[FAIL] Real mode failed with valid RDS. Stderr: {stderr}")
        sys.exit(1)
        
    # Check 4: No mode loads the large RDS during lightweight preflight
    with open("scripts/python/preflight_checker.py", "r") as f:
        code = f.read()
    
    if "readRDS" not in code and "read_rds" not in code:
        print("[PASS] Verified that the preflight checker does not load the RDS file via R.")
    else:
        print("[FAIL] Preflight checker code contains loading functions (readRDS/read_rds)!")
        sys.exit(1)
        
    print("=== All Preflight Tests Passed successfully! ===")

if __name__ == "__main__":
    main()
