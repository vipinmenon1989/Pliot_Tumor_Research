#!/usr/bin/env python3
"""
Project State Consistency Validator for MPNST Phase 1 Handoff
Verifies that README.md, PROGRESS.md, CHANGELOG.md, and milestone reports
consistently reflect the M0-M9 complete, Phase 1 frozen status.
"""

import os
import sys
import re

def log_success(msg):
    print(f"\033[92m[PASS]\033[0m {msg}")

def log_failure(msg):
    print(f"\033[91m[FAIL]\033[0m {msg}")

def file_contains(path, pattern, flags=re.IGNORECASE):
    if not os.path.exists(path):
        return False
    with open(path, 'r', encoding='utf-8') as f:
        content = f.read()
    return bool(re.search(pattern, content, flags))

def main():
    print("=== Starting Project State Consistency Validator ===")
    errors = 0

    # Files to check
    files = {
        "PROJECT.md": "PROJECT.md",
        "README.md": "README.md",
        "PROGRESS.md": "PROGRESS.md",
        "CHANGELOG.md": "CHANGELOG.md"
    }

    # Verify basic file existence
    for name, path in files.items():
        if os.path.exists(path):
            log_success(f"File '{name}' exists.")
        else:
            log_failure(f"File '{name}' is missing.")
            errors += 1

    # 1. Verify milestone reports existence
    for m in range(10):
        rep_path = f"reports/milestones/M{m}_REPORT.md"
        if os.path.exists(rep_path):
            log_success(f"Milestone report '{rep_path}' exists.")
        else:
            log_failure(f"Milestone report '{rep_path}' is missing.")
            errors += 1

    # 2. Verify PROGRESS.md consistency
    # PROGRESS.md must contain M0-M9 complete indicators
    for m in range(10):
        pattern = rf"(\*\*M{m}\*\*\s*\|.*Completed|###\s*M{m}.*Status\*\*:\s*Completed)"
        # We read the file content as a single multiline block to allow dotall-like matching if needed
        with open("PROGRESS.md", "r", encoding="utf-8") as pf:
            progress_content = pf.read()
        # Find if either the table row or status log matches
        match_table = re.search(rf"\*\*M{m}\*\*\s*\|.*Completed", progress_content, re.IGNORECASE)
        match_log = re.search(rf"###\s*M{m}.*?\-\s*\*\*Status\*\*:\s*Completed", progress_content, re.IGNORECASE | re.DOTALL)
        if match_table or match_log:
            log_success(f"PROGRESS.md shows M{m} as Completed.")
        else:
            log_failure(f"PROGRESS.md does not show M{m} as Completed.")
            errors += 1
            
    # Check Phase 1 Frozen & Phase 2 Not Started in PROGRESS
    if file_contains("PROGRESS.md", r"PHASE\s*1\s*=\s*FROZEN"):
        log_success("PROGRESS.md shows 'PHASE 1 = FROZEN'.")
    else:
        log_failure("PROGRESS.md is missing 'PHASE 1 = FROZEN'.")
        errors += 1

    if file_contains("PROGRESS.md", r"PHASE\s*2\s*=\s*NOT\s*STARTED"):
        log_success("PROGRESS.md shows 'PHASE 2 = NOT STARTED'.")
    else:
        log_failure("PROGRESS.md is missing 'PHASE 2 = NOT STARTED'.")
        errors += 1

    # 3. Verify CHANGELOG.md consistency
    # CHANGELOG.md must contain entries for all milestones [M0] ... [M9]
    for m in range(10):
        pattern = rf"\[M{m}\]"
        if file_contains("CHANGELOG.md", pattern):
            log_success(f"CHANGELOG.md represents Milestone M{m}.")
        else:
            log_failure(f"CHANGELOG.md is missing entries for [M{m}].")
            errors += 1

    # 4. Verify README.md consistency
    # README must consistently show M0-M9 complete, Phase 1 frozen, Phase 2 not started
    if file_contains("README.md", r"Completed:\s*Milestones\s*M0–M9|Completed\s*M0–M9"):
        log_success("README.md shows Milestones M0-M9 as completed.")
    else:
        log_failure("README.md does not show 'Completed: Milestones M0-M9'.")
        errors += 1

    if file_contains("README.md", r"PHASE\s*1\s*=\s*FROZEN"):
        log_success("README.md shows 'PHASE 1 = FROZEN'.")
    else:
        log_failure("README.md is missing 'PHASE 1 = FROZEN'.")
        errors += 1

    if file_contains("README.md", r"PHASE\s*2\s*=\s*NOT\s*STARTED"):
        log_success("README.md shows 'PHASE 2 = NOT STARTED'.")
    else:
        log_failure("README.md is missing 'PHASE 2 = NOT STARTED'.")
        errors += 1

    # Check for M9_REPORT.md in README
    if file_contains("README.md", r"M9_REPORT.md"):
        log_success("README.md references M9_REPORT.md.")
    else:
        log_failure("README.md is missing reference to M9_REPORT.md.")
        errors += 1

    # README must contain no stale current-status phrases
    stale_phrases = [
        r"Through\s+Milestone\s+[5-8]\b",
        r"M[5-8]\s+next\b",
        r"M9\s+next\b",
        r"M[8-9]\s+not\s+started\b",
        r"Completed:\s*M0–M[5-8]\b"
    ]
    # Scan README.md for active stale phrases (excluding lines with "historical" or "previous")
    if os.path.exists("README.md"):
        with open("README.md", 'r', encoding='utf-8') as f:
            readme_lines = f.readlines()
        for idx, line in enumerate(readme_lines):
            # Skip historical explanation lines
            if "historical" in line.lower() or "previous" in line.lower() or "checkpoint" in line.lower():
                continue
            for sp in stale_phrases:
                if re.search(sp, line, re.IGNORECASE):
                    log_failure(f"README.md line {idx+1} contains stale phrase matching '{sp}': '{line.strip()}'")
                    errors += 1

    # 5. Verify required milestone outputs (manifest and handoff reports)
    required_outputs = [
        "results/phase1_manifest.json",
        "reports/PHASE1_HANDOFF.md",
        "reports/audits/M0_M8_RECONCILIATION.md",
        "reports/audits/M0_M8_RECONCILIATION.tsv"
    ]
    for ro in required_outputs:
        if os.path.exists(ro):
            log_success(f"Final target output '{ro}' exists.")
        else:
            log_failure(f"Final target output '{ro}' is missing.")
            errors += 1

    print("=== Project State Consistency Validation Summary ===")
    if errors == 0:
        print("\033[92mAll project state consistency checks PASSED.\033[0m")
        sys.exit(0)
    else:
        print(f"\033[91mFailed {errors} consistency check(s).\033[0m")
        sys.exit(1)

if __name__ == "__main__":
    main()
