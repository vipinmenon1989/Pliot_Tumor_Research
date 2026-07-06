#!/bin/bash
# Reusable Project Verification Utility for MPNST Phase 1 Workflow
#
# Purpose:
#   This utility automates post-milestone sanity checks including Git status,
#   whitespace compliance, directory structure auditing, and repository safety
#   checks (verifying that the real processed_mpnst.rds dataset is not referenced
#   or loaded directly via readRDS in critical workflow files).
#
# Usage:
#   ./scripts/shell/verify_milestone.sh

set -euo pipefail

# Ensure working directory is the project root
PROJECT_ROOT="/local/projects-t3/lilab/vmenon/Pilot_tumor"
cd "${PROJECT_ROOT}"

echo "=== 1. Git Status ==="
git status

echo -e "\n=== 2. Git Cached Diff Stats ==="
git diff --cached --stat

echo -e "\n=== 3. Git Whitespace and Check Checks ==="
git diff --cached --check

echo -e "\n=== 4. Repository Directory Structure Audit ==="
find . -maxdepth 4 -type f -not -path '*/.git/*' -not -path '*/.snakemake/*' | sort

echo -e "\n=== 5. Environment Specifications Listing ==="
find workflow/envs -maxdepth 1 -type f -print -exec ls -lh {} \;

echo -e "\n=== 6. Safety Check: Real RDS References ==="
# Search for occurrences of processed_mpnst.rds
RDS_GREP=$(grep -rn "processed_mpnst.rds" Snakefile workflow scripts tests config .github 2>/dev/null || true)
if [ -n "${RDS_GREP}" ]; then
  echo "Found references to processed_mpnst.rds:"
  echo "${RDS_GREP}"
else
  echo "No direct code references to processed_mpnst.rds (safe)."
fi

echo -e "\n=== 7. Safety Check: readRDS Occurrences ==="
# Search for occurrences of readRDS in critical workflow/code files
READRDS_GREP=$(grep -rn "readRDS" Snakefile workflow scripts tests 2>/dev/null || true)
if [ -n "${READRDS_GREP}" ]; then
  echo "Found occurrences of readRDS:"
  echo "${READRDS_GREP}"
else
  echo "No occurrences of readRDS in workflow or analysis scripts (safe)."
fi

echo -e "\n=== 8. Absolute User Path References ==="
# Search for user-specific path leaks
USERPATH_GREP=$(grep -rn "/local/projects-t3/lilab/vmenon" Snakefile workflow scripts tests config .github 2>/dev/null || true)
if [ -n "${USERPATH_GREP}" ]; then
  echo "Found user-specific absolute path references:"
  echo "${USERPATH_GREP}"
else
  echo "No user-specific absolute path references found."
fi
