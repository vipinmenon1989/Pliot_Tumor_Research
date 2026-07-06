# Python script to compute Conda environment package diff
# Generates machine-readable JSON and formatted Markdown reports.

import json
import os
import re

def parse_conda_list(file_path):
    packages = {}
    if not os.path.exists(file_path):
        return packages
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = re.split(r'\s+', line.strip())
            if len(parts) >= 2:
                name = parts[0]
                version = parts[1]
                build = parts[2] if len(parts) >= 3 else ""
                channel = parts[3] if len(parts) >= 4 else ""
                packages[name] = {
                    "version": version,
                    "build": build,
                    "channel": channel
                }
    return packages

def compare_envs(before_path, after_path):
    before = parse_conda_list(before_path)
    after = parse_conda_list(after_path)

    added = {}
    removed = {}
    changed = {}

    for name, info in after.items():
        if name not in before:
            added[name] = info
        elif before[name]["version"] != info["version"]:
            changed[name] = {
                "before": before[name],
                "after": info
            }

    for name, info in before.items():
        if name not in after:
            removed[name] = info

    # Classify changed as upgrade/downgrade based on alphanumeric version check
    upgraded = {}
    downgraded = {}
    for name, diff in changed.items():
        v1 = diff["before"]["version"]
        v2 = diff["after"]["version"]
        # Simple alphanumeric version compare as fallback
        if v2 > v1:
            upgraded[name] = diff
        else:
            downgraded[name] = diff

    # Check critical R and single-cell packages
    # Conda packages can be named directly or prefixed with r- (e.g. r-seurat, r-base)
    critical_names = [
        "R", "Seurat", "SeuratObject", "Matrix", "sctransform",
        "SingleCellExperiment", "SummarizedExperiment", "BiocGenerics",
        "future", "future.apply", "ggplot2", "patchwork", "data.table",
        "harmony", "DoubletFinder", "scDblFinder"
    ]

    critical_changes = {}
    for name in critical_names:
        conda_variants = [name, name.lower(), f"r-{name.lower()}", f"bioconductor-{name.lower()}"]
        if name == "R":
            conda_variants.append("r-base")
        found_in_diff = False
        for var in conda_variants:
            if var in added:
                critical_changes[name] = {"status": "added", "details": added[var]}
                found_in_diff = True
            elif var in removed:
                critical_changes[name] = {"status": "removed", "details": removed[var]}
                found_in_diff = True
            elif var in changed:
                critical_changes[name] = {"status": "changed", "details": changed[var]}
                found_in_diff = True
        if not found_in_diff:
            # Report current version if same
            current_ver = "NOT_INSTALLED"
            for var in conda_variants:
                if var in after:
                    current_ver = after[var]["version"]
                    break
            critical_changes[name] = {"status": "unchanged", "version": current_ver}

    return {
        "added": added,
        "removed": removed,
        "upgraded": upgraded,
        "downgraded": downgraded,
        "critical_package_status": critical_changes
    }

def main():
    before_file = "workflow/envs/R_env_list.txt"
    after_file = "workflow/envs/R_env_final_list.txt"

    diff = compare_envs(before_file, after_file)

    # Save JSON report
    os.makedirs("reports/audits", exist_ok=True)
    with open("reports/audits/conda_env_diff.json", "w") as f:
        json.dump(diff, f, indent=2)

    # Generate Markdown report
    md_lines = [
        "# Conda Environment Comparison Report",
        "**Before vs After Snakemake Installation in R_env**",
        "",
        "## 1. Summary of Changes",
        f"- **Packages Added**: {len(diff['added'])}",
        f"- **Packages Removed**: {len(diff['removed'])}",
        f"- **Packages Upgraded**: {len(diff['upgraded'])}",
        f"- **Packages Downgraded**: {len(diff['downgraded'])}",
        "",
        "## 2. Critical Package Integrity Verification",
        "The following table lists the status of critical R and single-cell packages:",
        "",
        "| Package | Status | Version | Change details |",
        "| --- | --- | --- | --- |"
    ]

    for pkg, status_info in diff["critical_package_status"].items():
        status = status_info["status"]
        if status == "unchanged":
            md_lines.append(f"| {pkg} | Unchanged | {status_info['version']} | None |")
        elif status == "changed":
            v_before = status_info["details"]["before"]["version"]
            v_after = status_info["details"]["after"]["version"]
            md_lines.append(f"| {pkg} | CHANGED | {v_after} | Upgraded/downgraded from {v_before} |")
        elif status == "added":
            v = status_info["details"]["version"]
            md_lines.append(f"| {pkg} | ADDED | {v} | Package newly installed |")
        elif status == "removed":
            md_lines.append(f"| {pkg} | REMOVED | N/A | Package uninstalled |")

    md_lines.extend([
        "",
        "## 3. Detailed Package Transactions",
        ""
    ])

    if diff["added"]:
        md_lines.extend([
            "### Packages Added",
            "| Package | Version | Build | Channel |",
            "| --- | --- | --- | --- |"
        ])
        for name, info in sorted(diff["added"].items()):
            md_lines.append(f"| {name} | {info['version']} | {info['build']} | {info['channel']} |")
        md_lines.append("")

    if diff["upgraded"] or diff["downgraded"]:
        md_lines.extend([
            "### Packages Updated",
            "| Package | Before Version | After Version | Channel |",
            "| --- | --- | --- | --- |"
        ])
        for name, info in sorted({**diff["upgraded"], **diff["downgraded"]}.items()):
            md_lines.append(f"| {name} | {info['before']['version']} | {info['after']['version']} | {info['after']['channel']} |")
        md_lines.append("")

    if diff["removed"]:
        md_lines.extend([
            "### Packages Removed",
            "| Package | Version | Build | Channel |",
            "| --- | --- | --- | --- |"
        ])
        for name, info in sorted(diff["removed"].items()):
            md_lines.append(f"| {name} | {info['version']} | {info['build']} | {info['channel']} |")
        md_lines.append("")

    with open("reports/audits/conda_env_diff.md", "w") as f:
        f.write("\n".join(md_lines) + "\n")

    print("Diff report successfully generated.")

if __name__ == "__main__":
    main()
