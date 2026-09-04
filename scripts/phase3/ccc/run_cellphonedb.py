#!/usr/bin/env python
"""Phase 3 / Milestone M20 - CellPhoneDB execution, sample-aware.

CellPhoneDB is used as a third INDEPENDENT permutation-based ligand-receptor framework.
It runs in an isolated conda environment (cpdb_env) deliberately: its numpy/pandas pins
would destabilise R_env, which carries the frozen Phase 2 stack. R_env is left untouched.

Inputs are the per-sample sparse MTX directories written by M19 (results/phase3/ccc/
cellphonedb_inputs/<sample>/), each with matrix.mtx, barcodes.tsv, features.tsv and meta.tsv.

CellPhoneDB p-values are NOT numerically comparable to a LIANA rank or a CellChat
probability. Only a boolean support flag enters the M21 concordance model.
"""
import argparse, json, os, subprocess, sys, time, hashlib
from pathlib import Path

def log(msg):
    print(f"[{time.strftime('%Y-%m-%dT%H:%M:%S%z')}] [INFO] [Stage: phase3_m20_cellphonedb] {msg}", flush=True)

def warn(msg, store):
    store.append(msg)
    print(f"[{time.strftime('%Y-%m-%dT%H:%M:%S%z')}] [WARNING] [Stage: phase3_m20_cellphonedb] {msg}",
          file=sys.stderr, flush=True)

def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", default="results/phase3/ccc/cellphonedb_inputs")
    ap.add_argument("--out-dir", default="results/phase3/ccc/by_method")
    ap.add_argument("--db-dir", default="external/cellphonedb")
    ap.add_argument("--iterations", type=int, default=1000)
    ap.add_argument("--threshold", type=float, default=0.10)
    ap.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()

    warns = []
    t0 = time.time()
    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    dbd = Path(a.db_dir); dbd.mkdir(parents=True, exist_ok=True)

    import cellphonedb
    import pandas as pd, numpy as np, anndata as ad, scipy.io as sio
    log(f"cellphonedb import OK | pandas {pd.__version__} | numpy {np.__version__} | anndata {ad.__version__}")

    # ---- database ----
    from cellphonedb.utils import db_utils
    existing = sorted(dbd.rglob("*.zip"))
    if not existing:
        log("Downloading the CellPhoneDB database ...")
        # CellPhoneDB 5.0.1 signature is download_database(target_dir, cpdb_version)
        try:
            db_utils.download_database(str(dbd), cpdb_version="v5.0.0")
        except Exception as e:
            warn(f"download_database(cpdb_version='v5.0.0') failed: {e}; trying 'latest'", warns)
            db_utils.download_database(str(dbd), cpdb_version="latest")
        existing = sorted(dbd.rglob("*.zip"))
    if not existing:
        print("FATAL: no CellPhoneDB database available", file=sys.stderr); sys.exit(1)
    cpdb_file = str(existing[-1])
    log(f"CellPhoneDB database: {cpdb_file} (md5 {md5(cpdb_file)})")

    from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

    samples = sorted([p.name for p in Path(a.input_dir).iterdir() if p.is_dir()])
    log(f"Samples found: {samples}")
    all_rows, per_sample_counts, run_records = [], {}, {}

    for s in samples:
        sd = Path(a.input_dir) / s
        log(f"--- sample {s} ---")
        # MTX -> AnnData (cells x genes). Written by M19 as genes x cells.
        m = sio.mmread(sd / "matrix.mtx").tocsr()
        genes = [l.strip() for l in open(sd / "features.tsv")]
        cells = [l.strip() for l in open(sd / "barcodes.tsv")]
        meta = pd.read_csv(sd / "meta.tsv", sep="\t")
        adata = ad.AnnData(X=m.T.tocsr())
        adata.var_names = genes
        adata.obs_names = cells
        h5 = sd / f"{s}_counts.h5ad"
        adata.write_h5ad(h5)
        meta_p = sd / f"{s}_meta.tsv"
        meta.to_csv(meta_p, sep="\t", index=False)
        log(f"  {adata.n_obs} cells x {adata.n_vars} genes | {meta['cell_type'].nunique()} populations")

        odir = out / "cellphonedb_out" / s
        odir.mkdir(parents=True, exist_ok=True)
        try:
            res = cpdb_statistical_analysis_method.call(
                cpdb_file_path=cpdb_file,
                meta_file_path=str(meta_p),
                counts_file_path=str(h5),
                counts_data="hgnc_symbol",
                output_path=str(odir),
                iterations=a.iterations,
                threshold=a.threshold,
                threads=a.threads,
                debug_seed=a.seed,
                pvalue=0.05,
                output_suffix=s,
            )
        except Exception as e:
            warn(f"CellPhoneDB failed for {s}: {type(e).__name__}: {e}", warns)
            continue

        # relevant_interactions / pvalues -> long form
        try:
            pv = res["pvalues"] if isinstance(res, dict) and "pvalues" in res else None
            mn = res["means"] if isinstance(res, dict) and "means" in res else None
            if pv is None:
                cand = sorted(odir.glob("statistical_analysis_pvalues*.txt"))
                pv = pd.read_csv(cand[-1], sep="\t") if cand else None
                cand = sorted(odir.glob("statistical_analysis_means*.txt"))
                mn = pd.read_csv(cand[-1], sep="\t") if cand else None
            if pv is None:
                warn(f"No pvalues table for {s}", warns); continue
            idcols = [c for c in pv.columns if "|" not in c]
            paircols = [c for c in pv.columns if "|" in c]
            pl = pv.melt(id_vars=idcols, value_vars=paircols,
                         var_name="pair", value_name="pvalue")
            ml = mn.melt(id_vars=[c for c in mn.columns if "|" not in c],
                         value_vars=[c for c in mn.columns if "|" in c],
                         var_name="pair", value_name="mean_expr")
            key = ["interacting_pair", "pair"]
            df = pl.merge(ml[key + ["mean_expr"]], on=key, how="left")
            df[["source", "target"]] = df["pair"].str.split("|", expand=True)
            df["sample_id"] = s
            df["cellphonedb_supported"] = (df["pvalue"] < 0.05) & (df["mean_expr"] > 0)
            all_rows.append(df)
            per_sample_counts[s] = int(df["cellphonedb_supported"].sum())
            run_records[s] = {"cells": int(adata.n_obs), "genes": int(adata.n_vars),
                              "rows": int(len(df)), "supported": per_sample_counts[s]}
            log(f"  {len(df)} pair rows | {per_sample_counts[s]} supported (p<0.05 & mean>0)")
        except Exception as e:
            warn(f"Post-processing failed for {s}: {type(e).__name__}: {e}", warns)

    if not all_rows:
        print("FATAL: CellPhoneDB produced no usable results for any sample", file=sys.stderr)
        json.dump({"milestone": "M20", "method": "CellPhoneDB", "status": "FAILED",
                   "warnings": warns}, open(out / "m20_cellphonedb_record.json", "w"), indent=2)
        sys.exit(1)

    import pandas as pd
    full = pd.concat(all_rows, ignore_index=True)
    full.to_csv(out / "all_CellPhoneDB_interactions.tsv", sep="\t", index=False)
    sup = full[full["cellphonedb_supported"]]
    sup.to_csv(out / "CellPhoneDB_supported_interactions.tsv", sep="\t", index=False)
    log(f"Pooled rows: {len(full)} | supported: {len(sup)}")
    log(f"Supported per sample: {per_sample_counts}")
    tum = sup[(sup['source'] == 'MPNST-Tumor') | (sup['target'] == 'MPNST-Tumor')]
    log(f"Tumour-involving supported: {len(tum)}")

    rec = {
        "milestone": "M20", "phase": "phase3", "method": "CellPhoneDB",
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "script": "scripts/phase3/ccc/run_cellphonedb.py",
        "environment": "isolated conda env cpdb_env (R_env left untouched; cellphonedb's numpy/pandas pins would destabilise the frozen Phase 2 stack)",
        "parameters": {"database_file": cpdb_file, "database_md5": md5(cpdb_file),
                       "analysis_mode": "statistical_analysis", "iterations": a.iterations,
                       "threshold": a.threshold, "pvalue_cut": 0.05, "counts_data": "hgnc_symbol",
                       "threads": a.threads, "seed": a.seed,
                       "expression_basis": "RNA assay, joined, LogNormalize (prepared in M19)",
                       "support_rule": "pvalue < 0.05 AND mean_expr > 0"},
        "comparability_caveat": "CellPhoneDB p-values are not numerically comparable to a LIANA rank or a CellChat probability. Only a boolean support flag enters the concordance model.",
        "results": {"samples_run": list(run_records), "per_sample": run_records,
                    "n_rows_pooled": int(len(full)), "n_supported": int(len(sup)),
                    "tumour_involving_supported": int(len(tum))},
        "warnings": warns,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID", "NOT_IN_SLURM"),
        "elapsed_seconds": time.time() - t0,
    }
    json.dump(rec, open(out / "m20_cellphonedb_record.json", "w"), indent=2)
    log(f"CellPhoneDB COMPLETE in {(time.time()-t0)/60:.1f} min. Warnings: {len(warns)}")

if __name__ == "__main__":
    main()
