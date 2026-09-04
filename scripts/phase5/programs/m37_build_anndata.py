#!/usr/bin/env python
"""Phase 5 - M37b - turn the R-exported Matrix Market counts into the AnnData
object cNMF consumes.

Runs in the ISOLATED p5_cnmf_env. R_env is never touched by this step.
The matrix written by R is genes x cells; cNMF wants cells x genes with raw
integer counts in .X, so it is transposed here and nowhere else.
"""
import sys, os, json
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import anndata as ad

run_dir = sys.argv[1]
out_h5ad = sys.argv[2]
cov_file = sys.argv[3] if len(sys.argv) > 3 else None

X = sio.mmread(os.path.join(run_dir, "counts.mtx")).tocsr()      # genes x cells
genes = [l.strip() for l in open(os.path.join(run_dir, "genes.txt"))]
cells = [l.strip() for l in open(os.path.join(run_dir, "cells.txt"))]
assert X.shape == (len(genes), len(cells)), (X.shape, len(genes), len(cells))

X = X.T.tocsr()                                                   # cells x genes
assert np.allclose(X.data, np.round(X.data)), "counts must be integers"
X = X.astype(np.float32)

obs = pd.DataFrame(index=pd.Index(cells, name="cell_id"))
if cov_file and os.path.exists(cov_file):
    cov = pd.read_csv(cov_file, sep="\t").set_index("cell_id")
    obs = obs.join(cov, how="left")
var = pd.DataFrame(index=pd.Index(genes, name="gene"))

adata = ad.AnnData(X=X, obs=obs, var=var)
adata.write_h5ad(out_h5ad, compression="gzip")

print(json.dumps({
    "run_dir": run_dir, "h5ad": out_h5ad,
    "n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars),
    "nonzero_frac": float(X.nnz / (X.shape[0] * X.shape[1])),
    "total_counts": float(X.sum()),
    "min_cell_counts": float(np.asarray(X.sum(axis=1)).min()),
    "patients": (obs["sample_id"].value_counts().to_dict()
                 if "sample_id" in obs.columns else None),
}, indent=2))
