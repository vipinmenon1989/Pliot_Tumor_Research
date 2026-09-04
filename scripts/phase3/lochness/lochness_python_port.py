#!/usr/bin/env python
"""LochNESS computed with the perturb-seq pipeline's FORMULATION, for cross-implementation
comparison only (Phase 3 M24, spec SS36).

Algorithm transcribed from the read-only reference
  perturbseq-pipeline/src/perturbseq_pipeline/lochness.py
whose module docstring defines
      lochNESS(cell, g) = local_fraction(g) / overall_fraction(g) - 1
and whose implementation choices are:
  * vectorised sparse matrix-vector product over a binarised kNN adjacency
  * the local fraction divided by the ACTUAL neighbour count per row (not the requested k)
  * NO same-sample exclusion
  * k taken from configuration rather than from round(0.5*sqrt(N))

This script is a faithful re-expression of that formulation. It does NOT import the
perturb-seq package, so Pilot_MPNST has no runtime dependency on it.

Two modes:
  --mode perturbseq  : the formulation as written there (no exclusion, actual-count denominator)
  --mode mmca        : same code path but with same-sample exclusion and a fixed-k denominator,
                       i.e. mathematically equivalent inputs to the R MMCA implementation
"""
import argparse, json, sys
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.neighbors import NearestNeighbors

ap = argparse.ArgumentParser()
ap.add_argument("--embedding", required=True)   # TSV: cell + dims
ap.add_argument("--meta", required=True)        # TSV: cell, sample_id, context
ap.add_argument("--out", required=True)
ap.add_argument("--k", type=int, required=True)
ap.add_argument("--target", default="high")
ap.add_argument("--mode", choices=["perturbseq", "mmca"], default="perturbseq")
a = ap.parse_args()

emb = pd.read_csv(a.embedding, sep="\t", index_col=0)
meta = pd.read_csv(a.meta, sep="\t", index_col=0).loc[emb.index]
X = emb.to_numpy(dtype=np.float64)
samples = meta["sample_id"].to_numpy()
context = meta["context"].to_numpy()
n = X.shape[0]
is_target = (context == a.target)

if a.mode == "perturbseq":
    # perturbseq formulation: one global kNN graph, self excluded, NO same-sample exclusion,
    # local fraction divided by the ACTUAL neighbour count, global fraction over ALL cells.
    nn = NearestNeighbors(n_neighbors=a.k + 1, algorithm="kd_tree").fit(X)
    _, idx = nn.kneighbors(X)
    idx = idx[:, 1:]                                   # drop self
    rows = np.repeat(np.arange(n), idx.shape[1])
    adj = sparse.csr_matrix((np.ones(rows.size, dtype=np.uint8), (rows, idx.ravel())), shape=(n, n))
    neighbour_counts = np.asarray(adj.sum(axis=1)).ravel().astype(float)
    neighbour_counts[neighbour_counts == 0] = np.nan
    local = (adj @ is_target.astype(float)) / neighbour_counts
    overall = is_target.mean()
    score = local / overall - 1.0
else:
    # mmca-equivalent: query = one sample, reference = all OTHER samples (disjoint sets),
    # local fraction divided by the requested k, global fraction over the REFERENCE set.
    score = np.full(n, np.nan)
    for s in np.unique(samples):
        qi = np.where(samples == s)[0]
        ri = np.where(samples != s)[0]
        if qi.size == 0 or ri.size <= a.k:
            continue
        nn = NearestNeighbors(n_neighbors=a.k, algorithm="kd_tree").fit(X[ri])
        _, idx = nn.kneighbors(X[qi])
        n_target = is_target[ri][idx].sum(axis=1).astype(float)
        overall = is_target[ri].mean()
        if overall <= 0 or overall >= 1:
            continue
        score[qi] = (n_target / a.k) / overall - 1.0

out = pd.DataFrame({"cell": emb.index, "lochness_python": score,
                    "sample_id": samples, "context": context})
out.to_csv(a.out, sep="\t", index=False)
print(json.dumps({"mode": a.mode, "k": a.k, "n_cells": int(n),
                  "n_scored": int(np.isfinite(score).sum()),
                  "mean": float(np.nanmean(score)), "sd": float(np.nanstd(score))}), flush=True)
