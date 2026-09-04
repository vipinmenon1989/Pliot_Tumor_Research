#!/usr/bin/env python
"""Phase 5 - M37d - measure every cNMF solution on the K grid and select K.

The K-selection rule and the recurrence framework below are DECLARED HERE and
reference only measured properties of the factorization (stability, error,
program duplication, dead programs) plus patient bookkeeping. Neither rule
looks at a gene name, a pathway or a biological label, so K cannot be chosen to
produce a pleasing biological answer.

K-SELECTION RULE (pre-declared)
    K* = the LARGEST K on the grid satisfying all of:
      (i)   max pairwise cosine similarity between consensus spectra <= 0.75
            (no duplicated program),
      (ii)  every program is the dominant program for >= 1% of malignant cells
            (no dead program),
      (iii) cNMF stability >= the median stability of the K values that satisfy
            (i) and (ii).

RECURRENCE FRAMEWORK (pre-declared, and deliberately the Phase 4 shape so the
result is directly comparable to Phase 4's "0 of 8 states recurrent")
    per-cell usages are normalized to sum to 1 -> relative program usage
    a cell is "program-active" for p if relative usage(p) >= ACT (primary 0.20)
    a patient CARRIES p if it has >= MINCELL program-active cells (primary 10)
      AND >= MINFRAC of its malignant cells are program-active (primary 0.05)
    recurrent      >= 3 carrying patients
    shared-limited    2 carrying patients
    patient-private   1 carrying patient
    uncertain         0 carrying patients
"""
import os, sys, json, glob, itertools
import numpy as np
import pandas as pd

OUTDIR = "results/phase5/programs/cnmf"
TAB    = "results/phase5/tables"
VAL    = "results/phase5/validation"
COV    = "results/phase5/programs/cnmf_input/malignant_cell_covariates.tsv"
DT     = "0_1"                      # local-density-threshold 0.10, declared
KGRID  = list(range(4, 16))
MAX_COS, MIN_DOM_SHARE = 0.75, 0.01
ACT, MINCELL, MINFRAC, MINPAT = 0.20, 10, 0.05, 3
ACT_GRID   = [0.10, 0.15, 0.20, 0.25, 0.30]
FRAC_GRID  = [0.025, 0.05, 0.10]

os.makedirs(TAB, exist_ok=True); os.makedirs(VAL, exist_ok=True)
cov = pd.read_csv(COV, sep="\t").set_index("cell_id")


def paths(run, k, dt=DT):
    b = os.path.join(OUTDIR, run)
    return dict(
        spectra=os.path.join(b, f"{run}.spectra.k_{k}.dt_{dt}.consensus.txt"),
        score=os.path.join(b, f"{run}.gene_spectra_score.k_{k}.dt_{dt}.txt"),
        tpm=os.path.join(b, f"{run}.gene_spectra_tpm.k_{k}.dt_{dt}.txt"),
        usage=os.path.join(b, f"{run}.usages.k_{k}.dt_{dt}.consensus.txt"),
        stats=os.path.join(b, f"{run}.stats.k_{k}.dt_{dt}.df.npz"))


def read_df(p):
    return pd.read_csv(p, sep="\t", index_col=0)


def k_selection_stats(run):
    """cNMF 1.7.1 stores a DataFrame as data/index/columns inside the .npz."""
    f = os.path.join(OUTDIR, run, f"{run}.k_selection_stats.df.npz")
    d = np.load(f, allow_pickle=True)
    return pd.DataFrame(d["data"], index=d["index"],
                        columns=[str(c) for c in d["columns"]])


def cosine_matrix(S):
    X = S.values.astype(float)
    X = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-12)
    return X @ X.T


def recurrence(rel, patients, act, mincell, minfrac, minpat=MINPAT):
    """rel: cells x K relative usage; returns per-program carriage bookkeeping."""
    out = []
    npat_cells = patients.value_counts()
    for p in rel.columns:
        active = rel[p] >= act
        n_act = int(active.sum())
        per = {}
        carry = []
        for s in sorted(npat_cells.index):
            m = (patients == s)
            n_s = int((active & m).sum())
            f_s = n_s / int(m.sum())
            per[s] = dict(n_active=n_s, frac_of_patient=f_s)
            if n_s >= mincell and f_s >= minfrac:
                carry.append(s)
        dom_frac = (0.0 if n_act == 0 else
                    float(patients[active].value_counts().iloc[0]) / n_act)
        dom_pat = (None if n_act == 0 else
                   str(patients[active].value_counts().index[0]))
        ncar = len(carry)
        cls = ("recurrent" if ncar >= minpat else
               "shared-limited" if ncar == 2 else
               "patient-private" if ncar == 1 else "uncertain")
        out.append(dict(program=p, n_active_cells=n_act,
                        carrying_patients=",".join(carry), n_carrying=ncar,
                        dominant_patient=dom_pat,
                        dominant_patient_fraction=dom_frac,
                        recurrence=cls, per_patient=per))
    return out


rows, rec_rows, sens_rows = [], [], []
for run in ["primary", "balanced", "nocc"]:
    try:
        ks = k_selection_stats(run)
    except Exception as e:                                   # noqa: BLE001
        print(f"[warn] {run}: k_selection_stats unreadable ({e})")
        ks = pd.DataFrame()
    for k in KGRID:
        p = paths(run, k)
        if not os.path.exists(p["usage"]):
            print(f"[skip] {run} K={k}: no consensus yet")
            continue
        S = read_df(p["spectra"])                # K x genes
        U = read_df(p["usage"])                  # cells x K
        U.columns = [f"P{c}" for c in range(1, U.shape[1] + 1)]
        S.index = U.columns
        rel = U.div(U.sum(axis=1), axis=0)
        pats = cov.loc[U.index, "sample_id"] if run != "balanced" else \
            cov.loc[U.index, "sample_id"]

        C = cosine_matrix(S)
        off = C[np.triu_indices_from(C, k=1)]
        dom = rel.idxmax(axis=1).value_counts().reindex(U.columns).fillna(0)
        dom_share = dom / len(U)
        # top-50 gene overlap between programs (duplication, gene-level view)
        top50 = {c: set(S.loc[c].sort_values(ascending=False).index[:50])
                 for c in U.columns}
        jac = [len(top50[a] & top50[b]) / len(top50[a] | top50[b])
               for a, b in itertools.combinations(U.columns, 2)]

        st = np.nan; pe = np.nan
        if len(ks) and "k" in ks.columns:
            hit = ks[ks["k"] == k]
            if len(hit):
                for cand in ("silhouette", "stability"):
                    if cand in ks.columns:
                        st = float(hit[cand].iloc[0]); break
                if "prediction_error" in ks.columns:
                    pe = float(hit["prediction_error"].iloc[0])

        r = recurrence(rel, pats, ACT, MINCELL, MINFRAC)
        n_rec = sum(x["recurrence"] == "recurrent" for x in r)
        rows.append(dict(
            run=run, K=k, stability=st, prediction_error=pe,
            max_program_cosine=float(off.max()),
            mean_program_cosine=float(off.mean()),
            max_top50_jaccard=float(np.max(jac)),
            min_dominant_share=float(dom_share.min()),
            n_dead_programs=int((dom_share < MIN_DOM_SHARE).sum()),
            no_duplicate=bool(off.max() <= MAX_COS),
            no_dead=bool(dom_share.min() >= MIN_DOM_SHARE),
            n_recurrent_at_primary_threshold=n_rec,
            n_cells=int(len(U)), n_genes=int(S.shape[1])))
        for x in r:
            rec_rows.append(dict(run=run, K=k, **{q: x[q] for q in
                            ("program", "n_active_cells", "carrying_patients",
                             "n_carrying", "dominant_patient",
                             "dominant_patient_fraction", "recurrence")},
                            **{f"n_active_{s}": v["n_active"]
                               for s, v in x["per_patient"].items()},
                            **{f"frac_{s}": v["frac_of_patient"]
                               for s, v in x["per_patient"].items()}))
        if run == "primary":
            for a in ACT_GRID:
                for fr in FRAC_GRID:
                    rr = recurrence(rel, pats, a, MINCELL, fr)
                    sens_rows.append(dict(
                        run=run, K=k, activity_threshold=a, min_frac=fr,
                        min_cells=MINCELL,
                        n_recurrent=sum(x["recurrence"] == "recurrent" for x in rr),
                        n_shared_limited=sum(x["recurrence"] == "shared-limited" for x in rr),
                        n_patient_private=sum(x["recurrence"] == "patient-private" for x in rr),
                        n_uncertain=sum(x["recurrence"] == "uncertain" for x in rr)))

sel = pd.DataFrame(rows).sort_values(["run", "K"])
sel.to_csv(os.path.join(TAB, "PROGRAM_K_SELECTION.tsv"), sep="\t", index=False)
pd.DataFrame(rec_rows).to_csv(
    os.path.join(TAB, "PROGRAM_RECURRENCE_BY_K.tsv"), sep="\t", index=False)
pd.DataFrame(sens_rows).to_csv(
    os.path.join(TAB, "PROGRAM_RECURRENCE_THRESHOLD_SENSITIVITY.tsv"),
    sep="\t", index=False)

# ---- apply the declared K-selection rule to the PRIMARY run -----------------
pri = sel[sel.run == "primary"].copy()
elig = pri[pri.no_duplicate & pri.no_dead]
decision = {}
if len(elig):
    med = float(elig.stability.median()) if elig.stability.notna().any() else -np.inf
    ok = elig[(elig.stability >= med) | elig.stability.isna()]
    kstar = int(ok.K.max()) if len(ok) else int(elig.K.max())
    decision = dict(rule="largest K with max program cosine <= 0.75, no dead "
                         "program (>= 1% dominant share), and stability >= the "
                         "median stability of the K values satisfying both",
                    eligible_K=[int(x) for x in elig.K],
                    median_stability_of_eligible=med,
                    selected_K=kstar)
else:
    decision = dict(rule="no K satisfied the declared constraints", selected_K=None,
                    eligible_K=[])
print(json.dumps(decision, indent=2))
with open(os.path.join(VAL, "m37d_k_selection_decision.json"), "w") as fh:
    json.dump(dict(constants=dict(dt=DT, max_cosine=MAX_COS,
                                  min_dominant_share=MIN_DOM_SHARE,
                                  activity_threshold=ACT, min_cells=MINCELL,
                                  min_frac=MINFRAC, min_patients=MINPAT,
                                  activity_grid=ACT_GRID, frac_grid=FRAC_GRID),
                   decision=decision), fh, indent=2)
print(sel.to_string(index=False))
