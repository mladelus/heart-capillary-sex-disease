"""Shared code for the endothelial-state analyses (plan section M).
Reads the per-donor count files saved by R script 02 (dl/*.rds) and the cell table of script 03 (ec_subtypes.rds),
and reproduces the donor-level model and meta-analysis of the R pipeline (sex_fit and meta_one in 00_setup.R)."""
import os, sys, glob, pickle
import numpy as np, pandas as pd
from scipy import sparse, stats, optimize
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "figures")); sys.path.insert(0, "/home/claude/data")
from rds import read_rds
P1 = os.path.expanduser(os.environ.get("HCAP_DATA", "~/Documents/heart_capillary_project"))
P2 = os.path.expanduser(os.environ.get("HCAP2_DATA", "~/Documents/heart_capillary_disease"))
DIRS = {"healthy": P1, "disease": P2}
OUT = os.path.join(HERE, "out"); os.makedirs(OUT, exist_ok=True)
MIN_CAP, MIN_PER_SEX, MIN_POS = 30, 3, 5
ENDOTHELIN = ["EDN1", "EDN2", "EDN3", "ECE1", "ECE2", "EDNRA", "EDNRB"]
VESSEL_MARKERS = ["CA4", "RGCC", "BTNL9", "GJA5", "HEY1", "SEMA3G", "GJA4", "DKK2", "FBLN5", "ACKR1", "NR2F2", "CPE"]
FLOW = ["KLF2", "KLF4", "NOS3", "THBD"]


def samples():
    """The 138 donors with stratum, cohort, sex, age and assay (from the saved results of script 18)."""
    d = read_rds(os.path.join(P2, "stage3_endothelin.rds"))["donors"]
    d = d[["sample", "stratum", "cohort", "sex", "age", "assay"]].copy()
    for c in ("stratum", "cohort", "sex", "assay"): d[c] = d[c].astype(str)
    return d.set_index("sample")


def donors(S):
    """Yield (sample, cohort, genes, counts[genes x cells, CSC], cell table aligned to the columns) for each included donor.
    Uses the .npz cache written by cache_cells.py."""
    cache = os.path.join(OUT, "cells")
    for co, folder in DIRS.items():
        ec = read_rds(os.path.join(folder, "ec_subtypes.rds"))
        for c in ("dataset_id", "donor_id", "soma_joinid", "cell_type", "subtype"): ec[c] = ec[c].astype(str)
        ec["key"] = ec.dataset_id.str[:8] + "|" + ec.donor_id
        groups = {k: g.set_index("soma_joinid") for k, g in ec.groupby("key")}
        for f in sorted(glob.glob(os.path.join(cache, co + "__*.npz"))):
            z = np.load(f, allow_pickle=False); s = str(z["dataset_id"])[:8] + "|" + str(z["donor_id"])
            if s not in S.index or S.loc[s, "cohort"] != co or s not in groups: continue
            genes, cells = z["genes"], z["cells"]
            M = sparse.csc_matrix((z["x"].astype(float), z["i"], z["p"]), shape=tuple(int(v) for v in z["dim"]))
            a = groups[s]; keep = np.isin(cells, a.index.values)
            M = M[:, keep]; tab = a.loc[cells[keep]].reset_index()
            tab["umi"] = np.asarray(M.sum(axis=0)).ravel()
            M.r = M.tocsr()                                   # row-oriented copy for fast gene look-ups
            yield s, co, genes, M, tab


def row(M, genes, g):
    i = np.where(genes == g)[0]
    return np.asarray(M.r[i[0], :].todense()).ravel() if len(i) else np.zeros(M.shape[1])


def rows_sum(M, genes, gs):
    i = np.where(np.isin(genes, gs))[0]
    return np.asarray(M.r[i, :].sum(axis=0)).ravel() if len(i) else np.zeros(M.shape[1])


# ---------------- statistics, as in 00_setup.R ----------------
def sex_fit(d, y, extra=()):
    """Female - male in SD units of y, adjusted for age (and assay when it varies) within one stratum."""
    d = d[np.isfinite(d[y]) & d.age.notna()]
    nf, nm = int((d.sex == "female").sum()), int((d.sex == "male").sum())
    if nf < MIN_PER_SEX or nm < MIN_PER_SEX or d[y].std() == 0: return None
    yy = (d[y] - d[y].mean()) / d[y].std()
    X = [np.ones(len(d)), (d.sex == "female").astype(float).values, d.age.values.astype(float)] + [d[e].values.astype(float) for e in extra]
    if d.assay.nunique() > 1:
        for lv in sorted(d.assay.unique())[1:]: X.append((d.assay == lv).astype(float).values)
    X = np.column_stack(X); n, k = X.shape
    if n <= k: return None
    XtXi = np.linalg.pinv(X.T @ X); b = XtXi @ X.T @ yy.values; res = yy.values - X @ b
    s2 = res @ res / (n - k); se = np.sqrt(s2 * XtXi[1, 1]); t = b[1] / se
    return dict(est=b[1], se=se, p=2 * stats.t.sf(abs(t), n - k), n_w=nf, n_m=nm)


def meta_one(r):
    """Random-effects meta-analysis (REML; fixed effect with fewer than three strata), as metafor::rma."""
    r = r[np.isfinite(r.est) & np.isfinite(r.se)]
    if not len(r): return None
    y, v = r.est.values, r.se.values ** 2; k = len(y)
    tau2 = 0.0
    if k >= 3:
        def nll(t2):
            w = 1 / (v + t2); mu = (w * y).sum() / w.sum()
            return 0.5 * (np.log(v + t2).sum() + np.log(w.sum()) + (w * (y - mu) ** 2).sum())
        o = optimize.minimize_scalar(nll, bounds=(0, 100 * max(v.max(), np.var(y) + 1e-9)), method="bounded", options={"xatol": 1e-10})
        tau2 = o.x if nll(o.x) < nll(0.0) - 1e-12 else 0.0
        if tau2 < 1e-8: tau2 = 0.0
    w = 1 / (v + tau2); mu = (w * y).sum() / w.sum(); se = np.sqrt(1 / w.sum()); z = mu / se
    w0 = 1 / v; Q = (w0 * (y - (w0 * y).sum() / w0.sum()) ** 2).sum()
    vt = (k - 1) * w0.sum() / (w0.sum() ** 2 - (w0 ** 2).sum()) if k > 1 else np.nan
    return dict(k=k, n_w=int(r.n_w.sum()), n_m=int(r.n_m.sum()), est=mu, lo=mu - 1.96 * se, hi=mu + 1.96 * se, p=2 * stats.norm.sf(abs(z)),
                I2=100 * tau2 / (tau2 + vt) if k >= 3 else np.nan, tau2=tau2, same_dir=int((np.sign(y) == np.sign(mu)).sum()))


def pooled(don, y, extra=(), by="stratum"):
    """Per-stratum estimates and the pooled estimate for outcome y."""
    per = []
    for s, d in don.groupby(by):
        r = sex_fit(d, y, extra)
        if r: per.append(dict(r, stratum=s, cohort=d.cohort.iloc[0]))
    per = pd.DataFrame(per)
    out = {"per": per}
    if len(per):
        out["both"] = meta_one(per)
        for co in ("healthy", "disease"): out[co] = meta_one(per[per.cohort == co])
    return out


def short(x):
    x = str(x)
    for k, v in (("dilated", "DCM"), ("arrhythmogenic", "ACM"), ("myocardial infarction", "Infarction"), ("myocarditis", "Myocarditis")):
        if x.startswith(k): return v
    return "Healthy 1" if x == "Heart" else "Healthy 2"
