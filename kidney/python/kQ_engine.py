"""Pseudobulk sex-effect engine for the kidney (plan section Q). OLS on log2 CPM within stratum; exact t-tests."""
from k_common import *
POS = pd.read_csv('/home/claude/proj/p2/gene_positions.csv'); CHR = POS.drop_duplicates('gene_name').set_index('gene_name').seq_name.astype(str)
def chr_class(sym):
    c = CHR.reindex(sym).values; return np.where(c == 'X', 'X', np.where(c == 'Y', 'Y', np.where(pd.Series(c).isin([str(i) for i in range(1, 23)]).values, 'autosome', 'unknown')))
HARM = {'snRNA': {'PT': 'PT', 'TAL': 'TAL', 'PC': 'PC', 'IC': 'IC', 'EC': 'EC', 'DCT': 'DCT/CNT', 'CNT': 'DCT/CNT', 'FIB': 'stroma', 'VSM/P': 'stroma', 'IMM': 'immune', 'POD': 'POD/PEC', 'PEC': 'POD/PEC', 'DTL': 'thin limb', 'ATL': 'thin limb'},
        'scRNA': {'PT': 'PT', 'TAL': 'TAL', 'PC': 'PC', 'IC': 'IC', 'EC': 'EC', 'DCT/CNT': 'DCT/CNT', 'VSMC/MC/FIB': 'stroma', 'T': 'immune', 'MYL': 'immune', 'B': 'immune', 'POD/PEC': 'POD/PEC', 'DTL/ATL': 'thin limb'}}
_cache = {}
def dataset(ds):
    if ds not in _cache:
        d = load(ds); o = d['obs']; d['meta'] = donor_meta(o); d['meta'] = d['meta'][(d['meta'].age >= 18) & d['meta'].sex.isin(['female', 'male'])]
        _cache[ds] = d
    return _cache[ds]
def pseudobulk(ds, key, groups):
    """Counts (genes x donors, summed over the labels in `groups`) and cells per donor. key: 'subclassl1', 'subclassl2' or 'celltype'."""
    d = dataset(ds); M = d['pb_' + key]; cols = pd.Series(d['pb_' + key + '_cols']).str.split('||', regex=False, expand=True); cols.columns = ['donor', 'lab']
    col = {'subclassl1': 'subclass.l1', 'subclassl2': 'subclass.l2', 'celltype': 'cell_type'}[key]
    o = d['obs']; ncell = o[o[col].isin(groups)].groupby('donor_id').size()
    sel = cols.lab.isin(groups).values; donors = sorted(set(cols.donor[sel]) & set(d['meta'].index))
    ind = sp.csr_matrix((np.ones(sel.sum()), (np.where(sel)[0], pd.Categorical(cols.donor[sel], categories=sorted(set(cols.donor[sel]))).codes)), shape=(M.shape[1], len(set(cols.donor[sel]))))
    P = pd.DataFrame((M @ ind).toarray(), index=d['sym'], columns=sorted(set(cols.donor[sel])))[donors]
    P = P.groupby(level=0).sum()                                    # sum gene ids that share a symbol
    return P, ncell.reindex(donors).fillna(0).astype(int)
def logcpm(P, keep=None):
    lib = P.sum(axis=0); L = np.log2((P + 1) / (lib + 2) * 1e6 * 1.0); return L if keep is None else L.loc[keep]
def expressed(P, sex):
    lib = P.sum(axis=0); cut = 10 / np.median(lib) * 1e6; cpm = P / lib * 1e6; return P.index[(cpm > cut).sum(axis=1) >= min((sex == 'female').sum(), (sex == 'male').sum())]
def ols(L, X, j):
    """Per-gene OLS; returns estimate, se, t, p and df for column j of X. L: genes x samples."""
    Y = L.values; n, k = X.shape; XtXi = np.linalg.inv(X.T @ X); B = Y @ X @ XtXi; R = Y - B @ X.T; s2 = (R ** 2).sum(axis=1) / (n - k); se = np.sqrt(s2 * XtXi[j, j]); t = B[:, j] / se
    return pd.DataFrame(dict(gene=L.index, est=B[:, j], se=se, t=t, p=2 * stats.t.sf(np.abs(t), n - k), df=n - k))
def bh(p):
    p = np.asarray(p); o = np.argsort(p); r = np.empty(len(p)); r[o] = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]; return np.minimum(r, 1)
def sex_de(P, meta, sex=None):
    """Female - male (log2), adjusted for age, in one group of donors. `sex` may be a permuted label vector."""
    m = meta.loc[P.columns]; sx = (m.sex if sex is None else pd.Series(sex, index=m.index))
    keep = expressed(P, m.sex); L = logcpm(P, keep); X = np.column_stack([np.ones(len(m)), (sx == 'female').astype(float), m.age.values])
    r = ols(L, X, 1); r['FDR'] = bh(r.p.values); r['chr_class'] = chr_class(r.gene.values); return r
