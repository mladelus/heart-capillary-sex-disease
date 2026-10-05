"""Convert the K2 exports (rds) into scipy sparse matrices and pandas tables."""
import sys, os, pickle, numpy as np, pandas as pd, scipy.sparse as sp
sys.path.insert(0, '/home/claude/data'); from rds import read_rds
U = '/mnt/user-data/uploads/heart_capillary_disease/kidney/export/'; C = '/home/claude/kidney/cache/'
def mat(path):
    m = read_rds(path)
    if not isinstance(m, dict): raise SystemExit(f'unexpected object {type(m)} in {path}')
    keys = {k.lower(): k for k in m}; i = np.asarray(m[keys['i']]); p = np.asarray(m[keys['p']]); x = np.asarray(m[keys['x']], dtype=np.float32); dim = [int(v) for v in m[keys['dim']]]
    dn = m[keys['dimnames']]; M = sp.csc_matrix((x, i, p), shape=dim)
    return M, list(dn[0]), list(dn[1])
for ds in ('snRNA', 'scRNA'):
    obs = read_rds(U + ds + '_obs.rds'); obs = pd.DataFrame(obs) if not isinstance(obs, pd.DataFrame) else obs
    genes = read_rds(U + ds + '_genes.rds'); genes = pd.DataFrame(genes) if not isinstance(genes, pd.DataFrame) else genes
    print(ds, 'obs', obs.shape, 'genes', genes.shape, list(obs.columns)[:60])
    out = dict(obs=obs, genes=genes)
    a, g, c1 = mat(U + ds + '_EC_counts_1.rds'); b, g2, c2 = mat(U + ds + '_EC_counts_2.rds'); assert g == g2
    out['ec'] = sp.hstack([a, b]).tocsr(); out['ec_cells'] = c1 + c2; out['gene_ids'] = g
    print('  EC', out['ec'].shape, 'nnz', out['ec'].nnz)
    for lab in ('subclassl1', 'subclassl2', 'celltype'):
        f = U + f'{ds}_pb_{lab}.rds'
        if os.path.exists(f):
            M, g3, cols = mat(f); assert g3 == g; out['pb_' + lab] = M.tocsc(); out['pb_' + lab + '_cols'] = cols; print('  pb', lab, M.shape, 'nnz', M.nnz)
    pickle.dump(out, open(C + ds + '.pkl', 'wb'), protocol=4)
