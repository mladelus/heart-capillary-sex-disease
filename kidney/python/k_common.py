"""Shared loading for the kidney analyses (plan sections P and Q)."""
import sys, os, pickle, re, numpy as np, pandas as pd, scipy.sparse as sp
from scipy import stats
sys.path.insert(0, '/home/claude/state'); os.environ.setdefault('HCAP_DATA', '/home/claude/proj/p1'); os.environ.setdefault('HCAP2_DATA', '/home/claude/proj/p2')
from common import sex_fit, meta_one, MIN_CAP, MIN_PER_SEX, MIN_POS, FLOW
C = '/home/claude/kidney/cache/'; OUTK = '/home/claude/kidney/out'; os.makedirs(OUTK + '/tables', exist_ok=True)
DEC = {'third': 25, 'fourth': 35, 'fifth': 45, 'sixth': 55, 'seventh': 65, 'eighth': 75, 'ninth': 85, 'second': 15}
AMBIENT = ['LRP2', 'CUBN', 'SLC34A1', 'UMOD']
ECLAB = {'EC-PTC': 'peritubular capillary', 'EC-GC': 'glomerular capillary', 'EC-AEA': 'arteriole', 'EC-DVR': 'descending vasa recta', 'EC-AVR': 'ascending vasa recta'}
def age_of(s):
    s = str(s); m = re.match(r'^(\d+)-year', s)
    if m: return float(m.group(1))
    for k, v in DEC.items():
        if s.startswith(k + ' decade'): return float(v)
    return np.nan
def load(ds):
    d = pickle.load(open(C + ds + '.pkl', 'rb')); o = d['obs'].copy()
    o['dataset'] = ds; o['age'] = o.development_stage.map(age_of); o['label'] = o['subclass.l2'] if ds == 'snRNA' else o['author_cell_type']
    o['dis'] = o.disease.map({'normal': 'healthy', 'chronic kidney disease': 'CKD', 'acute kidney failure': 'AKI'})
    d['obs'] = o; d['sym'] = np.array(d['genes'].feature_name.astype(str)); return d
def donor_meta(o):
    cols = ['donor_id', 'sex', 'age', 'dis', 'disease_category', 'eGFR', 'diabetes_history', 'hypertension', 'dataset']
    return o.groupby('donor_id').agg({c: 'first' for c in cols[1:]}).assign(cells=o.groupby('donor_id').size())
def gsum(M, sym, names):
    idx = np.where(np.isin(sym, list(names)))[0]
    return np.asarray(M[idx, :].sum(axis=0)).ravel() if len(idx) else np.zeros(M.shape[1])
