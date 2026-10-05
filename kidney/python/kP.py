"""Plan section P: the endothelin-1-producing, low-flow-response endothelial state in kidney. Everything as fixed in the plan."""
from k_common import *
import warnings; warnings.filterwarnings('ignore')
pd.set_option('display.width', 250, 'display.max_columns', 40, 'display.max_rows', 200, 'display.float_format', lambda v: '%.3g' % v)
SIG = pickle.load(open('/home/claude/state/out/signatures.pkl', 'rb'))['SIG']['all']
cells = []
for ds in ('snRNA', 'scRNA'):
    d = load(ds); o = d['obs'].set_index('_index').loc[d['ec_cells']]; M = d['ec']; sym = d['sym']
    E = pd.DataFrame(dict(cell=d['ec_cells'], dataset=ds, donor=o.donor_id.values, sex=o.sex.values, age=o.age.values, dis=o.dis.values, label=o.label.values,
                          cat=o.disease_category.values, umi=np.asarray(M.sum(axis=0)).ravel()))
    E['EDN1'] = gsum(M, sym, ['EDN1']); E['flow'] = gsum(M, sym, FLOW); E['up'] = gsum(M, sym, SIG['up']); E['down'] = gsum(M, sym, SIG['down']); E['amb'] = gsum(M, sym, AMBIENT)
    print(ds, 'endothelial cells', len(E), '| signature genes found:', int(np.isin(SIG['up'], sym).sum()), 'of', len(SIG['up']), 'up;', int(np.isin(SIG['down'], sym).sum()), 'of', len(SIG['down']), 'down')
    cells.append(E)
E = pd.concat(cells, ignore_index=True); E = E[E.label.isin(ECLAB)].copy(); E['pop'] = E.label.map(ECLAB)
E = E[(E.age >= 18) & E.sex.isin(['female', 'male'])]
# a donor present in both datasets is used once, in the dataset with more endothelial cells
n = E.groupby(['donor', 'dataset']).size().unstack(fill_value=0); keep = n.idxmax(axis=1)
both = (n > 0).sum(axis=1) == 2; print('donors in both datasets:', int(both.sum()), '-> kept in', keep[both].value_counts().to_dict())
E = E[E.dataset.values == keep.loc[E.donor].values].copy()
E['pos'] = E.EDN1 > 0; E['score'] = np.log2((E.up + 1) / (E.down + 1)); E['flow_score'] = np.log1p(1e4 * E.flow / E.umi); E['amb_score'] = np.log1p(1e4 * E.amb / E.umi)
E['stratum'] = E.dis + ', ' + E.dataset; E['assay'] = E.dataset
PRIM = 'peritubular capillary'
thr = E[E['pop'] == PRIM].groupby('dataset').score.quantile(0.90); E['hi'] = E.score > E.dataset.map(thr)
print('\ncells by population and dataset:\n', pd.crosstab(E['pop'], E.dataset).to_string())
print('\nEDN1-positive (%), by population and dataset:\n', (100 * E.groupby(['pop', 'dataset']).pos.mean()).unstack().round(2).to_string())
print('median UMI:', E.groupby('dataset').umi.median().to_dict())

def within(mask, seed=2033):
    """Per donor: AUC of the state score (EDN1+ vs EDN1-) within depth deciles; flow-response and ambient difference vs depth-matched negatives."""
    rng = np.random.default_rng(seed); rows = []
    for s, d in E[mask].groupby('donor'):
        if d.pos.sum() < MIN_POS: continue
        dec = np.ceil(stats.rankdata(d.umi, method='ordinal') * 10 / len(d)); num = den = 0; fp = []; fn = []; ap = []; an = []
        for k in np.unique(dec):
            a = d[(dec == k) & d.pos.values]; b = d[(dec == k) & ~d.pos.values]
            if len(a) and len(b):
                num += stats.mannwhitneyu(a.score, b.score).statistic; den += len(a) * len(b)
                bs = b.sample(min(len(b), 3 * len(a)), random_state=int(rng.integers(1e9)))
                fp.extend(a.flow_score); fn.extend(bs.flow_score); ap.extend(a.amb_score); an.extend(bs.amb_score)
        if den: rows.append(dict(donor=s, dataset=d.dataset.iloc[0], dis=d.dis.iloc[0], sex=d.sex.iloc[0], auc=num / den, flow_diff=np.mean(fp) - np.mean(fn), amb_diff=np.mean(ap) - np.mean(an), n_pos=int(d.pos.sum()), n=len(d)))
    return pd.DataFrame(rows)
def report(lab, w):
    if len(w) < 3: print(f'  {lab:58s} too few donors ({len(w)})'); return None
    t1 = stats.ttest_1samp(w.auc, 0.5); t2 = stats.ttest_1samp(w.flow_diff, 0)
    print(f'  {lab:58s} donors {len(w):3d} | AUC {w.auc.mean():.3f} (p = {t1.pvalue:.2e}) | flow diff {w.flow_diff.mean():+.4f}, lower in {int((w.flow_diff < 0).sum())}/{len(w)} (p = {t2.pvalue:.2e})')
    return dict(analysis=lab, donors=len(w), auc=w.auc.mean(), auc_lo=w.auc.mean() - 1.96 * w.auc.sem(), auc_hi=w.auc.mean() + 1.96 * w.auc.sem(), p_auc=t1.pvalue, flow_diff=w.flow_diff.mean(), flow_lower=int((w.flow_diff < 0).sum()), p_flow=t2.pvalue)
rows = []
print('\n== PRIMARY (plan P; alpha 0.025 each): peritubular capillary endothelium ==')
W = within(E['pop'] == PRIM); rows.append(report('PRIMARY peritubular capillary', W))
print('\n== Ambient control ==')
ta = stats.ttest_1samp(W.amb_diff, 0); print(f'  tubular ambient score, EDN1+ minus depth-matched EDN1-: {W.amb_diff.mean():+.4f} (p = {ta.pvalue:.3g}); cells with any tubular marker count: {100 * (E[E["pop"] == PRIM].amb > 0).mean():.1f}%')
W2 = within((E['pop'] == PRIM) & (E.amb == 0)); rows.append(report('PRIMARY after excluding cells with any tubular marker', W2))
print('\n== Secondary: by dataset, disease group and population ==')
for ds in ('snRNA', 'scRNA'): rows.append(report(f'peritubular, {ds}', W[W.dataset == ds]))
for g in ('healthy', 'CKD', 'AKI'): rows.append(report(f'peritubular, {g}', W[W.dis == g]))
WP = {PRIM: W}
for p in ECLAB.values():
    if p != PRIM: WP[p] = within(E['pop'] == p); rows.append(report(p, WP[p]))
Wall = within(E['pop'].notna()); rows.append(report('all blood endothelial cells', Wall))

# ---- sex (secondary, exploratory) ----
def donor_table(mask):
    g = E[mask].groupby('donor'); d = pd.DataFrame(dict(n=g.size(), n_hi=g.hi.sum(), n_pos=g.pos.sum(), log_umi=g.umi.apply(lambda v: np.median(np.log(v)))))
    d = d[d.n >= MIN_CAP]; d['abundance'] = np.log((d.n_hi + 0.5) / (d.n - d.n_hi + 0.5)); d['edn1'] = np.log((d.n_pos + 0.5) / (d.n - d.n_pos + 0.5))
    d['pct_hi'] = 100 * d.n_hi / d.n; d['pct_pos'] = 100 * d.n_pos / d.n
    m = E[mask].groupby('donor')[['sex', 'age', 'dis', 'dataset', 'stratum', 'assay', 'cat']].first(); return d.join(m).reset_index()
def sextest(D, y, extra=()):
    per = pd.DataFrame([dict(r, stratum=s) for s, q in D.groupby('stratum') if (r := sex_fit(q, y, extra))]); return per, (meta_one(per) if len(per) else None)
def line(lab, res):
    per, b = res
    if b is None: print(f'  {lab:58s} not estimable'); return None
    print(f"  {lab:58s} {b['est']:+.3f} ({b['lo']:+.3f} to {b['hi']:+.3f})  p = {b['p']:.4f}  strata {b['k']}  F/M = {b['n_w']}/{b['n_m']}"); return dict(analysis=lab, **b)
srows = []; D0 = donor_table(E['pop'] == PRIM)
print('\n== Sex (secondary, exploratory): female - male, SD units, age-adjusted, within dataset x disease, pooled ==')
print(D0.groupby(['stratum', 'sex']).size().unstack(fill_value=0).to_string())
print('median % state-high:', D0.groupby(['stratum', 'sex']).pct_hi.median().round(1).to_dict()); print('median % EDN1-positive:', D0.groupby(['stratum', 'sex']).pct_pos.median().round(2).to_dict())
sprim = sextest(D0, 'abundance'); srows.append(line('abundance of state-high peritubular cells', sprim)); print(sprim[0].round(3).to_string(index=False))
sedn = sextest(D0, 'edn1'); srows.append(line('share of EDN1-positive peritubular cells', sedn)); print(sedn[0].round(3).to_string(index=False))
srows.append(line('abundance + median log UMI', sextest(D0, 'abundance', ('log_umi',))))
DP = {PRIM: D0}
for p in ECLAB.values():
    if p != PRIM:
        DP[p] = donor_table(E['pop'] == p); srows.append(line(f'abundance, {p}', sextest(DP[p], 'abundance'))); srows.append(line(f'EDN1-positive share, {p}', sextest(DP[p], 'edn1')))
# sex x CKD interaction within dataset (healthy and CKD donors)
def inter(D, y):
    out = []
    for ds, q in D[D.dis.isin(['healthy', 'CKD'])].groupby('dataset'):
        q = q[np.isfinite(q[y]) & q.age.notna()]
        if q.groupby(['sex', 'dis']).size().reindex(pd.MultiIndex.from_product([['female', 'male'], ['healthy', 'CKD']]), fill_value=0).min() < 3: continue
        yy = ((q[y] - q[y].mean()) / q[y].std()).values; f = (q.sex == 'female').astype(float).values; c = (q.dis == 'CKD').astype(float).values
        X = np.column_stack([np.ones(len(q)), f, c, f * c, q.age.values]); nn, k = X.shape; XtXi = np.linalg.inv(X.T @ X); b = XtXi @ X.T @ yy; r = yy - X @ b; se = np.sqrt(r @ r / (nn - k) * np.diag(XtXi))
        out.append(dict(stratum=ds, est=b[3], se=se[3], p=2 * stats.t.sf(abs(b[3] / se[3]), nn - k), n_w=int(f.sum()), n_m=int((1 - f).sum()), sex_in_healthy=b[1], sex_in_CKD=b[1] + b[3]))
    per = pd.DataFrame(out); return per, (meta_one(per) if len(per) else None)
print('\n== Sex x CKD interaction (difference of the female - male difference, CKD minus healthy) ==')
for p in (PRIM, 'glomerular capillary'):
    for y in ('abundance', 'edn1'):
        r = inter(DP[p], y); srows.append(line(f'sex x CKD, {y}, {p}', r)); print(r[0].round(3).to_string(index=False) if len(r[0]) else '')
res = pd.DataFrame([r for r in rows if r]); sres = pd.DataFrame([r for r in srows if r])
res.to_csv(OUTK + '/tables/S3P_state_in_kidney.csv', index=False); sres.to_csv(OUTK + '/tables/S3P_state_sex.csv', index=False); D0.to_csv(OUTK + '/tables/S3P_donors_peritubular.csv', index=False)
pickle.dump(dict(table=res, sex=sres, W=WP, Wall=Wall, W_noamb=W2, donors=DP, thr=thr, E=E[['donor', 'dataset', 'dis', 'sex', 'pop', 'score', 'hi', 'pos', 'umi', 'flow_score', 'amb']], pos_share=E.groupby(['pop', 'dataset']).pos.mean()), open(OUTK + '/P.pkl', 'wb'))
