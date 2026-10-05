"""Plan section Q (exploratory parts Q2, Q3, Q5, Q6, Q7). Q1 and Q4 need Hallmark gene sets and limma (R script K3)."""
from kQ_engine import *
import warnings; warnings.filterwarnings('ignore')
pd.set_option('display.width', 250, 'display.max_columns', 40, 'display.max_rows', 300, 'display.float_format', lambda v: '%.3g' % v)
rng = np.random.default_rng(2051)
sn = dataset('snRNA'); sc = dataset('scRNA'); META = {'snRNA': sn['meta'], 'scRNA': sc['meta'][~sc['meta'].index.isin(sn['meta'].index)]}   # single-cell donors not in the single-nucleus set
TYPES = ['PT', 'TAL', 'DCT/CNT', 'PC', 'IC', 'EC', 'stroma', 'immune']
PB = {}
for ct in TYPES + ['POD/PEC', 'thin limb']:
    for ds in ('snRNA', 'scRNA'):
        P, nc = pseudobulk(ds, 'subclassl1', [k for k, v in HARM[ds].items() if v == ct]); ok = nc[nc >= 30].index.intersection(META[ds].index); PB[(ct, ds)] = P[ok]
# endothelial subtypes from the cell-level matrices
ECSUB = {}
for ds in ('snRNA', 'scRNA'):
    d = dataset(ds); o = d['obs'].set_index('_index').loc[d['ec_cells']]; lab = o.label.values; don = o.donor_id.values
    for l, name in ECLAB.items():
        sel = np.where((lab == l) & np.isin(don, META[ds].index))[0]
        if not len(sel): continue
        cat = pd.Categorical(don[sel]); ind = sp.csr_matrix((np.ones(len(sel)), (sel, cat.codes)), shape=(len(lab), len(cat.categories)))
        P = pd.DataFrame((d['ec'] @ ind).toarray(), index=d['sym'], columns=list(cat.categories)).groupby(level=0).sum(); nc = pd.Series(cat.value_counts()); ECSUB[(name, ds)] = P[nc[nc >= 30].index]
def strata(P, ds, groups=('healthy', 'CKD', 'AKI')):
    m = META[ds].loc[P.columns]
    for g in groups:
        q = m[m.dis == g]
        if (q.sex == 'female').sum() >= MIN_PER_SEX and (q.sex == 'male').sum() >= MIN_PER_SEX: yield g, q

# ================= Q2: where sex differences lie (healthy donors; single-nucleus vs independent single-cell donors) =================
def q2(perm=False, donors=None):
    out = {}
    for ct in TYPES:
        de = {}
        for ds in ('snRNA', 'scRNA'):
            P = PB[(ct, ds)]; m = META[ds].loc[P.columns]; m = m[m.dis == 'healthy']
            if donors is not None: m = m[m.index.isin(donors[ds])]
            if (m.sex == 'female').sum() < MIN_PER_SEX or (m.sex == 'male').sum() < MIN_PER_SEX: de = None; break
            sx = rng.permutation(m.sex.values) if perm else None
            de[ds] = sex_de(P[m.index], META[ds], sx).set_index('gene'); de[ds].attrs['n'] = (int((m.sex == 'female').sum()), int((m.sex == 'male').sum()))
        if de is None: continue
        j = de['snRNA'].join(de['scRNA'], lsuffix='_a', rsuffix='_b', how='inner'); a = j[j.chr_class_a == 'autosome']
        rep = ((a.FDR_a < 0.05) & (a.p_b < 0.05) & (np.sign(a.t_a) == np.sign(a.t_b))) | ((a.FDR_b < 0.05) & (a.p_a < 0.05) & (np.sign(a.t_a) == np.sign(a.t_b)))
        xy = j[j.chr_class_a.isin(['X', 'Y'])]
        out[ct] = dict(cell_type=ct, sn_F=de['snRNA'].attrs['n'][0], sn_M=de['snRNA'].attrs['n'][1], sc_F=de['scRNA'].attrs['n'][0], sc_M=de['scRNA'].attrs['n'][1], autosomal_genes=len(a),
                       rho=stats.spearmanr(a.t_a, a.t_b)[0], FDR05_sn=int((a.FDR_a < 0.05).sum()), FDR05_sc=int((a.FDR_b < 0.05).sum()), replicated=int(rep.sum()),
                       XY_both=int(((xy.FDR_a < 0.05) & (xy.FDR_b < 0.05) & (np.sign(xy.t_a) == np.sign(xy.t_b))).sum()), genes=list(a.index[rep]), table=j)
    return out
print('== Q2. Replicated autosomal sex-biased genes by kidney cell type (healthy donors) ==')
R2 = q2(); t2 = pd.DataFrame([{k: v for k, v in r.items() if k not in ('genes', 'table')} for r in R2.values()])
NPERM = 30; null = pd.DataFrame([{ct: (r['replicated'], r['rho']) for ct, r in q2(perm=True).items()} for _ in range(NPERM)])
t2['null_replicated_max'] = [max(v[0] for v in null[ct]) for ct in t2.cell_type]; t2['null_replicated_median'] = [np.median([v[0] for v in null[ct]]) for ct in t2.cell_type]
t2['null_rho_max'] = [max(v[1] for v in null[ct]) for ct in t2.cell_type]; print(t2.sort_values('replicated', ascending=False).to_string(index=False))
for ct, r in R2.items():
    if r['genes']: j = r['table'].loc[r['genes']]; print(f"  {ct}: " + ', '.join(f"{g}({'F' if j.loc[g, 't_a'] > 0 else 'M'})" for g in j.sort_values('p_a').index[:40]))
# equal donors: those with every core cell type
core = ['PT', 'TAL', 'PC', 'IC', 'EC', 'stroma']; eq = {ds: set.intersection(*[set(PB[(ct, ds)].columns) for ct in core]) for ds in ('snRNA', 'scRNA')}
R2e = q2(donors=eq); t2e = pd.DataFrame([{k: v for k, v in r.items() if k not in ('genes', 'table')} for r in R2e.values()])
nulle = pd.DataFrame([{ct: r['replicated'] for ct, r in q2(perm=True, donors=eq).items()} for _ in range(NPERM)])
t2e['null_replicated_max'] = [nulle[ct].max() for ct in t2e.cell_type]; print('\n-- the same on donors that have all six core cell types --'); print(t2e[['cell_type', 'sn_F', 'sn_M', 'sc_F', 'sc_M', 'rho', 'replicated', 'null_replicated_max', 'XY_both']].to_string(index=False))

# ================= kidney-wide per-gene sex effect (all strata pooled), used by Q5 and Q7 =================
def pooled_de(getP, perm=False):
    per = []
    for ds in ('snRNA', 'scRNA'):
        P = getP(ds)
        if P is None: continue
        for g, q in strata(P, ds):
            sx = rng.permutation(q.sex.values) if perm else None
            r = sex_de(P[q.index], META[ds], sx); r['w'] = 1 / r.se ** 2; per.append(r[['gene', 'est', 'w']])
    if not per: return None
    a = pd.concat(per); a['we'] = a.est * a.w; s = a.groupby('gene').agg(we=('we', 'sum'), w=('w', 'sum'), k=('est', 'size')); s = s[s.k == s.k.max()]
    return (s.we / s.w) / np.sqrt(1 / s.w)                               # z per gene

# ================= Q5: X-Y dosage in every kidney cell type =================
XY = [('KDM6A', 'UTY'), ('KDM5C', 'KDM5D'), ('USP9X', 'USP9Y'), ('DDX3X', 'DDX3Y'), ('EIF1AX', 'EIF1AY'), ('ZFX', 'ZFY'), ('RPS4X', 'RPS4Y1'), ('NLGN4X', 'NLGN4Y')]
def donor_fit(tab, y):
    per = pd.DataFrame([dict(r, stratum=s) for s, q in tab.groupby('stratum') if (r := sex_fit(q, y))]); return meta_one(per) if len(per) else None
rows5 = []; ALLP = dict(PB); ALLP.update(ECSUB)
for (ct, ds) in list(ALLP):
    pass
for ct in TYPES + list(ECLAB.values()):
    tabs = []
    for ds in ('snRNA', 'scRNA'):
        P = ALLP.get((ct, ds))
        if P is None or P.shape[1] < 6: continue
        cpm = P / P.sum(axis=0) * 1e6; m = META[ds].loc[P.columns]; t = pd.DataFrame(dict(sex=m.sex, age=m.age, stratum=m.dis + ', ' + ds, assay=ds))
        for x, y in XY:
            if x in cpm.index:
                t['X_' + x] = np.log2(cpm.loc[x] + 1).values; t['XY_' + x] = np.log2(cpm.loc[x] + (cpm.loc[y] if y in cpm.index else 0) + 1).values
                t['Yshare_' + x] = (cpm.loc[y] / (cpm.loc[x] + cpm.loc[y] + 1e-9)).values if y in cpm.index else np.nan
        tabs.append(t)
    if not tabs: continue
    T = pd.concat(tabs)
    for x, y in XY:
        if 'X_' + x not in T: continue
        # raw log2 differences (not SD units): fit on the unscaled outcome
        def raw(col):
            per = []
            for s, q in T.groupby('stratum'):
                q = q[np.isfinite(q[col])]; nf, nm = (q.sex == 'female').sum(), (q.sex == 'male').sum()
                if nf < 3 or nm < 3: continue
                X = np.column_stack([np.ones(len(q)), (q.sex == 'female').astype(float), q.age]); b, res, *_ = np.linalg.lstsq(X, q[col].values, rcond=None); r = q[col].values - X @ b
                se = np.sqrt(r @ r / (len(q) - 3) * np.linalg.inv(X.T @ X)[1, 1]); per.append(dict(est=b[1], se=se, n_w=nf, n_m=nm))
            return meta_one(pd.DataFrame(per)) if per else None
        a, b = raw('X_' + x), raw('XY_' + x)
        if a and b: rows5.append(dict(cell_type=ct, pair=f'{x}/{y}', X_only=a['est'], X_p=a['p'], X_plus_Y=b['est'], lo=b['lo'], hi=b['hi'], p=b['p'], strata=b['k'], F=b['n_w'], M=b['n_m'], Y_share_in_males=T.loc[T.sex == 'male', 'Yshare_' + x].median()))
t5 = pd.DataFrame(rows5); print('\n== Q5. X-Y pairs: female - male (log2 CPM), X copy alone and X + Y combined, pooled across strata ==')
print(t5.pivot(index='pair', columns='cell_type', values='X_plus_Y').round(2).to_string()); print('\nY share of the pair in males (median):'); print(t5.pivot(index='pair', columns='cell_type', values='Y_share_in_males').round(2).to_string())
print('\nendothelium and proximal tubule in detail:'); print(t5[t5.cell_type.isin(['peritubular capillary', 'glomerular capillary', 'EC', 'PT'])].sort_values(['pair', 'cell_type']).to_string(index=False))

# ================= Q3: endothelial programs with equivalence bounds =================
PROG = dict(NO_eNOS=['NOS3', 'KLF2', 'KLF4', 'CAV1', 'GCH1', 'SLC7A1', 'DDAH1', 'DDAH2'], Endothelin_ACE=['EDN1', 'ECE1', 'EDNRB', 'ACE'], Prostacyclin=['PTGS1', 'PTGS2', 'PTGIS', 'PLA2G4A'],
            Barrier=['CLDN5', 'CDH5', 'ESAM', 'OCLN', 'TJP1', 'JAM2', 'PECAM1'], FA_transport=['CD36', 'FABP4', 'FABP5', 'LPL', 'GPIHBP1'], Angiogenic_tip=['ESM1', 'APLN', 'DLL4', 'ANGPT2', 'PGF', 'KCNE3'],
            X_escape_control=['KDM6A', 'KDM5C', 'DDX3X', 'EIF1AX', 'ZFX', 'USP9X', 'JPX'])
rows3 = []
for pop in ('peritubular capillary', 'glomerular capillary'):
    tabs = []
    for ds in ('snRNA', 'scRNA'):
        P = ECSUB.get((pop, ds))
        if P is None: continue
        L = logcpm(P); m = META[ds].loc[P.columns]; t = pd.DataFrame(dict(sex=m.sex, age=m.age, stratum=m.dis + ', ' + ds, assay=ds))
        for k, g in PROG.items():
            g = [x for x in g if x in L.index]; z = L.loc[g].sub(L.loc[g].mean(axis=1), axis=0).div(L.loc[g].std(axis=1).replace(0, np.nan), axis=0); t[k] = z.mean(axis=0).values
        tabs.append(t)
    T = pd.concat(tabs)
    for k in PROG:
        per = pd.DataFrame([dict(r, stratum=s) for s, q in T.groupby('stratum') if (r := sex_fit(q, k))]); b = meta_one(per)
        if b is None: continue
        se = (b['hi'] - b['lo']) / 3.92; ptost = max(stats.norm.sf((b['est'] + 0.8) / se), stats.norm.cdf((b['est'] - 0.8) / se)); rows3.append(dict(population=pop, program=k, est=b['est'], lo=b['lo'], hi=b['hi'], p=b['p'], p_equiv_0_8=ptost, strata=b['k'], F=b['n_w'], M=b['n_m'], I2=b['I2']))
t3 = pd.DataFrame(rows3)
def holm(p):
    p = np.asarray(p); o = np.argsort(p); a = np.maximum.accumulate((len(p) - np.arange(len(p))) * p[o]); r = np.empty(len(p)); r[o] = np.minimum(a, 1); return r
for pop, q in t3.groupby('population'):
    i = q.index[q.program != 'X_escape_control']; t3.loc[i, 'p_holm'] = holm(t3.loc[i, 'p'].values)
print('\n== Q3. Endothelial programs, female - male (SD units; age-adjusted; pooled across dataset x disease strata) =='); print(t3.to_string(index=False))

# ================= Q6: injury-state composition by sex =================
STATES = {'snRNA': {'adaptive PT (aPT) among PT': (['aPT'], 'PT'), 'degenerative PT (dPT) among PT': (['dPT'], 'PT'), 'adaptive TAL among TAL': (['aTAL1', 'aTAL2'], 'TAL'),
                    'activated fibroblasts (aFIB, MYOF) among fibroblasts': (['aFIB', 'MYOF'], 'FIB'), 'degenerative states among all cells': (None, None)},
          'scRNA': {'adaptive PT (aPT) among PT': (['aPT'], 'PT'), 'degenerative PT (dPT) among PT': (['dPT'], 'PT'), 'adaptive TAL among TAL': (['aTAL1', 'aTAL2'], 'TAL'), 'degenerative states among all cells': (None, None)}}
tabs = []
for ds in ('snRNA', 'scRNA'):
    o = dataset(ds)['obs']; o = o[o.donor_id.isin(META[ds].index)]
    for name, (labs, parent) in STATES[ds].items():
        par = o if parent is None else o[o['subclass.l1'] == parent]; hit = par.label.str.match(r'^d[A-Z]') if labs is None else par.label.isin(labs)
        g = pd.DataFrame(dict(n=par.groupby('donor_id').size(), k=hit.groupby(par.donor_id).sum())); g = g[g.n >= 100]; m = META[ds].loc[g.index]
        tabs.append(pd.DataFrame(dict(state=name, y=np.log((g.k + 0.5) / (g.n - g.k + 0.5)), pct=100 * g.k / g.n, sex=m.sex, age=m.age, dis=m.dis, stratum=m.dis + ', ' + ds, assay=ds, dataset=ds)))
T6 = pd.concat(tabs); rows6 = []
for st, q in T6.groupby('state'):
    per = pd.DataFrame([dict(r, stratum=s) for s, d in q.groupby('stratum') if (r := sex_fit(d, 'y'))]); b = meta_one(per)
    row = dict(state=st, sex_est=b['est'], lo=b['lo'], hi=b['hi'], p=b['p'], strata=b['k'], F=b['n_w'], M=b['n_m'])
    inter = []
    for ds, d in q[q.dis.isin(['healthy', 'CKD'])].groupby('dataset'):
        d = d[d.age.notna()]
        if d.groupby(['sex', 'dis']).size().min() < 3 or d.groupby(['sex', 'dis']).ngroups < 4: continue
        yy = ((d.y - d.y.mean()) / d.y.std()).values; f = (d.sex == 'female').astype(float).values; c = (d.dis == 'CKD').astype(float).values; X = np.column_stack([np.ones(len(d)), f, c, f * c, d.age.values])
        XtXi = np.linalg.inv(X.T @ X); bb = XtXi @ X.T @ yy; r = yy - X @ bb; se = np.sqrt(r @ r / (len(d) - 5) * np.diag(XtXi)); inter.append(dict(est=bb[3], se=se[3], n_w=int(f.sum()), n_m=int((1 - f).sum()), ckd=bb[2]))
    if inter: ib = meta_one(pd.DataFrame(inter)); row.update(sexCKD_est=ib['est'], sexCKD_p=ib['p'])
    med = q.groupby(['dis', 'sex']).pct.median(); row.update({f'median_pct_{a}_{s}': med.get((a, s), np.nan) for a in ('healthy', 'CKD') for s in ('female', 'male')}); rows6.append(row)
t6 = pd.DataFrame(rows6); print('\n== Q6. Injury-state composition: female - male (SD of logit fraction; age-adjusted; pooled) and sex x CKD =='); print(t6.to_string(index=False))

# ================= Q7: kidney - heart links =================
HD = pickle.load(open('/home/claude/data/D.pkl', 'rb')); hde = HD['stage3_de']['de']
def heart_z(ct): h = hde[ct]; h = h[h.chr_class == 'autosome']; return h.set_index('gene').z_comb
pairs = [('peritubular capillary', 'capillary'), ('glomerular capillary', 'capillary'), ('EC', 'capillary'), ('PT', 'cardiomyocyte'), ('PT', 'capillary'), ('stroma', 'fibroblast'), ('immune', 'myeloid')]
getP = lambda k: (lambda ds: ALLP.get((k, ds)))
KZ = {k: pooled_de(getP(k)) for k in set(p[0] for p in pairs) | set(TYPES)}
rows7 = []; NP7 = 200
for kct, hct in pairs:
    hz = heart_z(hct); kz = KZ[kct]; g = kz.index.intersection(hz.index); g = g[chr_class(g.values) == 'autosome']; rho = stats.spearmanr(kz[g], hz[g])[0]
    nullr = []
    for _ in range(NP7):
        pz = pooled_de(getP(kct), perm=True); gg = pz.index.intersection(g); nullr.append(stats.spearmanr(pz[gg], hz[gg])[0])
    rows7.append(dict(kidney=kct, heart=hct, genes=len(g), rho=rho, null_sd=np.std(nullr), p_perm=(1 + np.sum(np.abs(nullr) >= abs(rho))) / (NP7 + 1)))
t7a = pd.DataFrame(rows7); print('\n== Q7a. Concordance of autosomal sex effects, kidney vs heart (Spearman rho of per-gene z; permutation of sex within kidney strata) =='); print(t7a.to_string(index=False))
G55 = HD['stage3_cm_pathways']['genes']; G55 = G55[G55.consistent.astype(bool)] if 'consistent' in G55 else G55; hdir = np.sign(G55.set_index('gene').z_comb)
rows = []
for ct in TYPES + ['peritubular capillary', 'glomerular capillary']:
    kz = KZ.get(ct) if ct in KZ else pooled_de(getP(ct))
    if kz is None: continue
    g = hdir.index.intersection(kz.index); same = int((np.sign(kz[g]) == hdir[g]).sum()); rows.append(dict(kidney_cell_type=ct, genes_found=len(g), same_direction=same, pct=100 * same / len(g), p_binom=stats.binomtest(same, len(g), 0.5).pvalue,
                                                                                                           nominal_same=int(((np.sign(kz[g]) == hdir[g]) & (np.abs(kz[g]) > 1.96)).sum()), nominal_opposite=int(((np.sign(kz[g]) != hdir[g]) & (np.abs(kz[g]) > 1.96)).sum())))
t7b = pd.DataFrame(rows); print('\n== Q7b. The 55 cardiomyocyte sex genes in kidney cell types =='); print(t7b.to_string(index=False))
hx = HD['combined']['xy_comb']; hx = hx[hx.measure == 'X_plus_Y'].set_index('pair') if 'measure' in hx else None
if hx is not None:
    k = t5[t5.cell_type == 'peritubular capillary'].set_index('pair').X_plus_Y; j = pd.DataFrame(dict(heart_capillary=hx.est_comb, kidney_peritubular=k, kidney_glomerular=t5[t5.cell_type == 'glomerular capillary'].set_index('pair').X_plus_Y, kidney_PT=t5[t5.cell_type == 'PT'].set_index('pair').X_plus_Y)).dropna(subset=['heart_capillary', 'kidney_peritubular'])
    print('\n== Q7c. X + Y dosage difference (female - male, log2), heart capillary vs kidney =='); print(j.round(2).to_string()); print('Pearson r heart capillary vs kidney peritubular:', round(stats.pearsonr(j.heart_capillary, j.kidney_peritubular)[0], 2), '| n pairs', len(j))
else: j = None
for nm, t in (('S3Q2_celltype_ranking', t2), ('S3Q2_equal_donors', t2e), ('S3Q5_xy_dosage', t5), ('S3Q3_programs', t3), ('S3Q6_injury_states', t6), ('S3Q7a_concordance', t7a), ('S3Q7b_cm_genes', t7b)): t.to_csv(f'{OUTK}/tables/{nm}.csv', index=False)
pickle.dump(dict(q2=t2, q2e=t2e, q2_genes={k: v['genes'] for k, v in R2.items()}, q2_tables={k: v['table'] for k, v in R2.items()}, q5=t5, q3=t3, q6=t6, q6_donors=T6, q7a=t7a, q7b=t7b, q7c=j, KZ=KZ), open(OUTK + '/Q.pkl', 'wb'))
