"""M2-M4: abundance and intensity of the EDN1-independent state by sex; validation; continuum; strict capillary definitions.
Reads out/cells.pkl (from m1_signature.py). Writes out/state_results.pkl and tables/S3M_*.csv."""
from common import *
S = samples(); C = pd.read_pickle(os.path.join(OUT, "cells.pkl"))
C["is_cap"] = C.subtype == "capillary"; C["pos"] = C.EDN1 > 0
TAB = os.path.join(OUT, "tables"); os.makedirs(TAB, exist_ok=True)
R = {}


def thresholds(score, q):
    return C[C.is_cap].groupby("dataset_id")[score].quantile(q)


def donor_table(sub, score="score_cf", q=0.90, min_cells=MIN_CAP):
    """One row per donor for the cells in `sub` (boolean mask): abundance, intensity, mean score, depth."""
    thr = C.dataset_id.map(thresholds(score, q)); hi = (C[score] > thr) & sub
    g = C[sub].assign(hi=hi[sub]).groupby("sample")
    d = pd.DataFrame(dict(n=g.size(), n_hi=g.hi.sum(), mean_score=g[score].mean(), log_umi=g.umi.apply(lambda v: np.median(np.log(v))), pct_pos=100 * g.pos.mean()))
    d["intensity"] = C[hi].groupby("sample")[score].mean()
    d = d[d.n >= min_cells]
    d["abundance"] = np.log((d.n_hi + 0.5) / (d.n - d.n_hi + 0.5)); d["pct_hi"] = 100 * d.n_hi / d.n
    return d.join(S, how="inner").reset_index().rename(columns={"index": "sample"})


def line(label, r, extra=""):
    b = r.get("both")
    if b is None: print(f"  {label:58s} not estimable"); return None
    print(f"  {label:58s} {b['est']:+.3f} ({b['lo']:+.3f} to {b['hi']:+.3f})  p = {b['p']:.4f}  k = {b['k']}  F/M = {b['n_w']}/{b['n_m']}  same direction {b['same_dir']}/{b['k']}  I2 = {b['I2']:.0f}% {extra}")
    return dict(analysis=label, **b)


rows = []
print("\n== M2. State abundance and intensity, female - male (SD units) ==")
D0 = donor_table(C.is_cap); R["donors"] = D0
print("  realised % of capillary cells called state-high, by cohort:", (100 * D0.groupby("cohort").n_hi.sum() / D0.groupby("cohort").n.sum()).round(1).to_dict())
r = pooled(D0, "abundance"); R["primary"] = r; rows.append(line("PRIMARY abundance (cross-fitted signature, 90th pct)", r))
print(r["per"].assign(stratum=r["per"].stratum.map(short)).round(3).to_string(index=False))
for co in ("healthy", "disease"): print(f"    {co}: {r[co]['est']:+.3f} ({r[co]['lo']:+.3f} to {r[co]['hi']:+.3f}) p = {r[co]['p']:.4f}")
rows.append(line("abundance + median log UMI", pooled(D0, "abundance", ("log_umi",))))
rows.append(line("intensity (mean score of state-high cells)", pooled(D0, "intensity")))
rows.append(line("intensity + median log UMI", pooled(D0, "intensity", ("log_umi",))))
rows.append(line("mean score of all capillary cells", pooled(D0, "mean_score")))
rows.append(line("abundance, 80th percentile", pooled(donor_table(C.is_cap, q=0.80), "abundance")))
rows.append(line("abundance, 95th percentile", pooled(donor_table(C.is_cap, q=0.95), "abundance")))
Dall = donor_table(C.is_cap, score="score_all")
rows.append(line("abundance, all-donor signature (not cross-fitted)", pooled(Dall, "abundance")))
rows.append(line("abundance, all-donor signature + median log UMI", pooled(Dall, "abundance", ("log_umi",))))
rows.append(line("for comparison: share of EDN1-positive cells", pooled(D0.assign(lp=np.log((D0.pct_pos / 100 * D0.n + 0.5) / (D0.n - D0.pct_pos / 100 * D0.n + 0.5))), "lp")))
print("  median % state-high per donor:", D0.groupby(["cohort", "sex"]).pct_hi.median().round(2).to_dict())

# excluding EDN1-positive cells altogether: is the state more abundant among cells with no EDN1 transcript?
rows.append(line("abundance among EDN1-negative capillary cells only", pooled(donor_table(C.is_cap & ~C.pos), "abundance")))

print("\n== Validation: does the score mark EDN1-positive cells in held-out donors (no EDN1 in the score)? ==")
auc = []
for s, d in C[C.is_cap].groupby("sample"):
    if d.pos.sum() < MIN_POS: continue
    u = stats.mannwhitneyu(d.score_cf[d.pos], d.score_cf[~d.pos]).statistic / (d.pos.sum() * (~d.pos).sum())
    ud = stats.mannwhitneyu(d.umi[d.pos], d.umi[~d.pos]).statistic / (d.pos.sum() * (~d.pos).sum())
    # depth-stratified AUC: compare within depth deciles
    dec = np.ceil(stats.rankdata(d.umi, method="ordinal") * 10 / len(d)); num = den = 0
    for k in np.unique(dec):
        a = d.score_cf[(dec == k) & d.pos.values]; b = d.score_cf[(dec == k) & ~d.pos.values]
        if len(a) and len(b): num += stats.mannwhitneyu(a, b).statistic; den += len(a) * len(b)
    auc.append(dict(sample=s, auc=u, auc_depth=ud, auc_within_depth=num / den if den else np.nan, pct_hi_pos=100 * (d.score_cf[d.pos] > thresholds("score_cf", 0.9)[d.dataset_id.iloc[0]]).mean()))
auc = pd.DataFrame(auc).set_index("sample").join(S); R["auc"] = auc
for c in ("auc", "auc_within_depth", "auc_depth"):
    t = stats.ttest_1samp(auc[c].dropna(), 0.5); print(f"  {c:18s} mean {auc[c].mean():.3f} (SD {auc[c].std():.3f}), {len(auc[c].dropna())} donors, p vs 0.5 = {t.pvalue:.2e}; by cohort {auc.groupby('cohort')[c].mean().round(3).to_dict()}")
print(f"  EDN1-positive capillary cells that are state-high: mean {auc.pct_hi_pos.mean():.1f}% (10% expected if unrelated)")

print("\n== M3. Continuum ==")
ca = C[C.subtype.isin(["capillary", "arterial"])].copy()
ca["bin"] = ca.groupby("dataset_id").art_pos.transform(lambda v: np.ceil(stats.rankdata(v, method="ordinal") * 10 / len(v))).astype(int)
g = ca.groupby(["sample", "bin"]); B = pd.DataFrame(dict(n=g.size(), score=g.score_cf.mean(), flow=g.flow_score.mean(), pct_pos=100 * g.pos.mean(), pct_cap=100 * g.is_cap.mean())).reset_index().join(S, on="sample")
B = B[B.n >= 5]; R["bins"] = B
print("  by bin of arterial position (1 = most capillary-like, 10 = most arterial-like), mean of donor means:")
print(B.groupby("bin").agg(donors=("sample", "nunique"), pct_capillary=("pct_cap", "mean"), state_score=("score", "mean"), flow=("flow", "mean"), pct_EDN1_pos=("pct_pos", "mean")).round(3).to_string())
occ = ca.groupby(["sample", "bin"]).size().unstack(fill_value=0); occ = occ.div(occ.sum(axis=1), axis=0)
occd = pd.DataFrame(dict(top3=np.log((occ[[8, 9, 10]].sum(axis=1) + 0.01) / (1 - occ[[8, 9, 10]].sum(axis=1) + 0.01)), mean_pos=ca.groupby("sample").art_pos.mean(), n=ca.groupby("sample").size())).join(S).reset_index()
rows.append(line("M3 occupancy of the 3 most arterial-like bins (cap + art cells)", pooled(occd, "top3")))
rows.append(line("M3 mean arterial position (cap + art cells)", pooled(occd, "mean_pos")))
cor = []
for s, d in C[C.is_cap].groupby("sample"):
    if len(d) < MIN_CAP: continue
    cor.append(dict(sample=s, rho_art=stats.spearmanr(d.score_cf, d.art_pos).statistic, rho_flow=stats.spearmanr(d.score_cf, d.flow_score).statistic,
                    rho_art_all=np.nan))
cor = pd.DataFrame(cor).set_index("sample").join(S); R["cor"] = cor
for c, lab in (("rho_art", "state score vs arterial position"), ("rho_flow", "state score vs flow-response score")):
    v = cor[c].dropna(); t = stats.ttest_1samp(np.arctanh(v), 0)
    print(f"  within-donor Spearman, capillary cells: {lab:38s} mean rho {v.mean():+.3f}, {int((v > 0).sum())}/{len(v)} donors positive, p = {t.pvalue:.2e}")
h, edges = np.histogram(C.score_cf[C.is_cap], bins=30); R["hist"] = (h, edges)
print("  distribution of the state score among capillary cells (30 bins, counts):", h.tolist())

print("\n== M4. Capillary cells or arteriolar cells? abundance of state-high cells, female - male ==")
margin = C.z_capillary - np.maximum(C.z_arterial, C.z_venous)
strict = C.is_cap & (margin >= 1); atlas = C.is_cap & (C.cell_type == "capillary endothelial cell")
print("  strict capillary cells:", int(strict.sum()), f"({100 * strict.sum() / C.is_cap.sum():.0f}% of capillary-assigned); atlas-concordant:", int(atlas.sum()), f"({100 * atlas.sum() / C.is_cap.sum():.0f}%)")
print("  atlas capillary label available, cells by stratum:", C[atlas].join(S[["stratum"]], on="sample").stratum.map(short).value_counts().to_dict())
third = C[C.is_cap].groupby("dataset_id").art_pos.transform(lambda v: np.ceil(stats.rankdata(v, method="ordinal") * 3 / len(v))).reindex(C.index)
M4 = [("(i) capillary, original rule", C.is_cap), ("(ii) strict capillary (margin >= 1)", strict), ("(iii) atlas-concordant capillary", atlas),
      ("(iv) capillary, capillary-like third", C.is_cap & (third == 1)), ("(iv) capillary, middle third", C.is_cap & (third == 2)), ("(iv) capillary, arterial-like third", C.is_cap & (third == 3)),
      ("(v) arterial-assigned cells", C.subtype == "arterial"), ("(v) venous-assigned cells", C.subtype == "venous")]
R["m4"] = {}
for lab, mask in M4:
    d = donor_table(mask); r = pooled(d, "abundance"); R["m4"][lab] = r
    x = line(lab, r, extra=f"| {100 * d.n_hi.sum() / d.n.sum():.1f}% state-high")
    if x: rows.append(dict(x, pct_state_high=100 * d.n_hi.sum() / d.n.sum(), cells=int(d.n.sum())))
    r2 = pooled(d.assign(lp=np.log((d.pct_pos / 100 * d.n + 0.5) / (d.n - d.pct_pos / 100 * d.n + 0.5))), "lp"); R["m4"][lab + " EDN1"] = r2
    x = line("     share of EDN1-positive cells, same cells", r2, extra=f"| {(C[mask].pos.mean() * 100):.1f}% EDN1-positive")
    if x: rows.append(x)
res = pd.DataFrame([x for x in rows if x]); R["table"] = res
res.to_csv(os.path.join(TAB, "S3M_state_sex_differences.csv"), index=False)
D0.to_csv(os.path.join(TAB, "S3M_state_donors.csv"), index=False); B.to_csv(os.path.join(TAB, "S3M_continuum_bins.csv"), index=False); auc.to_csv(os.path.join(TAB, "S3M_validation_auc.csv"))
pickle.dump(R, open(os.path.join(OUT, "state_results.pkl"), "wb"))
