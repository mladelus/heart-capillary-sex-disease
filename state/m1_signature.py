"""M1: learn the signature of EDN1-positive capillary cells without reference to sex (cross-fitted by cohort),
then M2 scoring: per-cell sums of signature genes for every blood endothelial cell. Writes out/signatures.pkl and out/cells.pkl."""
from common import *
rng_seed = 2032
S = samples()
de = read_rds(os.path.join(P2, "stage3_de.rds"))["de"]
chrom = pd.concat([d[["gene", "chr_class"]] for d in de.values()]).drop_duplicates("gene").set_index("gene").chr_class.astype(str)

# ---- pass A: positive and matched-negative pseudobulks per donor ----
pos, neg, info = {}, {}, []
for k, (s, co, genes, M, tab) in enumerate(donors(S)):
    cap = np.where((tab.subtype == "capillary").values)[0]
    if len(cap) < MIN_CAP or "EDN1" not in genes: continue
    e = row(M, genes, "EDN1")[cap]; umi = tab.umi.values[cap]
    ip = np.where(e > 0)[0]
    if len(ip) < MIN_POS: continue
    dec = np.ceil(stats.rankdata(umi, method="ordinal") * 10 / len(umi)).astype(int)        # within-donor depth deciles
    rng = np.random.default_rng(rng_seed + k); ineg = []
    for dd in np.unique(dec[ip]):
        cand = np.where((e == 0) & (dec == dd))[0]; n = min(len(cand), 3 * int((dec[ip] == dd).sum()))
        if n: ineg.extend(rng.choice(cand, n, replace=False))
    ineg = np.array(ineg, dtype=int)
    if len(ineg) < MIN_POS: continue
    pos[s] = pd.Series(np.asarray(M[:, cap[ip]].sum(axis=1)).ravel(), index=genes)
    neg[s] = pd.Series(np.asarray(M[:, cap[ineg]].sum(axis=1)).ravel(), index=genes)
    info.append(dict(sample=s, cohort=co, n_cap=len(cap), n_pos=len(ip), n_neg=len(ineg), depth_pos=np.median(umi[ip]), depth_neg=np.median(umi[ineg])))
info = pd.DataFrame(info).set_index("sample")
Pm = pd.DataFrame(pos).fillna(0); Nm = pd.DataFrame(neg).fillna(0).reindex(Pm.index).fillna(0)
print("donors with >= 5 EDN1-positive capillary cells:", len(info), info.groupby("cohort").size().to_dict(),
      "| median depth pos/neg:", info.depth_pos.median(), info.depth_neg.median())

EXCLUDE = set(ENDOTHELIN) | set(VESSEL_MARKERS) | set(FLOW)


def learn(ids, label):
    P, N = Pm[ids].drop(index="EDN1"), Nm[ids].drop(index="EDN1")
    keep = ((P + N) >= 10).mean(axis=1) >= 0.5
    P, N = P[keep], N[keep]
    d = np.log2((P + 0.5) / (P.sum() + 1) * 1e6) - np.log2((N + 0.5) / (N.sum() + 1) * 1e6)
    t, p = stats.ttest_1samp(d.values, 0, axis=1)
    o = np.argsort(p); fdr = np.empty(len(p)); fdr[o] = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    T = pd.DataFrame(dict(gene=d.index, mean_log2_diff=d.mean(axis=1).values, t=t, p=p, FDR=np.minimum(fdr, 1), donors=len(ids), mean_log2cpm=np.log2(((P + N).sum(axis=1) + 0.5) / ((P + N).sum().sum() + 1) * 1e6).values))
    T["chr_class"] = T.gene.map(chrom).fillna("unknown")
    T["excluded"] = T.gene.isin(EXCLUDE) | T.chr_class.isin(["X", "Y"])
    sig = T[(T.FDR < 0.05) & ~T.excluded]
    up, down = sorted(sig.gene[sig.mean_log2_diff > 0]), sorted(sig.gene[sig.mean_log2_diff < 0])
    print(f"signature learned in {label}: {len(ids)} donors, {len(T)} genes tested, {int((T.FDR < 0.05).sum())} at FDR < 0.05; after exclusions up {len(up)}, down {len(down)}")
    return dict(table=T.sort_values("p"), up=up, down=down, donors=list(ids))


SIG = {"healthy": learn(info.index[info.cohort == "healthy"], "healthy"), "disease": learn(info.index[info.cohort == "disease"], "disease"), "all": learn(info.index, "all donors")}
ov = lambda a, b: len(set(a) & set(b))
print("overlap healthy vs disease signatures: up", ov(SIG["healthy"]["up"], SIG["disease"]["up"]), "down", ov(SIG["healthy"]["down"], SIG["disease"]["down"]))
old = read_rds(os.path.join(P2, "stage3_edn1_cells.rds"))["de"]; oldsig = set(old.gene[old["adj.P.Val"] < 0.05])
newsig = set(SIG["all"]["table"].query("FDR < 0.05").gene)
print("all-donor signature vs script 19 (limma): script 19", len(oldsig), "| here", len(newsig), "| shared", len(oldsig & newsig))
pickle.dump(dict(SIG=SIG, info=info), open(os.path.join(OUT, "signatures.pkl"), "wb"))

# ---- pass B: per-cell sums ----
cells = []
for s, co, genes, M, tab in donors(S):
    t = tab[["soma_joinid", "cell_type", "subtype", "z_capillary", "z_arterial", "z_venous", "umi", "dataset_id"]].copy()
    t["sample"] = s; t["cohort"] = co
    t["EDN1"] = row(M, genes, "EDN1"); t["EDNRB"] = row(M, genes, "EDNRB"); t["flow"] = rows_sum(M, genes, FLOW)
    for k in SIG:
        t["up_" + k] = rows_sum(M, genes, SIG[k]["up"]); t["down_" + k] = rows_sum(M, genes, SIG[k]["down"])
    cells.append(t)
cells = pd.concat(cells, ignore_index=True)
other = np.where(cells.cohort == "healthy", "disease", "healthy")
cells["score_cf"] = np.where(other == "healthy", np.log2((cells.up_healthy + 1) / (cells.down_healthy + 1)), np.log2((cells.up_disease + 1) / (cells.down_disease + 1)))
cells["score_all"] = np.log2((cells.up_all + 1) / (cells.down_all + 1))
cells["flow_score"] = np.log1p(1e4 * cells.flow / cells.umi)
cells["art_pos"] = cells.z_arterial - cells.z_capillary
cells.to_pickle(os.path.join(OUT, "cells.pkl"))
print("cells:", len(cells), cells.groupby(["cohort", "subtype"]).size().to_dict())
