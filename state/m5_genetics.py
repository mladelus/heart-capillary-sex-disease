"""M5: are the genes of the EDN1-associated state mapped to vascular GWAS associations more often than the rest of the
capillary transcriptome? GWAS Catalog (file already downloaded by R script 16); logistic regression with log gene length
and mean capillary expression as covariates. Needs gene_positions.csv (gene_name, seq_name, gene_seq_start, gene_seq_end)."""
from common import *
import re
SIG = pickle.load(open(os.path.join(OUT, "signatures.pkl"), "rb"))["SIG"]["all"]
T = SIG["table"].copy(); T = T[~T.excluded]
T["state"] = T.gene.isin(SIG["up"] + SIG["down"]).astype(int); T["state_up"] = T.gene.isin(SIG["up"]).astype(int)
pos = pd.read_csv(os.path.join(P2, "gene_positions.csv")); pos = pos[pos.seq_name.astype(str).isin([str(i) for i in range(1, 23)])]
pos["len"] = pos.gene_seq_end - pos.gene_seq_start + 1; T = T.merge(pos.groupby("gene_name").len.max().rename("length"), left_on="gene", right_index=True, how="inner")
gw_file = glob.glob(os.path.join(P2, "gwas", "*.tsv"))[0]
G = pd.read_csv(gw_file, sep="\t", low_memory=False, usecols=["DISEASE/TRAIT", "MAPPED_TRAIT", "MAPPED_GENE", "PVALUE_MLOG", "SNPS"])
G = G[(G.PVALUE_MLOG >= -np.log10(5e-8)) & G.MAPPED_GENE.notna()]
G["trait"] = (G["DISEASE/TRAIT"].fillna("") + " | " + G["MAPPED_TRAIT"].fillna("")).str.lower()
GROUPS = {  # fixed before the analysis (plan M5)
    "Coronary artery disease or myocardial infarction": r"coronary artery disease|coronary heart disease|myocardial infarction|coronary atherosclerosis|ischemic heart|ischaemic heart",
    "Blood pressure or hypertension": r"blood pressure|hypertension|pulse pressure",
    "Migraine": r"migraine",
    "Arterial dissection or fibromuscular dysplasia": r"dissection|fibromuscular",
    "Heart failure": r"heart failure",
    "CONTROL any trait": r".",
    "CONTROL height": r"\bheight\b",
    "CONTROL educational attainment": r"educational attainment",
}


def genes_of(rx):
    m = G[G.trait.str.contains(rx, regex=True)]
    return set(g.strip() for s in m.MAPPED_GENE for g in re.split(r"[,;]| - ", str(s)) if g.strip())


def logit(y, X):
    """Logistic regression by iteratively reweighted least squares; returns coefficients and standard errors."""
    b = np.zeros(X.shape[1])
    for _ in range(100):
        eta = X @ b; mu = 1 / (1 + np.exp(-eta)); W = mu * (1 - mu) + 1e-12
        H = X.T @ (X * W[:, None]); step = np.linalg.solve(H, X.T @ (y - mu)); b = b + step
        if np.abs(step).max() < 1e-9: break
    return b, np.sqrt(np.diag(np.linalg.inv(H)))


rows = []
z = lambda v: (v - v.mean()) / v.std()
for name, rx in GROUPS.items():
    hit = T.gene.isin(genes_of(rx)).astype(float).values
    for setname in ("state", "state_up"):
        x = T[setname].values.astype(float)
        if hit.sum() < 5 or hit[x == 1].sum() == 0 or hit.sum() == len(hit):
            rows.append(dict(trait_group=name, gene_set=setname, set_genes=int(x.sum()), set_with_hit=int(hit[x == 1].sum()), pct_set=100 * hit[x == 1].mean(), pct_rest=100 * hit[x == 0].mean())); continue
        X = np.column_stack([np.ones(len(T)), x, z(np.log(T.length.values)), z(T.mean_log2cpm.values)])
        b, se = logit(hit, X); X0 = np.column_stack([np.ones(len(T)), x]); b0, se0 = logit(hit, X0)
        rows.append(dict(trait_group=name, gene_set=setname, set_genes=int(x.sum()), set_with_hit=int(hit[x == 1].sum()), pct_set=100 * hit[x == 1].mean(), pct_rest=100 * hit[x == 0].mean(),
                         OR_unadjusted=np.exp(b0[1]), p_unadjusted=2 * stats.norm.sf(abs(b0[1] / se0[1])), OR=np.exp(b[1]), OR_lo=np.exp(b[1] - 1.96 * se[1]), OR_hi=np.exp(b[1] + 1.96 * se[1]),
                         p=2 * stats.norm.sf(abs(b[1] / se[1])), OR_per_SD_log_length=np.exp(b[2]), OR_per_SD_expression=np.exp(b[3])))
res = pd.DataFrame(rows)
pd.set_option("display.width", 250, "display.max_columns", 30, "display.float_format", lambda v: "%.3g" % v)
print(f"genes: {len(T)} tested autosomal capillary genes with a position; state signature {int(T.state.sum())} (up {int(T.state_up.sum())})")
print(res.to_string(index=False))
res.to_csv(os.path.join(OUT, "tables", "S3M_state_genetics.csv"), index=False)
pickle.dump(dict(table=res, genes=T), open(os.path.join(OUT, "genetics.pkl"), "wb"))
