"""Figure 6: an endothelial state defined without EDN1 (plan section M). Reads the outputs of state/m1, m2 (and m5 if present)."""
import sys
from figlib import *
from figlib import _rect
from figs_main import dot, sex_violins, SEX_HANDLES, MEDIAN_H, SORDER, HI_, LO_
ST = os.environ.get("STATE_OUT", os.path.join(HERE, "..", "state", "out"))
R = pickle.load(open(os.path.join(ST, "state_results.pkl"), "rb"))
C = pd.read_pickle(os.path.join(ST, "cells.pkl"))
GEN = pickle.load(open(os.path.join(ST, "genetics.pkl"), "rb")) if os.path.exists(os.path.join(ST, "genetics.pkl")) else None


def short(x):
    x = str(x)
    for k, v in (("dilated", "DCM"), ("arrhythmogenic", "ACM"), ("myocardial infarction", "Infarction"), ("myocarditis", "Myocarditis")):
        if x.startswith(k): return v
    return "Healthy 1" if x == "Heart" else "Healthy 2"


def fig_state(name="Figure6_endothelial_state"):
    h_gen = 62 if GEN is not None else 0
    fig = new_fig(180, 199 + h_gen)
    U = D.get("umap")
    # A: map coloured by the state score
    ax = panel(fig, (0, 0, 62, 66), "A", "A state scored without EDN1", "Signature learned in the other cohort", l=4, r=2, b=11)
    if U is not None:
        key = C.cohort + "_" + C["sample"] + "_" + C.soma_joinid; mu0 = C[C.subtype == "capillary"].groupby("dataset_id").score_cf.agg(["mean", "std"])
        sc = pd.Series(((C.score_cf - C.dataset_id.map(mu0["mean"])) / C.dataset_id.map(mu0["std"])).values, index=key.values)
        V = U.assign(score=U.cell.map(sc)).dropna(subset=["score"]).sort_values("score")
        lo, hi = np.percentile(V.score, [2, 98])
        scm = ax.scatter(V.UMAP1, V.UMAP2, s=0.5, c=V.score, cmap=LinearSegmentedColormap.from_list("st", ["#E6E6E6", "#F2C4A0", "#E0672F", "#8B1A1A"]), vmin=lo, vmax=hi, linewidths=0, rasterized=True)
        ax.set_xlim(*np.percentile(U.UMAP1, [0.3, 99.7]) + np.array([-0.8, 0.8])); ax.set_ylim(*np.percentile(U.UMAP2, [0.3, 99.7]) + np.array([-0.8, 0.8]))
        cax = fig.add_axes(_rect(fig, 18, 60.5, 28, 1.8)); cb = fig.colorbar(scm, cax=cax, orientation="horizontal"); cb.outline.set_visible(False); cb.ax.tick_params(length=1.5, labelsize=5.6, pad=1); cb.set_label("State score (z within dataset)", fontsize=6.2, labelpad=1)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
    # B: validation AUC
    auc = R["auc"]
    ax = panel(fig, (62, 0, 50, 66), "B", "It finds EDN1 cells", "Held-out donors; 0.5 = chance", l=10, r=3, b=11)
    rng = np.random.default_rng(4); grid(ax, "y"); ax.axhline(0.5, color="#555555", linewidth=0.7, linestyle=(0, (3, 2)))
    for i, co in enumerate(("healthy", "disease")):
        v = auc.auc_within_depth[auc.cohort == co].dropna().values
        vp = ax.violinplot(v, positions=[i], widths=0.7, showextrema=False)
        for b in vp["bodies"]: b.set_facecolor(list(COH.values())[i]); b.set_alpha(0.2); b.set_edgecolor(list(COH.values())[i])
        ax.scatter(i + rng.uniform(-0.12, 0.12, len(v)), v, s=6, color=list(COH.values())[i], alpha=0.85, linewidths=0, zorder=3)
        ax.hlines(np.median(v), i - 0.25, i + 0.25, color=INK, linewidth=1.4, zorder=4)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Healthy", "Diseased"]); ax.set_xlim(-0.6, 1.6); ax.set_ylim(0.3, 1.0); ax.tick_params(axis="x", length=0)
    ax.set_ylabel("AUC, EDN1-positive vs negative cells", labelpad=1)
    ax.text(0.5, 0.03, f"mean AUC {auc.auc_within_depth.mean():.2f}; 83 donors,\ncompared within depth deciles", transform=ax.transAxes, ha="center", va="bottom", fontsize=6, color=GREY, linespacing=1.1)
    # C: continuum (scores standardized within dataset, because each cohort is scored with the other cohort's signature)
    ca = C[C.subtype.isin(["capillary", "arterial"])].copy()
    mu = C[C.subtype == "capillary"].groupby("dataset_id").score_cf.agg(["mean", "std"])
    ca["z"] = (ca.score_cf - ca.dataset_id.map(mu["mean"])) / ca.dataset_id.map(mu["std"])
    ca["bin"] = ca.groupby("dataset_id").art_pos.transform(lambda v: np.ceil(pd.Series(v).rank(method="first") * 10 / len(v))).astype(int)
    g = ca.groupby(["cohort", "sample", "bin"]); B = pd.DataFrame(dict(n=g.size(), score=g.z.mean(), flow=g.flow_score.mean(), pct_pos=100 * (g.EDN1.apply(lambda v: (v > 0).mean())))).reset_index(); B = B[B.n >= 5]
    axs = panel_grid(fig, (112, 0, 68, 66), "C", "Highest at the arterial end", "Cells ordered from capillary-like to arterial-like", ncol=1, nrow=3, l=13, r=2, t=10.5, b=14, vgap=2.5, sharey=False)
    for ax, (col, lab) in zip(axs, (("score", "State\nscore (z)"), ("pct_pos", "EDN1+\ncells (%)"), ("flow", "Flow\nresponse"))):
        grid(ax, "y")
        for co, c in (("healthy", COH["Healthy"]), ("disease", COH["Diseased"])):
            gg = B[B.cohort == co].groupby("bin")[col]; mm, se = gg.mean(), gg.sem()
            ax.fill_between(mm.index, mm - se, mm + se, color=c, alpha=0.18, linewidth=0); ax.plot(mm.index, mm, "-o", color=c, markersize=2.6, linewidth=1.1)
        ax.set_ylabel(lab, fontsize=6.3, labelpad=1, linespacing=1.0); ax.set_xlim(0.6, 10.4); ax.set_xticks(range(1, 11)); ax.tick_params(labelsize=6)
        if ax is not axs[-1]: ax.tick_params(labelbottom=False)
        ax.axvspan(7.5, 10.4, color="#B2182B", alpha=0.05, linewidth=0)
    axs[-1].set_xlabel("Arterial position (decile; 1 = capillary, 10 = arterial)", fontsize=6.6, labelpad=1.5)
    legend_at(fig, 150, 60.5, [dot(COH["Healthy"], "Healthy"), dot(COH["Diseased"], "Diseased")], columnspacing=0.9)
    # D: abundance per donor
    D0 = R["donors"]; D0 = D0.assign(sex=D0.sex.astype(str)); pr = R["primary"]["both"]
    ax = panel(fig, (0, 68, 56, 64), "D", "State-high capillary cells", "Share per donor; top 10% of each dataset", l=10, r=2, b=11)
    g = [(D0.pct_hi[(D0.cohort == c) & (D0.sex == "female")], D0.pct_hi[(D0.cohort == c) & (D0.sex == "male")]) for c in ("healthy", "disease")]
    grid(ax, "y"); sex_violins(ax, g, ["Healthy", "Diseased"]); ax.set_ylabel("State-high capillary cells (%)", labelpad=1); ax.set_ylim(-1, min(60, D0.pct_hi.max() * 1.18))
    ax.text(0.97, 0.97, f"pooled {pr['est']:+.2f} SD\np = {pr['p']:.3f}", transform=ax.transAxes, ha="right", va="top", fontsize=6.2, fontweight="bold", color="#555555", linespacing=1.05)
    legend_at(fig, 31, 126.5, SEX_HANDLES + [MEDIAN_H], columnspacing=1.0)
    # E: by stratum
    per = R["primary"]["per"].assign(S=lambda d: d.stratum.map(short)).set_index("S").loc[SORDER]
    ax = panel(fig, (56, 68, 56, 64), "E", "Primary test: abundance", "Female minus male; cross-fitted signature", l=17, r=3, b=11)
    zero(ax); grid(ax, "x")
    for i, (s, rw) in enumerate(per.iterrows()): forest(ax, len(per) - i, rw.est, rw.est - 1.96 * rw.se, rw.est + 1.96 * rw.se, color=STRATA[s], ms=2.5 + 0.09 * (rw.n_w + rw.n_m), lw=1.2)
    ax.axhspan(-0.55, 0.55, color="#F1F1F1", linewidth=0, zorder=0); ax.hlines(0, pr["lo"], pr["hi"], color=INK, linewidth=1.8, zorder=4)
    ax.plot(pr["est"], 0, "D", markersize=6, markerfacecolor="white", markeredgecolor=INK, markeredgewidth=1.2, zorder=5)
    ax.set_yticks(range(len(per) + 1)); ax.set_yticklabels(["Pooled"] + SORDER[::-1]); ax.get_yticklabels()[0].set_fontweight("bold"); ax.set_ylim(-0.6, len(per) + 0.5)
    ax.set_xlim(-2.2, 3.3); ax.set_xlabel("SD units, 95% CI", labelpad=1.5); ax.tick_params(axis="y", length=0)
    # F: robustness and capillary definitions
    T = R["table"].set_index("analysis")
    rows = [("Primary (cross-fitted)", "PRIMARY abundance (cross-fitted signature, 90th pct)", INK), ("+ sequencing depth", "abundance + median log UMI", "#555555"),
            ("EDN1-negative cells only", "abundance among EDN1-negative capillary cells only", "#555555"), ("Strict capillary cells", "(ii) strict capillary (margin >= 1)", CELL["capillary"]),
            ("Atlas-labelled capillary*", "(iii) atlas-concordant capillary", CELL["capillary"]), ("Capillary-like third", "(iv) capillary, capillary-like third", CELL["capillary"]),
            ("Middle third", "(iv) capillary, middle third", CELL["capillary"]), ("Arterial-like third", "(iv) capillary, arterial-like third", CELL["capillary"]),
            ("Arterial cells", "(v) arterial-assigned cells", "#B2182B"), ("Venous cells", "(v) venous-assigned cells", "#5E3C99"),
            ("Intensity of the state", "intensity (mean score of state-high cells)", "#1B9E77"), ("Signature from all donors†", "abundance, all-donor signature (not cross-fitted)", "#9A9A9A")]
    ax = panel(fig, (112, 68, 68, 64), "F", "The same test, other definitions", "Abundance of state-high cells, female minus male", l=30, r=3, b=11)
    zero(ax); grid(ax, "x"); n = len(rows)
    for i, (lab, key, c) in enumerate(rows):
        rw = T.loc[key]; forest(ax, n - 1 - i, rw.est, rw.lo, rw.hi, color=c, ms=3.8, lw=1.2, marker="s" if "Intensity" in lab else "o")
    ax.set_yticks(range(n)); ax.set_yticklabels([r[0] for r in rows][::-1]); ax.get_yticklabels()[-1].set_fontweight("bold"); ax.set_ylim(-0.6, n - 0.4)
    ax.set_xlabel("SD units, 95% CI", labelpad=1.5); ax.tick_params(axis="y", length=0); ax.axhline(1.5, color="#BBBBBB", linewidth=0.5, linestyle=(0, (2, 2)))
    # G: signature genes (cross-cohort agreement)
    SIG = pickle.load(open(os.path.join(ST, "signatures.pkl"), "rb"))["SIG"]
    th, td = SIG["healthy"]["table"].set_index("gene"), SIG["disease"]["table"].set_index("gene"); j = th.join(td, lsuffix="_h", rsuffix="_d", how="inner"); j = j[~j.excluded_h]
    ax = panel(fig, (0, 134, 66, 62), "G", "The same genes in both cohorts", f"EDN1-positive vs matched cells; {len(j):,} genes", l=10, r=2, b=10)
    ax.axhline(0, color="#555555", linewidth=0.5); ax.axvline(0, color="#555555", linewidth=0.5)
    both = (j.FDR_h < 0.05) & (j.FDR_d < 0.05) & (np.sign(j.t_h) == np.sign(j.t_d))
    ax.scatter(j.t_h[~both], j.t_d[~both], s=2, color="#C9C9C9", linewidths=0, alpha=0.7, rasterized=True)
    ax.scatter(j.t_h[both], j.t_d[both], s=11, c=np.where(j.t_h[both] > 0, HI_, LO_), edgecolors="white", linewidths=0.3, zorder=3)
    rho = stats_spearman(j.t_h, j.t_d); ax.text(0.03, 0.97, f"Spearman rho = {rho:.2f}", transform=ax.transAxes, va="top", fontsize=6.2, color=GREY)
    lim = max(abs(j.t_h).max(), abs(j.t_d).max()) * 1.05; ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("Healthy hearts (t)", labelpad=1.5); ax.set_ylabel("Diseased hearts (t)", labelpad=1)
    top = j[both].assign(a=lambda d: d.t_h.abs() + d.t_d.abs()).sort_values("a", ascending=False).head(14)
    repel(ax, top.t_h.values, top.t_d.values, top.index.values, fs=5.8, fontstyle="italic")
    # H: what the signature contains (top genes heat strip)
    ta = SIG["all"]["table"]; ta = ta[~ta.excluded & (ta.FDR < 0.05)]
    up = ta[ta.mean_log2_diff > 0].sort_values("t", ascending=False).head(16); dn = ta[ta.mean_log2_diff < 0].sort_values("t").head(16)
    ax = panel(fig, (66, 134, 114, 62), "H", "What defines the state", "Largest differences, EDN1-positive vs matched cells (all donors; EDN1 and marker genes removed)", l=11, r=2, b=18)
    allg = pd.concat([up, dn.iloc[::-1]]); x = np.arange(len(allg))
    ax.bar(x, allg.mean_log2_diff, color=np.where(allg.mean_log2_diff > 0, HI_, LO_), width=0.72, zorder=2); ax.axhline(0, color="#333333", linewidth=0.6)
    ax.set_xticks(x); ax.set_xticklabels(allg.gene, rotation=60, ha="right", rotation_mode="anchor", fontsize=6, fontstyle="italic"); ax.set_xlim(-0.7, len(allg) - 0.3)
    ax.set_ylabel("log2 difference", labelpad=1); grid(ax, "y"); ax.tick_params(axis="x", length=0)
    ax.text(0.01, 0.04, "higher in the state", transform=ax.transAxes, fontsize=6, color=HI_); ax.text(0.99, 0.94, "lower in the state", transform=ax.transAxes, fontsize=6, color=LO_, ha="right")
    fig.text(0.995, 0.003 + (h_gen / fig._H if GEN is not None else 0), "* healthy cohort only (the disease atlases do not label capillaries).  † not cross-fitted: learned and scored in the same donors.", fontsize=5.6, color=GREY, ha="right", va="bottom")
    if GEN is not None:
        G = GEN["table"]; G = G[G.gene_set == "state"].set_index("trait_group")
        order = ["Coronary artery disease or myocardial infarction", "Blood pressure or hypertension", "Migraine", "Arterial dissection or fibromuscular dysplasia", "Heart failure", "CONTROL any trait", "CONTROL height", "CONTROL educational attainment"]
        lab = ["Coronary disease / infarction", "Blood pressure", "Migraine", "Dissection / fibromuscular dysplasia", "Heart failure", "Any trait (control)", "Height (control)", "Education (control)"]
        ax = panel(fig, (0, 203, 92, 58), "I", "Genetics matched to the state", "State genes vs the other capillary genes; adjusted for gene length and expression", l=42, r=12, b=10)
        ax.axvline(1, color="#333333", linewidth=0.7); grid(ax, "x")
        for i, (k, lb) in enumerate(zip(order, lab)):
            if k not in G.index or not np.isfinite(G.loc[k].get("OR", np.nan)): continue
            rw = G.loc[k]; c = "#7B3294" if i < 5 else "#9A9A9A"; y = len(order) - 1 - i
            ax.hlines(y, rw.OR_lo, rw.OR_hi, color=c, linewidth=1.3); ax.plot(rw.OR, y, "o", color=c, markersize=4.5)
            ax.text(1.02, y, f"p = {rw.p:.2f}", transform=ax.get_yaxis_transform(), va="center", fontsize=5.8, color=GREY)
        ax.set_xscale("log"); ax.set_yticks(range(len(order))); ax.set_yticklabels(lab[::-1]); ax.set_ylim(-0.6, len(order) - 0.4); ax.tick_params(axis="y", length=0)
        ax.set_xlabel("Odds ratio (95% CI) of being mapped to an association", labelpad=1.5)
        from matplotlib.ticker import ScalarFormatter, NullFormatter
        ax.set_xticks([0.5, 1, 2, 4]); ax.xaxis.set_major_formatter(ScalarFormatter()); ax.xaxis.set_minor_formatter(NullFormatter())
        ax = panel(fig, (92, 203, 88, 58), "J", "Share of genes with an association", "State genes and the other capillary genes", l=36, r=4, b=14)
        y = np.arange(len(order))[::-1]; grid(ax, "x")
        for yi, k in zip(y, order):
            if k not in G.index: continue
            rw = G.loc[k]; ax.hlines(yi, rw.pct_rest, rw.pct_set, color="#BBBBBB", linewidth=1.2); ax.plot(rw.pct_rest, yi, "o", color="#8C8C8C", markersize=4.5); ax.plot(rw.pct_set, yi, "o", color="#7B3294", markersize=5)
        ax.set_yticks(range(len(order))); ax.set_yticklabels(lab[::-1]); ax.set_ylim(-0.6, len(order) - 0.4); ax.tick_params(axis="y", length=0); ax.set_xlabel("Genes mapped to an association (%)", labelpad=1.5)
        legend_at(fig, 150, 255.5, [dot("#7B3294", "State genes"), dot("#8C8C8C", "Other capillary genes")], columnspacing=0.9)
    finish(fig, name)


def stats_spearman(a, b):
    ra, rb = pd.Series(a).rank(), pd.Series(b).rank(); return float(np.corrcoef(ra, rb)[0, 1])


if __name__ == "__main__":
    fig_state()
