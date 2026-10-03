"""Main figures of the combined manuscript, drawn from the saved result files.
They display results; no new statistic is computed except simple per-donor summaries stated in the legends.
Usage: python3 figs_main.py [figure names]"""
import sys
from figlib import *
from figlib import _rect

K = D["combined"]; A1 = D["aim1"]; A2 = D["p2_primary"]; A4 = D["aim4"]; S3 = D["stage3_de"]; GX = D["stage3_gtex"]
EN = D["stage3_endothelin"]; EC = D["stage3_edn1_cells"]
SORDER = list(STRATA)
HI_, LO_ = "#C5283D", "#2166AC"
SHORT = {"NO_eNOS": "Nitric oxide /\neNOS", "Endothelin_ACE": "Endothelin /\nACE", "Barrier": "Barrier /\njunction", "FA_transport": "Fatty-acid\ntransport",
         "Angiogenic_tip": "Angiogenic\ntip", "X_escape": "X-escape genes"}


def dot(c, label, m="o", ms=4.5): return Line2D([], [], marker=m, color=c, linestyle="none", markersize=ms, label=label)


SEX_HANDLES = [dot(FEM, "Female"), dot(MAL, "Male")]
MEDIAN_H = Line2D([], [], color=INK, linewidth=1.3, label="Median")


def sex_violins(ax, groups, labels, width=0.36, s=4.5):
    """groups: list of (female values, male values) pairs."""
    rng = np.random.default_rng(1)
    for i, (f, m) in enumerate(groups):
        for v, off, col in ((f, -0.21, FEM), (m, 0.21, MAL)):
            v = np.asarray(v, float); v = v[np.isfinite(v)]
            if len(v) > 3:
                vp = ax.violinplot(v, positions=[i + off], widths=width, showextrema=False)
                for b in vp["bodies"]:
                    b.set_facecolor(col); b.set_alpha(0.22); b.set_edgecolor(col); b.set_linewidth(0.6)
            ax.scatter(i + off + rng.uniform(-0.07, 0.07, len(v)), v, s=s, color=col, alpha=0.8, linewidths=0, zorder=3)
            ax.hlines(np.median(v), i + off - 0.15, i + off + 0.15, color=INK, linewidth=1.3, zorder=4)
    ax.set_xticks(range(len(groups))); ax.set_xticklabels(labels); ax.set_xlim(-0.55, len(groups) - 0.45)
    ax.tick_params(axis="x", length=0)


# ============================ Figure 2: where the differences lie ============================
def fig_cm(name="Figure2_cardiomyocytes"):
    fig = new_fig(180, 178)
    ax = panel(fig, (0, 0, 60, 76), "A", "Which cells differ by sex?", "Autosomal genes that differ in both\ncohorts, same donors for every cell type", l=22, r=6, b=12)
    F = S3["fair"].sort_values("replicated_equal_donors"); y = np.arange(len(F))
    ax.barh(y, F.replicated_equal_donors, color=[CELL[c] for c in F.cell_type], height=0.66, zorder=2)
    ax.scatter(F.null_replicated_max, y, marker="|", s=170, color=INK, linewidths=1.6, zorder=3)
    for yi, (v, nmax) in enumerate(zip(F.replicated_equal_donors, F.null_replicated_max)):
        ax.text(max(v, nmax) + 1.3, yi, str(int(v)), va="center", fontsize=7, fontweight="bold" if v > nmax else "normal")
    ax.set_yticks(y); ax.set_yticklabels(F.cell_type); ax.set_xlim(0, 40); grid(ax, "x"); ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Replicated genes", labelpad=1.5)
    ax.text(39.5, 0.1, "black tick = largest count\nin 20 runs with shuffled\nsex labels", ha="right", va="bottom", fontsize=5.8, color=GREY, linespacing=1.1)

    def zz(box, letter, ct, ttl):
        d = S3["de"][ct]; d = d[d.chr_class == "autosome"].copy(); d["rep"] = (d.rep_h2d.astype(bool) | d.rep_d2h.astype(bool))
        ax = panel(fig, box, letter, ttl, f"{int(d.rep.sum())} of {len(d):,} autosomal genes replicated", l=9, r=2, b=12)
        ax.axhline(0, color="#555555", linewidth=0.5); ax.axvline(0, color="#555555", linewidth=0.5)
        n = d[~d.rep]; ax.scatter(n.z_h, n.z_d, s=1.2, color="#C9C9C9", alpha=0.55, linewidths=0, rasterized=True)
        rr = d[d.rep]; ax.scatter(rr.z_h, rr.z_d, s=13, c=np.where(rr.z_comb > 0, FEM, MAL), edgecolors="white", linewidths=0.4, zorder=3)
        ax.set_xlim(-9.5, 9.5); ax.set_ylim(-9.5, 9.5); ax.set_xticks([-8, -4, 0, 4, 8]); ax.set_yticks([-8, -4, 0, 4, 8])
        ax.set_xlabel("Sex effect, healthy hearts (z)", labelpad=1.5); ax.set_ylabel("Sex effect, diseased hearts (z)", labelpad=0)
        top = rr.reindex(rr.z_comb.abs().sort_values(ascending=False).index).head(9)
        repel(ax, top.z_h.values, top.z_d.values, top.gene.values, fs=5.8, fontstyle="italic")
    zz((60, 0, 60, 76), "B", "cardiomyocyte", "Cardiomyocytes")
    zz((120, 0, 60, 76), "C", "capillary", "Capillary endothelium")
    legend_at(fig, 120, 71.5, [dot(FEM, "Higher in females"), dot(MAL, "Higher in males"), dot("#C9C9C9", "Not replicated", ms=3)])

    G = GX["cm_genes"]; G = G[G.logFC_gtex.notna()].copy()
    ax = panel(fig, (0, 80, 98, 98), "D", "Tested in a third cohort: GTEx, 432 hearts",
               f"{int(G.same_direction.sum())} of {len(G)} cardiomyocyte genes in the same direction in bulk left ventricle; {int(G.validated.sum())} with p < 0.05", l=11, r=3, b=14)
    ax.axhline(0, color="#555555", linewidth=0.5); ax.axvline(0, color="#555555", linewidth=0.5); grid(ax)
    ax.fill_between([0, 7], 0, 1.15, color=FEM, alpha=0.05, linewidth=0); ax.fill_between([-11.5, 0], -2.7, 0, color=MAL, alpha=0.05, linewidth=0)
    col = np.where(G.z_comb > 0, FEM, MAL); val = G.validated.astype(bool).values
    ax.scatter(G.z_comb[val], G.logFC_gtex[val], s=20, c=col[val], edgecolors="white", linewidths=0.4, zorder=3)
    ax.scatter(G.z_comb[~val], G.logFC_gtex[~val], s=16, facecolors="white", edgecolors=col[~val], linewidths=0.9, zorder=3)
    ax.set_xlim(-11.5, 7); ax.set_ylim(-2.7, 1.15)
    ax.set_xlabel("Single-cell cardiomyocytes, both cohorts (z)", labelpad=1.5); ax.set_ylabel("GTEx bulk left ventricle, female minus male (log2)", labelpad=1)
    lab = pd.concat([G.sort_values("p_gtex").head(16), G[G.gene.isin(["FRMD5"])]]).drop_duplicates("gene")
    repel(ax, lab.z_comb.values, lab.logFC_gtex.values, lab.gene.values, fs=6, fontstyle="italic")
    legend_at(fig, 52, 172.5, [dot("#555555", "p < 0.05 in GTEx"), Line2D([], [], marker="o", markerfacecolor="white", markeredgecolor="#555555", linestyle="none", markersize=4.5, label="Not significant in GTEx")])

    top = G.reindex(G.z_comb.abs().sort_values(ascending=False).index).head(26).sort_values("z_comb", ascending=False)
    M = top.set_index("gene")[["logFC_healthy", "logFC_disease", "logFC_gtex"]]; M.columns = ["Healthy\nsingle-cell", "Diseased\nsingle-cell", "GTEx\nbulk"]
    ax = panel(fig, (98, 80, 82, 98), "E", "Top genes in three cohorts", "Female minus male (log2)", l=17, r=19, b=9)
    sm = heat(ax, M, vlim=1.6, fs=5.6); ax.set_yticklabels(M.index, fontstyle="italic", fontsize=6.3)
    cbar_at(fig, sm, (164.5, 100, 2.6, 45), "log2", vertical=True, ticks=[-1.5, -1, -0.5, 0, 0.5, 1, 1.5])
    finish(fig, name)


# ============================ Figure 3: X-Y dosage ============================
def fig_xy(name="Figure3_XY_dosage"):
    fig = new_fig(180, 196)
    order = ["KDM5C/KDM5D", "KDM6A/UTY", "USP9X/USP9Y", "ZFX/ZFY", "DDX3X/DDX3Y", "EIF1AX/EIF1AY", "RPS4X/RPS4Y1", "NLGN4X/NLGN4Y"]
    groups = ["capillary", "arterial", "pericyte", "smooth muscle", "fibroblast", "myeloid", "cardiomyocyte"]
    B = A4["by_group"]

    def hm(box, letter, meas, ttl, sub, r):
        d = B[B.measure == meas]
        M = d.pivot(index="pair", columns="group", values="est").loc[order, groups]
        St = d.assign(s=np.where(d.FDR < 0.05, "*", "")).pivot(index="pair", columns="group", values="s").loc[order, groups]
        ax = panel(fig, box, letter, ttl, sub, l=24, r=r, b=15)
        sm = heat(ax, M, vlim=2.0, stars=St, fs=5.9); ax.set_xticklabels(groups, rotation=32, ha="right", rotation_mode="anchor"); ax.set_yticklabels(order, fontstyle="italic")
        return sm
    hm((0, 0, 86, 66), "A", "X_plus_Y", "X + Y copies together: less in females", "Healthy hearts; female minus male (log2); * FDR < 0.05", 1)
    sm = hm((86, 0, 94, 66), "B", "X_only", "X copy alone: more in females", "Escape from X inactivation does not make up for the Y copy", 14)
    cbar_at(fig, sm, (169, 14, 2.4, 34), "log2", vertical=True, ticks=[-2, -1, 0, 1, 2])

    xy = K["xy"].copy(); xy["Cohort"] = xy.cohort.map({"healthy": "Healthy", "disease": "Diseased"})
    axs = panel_grid(fig, (0, 68, 116, 64), "C", "Capillary endothelium: replicated in diseased hearts", "Female minus male (log2), 95% CI", ncol=2, l=24, r=2, b=14, gap=4, t=14, sharey=False)
    for ax, meas, tt in zip(axs, ["X_plus_Y", "X_only"], ["X + Y copies together", "X copy alone"]):
        zero(ax); grid(ax, "x")
        for i, p in enumerate(order):
            y = len(order) - 1 - i
            for co, off in (("Healthy", 0.17), ("Diseased", -0.17)):
                rw = xy[(xy.pair == p) & (xy.measure == meas) & (xy.Cohort == co)]
                if len(rw): forest(ax, y + off, rw.est.values[0], rw.lo.values[0], rw.hi.values[0], color=COH[co], ms=3.6, lw=1.1)
        ax.set_yticks(range(len(order))); ax.set_ylim(-0.6, len(order) - 0.4)
        ax.set_title(tt, fontsize=7.2, fontweight="bold", pad=3); ax.tick_params(axis="y", length=0)
    axs[0].set_yticklabels(order[::-1], fontstyle="italic"); axs[1].set_yticklabels([])
    axs[0].set_xlim(-4.6, 0.6); axs[1].set_xlim(-1.4, 1.5)
    legend_at(fig, 70, 127, [dot(COH["Healthy"], "Healthy (45 donors)"), dot(COH["Diseased"], "Diseased (93 donors)")])

    ab = GX["age_bands"]
    xl = {"XplusY_KDM5C_KDM5D": "KDM5C + KDM5D", "XplusY_KDM6A_UTY": "KDM6A + UTY", "XplusY_USP9X_USP9Y": "USP9X + USP9Y", "XplusY_ZFX_ZFY": "ZFX + ZFY", "XplusY_DDX3X_DDX3Y": "DDX3X + DDX3Y", "X_escape": "X-escape genes"}
    ax = panel(fig, (116, 68, 64, 64), "D", "Tested in bulk tissue", "GTEx, 432 donors; ZFX + ZFY is the exception", l=25, r=3, b=14)
    zero(ax); grid(ax, "x")
    for i, (kk, lb) in enumerate(xl.items()):
        y = len(xl) - 1 - i
        for band, off in (("20-49", 0.17), ("50-79", -0.17)):
            rw = ab[(ab.outcome == kk) & (ab.contrast == f"female - male, age {band}")]
            forest(ax, y + off, rw.est.values[0], rw.lo.values[0], rw.hi.values[0], color=AGE[band + " years"], ms=3.6, lw=1.1)
    ax.set_yticks(range(len(xl))); ax.set_yticklabels(list(xl.values())[::-1], fontstyle="italic"); ax.get_yticklabels()[0].set_fontstyle("normal")
    ax.set_ylim(-0.6, len(xl) - 0.4); ax.set_xlabel("Female minus male (SD units)", labelpad=1.5); ax.tick_params(axis="y", length=0)
    legend_at(fig, 150, 127, [dot(AGE["20-49 years"], "20–49 years"), dot(AGE["50-79 years"], "50–79 years")])

    pl = A4["pair_long"]; pl = pl[(pl.sex.astype(str) == "male") & np.isfinite(pl.Y_share)]
    pairs = ["KDM5C/KDM5D", "KDM6A/UTY", "USP9X/USP9Y", "DDX3X/DDX3Y"]
    axs = panel_grid(fig, (0, 134, 180, 62), "E", "In males, endothelial cells draw less from the Y copy than muscle cells", "Y share of X + Y transcripts in each male donor; box = median and quartiles", ncol=4, l=11, r=2, b=14, gap=3, t=14)
    rng = np.random.default_rng(2)
    for ax, p in zip(axs, pairs):
        d = pl[pl.pair == p]
        for i, g in enumerate(groups):
            v = 100 * d.Y_share[d.group == g].values; c = CELL[g]
            ax.boxplot(v, positions=[i], widths=0.62, patch_artist=True, showfliers=False, medianprops=dict(color=INK, linewidth=1.2),
                       boxprops=dict(facecolor=c, alpha=0.35, edgecolor=c, linewidth=0.8), whiskerprops=dict(color=c, linewidth=0.8), capprops=dict(linewidth=0))
            ax.scatter(i + rng.uniform(-0.18, 0.18, len(v)), v, s=3.2, color=c, alpha=0.85, linewidths=0, zorder=3)
        ax.axvspan(-0.5, 1.5, color="#D55E00", alpha=0.06, linewidth=0)
        ax.set_xticks(range(len(groups))); ax.set_xticklabels(groups, rotation=40, ha="right", rotation_mode="anchor"); ax.set_xlim(-0.6, len(groups) - 0.4)
        ax.set_title(p, fontsize=7.2, fontweight="bold", fontstyle="italic", pad=3); grid(ax, "y"); ax.set_ylim(0, 100); ax.tick_params(axis="x", length=0)
    axs[0].set_ylabel("Y share (%)", labelpad=1); axs[0].text(0.5, 3, "endothelium", ha="center", va="bottom", fontsize=5.8, color="#B34E00")
    finish(fig, name)


# ============================ Figure 4: capillary programs ============================
def fig_capillary(name="Figure4_capillary_programs"):
    per = K["per"].copy(); per["S"] = per.stratum.map(stratum_name); comb = K["combined"].set_index("program")
    order = ["NO_eNOS", "Endothelin_ACE", "Barrier", "FA_transport", "Angiogenic_tip", "mural_per_capillary", "X_escape"]
    fig = new_fig(180, 168)
    ax = panel(fig, (0, 0, 90, 78), "A", "Capillary programs: 138 hearts pooled", "Circles: datasets and diseases (size = donors). Diamond: pooled, 95% CI", l=27, r=3, b=9)
    grid(ax, "x"); ax.axvspan(-0.8, 0.8, color="#7B3294", alpha=0.09, linewidth=0); zero(ax)
    rng = np.random.default_rng(3)
    for i, p in enumerate(order):
        y = len(order) - 1 - i; d = per[per.program == p]
        if p == "X_escape": ax.axhspan(y - 0.5, y + 0.5, color="#FFF4D6", linewidth=0, zorder=0)
        ax.scatter(d.est, y + rng.uniform(-0.24, 0.24, len(d)), s=6 + 2.2 * (d.n_w + d.n_m), c=d.S.map(STRATA), alpha=0.9, edgecolors="white", linewidths=0.5, zorder=3)
        c = comb.loc[p]; ax.hlines(y, c.lo, c.hi, color=INK, linewidth=1.6, zorder=4)
        ax.plot(c.est, y, "D", markersize=5.5, markerfacecolor="white", markeredgecolor=INK, markeredgewidth=1.3, zorder=5)
    ax.set_yticks(range(len(order))); ax.set_yticklabels([PROG[p] for p in order][::-1]); ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_xlabel("Female minus male (SD units)", labelpad=1.5); ax.set_xlim(-2.4, 2.6); ax.tick_params(axis="y", length=0)
    ax = panel(fig, (90, 0, 90, 78), "B", "Each dataset and disease", "Healthy 1 = harmonized cohort; Healthy 2 = Heart Cell Atlas", l=27, r=13, b=13)
    M = per.pivot(index="program", columns="S", values="est").loc[order, SORDER]; M.index = [PROG[p] for p in order]
    sm = heat(ax, M, vlim=2.2); ax.set_xticklabels(SORDER, rotation=30, ha="right", rotation_mode="anchor")
    cbar_at(fig, sm, (170, 17, 2.4, 36), "SD units", vertical=True)
    legend_at(fig, 90, 78.5, [dot(c, s, ms=5) for s, c in STRATA.items()] + [Patch(facecolor="#7B3294", alpha=0.15, label="Equivalence zone (±0.8 SD)")])
    progs = ["NO_eNOS", "Endothelin_ACE", "Barrier", "FA_transport", "Angiogenic_tip", "X_escape"]
    s1 = A1["scores"].assign(stratum=lambda d: d.dataset_id.astype(str), state="Healthy"); s2 = A2["scores"].assign(stratum=lambda d: d.stratum.astype(str), state="Diseased")
    don = pd.concat([s1[["stratum", "state", "sex"] + progs], s2[["stratum", "state", "sex"] + progs]]); don["sex"] = don.sex.astype(str)
    axs = panel_grid(fig, (0, 87, 146, 81), "C", "Each dot is one donor: females and males overlap", "Scores standardized within each dataset or disease", ncol=5, l=10, r=3, t=18, b=11, gap=2.5)
    axd = panel(fig, (146, 87, 34, 81), "D", "Positive control", "A true sex difference", l=8, r=2, t=18, b=11)
    for j, (ax, p) in enumerate(zip(axs + [axd], progs)):
        z = don.groupby("stratum")[p].transform(lambda v: (v - v.mean()) / v.std())
        g = [(z[(don.state == st) & (don.sex == "female")], z[(don.state == st) & (don.sex == "male")]) for st in ("Healthy", "Diseased")]
        ax.axhline(0, color="#BBBBBB", linewidth=0.6); grid(ax, "y"); sex_violins(ax, g, ["Healthy", "Diseased"])
        ax.set_title(SHORT[p], fontsize=7.2, fontweight="bold", color=INK if j < 5 else "#8A5A00", pad=3)
    axs[0].set_ylabel("Program score (SD units)", labelpad=1); axd.set_facecolor("#FFF9E8")
    legend_at(fig, 90, 162.5, SEX_HANDLES + [MEDIAN_H])
    finish(fig, name)


# ============================ Figure 5: endothelin-1 ============================
def umap_row(fig, U, y0, letters="ABC"):
    """Three views of the same endothelial-cell map."""
    rng = np.random.default_rng(0); U = U.iloc[rng.permutation(len(U))]
    SUB = {"capillary": "#F2A65A", "arterial": "#B2182B", "venous": "#5E3C99"}
    specs = [("Vessel type", f"{len(U):,} endothelial cells, {U['sample'].nunique()} donors", U.subtype.map(SUB).values, [dot(c, s.capitalize(), ms=4) for s, c in SUB.items()]),
             ("Female and male cells intermix", "Sex-chromosome genes left out of the map", np.where(U.sex.astype(str) == "female", FEM, MAL), SEX_HANDLES),
             ("EDN1-expressing cells", "Cells with at least one EDN1 transcript", None, [dot("#C5283D", "EDN1-positive", ms=4), dot("#D9D9D9", "EDN1-negative", ms=4)])]
    for j, (tt, sub, cols, hd) in enumerate(specs):
        ax = panel(fig, (j * 60, y0, 60, 62), letters[j], tt, sub, l=4, r=2, b=8)
        if cols is not None: ax.scatter(U.UMAP1, U.UMAP2, s=0.5, c=cols, linewidths=0, alpha=0.6, rasterized=True)
        else:
            pos = (U.EDN1_pos > 0).values
            ax.scatter(U.UMAP1[~pos], U.UMAP2[~pos], s=0.4, color="#D9D9D9", linewidths=0, alpha=0.6, rasterized=True)
            ax.scatter(U.UMAP1[pos], U.UMAP2[pos], s=1.3, color="#C5283D", linewidths=0, alpha=0.9, rasterized=True)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_xlim(*np.percentile(U.UMAP1, [0.3, 99.7]) + np.array([-0.8, 0.8])); ax.set_ylim(*np.percentile(U.UMAP2, [0.3, 99.7]) + np.array([-0.8, 0.8]))
        for s in ax.spines.values(): s.set_visible(False)
        for xy_ in ((0.16, 0.0), (0, 0.16)):
            ax.annotate("", xy=xy_, xytext=(0, 0), xycoords="axes fraction", arrowprops=dict(arrowstyle="-|>", color="#555555", lw=0.6, mutation_scale=5))
        ax.text(0.02, -0.01, "UMAP 1", transform=ax.transAxes, fontsize=5.5, va="top", color="#555555"); ax.text(-0.01, 0.02, "UMAP 2", transform=ax.transAxes, fontsize=5.5, rotation=90, ha="right", color="#555555")
        legend_at(fig, j * 60 + 33, y0 + 56.5, hd, columnspacing=0.9)


def schematic(fig, box, letter):
    """Working model, drawn as a cartoon (hypothesis, not data)."""
    x, y, w, h = box
    head(fig, x, y, letter, "Working model (a hypothesis to test)", "Similar capillary programs, but a larger share of endothelin-1-producing, low-flow-response cells in female capillaries")
    ax = fig.add_axes(_rect(fig, x + 2, y + 10, w - 4, h - 11)); ax.set_xlim(0, 176); ax.set_ylim(0, (h - 11) / (w - 4) * 176); ax.axis("off"); fig._boxes.append((letter, box, [ax]))
    T = ax.get_ylim()[1]
    cf = EN["cell_fractions"]; md = cf.assign(sex=cf.sex.astype(str)).groupby(["cohort", "sex"]).pct_EDN1_pos.median()
    for x0, sexname, col, posidx, l1, l2 in [(2, "Female", FEM, [1, 2, 4, 7], "Larger share of EDN1-expressing capillary cells", f"median {md[('healthy', 'female')]:.1f}% of cells (healthy), {md[('disease', 'female')]:.1f}% (diseased)"),
                                             (91, "Male", MAL, [1, 5], "Smaller share of EDN1-expressing capillary cells", f"median {md[('healthy', 'male')]:.1f}% of cells (healthy), {md[('disease', 'male')]:.1f}% (diseased)")]:
        ax.add_patch(FancyBboxPatch((x0, 1), 83, T - 2, boxstyle="round,pad=0,rounding_size=2.5", facecolor="white", edgecolor=col, linewidth=1.0))
        ax.text(x0 + 3, T - 3, sexname, fontsize=8, fontweight="bold", color=col, va="top")
        ax.text(x0 + 80, T - 3.4, "arteriole  →  capillary  →  venule", fontsize=5.8, color=GREY, va="top", ha="right")
        cw = 6.6; n = 11; yb = T - 31
        ax.add_patch(Rectangle((x0 + 5, yb + 4.2), n * cw - 0.7, 5.6, facecolor="#FBEDED", edgecolor="none"))
        for i in range(n):
            pos = i in posidx
            ax.add_patch(FancyBboxPatch((x0 + 5 + i * cw, yb + 9.8), cw - 0.7, 4.2, boxstyle="round,pad=0,rounding_size=1.1", facecolor="#E8743B" if pos else "#F7E3D6", edgecolor="#B5714A", linewidth=0.5))
            ax.add_patch(FancyBboxPatch((x0 + 5 + i * cw, yb), cw - 0.7, 4.2, boxstyle="round,pad=0,rounding_size=1.1", facecolor="#F7E3D6", edgecolor="#B5714A", linewidth=0.5))
            if pos:
                xm = x0 + 5 + i * cw + (cw - 0.7) / 2
                ax.annotate("", xy=(xm, yb + 17.6), xytext=(xm, yb + 14.2), arrowprops=dict(arrowstyle="-|>", color="#C0392B", lw=0.9, mutation_scale=6))
        ax.annotate("", xy=(x0 + 5 + n * cw - 2, yb + 7), xytext=(x0 + 6.5, yb + 7), arrowprops=dict(arrowstyle="-|>", color="#D98C8C", lw=0.8, mutation_scale=6))
        ax.text(x0 + 5 + n * cw / 2, yb + 7, "blood flow", fontsize=5.5, color="#B86A6A", ha="center", va="center", bbox=dict(facecolor="#FBEDED", edgecolor="none", pad=0.8))
        ax.add_patch(FancyBboxPatch((x0 + 9, yb + 18), 46, 3.4, boxstyle="round,pad=0,rounding_size=1.6", facecolor="#BFE3D5", edgecolor="#009E73", linewidth=0.6))
        ax.text(x0 + 32, yb + 19.7, "pericyte / smooth muscle (ET-A receptor)", fontsize=5.4, ha="center", va="center", color="#00664B")
        ax.text(x0 + 58, yb + 16.2, "endothelin-1", fontsize=5.8, color="#C0392B", fontweight="bold", va="center")
        ax.text(x0 + 5, yb - 3, l1, fontsize=6.8, fontweight="bold", color=INK, va="top"); ax.text(x0 + 5, yb - 7.2, l2, fontsize=6, color=GREY, va="top")
        ax.text(x0 + 5, yb - 10.8, "EDN1 cells (orange): low flow-response genes KLF2, KLF4, NOS3, THBD", fontsize=6, color=GREY, va="top")


def fig_edn1(name="Figure5_endothelin"):
    U = D.get("umap"); y0 = 64 if U is not None else 0
    fig = new_fig(180, 188 + y0)
    if U is not None: umap_row(fig, U, 0)
    let = iter("DEFGHIJ" if U is not None else "ABCDEFG")
    ps = EN["per_stratum"]; ps = ps[(ps.gene == "EDN1") & (ps.cell_type == "capillary")].set_index("group").loc[SORDER]
    pooled = EN["genes"]; pooled = pooled[(pooled.gene == "EDN1") & (pooled.cell_type == "capillary") & (pooled.cohort == "both")].iloc[0]
    ax = panel(fig, (0, y0, 62, 62), next(let), "EDN1 in capillary cells", "Positive in all six strata", l=17, r=3, b=10)
    zero(ax); grid(ax, "x")
    for i, (s, rw) in enumerate(ps.iterrows()):
        forest(ax, len(ps) - i, rw.est, rw.est - 1.96 * rw.se, rw.est + 1.96 * rw.se, color=STRATA[s], ms=2.5 + 0.09 * (rw.n_w + rw.n_m), lw=1.2)
    ax.axhspan(-0.55, 0.55, color="#FCE4EC", linewidth=0, zorder=0)
    ax.hlines(0, pooled.lo, pooled.hi, color=INK, linewidth=1.8, zorder=4); ax.plot(pooled.est, 0, "D", markersize=6, markerfacecolor=FEM, markeredgecolor=INK, markeredgewidth=1.1, zorder=5)
    ax.set_yticks(range(len(ps) + 1)); ax.set_yticklabels(["Pooled"] + SORDER[::-1]); ax.get_yticklabels()[0].set_fontweight("bold")
    ax.set_xlim(-2.2, 3.7); ax.set_ylim(-0.6, len(ps) + 0.5); ax.set_xlabel("Female minus male (log2), 95% CI", labelpad=1.5); ax.tick_params(axis="y", length=0)
    ax.text(3.6, 0, f"+{pooled.est:.2f}\np = {pooled.p:.3f}", ha="right", va="center", fontsize=6.2, fontweight="bold", color="#A01347", linespacing=1.05)

    cf = EN["cell_fractions"].copy(); cf["sex"] = cf.sex.astype(str); pc = EN["positive_cells"]; pb = pc[(pc.cohort == "both") & (pc.adjusted == "none")].iloc[0]
    ax = panel(fig, (62, y0, 52, 62), next(let), "A larger share of EDN1 cells", "Share of capillary cells with EDN1, per donor", l=10, r=2, b=10)
    g = [(cf.pct_EDN1_pos[(cf.cohort == c) & (cf.sex == "female")], cf.pct_EDN1_pos[(cf.cohort == c) & (cf.sex == "male")]) for c in ("healthy", "disease")]
    grid(ax, "y"); sex_violins(ax, g, ["Healthy", "Diseased"]); ax.set_ylabel("EDN1-positive capillary cells (%)", labelpad=1); ax.set_ylim(-0.4, 12.5)
    ax.text(0.97, 0.97, f"pooled +{pb.est:.2f} SD\np = {pb.p:.3f}", transform=ax.transAxes, ha="right", va="top", fontsize=6.2, fontweight="bold", color="#A01347", linespacing=1.05)
    legend_at(fig, 62 + 29, y0 + 56.5, SEX_HANDLES + [MEDIAN_H], columnspacing=1.0)

    dn = EC["donors"]; ft = EC["features"].set_index("feature")
    feats = [("capscore_diff", "Capillary identity", "capscore"), ("EDNRB", "Cells with EDNRB (ET-B)", "% of cells with EDNRB (percentage points)"), ("shear_KLF2_diff", "Flow-response program", "shear_KLF2"),
             ("axis_diff", "Venous-vs-arterial axis", "axis"), ("hypoxia_diff", "Hypoxia response", "hypoxia"), ("angiogenic_diff", "Angiogenic tip", "angiogenic"),
             ("inflammatory_diff", "Inflammatory activation", "inflammatory"), ("interferon_diff", "Interferon response", "interferon"), ("stress_diff", "Tissue stress", "stress")]
    ax = panel(fig, (114, y0, 66, 62), next(let), "A recognizable profile", "EDN1-positive vs depth-matched cells, 83 donors", l=29, r=3, b=10)
    zero(ax); grid(ax, "x")
    for i, (kk, lb, fk) in enumerate(feats):
        v = (dn.EDNRB_pct_pos - dn.EDNRB_pct_neg) if kk == "EDNRB" else dn[kk]; v = v[np.isfinite(v)]; n = len(v)
        dz = v.mean() / v.std(); se = np.sqrt(1 / n + dz ** 2 / (2 * n))
        col = (HI_ if dz > 0 else LO_) if ft.loc[fk, "FDR"] < 0.05 else "#9A9A9A"
        forest(ax, len(feats) - 1 - i, dz, dz - 1.96 * se, dz + 1.96 * se, color=col, ms=4.2, lw=1.3)
    ax.set_yticks(range(len(feats))); ax.set_yticklabels([f[1] for f in feats][::-1]); ax.set_ylim(-0.6, len(feats) + 0.75); ax.tick_params(axis="y", length=0)
    ax.set_xlim(-0.85, 0.65); ax.set_xlabel("Standardized paired difference", labelpad=1.5)
    ax.text(-0.83, len(feats) + 0.65, "lower in\nEDN1+ cells", ha="left", va="top", fontsize=5.8, color=LO_, linespacing=1.05); ax.text(0.63, len(feats) + 0.65, "higher in\nEDN1+ cells", ha="right", va="top", fontsize=5.8, color=HI_, linespacing=1.05)

    de = EC["de"].copy(); de["lp"] = -np.log10(de["P.Value"]); sig = (de["adj.P.Val"] < 0.05).values
    ax = panel(fig, (0, y0 + 64, 74, 64), next(let), "What EDN1-making cells express", f"{int(sig.sum())} of {len(de)} genes differ (FDR < 0.05); paired, 83 donors", l=9, r=2, b=10)
    grid(ax); ax.axvline(0, color="#555555", linewidth=0.5)
    ax.scatter(de.logFC[~sig], de.lp[~sig], s=3, color="#C9C9C9", linewidths=0, alpha=0.7)
    ax.scatter(de.logFC[sig], de.lp[sig], s=8, c=np.where(de.logFC[sig] > 0, HI_, LO_), linewidths=0.3, edgecolors="white", alpha=0.95, zorder=3)
    ax.set_xlim(-1.25, 1.6); ax.set_ylim(0, 24.5); ax.set_xlabel("EDN1-positive vs matched cells (log2)", labelpad=1.5); ax.set_ylabel("−log10 p", labelpad=1)
    want = ["SMAD6", "NEBL", "CRIM1", "ENG", "CXCL12", "MRTFB", "FUT8", "KDR", "NEURL1B", "AQP1", "ITGA6", "CLIC5", "CA4", "RAMP3", "NRP1", "LPL", "KLF2", "ADGRL2"]
    lb = de[de.gene.isin(want)].sort_values("lp", ascending=False)
    repel(ax, lb.logFC.values, lb.lp.values, lb.gene.values, fs=5.9, fontstyle="italic")

    gg = EN["genes"]; gg = gg[gg.cohort == "both"]
    genes = ["EDN1", "ECE1", "EDNRB", "EDNRA"]; cts = ["capillary", "arterial", "venous", "mural", "fibroblast", "cardiomyocyte"]
    M = gg.pivot(index="gene", columns="cell_type", values="est").reindex(index=genes, columns=cts)
    St = gg.assign(s=np.where(gg.p < 0.05, "*", "")).pivot(index="gene", columns="cell_type", values="s").reindex(index=genes, columns=cts).fillna("")
    ax = panel(fig, (74, y0 + 64, 52, 64), next(let), "The endothelin system", "Female minus male (log2); * p < 0.05", l=13, r=2, b=24)
    sm = heat(ax, M, vlim=0.7, stars=St, fmt="{:+.2f}", fs=5.6, na="–"); ax.set_xticklabels(cts, rotation=40, ha="right", rotation_mode="anchor")
    ax.set_yticklabels(["EDN1\npeptide", "ECE1\nenzyme", "EDNRB\nET-B", "EDNRA\nET-A"], fontsize=6.2, linespacing=1.0)
    cbar_at(fig, sm, (90, y0 + 120.5, 30, 2.2), "log2", ticks=[-0.6, 0, 0.6])

    adj = EN["adjusted"]; adj = adj[adj.cohort == "both"].set_index("adjusted_for"); ag = EN["age"]; ag = ag[ag.cohort == "both"].set_index("ages")
    gt = GX["endothelin"].set_index("gene").loc["EDN1"]; gse = abs(gt.logFC_gtex / gt.t_gtex)
    rows = [("No covariates", adj.loc["nothing extra"], INK), ("+ tissue stress", adj.loc["Stress"], "#555555"), ("+ ambient muscle RNA", adj.loc["CM_ambient"], "#555555"),
            ("+ hypoxia", adj.loc["Hypoxia"], "#555555"), ("+ flow response (KLF2)", adj.loc["Shear_KLF2"], "#555555"), ("+ all four technical", adj.loc["stress + ambient + hypoxia + inflammation"], "#555555"),
            ("Age < 50 years", ag.loc["< 50 years"], AGE["20-49 years"]), ("Age ≥ 50 years", ag.loc[">= 50 years"], AGE["50-79 years"])]
    ax = panel(fig, (126, y0 + 64, 54, 64), next(let), "How robust is it?", "Pooled estimates positive; not seen in bulk", l=27, r=2, b=10)
    zero(ax); grid(ax, "x"); n = len(rows) + 1
    for i, (lbx, rw, c) in enumerate(rows): forest(ax, n - 1 - i, rw.est, rw.lo, rw.hi, color=c, ms=3.8, lw=1.2)
    ax.axhspan(-0.5, 0.5, color="#EAF5F0", linewidth=0, zorder=0); forest(ax, 0, gt.logFC_gtex, gt.logFC_gtex - 1.96 * gse, gt.logFC_gtex + 1.96 * gse, color="#1B9E77", ms=4, lw=1.3, marker="s")
    ax.set_yticks(range(n)); ax.set_yticklabels(["GTEx bulk, 432 hearts"] + [r[0] for r in rows][::-1]); ax.get_yticklabels()[-1].set_fontweight("bold")
    ax.set_ylim(-0.6, n - 0.4); ax.set_xlim(-0.6, 1.75); ax.set_xlabel("Female − male (log2)", labelpad=1.5); ax.tick_params(axis="y", length=0)
    schematic(fig, (0, y0 + 130, 180, 58), next(let))
    finish(fig, name)


# ============================ Figure 6: GTEx by age ============================
def fig_gtex(name="Figure7_GTEx_age"):
    fig = new_fig(180, 88)
    m = GX["meta"]; ct = m.groupby(["AGE", "sex"], observed=True).size().unstack().loc[["20-29", "30-39", "40-49", "50-59", "60-69", "70-79"]]
    ax = panel(fig, (0, 0, 60, 88), "A", "GTEx adds the younger women", "Left ventricle, 138 females and 294 males", l=10, r=2, b=14)
    x = np.arange(len(ct)); ax.axvspan(-0.5, 2.5, color=AGE["20-49 years"], alpha=0.07, linewidth=0)
    for sx, off, c in (("female", -0.2, FEM), ("male", 0.2, MAL)):
        ax.bar(x + off, ct[sx], width=0.38, color=c, zorder=2)
        for xi, v in zip(x, ct[sx]): ax.text(xi + off, v + 1.5, str(int(v)), ha="center", va="bottom", fontsize=5.8)
    ax.set_xticks(x); ax.set_xticklabels([s.replace("-", "–") for s in ct.index]); ax.set_xlabel("Age at death (years)", labelpad=1.5); ax.set_ylabel("Donors", labelpad=1); grid(ax, "y"); ax.set_ylim(0, 122); ax.tick_params(axis="x", length=0)
    ax.text(1, 118, "20–49 years:\n41 females, 73 males", ha="center", va="top", fontsize=6, color="#A34700", linespacing=1.1)
    legend_at(fig, 33, 82.5, [Patch(facecolor=FEM, label="Female"), Patch(facecolor=MAL, label="Male")])

    def bands(ax, tab, labs):
        zero(ax); grid(ax, "x")
        for i, kk in enumerate(labs):
            y = len(labs) - 1 - i
            for band, off in (("20-49", 0.17), ("50-79", -0.17)):
                rw = tab[(tab.outcome == kk) & (tab.contrast == f"female - male, age {band}")].iloc[0]
                forest(ax, y + off, rw.est, rw.lo, rw.hi, color=AGE[band + " years"], ms=3.8, lw=1.2)
                if band == "20-49" and rw.p < 0.05: ax.text(rw.hi + 0.02, y + off, "p = %.3f" % rw.p, va="center", fontsize=5.6, color="#A34700")
        ax.set_yticks(range(len(labs))); ax.set_yticklabels(list(labs.values())[::-1]); ax.set_ylim(-0.6, len(labs) - 0.4); ax.tick_params(axis="y", length=0)
    gl = {"NO_eNOS": "Nitric oxide / eNOS", "Barrier": "Barrier / junction", "FA_transport": "Fatty-acid transport", "Prostacyclin": "Prostacyclin", "Endothelin_ACE": "Endothelin / ACE", "Angiogenic_tip": "Angiogenic tip", "cap_content": "Capillary content"}
    ax = panel(fig, (60, 0, 64, 88), "B", "Vessel programs before and after 50", "Three of six higher in women at 20–49 years (uncorrected)", l=24, r=3, b=14)
    bands(ax, GX["age_bands"], gl); ax.set_xlim(-0.45, 0.98); ax.set_xlabel("Female minus male (SD units)", labelpad=1.5)
    ax = panel(fig, (124, 0, 56, 88), "C", "Endothelin genes in bulk", "No significant difference at either age", l=14, r=3, b=14)
    ea = GX["endothelin_age"]; el = {"gene_EDN1": "EDN1", "gene_ECE1": "ECE1", "gene_EDNRB": "EDNRB", "gene_EDNRA": "EDNRA", "gene_ACE": "ACE"}
    bands(ax, ea, el); ax.set_yticklabels(list(el.values())[::-1], fontstyle="italic"); ax.set_xlabel("Female minus male (SD units)", labelpad=1.5)
    e2 = ea[ea.contrast != "difference (20-49 minus 50-79)"]; ax.set_xlim(e2.lo.min() - 0.05, e2.hi.max() + 0.05)
    legend_at(fig, 122, 82.5, [dot(AGE["20-49 years"], "20–49 years"), dot(AGE["50-79 years"], "50–79 years")])
    finish(fig, name)


FIGS = {"cm": fig_cm, "xy": fig_xy, "capillary": fig_capillary, "edn1": fig_edn1, "gtex": fig_gtex}
if __name__ == "__main__":
    for n in (sys.argv[1:] or FIGS):
        FIGS[n]()
