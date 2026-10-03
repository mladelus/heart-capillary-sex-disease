"""Figure 1 (study design) and the supplementary figures. Usage: python3 figs_extra.py [names]"""
import sys
from figlib import *
from figlib import _rect
from figs_main import dot, sex_violins, SEX_HANDLES, MEDIAN_H, SORDER, HI_, LO_, K, A1, A2, A4, S3, GX, EN, EC

U = D.get("umap")
if U is not None: U["S"] = U.stratum.map(stratum_name)


def umap_axes(ax, U):
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlim(*np.percentile(U.UMAP1, [0.3, 99.7]) + np.array([-0.8, 0.8])); ax.set_ylim(*np.percentile(U.UMAP2, [0.3, 99.7]) + np.array([-0.8, 0.8]))
    for s in ax.spines.values(): s.set_visible(False)


def box(ax, x, y, w, h, text, fc="#F4F4F4", ec="#888888", fs=6.4, bold_first=True, tc=INK, lw=0.8):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.6", facecolor=fc, edgecolor=ec, linewidth=lw))
    lines = text.split("\n")
    n = len(lines); lh = fs * 0.44
    for i, ln in enumerate(lines):
        ax.text(x + w / 2, y + h / 2 + (n - 1) / 2 * lh - i * lh, ln, ha="center", va="center", fontsize=fs + (0.6 if (i == 0 and bold_first) else 0), fontweight="bold" if (i == 0 and bold_first) else "normal", color=tc)


def arrow(ax, x1, y1, x2, y2, c="#777777"):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="-|>", color=c, lw=0.9, mutation_scale=7, shrinkA=0, shrinkB=0))


# ============================ Figure 1 ============================
def fig_design(name="Figure1_study_design"):
    fig = new_fig(180, 142)
    head(fig, 0, 0, "A", "Three open cohorts, one question: where do female and male hearts differ?", "All data are open access; analysis plans and code are public")
    ax = fig.add_axes(_rect(fig, 1, 10, 178, 62)); ax.set_xlim(0, 178); ax.set_ylim(0, 62); ax.axis("off"); fig._boxes.append(("A", (0, 0, 180, 73), [ax]))
    box(ax, 1, 43, 56, 18, "Healthy hearts, single-cell\n45 donors (23 F, 22 M), 20–75 years\n2 datasets · 516,224 ventricular cells\nCELLxGENE Census", fc="#EFEFEF", ec="#4D4D4D")
    box(ax, 61, 43, 56, 18, "Diseased hearts, single-cell\n93 donors (22 F, 71 M), 19–75 years\nDCM · ACM · infarction · myocarditis\n595,519 ventricular cells", fc="#FDEBDD", ec="#E66101")
    box(ax, 121, 43, 56, 18, "Bulk left ventricle (GTEx v8)\n432 donors (138 F, 294 M), 20–79 years\nincludes 41 women aged 20–49\nindependent validation", fc="#E3F2EC", ec="#1B9E77")
    box(ax, 21, 25, 96, 11, "One pipeline for 138 donors\none pseudobulk per donor and cell type · sex effect within each dataset or disease · random-effects meta-analysis", fc="white", ec="#555555", fs=6.0)
    arrow(ax, 29, 43, 45, 36.3); arrow(ax, 89, 43, 89, 36.3)
    qs = [("1  Which cells differ?", "discovery and replication\nacross cohorts, 8 cell types", CELL["cardiomyocyte"], "Fig. 2"),
          ("2  Does X–Y dose differ?", "X–Y gene pairs\nin 7 cell groups", "#7B3294", "Fig. 3"),
          ("3  Do capillaries differ?", "6 programs, equivalence\nbounds ±0.8 SD", CELL["capillary"], "Fig. 4"),
          ("4  Is endothelin-1 special?", "cells, a state without\nEDN1, genetics", "#C5283D", "Figs. 5, 6"),
          ("5  Does age matter?", "women before and\nafter 50 (GTEx)", "#5E3C99", "Fig. 7")]
    w = 33.2
    for i, (t, s, c, f) in enumerate(qs):
        x0 = 1 + i * (w + 2.5)
        ax.add_patch(FancyBboxPatch((x0, 1), w, 17, boxstyle="round,pad=0,rounding_size=1.6", facecolor="white", edgecolor=c, linewidth=1.1))
        ax.add_patch(FancyBboxPatch((x0, 13), w, 5, boxstyle="round,pad=0,rounding_size=1.6", facecolor=c, edgecolor=c, linewidth=1.1))
        ax.add_patch(Rectangle((x0, 13), w, 2, facecolor=c, edgecolor="none"))
        ax.text(x0 + w / 2, 15.5, t, ha="center", va="center", fontsize=6.5, fontweight="bold", color="white")
        ax.text(x0 + w / 2, 8.6, s, ha="center", va="center", fontsize=6, color=INK, linespacing=1.15)
        ax.text(x0 + w / 2, 3.2, f, ha="center", va="center", fontsize=5.8, color=c, fontweight="bold")
    for i in range(4): arrow(ax, 69, 25, 1 + i * (w + 2.5) + w / 2, 18.3)
    arrow(ax, 149, 43, 1 + 4 * (w + 2.5) + w / 2, 18.3, c="#1B9E77")
    ax.text(126, 30.5, "also tests 1, 2 and 4\nin bulk tissue", fontsize=5.8, color="#147A5B", ha="left", va="center", style="italic", linespacing=1.1)

    # B donors
    sc = pd.concat([A1["scores"].assign(S=lambda d: d.dataset_id.astype(str).map(stratum_name))[["S", "sex"]], A2["scores"].assign(S=lambda d: d.stratum.astype(str).map(stratum_name))[["S", "sex"]]])
    ct = sc.assign(sex=sc.sex.astype(str)).groupby(["S", "sex"]).size().unstack().loc[SORDER[::-1]]
    ax = panel(fig, (0, 76, 70, 66), "B", "138 donors in six strata", "Sex effects are estimated within each stratum", l=17, r=4, b=13)
    y = np.arange(len(ct))
    ax.barh(y + 0.19, ct.female, height=0.36, color=FEM, zorder=2); ax.barh(y - 0.19, ct.male, height=0.36, color=MAL, zorder=2)
    for yi, (f, m) in enumerate(zip(ct.female, ct.male)):
        ax.text(f + 0.7, yi + 0.19, str(int(f)), va="center", fontsize=6); ax.text(m + 0.7, yi - 0.19, str(int(m)), va="center", fontsize=6)
    ax.set_yticks(y); ax.set_yticklabels(ct.index); ax.set_xlim(0, 39); ax.set_xlabel("Donors", labelpad=1.5); grid(ax, "x"); ax.tick_params(axis="y", length=0)
    for lab, col in zip(ax.get_yticklabels(), [STRATA[s] for s in ct.index]): lab.set_color(col); lab.set_fontweight("bold")
    legend_at(fig, 40, 136.5, [Patch(facecolor=FEM, label="Female"), Patch(facecolor=MAL, label="Male")])
    # C / D UMAP
    if U is not None:
        rng = np.random.default_rng(0); V = U.iloc[rng.permutation(len(U))]
        ax = panel(fig, (70, 76, 55, 66), "C", "Endothelial cells of all donors", f"{len(V):,} cells, up to 300 per donor", l=3, r=1, b=13)
        ax.scatter(V.UMAP1, V.UMAP2, s=0.5, c=V.S.map(STRATA).values, linewidths=0, alpha=0.6, rasterized=True); umap_axes(ax, V)
        legend_at(fig, 97.5, 130.5, [dot(c, s, ms=3.8) for s, c in STRATA.items()], ncol=3, columnspacing=0.8, handletextpad=0.1)
        SUB = {"capillary": "#F2A65A", "arterial": "#B2182B", "venous": "#5E3C99"}
        ax = panel(fig, (125, 76, 55, 66), "D", "One rule assigns vessel type", "Capillary, arterial or venous", l=3, r=1, b=13)
        ax.scatter(V.UMAP1, V.UMAP2, s=0.5, c=V.subtype.map(SUB).values, linewidths=0, alpha=0.6, rasterized=True); umap_axes(ax, V)
        legend_at(fig, 152.5, 132.5, [dot(c, s.capitalize(), ms=3.8) for s, c in SUB.items()], columnspacing=0.8, handletextpad=0.1)
    finish(fig, name)


# ============================ S1: UMAP detail ============================
def fig_s_umap(name="FigureS1_endothelial_map"):
    fig = new_fig(180, 150)
    rng = np.random.default_rng(0); V = U.iloc[rng.permutation(len(U))]
    ax = panel(fig, (0, 0, 62, 70), "A", "Healthy and diseased cells", "Largely overlapping after alignment", l=3, r=1, b=9)
    ax.scatter(V.UMAP1, V.UMAP2, s=0.5, c=np.where(V.cohort == "healthy", COH["Healthy"], COH["Diseased"]), linewidths=0, alpha=0.55, rasterized=True); umap_axes(ax, V)
    legend_at(fig, 31, 64, [dot(COH["Healthy"], "Healthy"), dot(COH["Diseased"], "Diseased")])
    cl = U.groupby("cluster").agg(cells=("cell", "size"), cap=("subtype", lambda s: 100 * (s == "capillary").mean()), art=("subtype", lambda s: 100 * (s == "arterial").mean()),
                                  ven=("subtype", lambda s: 100 * (s == "venous").mean()), edn=("EDN1_pos", lambda s: 100 * s.mean()), fem=("sex", lambda s: 100 * (s == "female").mean())).sort_values("edn")
    cl = cl[cl.cells >= 300]
    axs = panel_grid(fig, (62, 0, 118, 70), "B", "EDN1-expressing cells are concentrated in a few clusters", "Unsupervised clusters with at least 300 cells; descriptive, cells are not independent observations", ncol=2, l=17, r=4, b=15, gap=5, sharey=False)
    y = np.arange(len(cl)); a = axs[0]
    a.barh(y, cl.cap, color="#F2A65A", height=0.7); a.barh(y, cl.art, left=cl.cap, color="#B2182B", height=0.7); a.barh(y, cl.ven, left=cl.cap + cl.art, color="#5E3C99", height=0.7)
    a.set_yticks(y); a.set_yticklabels([f"cluster {c}" for c in cl.index]); a.set_xlim(0, 100); a.set_xlabel("Vessel type of cells (%)", labelpad=1.5); a.tick_params(axis="y", length=0)
    b = axs[1]; b.barh(y, cl.edn, color="#C5283D", height=0.7, zorder=2); b.set_yticks(y); b.set_yticklabels([]); b.set_xlabel("EDN1-positive cells (%)", labelpad=1.5); grid(b, "x"); b.tick_params(axis="y", length=0)
    for yi, (v, n) in enumerate(zip(cl.edn, cl.cells)): b.text(v + 0.25, yi, f"{v:.1f}  (n = {n:,})", va="center", fontsize=5.8)
    b.set_xlim(0, 17)
    legend_at(fig, 100, 65, [Patch(facecolor="#F2A65A", label="Capillary"), Patch(facecolor="#B2182B", label="Arterial"), Patch(facecolor="#5E3C99", label="Venous")])
    genes = [("SEMA3G", "arterial marker"), ("CA4", "capillary marker"), ("ACKR1", "venous marker"), ("EDN1", "endothelin-1"), ("KLF2", "flow response"), ("SMAD6", "top gene of EDN1+ cells")]
    head(fig, 0, 73, "C", "Marker genes on the same map", "Log-normalized expression; grey = not detected")
    for j, (g, lab) in enumerate(genes):
        ax = fig.add_axes(_rect(fig, 2 + j * 29.8, 88, 28, 52)); fig._boxes.append(("C", (j * 30, 73, 30, 77), [ax]))
        v = V[g].values; pos = v > 0; o = np.argsort(v[pos])
        ax.scatter(V.UMAP1[~pos], V.UMAP2[~pos], s=0.3, color="#DDDDDD", linewidths=0, rasterized=True)
        scm = ax.scatter(V.UMAP1[pos].values[o], V.UMAP2[pos].values[o], s=0.7, c=v[pos][o], cmap=SEQ, vmin=0, vmax=np.percentile(v[pos], 98), linewidths=0, rasterized=True)
        umap_axes(ax, V); ax.set_title(g, fontsize=7.6, fontweight="bold", fontstyle="italic", pad=9); ax.text(0.5, 1.005, lab, transform=ax.transAxes, ha="center", va="bottom", fontsize=5.8, color=GREY)
        cax = fig.add_axes(_rect(fig, 2 + j * 29.8 + 6, 143, 16, 1.6)); cb = fig.colorbar(scm, cax=cax, orientation="horizontal"); cb.outline.set_visible(False); cb.ax.tick_params(length=1.5, labelsize=5.2, pad=1)
    finish(fig, name)


# ============================ S2: calibration ============================
def fig_s_calib(name="FigureS2_permutation_calibration"):
    fig = new_fig(180, 112)
    progs = ["NO_eNOS", "Endothelin_ACE", "Barrier", "FA_transport", "Angiogenic_tip", "X_escape"]
    for r, (A, lab, col, n) in enumerate([(A1, "Healthy hearts (45 donors)", COH["Healthy"], 45), (A2, "Diseased hearts (93 donors)", COH["Diseased"], 93)]):
        axs = panel_grid(fig, (0, r * 56, 180, 56), "AB"[r], lab + ": observed sex difference against 1,000 shuffles of the sex labels",
                         "Grey = estimates with shuffled labels; line = observed estimate; false positives = share of shuffles with p < 0.05 (expected 5%)", ncol=6, l=9, r=2, t=17, b=9, gap=2.5)
        pm = A["perm"]; mn = A["main"].set_index("program"); cb = A["calib"].set_index("program")
        for ax, p in zip(axs, progs):
            v = pm.est[pm.program.astype(str) == p].values; ob = mn.loc[p, "est"]
            hh = ax.hist(v, bins=30, color="#C9C9C9", edgecolor="white", linewidth=0.2); top = hh[0].max(); ax.set_ylim(0, top * 1.5)
            ax.vlines(ob, 0, top * 1.08, color=col if p != "X_escape" else "#C5283D", linewidth=1.6)
            ax.set_title(PROG[p].replace("\n", " ").replace(" (positive control)", "").replace(" / ", " /\n").replace("Fatty-acid ", "Fatty-acid\n").replace("Angiogenic ", "Angiogenic\n"), fontsize=6.8, fontweight="bold", pad=2)
            ax.set_xlim(-2.4, 2.4); ax.set_yticks([]); ax.spines["left"].set_visible(False)
            ax.text(0.03, 0.97, f"false positives {100 * cb.loc[p, 'false_pos_rate_at_0.05']:.1f}%\npermutation p = {cb.loc[p, 'p_perm']:.3f}", transform=ax.transAxes, va="top", fontsize=5.6, color=GREY, linespacing=1.1)
        axs[2].set_xlabel("Female minus male (SD units)", labelpad=1.5, x=1.05)
    finish(fig, name)


# ============================ S3: hormones, immune, states ============================
def fig_s_nulls(name="FigureS3_hormones_immune_states"):
    fig = new_fig(180, 150)
    H = D["stage3_hormone"]["expr"]; cts = ["capillary", "arterial", "venous", "mural", "fibroblast", "cardiomyocyte", "myeloid", "lymphoid"]; rec = ["ESR1", "ESR2", "GPER1", "PGR", "AR", "CYP19A1"]
    M = H.groupby(["gene", "cell_type"], observed=True).log2cpm.median().unstack().loc[rec, cts]
    ax = panel(fig, (0, 0, 92, 72), "A", "Sex-hormone receptors by cell type", "Median log2 CPM; ESR2 is likely inflated by reads from a neighbouring gene", l=14, r=12, b=15)
    sm = heat(ax, M, cmap=SEQ, fmt="{:.1f}", norm=Normalize(-1, 8)); ax.set_xticklabels(cts, rotation=32, ha="right", rotation_mode="anchor"); ax.set_yticklabels(rec, fontstyle="italic")
    cbar_at(fig, sm, (83, 16, 2.4, 34), "log2 CPM", vertical=True)
    I = D["stage3_immune"]; er = I["estrogen_oxphos"]; er = er[(er.outcome == "Estrogen_early") & (er.adjusted == "none") & (er.cohort == "both")].merge(I["esr1"], on="cell_type")
    ax = panel(fig, (92, 0, 88, 72), "B", "Estrogen response does not follow the receptor", "Each point is a cell type; Hallmark early estrogen response, 138 donors", l=16, r=3, b=10)
    ax.axhline(0, color="#333333", linewidth=0.7); grid(ax)
    for _, rw in er.iterrows():
        ax.vlines(rw.ESR1_log2cpm, rw.lo, rw.hi, color=CELL[rw.cell_type], linewidth=1.2); ax.plot(rw.ESR1_log2cpm, rw.est, "o", color=CELL[rw.cell_type], markersize=5)
    repel(ax, er.ESR1_log2cpm.values, er.est.values, er.cell_type.values, fs=6, dx=5, dy=5)
    ax.set_xlabel("ESR1 expression (median log2 CPM)", labelpad=1.5); ax.set_ylabel("Estrogen response,\nfemale minus male (SD)", labelpad=1)
    im = I["immune"]; im = im[(im.cohort == "both") & (im.adjusted == "none") & (im.outcome != "X_escape")].sort_values("est")
    ax = panel(fig, (0, 76, 96, 74), "C", "Immune-cell programs: no evidence of a sex difference", "Myeloid and lymphoid cells, 138 donors, 95% CI", l=44, r=3, b=10)
    zero(ax); grid(ax, "x"); ax.axvspan(-0.8, 0.8, color="#7B3294", alpha=0.07, linewidth=0)
    for i, (_, rw) in enumerate(im.iterrows()): forest(ax, i, rw.est, rw.lo, rw.hi, color=CELL[rw.cell_type], ms=4)
    ax.set_yticks(range(len(im))); ax.set_yticklabels([f"{c.capitalize()}: {o.replace('Interferon_I', 'type I interferon').replace('MHC_II', 'MHC class II').replace('_', ' ').lower().replace('mhc', 'MHC')}" for c, o in zip(im.cell_type, im.outcome)]); ax.set_ylim(-0.6, len(im) - 0.4)
    ax.set_xlim(-1.25, 1.25); ax.set_xlabel("Female minus male (SD units)", labelpad=1.5); ax.tick_params(axis="y", length=0)
    st = D["stage3_states_var"]["states"]; st = st[(st.cohort == "both") & (st.adjusted == "none")].sort_values("est")
    nm = {"hi_stress": "Stress-high cells", "hi_interferon": "Interferon-high cells", "hi_inflammatory": "Inflammation-high cells", "hi_angiogenic": "Angiogenic-high cells", "hi_proliferating": "Proliferating cells", "axis": "Arterial-to-venous position"}
    ax = panel(fig, (96, 76, 84, 74), "D", "Capillary cell states: no evidence of a difference", "Share of capillary cells in each state per donor, 95% CI", l=31, r=3, b=10)
    zero(ax); grid(ax, "x"); ax.axvspan(-0.8, 0.8, color="#7B3294", alpha=0.07, linewidth=0)
    for i, (_, rw) in enumerate(st.iterrows()): forest(ax, i, rw.est, rw.lo, rw.hi, color=CELL["capillary"], ms=4)
    ax.set_yticks(range(len(st))); ax.set_yticklabels([nm[o] for o in st.outcome]); ax.set_ylim(-0.6, len(st) - 0.4); ax.set_xlim(-1.25, 1.25)
    ax.set_xlabel("Female minus male (SD units)", labelpad=1.5); ax.tick_params(axis="y", length=0)
    finish(fig, name)


# ============================ S4: mural cells ============================
def fig_s_mural(name="FigureS4_mural_cells"):
    fig = new_fig(180, 78)
    per = K["per"]; d = per[per.program == "mural_per_capillary"].assign(S=lambda x: x.stratum.map(stratum_name)).set_index("S").loc[SORDER]; c = K["combined"].set_index("program").loc["mural_per_capillary"]
    ax = panel(fig, (0, 0, 80, 78), "A", "Mural cells per capillary cell", "Pericytes (healthy) or pericytes + smooth muscle (diseased)", l=17, r=3, b=10)
    zero(ax); grid(ax, "x"); ax.axvspan(-0.8, 0.8, color="#7B3294", alpha=0.07, linewidth=0)
    for i, (s, rw) in enumerate(d.iterrows()): forest(ax, len(d) - i, rw.est, rw.est - 1.96 * rw.se, rw.est + 1.96 * rw.se, color=STRATA[s], ms=2.5 + 0.09 * (rw.n_w + rw.n_m))
    ax.hlines(0, c.lo, c.hi, color=INK, linewidth=1.8); ax.plot(c.est, 0, "D", markersize=6, markerfacecolor="white", markeredgecolor=INK, markeredgewidth=1.2)
    ax.set_yticks(range(len(d) + 1)); ax.set_yticklabels(["Pooled"] + SORDER[::-1]); ax.get_yticklabels()[0].set_fontweight("bold"); ax.set_ylim(-0.6, len(d) + 0.5)
    ax.set_xlabel("Female minus male (SD units), 95% CI", labelpad=1.5); ax.tick_params(axis="y", length=0)
    ax.text(0.98, 0.04, f"pooled {c.est:+.2f}, p = {c.p:.3f}", transform=ax.transAxes, ha="right", fontsize=6.2, color=GREY)
    h = D["aim5"]["res"].set_index("program"); dd = D["p2_pericytes_xy"]["peri"].set_index("program")
    nm = {"Pericyte_identity": "Pericyte identity", "Contractile": "Contractile", "KATP_channel": "K-ATP channel", "Constrictor_receptors": "Constrictor receptors", "NO_cGMP_response": "NO–cGMP response", "X_escape": "X-escape genes (control)"}
    ax = panel(fig, (80, 0, 100, 78), "B", "Mural-cell programs", "Healthy: pericytes. Diseased: mural cells. 95% CI", l=31, r=3, b=14)
    zero(ax); grid(ax, "x"); ax.axvspan(-0.8, 0.8, color="#7B3294", alpha=0.07, linewidth=0)
    for i, kk in enumerate(nm):
        y = len(nm) - 1 - i
        if kk in h.index: forest(ax, y + 0.17, h.loc[kk, "est"], h.loc[kk, "lo"], h.loc[kk, "hi"], color=COH["Healthy"], ms=3.8)
        if kk in dd.index: forest(ax, y - 0.17, dd.loc[kk, "est"], dd.loc[kk, "lo"], dd.loc[kk, "hi"], color=COH["Diseased"], ms=3.8)
    ax.set_yticks(range(len(nm))); ax.set_yticklabels(list(nm.values())[::-1]); ax.set_ylim(-0.6, len(nm) - 0.4); ax.set_xlabel("Female minus male (SD units)", labelpad=1.5); ax.tick_params(axis="y", length=0)
    legend_at(fig, 145, 72.5, [dot(COH["Healthy"], "Healthy"), dot(COH["Diseased"], "Diseased")])
    finish(fig, name)


# ============================ S5: endothelin detail ============================
def fig_s_edn(name="FigureS5_endothelin_detail"):
    fig = new_fig(180, 150)
    ps = EN["per_stratum"]; ps = ps[ps.gene == "EDN1"]; gg = EN["genes"]; gg = gg[(gg.gene == "EDN1") & (gg.cohort == "both")].set_index("cell_type")
    cts = ["capillary", "arterial", "venous", "fibroblast", "cardiomyocyte"]
    axs = panel_grid(fig, (0, 0, 180, 74), "A", "EDN1 by cell type, in every dataset and disease", "Female minus male (log2), 95% CI; diamond = pooled. EDN1 is made mainly by endothelial cells", ncol=5, l=17, r=3, t=16, b=10, gap=3)
    for ax, ct in zip(axs, cts):
        d = ps[ps.cell_type == ct].set_index("group"); zero(ax); grid(ax, "x")
        for i, s in enumerate(SORDER):
            if s in d.index:
                rw = d.loc[s]; forest(ax, len(SORDER) - i, rw.est, rw.est - 1.96 * rw.se, rw.est + 1.96 * rw.se, color=STRATA[s], ms=2.5 + 0.08 * (rw.n_w + rw.n_m))
        c = gg.loc[ct]; ax.hlines(0, c.lo, c.hi, color=INK, linewidth=1.8); ax.plot(c.est, 0, "D", markersize=5.5, markerfacecolor=FEM if c.p < 0.05 else "white", markeredgecolor=INK, markeredgewidth=1.1)
        ax.set_title(ct, fontsize=7.4, fontweight="bold", color=CELL[ct], pad=3); ax.set_ylim(-0.6, len(SORDER) + 0.5); ax.set_xlim(-4.2, 5.6); ax.tick_params(axis="y", length=0)
        ax.text(0.97, 0.03, f"p = {c.p:.2f}", transform=ax.transAxes, ha="right", fontsize=5.8, color=GREY)
    axs[0].set_yticks(range(len(SORDER) + 1)); axs[0].set_yticklabels(["Pooled"] + SORDER[::-1]); axs[0].get_yticklabels()[0].set_fontweight("bold")
    axs[2].set_xlabel("Female minus male (log2)", labelpad=1.5)
    pc = EN["positive_cells"]
    ax = panel(fig, (0, 78, 62, 72), "B", "Share of EDN1-positive cells", "With and without adjustment for sequencing depth", l=17, r=3, b=17)
    zero(ax); grid(ax, "x")
    for i, co in enumerate(["healthy", "disease", "both"]):
        y = 2 - i
        for adjn, off, col in (("none", 0.17, "#C5283D"), ("+ sequencing depth per cell", -0.17, "#8C8C8C")):
            rw = pc[(pc.cohort == co) & (pc.adjusted == adjn)].iloc[0]; forest(ax, y + off, rw.est, rw.lo, rw.hi, color=col, ms=4)
    ax.set_yticks([0, 1, 2]); ax.set_yticklabels(["Pooled", "Diseased", "Healthy"]); ax.get_yticklabels()[0].set_fontweight("bold"); ax.set_ylim(-0.6, 2.6); ax.set_xlabel("Female minus male (SD units)", labelpad=1.5); ax.tick_params(axis="y", length=0)
    legend_at(fig, 36, 140.5, [dot("#C5283D", "Unadjusted"), dot("#8C8C8C", "Depth-adjusted")], columnspacing=0.9)
    tr = EN["tracks"]; tr = tr[tr.strata == 6].sort_values("slope_SD_per_SD")
    nm = {"Stress": "Tissue stress", "CM_ambient": "Ambient muscle RNA", "Hypoxia": "Hypoxia response", "Shear_KLF2": "Flow response (KLF2)", "Inflammation": "Inflammatory activation", "NO_eNOS": "Nitric oxide / eNOS",
          "ETB_receptor": "ET-B receptor (EDNRB)", "Estrogen_response": "Estrogen response", "age_z": "Age"}
    tr = tr[tr.feature.isin(nm)]
    ax = panel(fig, (62, 78, 60, 72), "C", "What tracks capillary EDN1?", "Slope across 138 donors, within stratum", l=27, r=6, b=14)
    zero(ax); grid(ax, "x")
    for i, (_, rw) in enumerate(tr.iterrows()):
        col = "#C5283D" if rw.p < 0.05 else "#9A9A9A"; ax.hlines(i, 0, rw.slope_SD_per_SD, color=col, linewidth=1.4); ax.plot(rw.slope_SD_per_SD, i, "o", color=col, markersize=4.5)
        ax.text(0.345, i, f"p = {rw.p:.2f}", va="center", ha="right", fontsize=5.6, color="#C5283D" if rw.p < 0.05 else GREY)
    ax.set_yticks(range(len(tr))); ax.set_yticklabels([nm[f] for f in tr.feature]); ax.set_ylim(-0.6, len(tr) - 0.4); ax.set_xlim(-0.25, 0.35); ax.set_xlabel("SD of EDN1 per SD of feature", labelpad=1.5); ax.tick_params(axis="y", length=0)
    dn = EC["donors"]; pcell = EC["per_cell"]
    ax = panel(fig, (122, 78, 58, 72), "D", "EDN1 per positive cell", "More cells, not more per cell", l=10, r=2, b=14)
    g = [(dn.EDN1_per_pos_cell[(dn.cohort == c) & (dn.sex.astype(str) == "female")], dn.EDN1_per_pos_cell[(dn.cohort == c) & (dn.sex.astype(str) == "male")]) for c in ("healthy", "disease")]
    grid(ax, "y"); sex_violins(ax, g, ["Healthy", "Diseased"]); ax.set_ylabel("EDN1 per positive cell (log CP10k)", labelpad=1)
    legend_at(fig, 154, 140.5, SEX_HANDLES + [MEDIAN_H], columnspacing=0.9)
    finish(fig, name)


# ============================ S6: genetics ============================
def fig_s_gen(name="FigureS6_genetics"):
    fig = new_fig(180, 96)
    W = D["stage3_gwas_wide"]["enrichment"]; e = W[W.trait_class == "cardiovascular"].iloc[::-1]
    nm = {"Cardiomyocyte sex-biased (55 consistent)": "Cardiomyocyte genes (55)", "Cardiomyocyte sex-biased (31 replicated)": "Cardiomyocyte genes (31)", "X-linked, higher in females in both cohorts": "X-escape genes (15)",
          "X members of X-Y pairs (18)": "X–Y pair genes (18)", "Capillary program genes": "Capillary program genes (34)", "Endothelin genes": "Endothelin genes (4)"}
    ax = panel(fig, (0, 0, 62, 96), "A", "No specific genetic enrichment", "Genes mapped to a cardiovascular GWAS hit:\nour sets vs other expressed genes", l=36, r=3, b=14, t=13.5)
    grid(ax, "x")
    for i, (_, rw) in enumerate(e.iterrows()):
        ax.hlines(i, rw.pct_background, rw.pct_set, color="#BBBBBB", linewidth=1.2); ax.plot(rw.pct_background, i, "o", color="#8C8C8C", markersize=4.5); ax.plot(rw.pct_set, i, "o", color="#7B3294", markersize=5)
        ax.text(84, i, f"p = {rw.p:.2f}", va="center", ha="right", fontsize=5.6, color=GREY)
    ax.set_yticks(range(len(e))); ax.set_yticklabels([nm[g] for g in e.gene_set]); ax.set_xlim(-3, 85); ax.set_xticks([0, 20, 40, 60]); ax.set_ylim(-0.6, len(e) - 0.4); ax.set_xlabel("Genes with a hit (%)", labelpad=1.5); ax.tick_params(axis="y", length=0)
    legend_at(fig, 40, 90.5, [dot("#7B3294", "Gene set"), dot("#8C8C8C", "Background")], columnspacing=0.9)
    L = D["stage3_gwas_wide"]["edn1_locus"]; t = L[L.snp == "rs9349379"].sort_values("p").drop_duplicates("mapped_trait")
    keep = {"coronary artery disorder": "Coronary artery disease", "myocardial infarction": "Myocardial infarction", "migraine disorder": "Migraine", "pulse pressure measurement": "Pulse pressure", "systolic blood pressure": "Systolic blood pressure",
            "spontaneous coronary artery dissection": "Spontaneous coronary dissection", "angina pectoris": "Angina", "coronary artery calcification": "Coronary calcification", "cervical artery dissection": "Cervical artery dissection", "fibromuscular dysplasia": "Fibromuscular dysplasia"}
    t = t[t.mapped_trait.isin(keep)].iloc[::-1]
    ax = panel(fig, (62, 0, 68, 96), "B", "One variant, many vessel diseases", "rs9349379 at the EDN1/PHACTR1 locus\n(GWAS Catalog); context, not our data", l=38, r=3, b=14, t=13.5)
    grid(ax, "x"); lp = -np.log10(t.p.values); wp = t.women_predominant.values.astype(bool)
    ax.barh(range(len(t)), lp, color=np.where(wp, FEM, "#9A9A9A"), height=0.66, zorder=2)
    ax.set_yticks(range(len(t))); ax.set_yticklabels([keep[m] for m in t.mapped_trait]); ax.set_xlabel("−log10 p", labelpad=1.5); ax.tick_params(axis="y", length=0); ax.set_ylim(-0.6, len(t) - 0.4)
    legend_at(fig, 100, 90.5, [Patch(facecolor=FEM, label="More common in women"), Patch(facecolor="#9A9A9A", label="Other")], columnspacing=0.9)
    q = EC["eqtl"]; tis = ["Artery_Tibial", "Artery_Aorta", "Artery_Coronary", "Heart_Atrial_Appendage", "Heart_Left_Ventricle"][::-1]
    ax = panel(fig, (130, 0, 50, 96), "C", "Linked to PHACTR1 in GTEx", "rs9349379 as an eQTL; weak for EDN1", l=20, r=3, b=14)
    grid(ax, "x")
    for i, ts in enumerate(tis):
        for gname, off, col in (("PHACTR1", 0.19, "#0072B2"), ("EDN1", -0.19, "#C5283D")):
            rw = q[(q.tissue == ts) & (q.gene == gname)]
            if len(rw): ax.barh(i + off, -np.log10(rw.p.values[0]), height=0.36, color=col, zorder=2)
    ax.axvline(-np.log10(0.05), color="#555555", linewidth=0.6, linestyle=(0, (3, 2))); ax.set_xscale("symlog", linthresh=3); ax.set_xlim(0, 60); ax.set_xticks([0, 1, 2, 3, 10, 40]); ax.set_xticklabels(["0", "1", "2", "3", "10", "40"])
    ax.set_yticks(range(len(tis))); ax.set_yticklabels([s.replace("Artery_", "Artery,\n").replace("Heart_Atrial_Appendage", "Atrial\nappendage").replace("Heart_Left_Ventricle", "Left\nventricle").lower().capitalize() for s in tis], linespacing=1.0)
    ax.set_xlabel("−log10 p", labelpad=1.5); ax.tick_params(axis="y", length=0)
    legend_at(fig, 158, 90.5, [Patch(facecolor="#0072B2", label="PHACTR1"), Patch(facecolor="#C5283D", label="EDN1")], columnspacing=0.9)
    finish(fig, name)


FIGS = {"design": fig_design, "s_umap": fig_s_umap, "s_calib": fig_s_calib, "s_nulls": fig_s_nulls, "s_mural": fig_s_mural, "s_edn": fig_s_edn, "s_gen": fig_s_gen}
if __name__ == "__main__":
    for n in (sys.argv[1:] or FIGS):
        FIGS[n]()
