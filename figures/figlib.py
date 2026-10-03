"""Shared style and helpers for the manuscript figures (matplotlib)."""
import pickle, os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, FancyBboxPatch, Rectangle
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
D = pickle.load(open(os.path.join(HERE, "D.pkl"), "rb"))
OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True)

FEM, MAL = "#D81B60", "#1E88E5"
INK, GREY, LIGHT = "#1A1A1A", "#6B6B6B", "#E9E9E9"
DIV = LinearSegmentedColormap.from_list("div", ["#0B3C7D", "#2166AC", "#8FB8DE", "#FFFFFF", "#F2A08B", "#C5283D", "#7A0F22"])
SEQ = LinearSegmentedColormap.from_list("seq", ["#F6F3FB", "#C9B8E8", "#8E6BBF", "#4A2A7F"])
STRATA = {"Healthy 1": "#4D4D4D", "Healthy 2": "#A6A6A6", "DCM": "#0072B2", "ACM": "#CC79A7", "Infarction": "#E69F00", "Myocarditis": "#009E73"}
CELL = {"capillary": "#D55E00", "arterial": "#E69F00", "venous": "#CC79A7", "mural": "#009E73", "pericyte": "#009E73",
        "fibroblast": "#56B4E9", "cardiomyocyte": "#0072B2", "myeloid": "#8C6BB1", "lymphoid": "#666666", "smooth muscle": "#A6761D"}
COH = {"Healthy": "#4D4D4D", "Diseased": "#E66101"}
AGE = {"20-49 years": "#E66101", "50-79 years": "#5E3C99"}
PROG = {"NO_eNOS": "Nitric oxide / eNOS", "Endothelin_ACE": "Endothelin / ACE", "Barrier": "Barrier / junction",
        "FA_transport": "Fatty-acid transport", "Angiogenic_tip": "Angiogenic tip",
        "mural_per_capillary": "Mural cells per\ncapillary cell", "X_escape": "X-escape genes\n(positive control)"}
MM = 1 / 25.4

plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 7, "axes.titlesize": 8.5, "axes.labelsize": 7.5,
    "xtick.labelsize": 6.8, "ytick.labelsize": 6.8, "legend.fontsize": 6.8, "axes.linewidth": 0.6,
    "axes.edgecolor": "#333333", "axes.labelcolor": INK, "text.color": INK, "xtick.color": "#333333", "ytick.color": "#333333",
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "legend.handletextpad": 0.4,
    "legend.columnspacing": 1.0, "legend.borderaxespad": 0.2, "pdf.fonttype": 42, "ps.fonttype": 42,
    "savefig.dpi": 300, "figure.dpi": 150, "axes.axisbelow": True, "hatch.linewidth": 0.5,
    "figure.constrained_layout.h_pad": 0.05, "figure.constrained_layout.w_pad": 0.05,
})


def stratum_name(x):
    x = str(x)
    if x.startswith("dilated"): return "DCM"
    if x.startswith("arrhythmogenic"): return "ACM"
    if x.startswith("myocardial infarction"): return "Infarction"
    if x.startswith("myocarditis"): return "Myocarditis"
    if x in ("Heart", "364bd0c7-f7fd-48ed-99c1-ae26872b1042") or x.startswith("364bd0c7"): return "Healthy 1"
    return "Healthy 2"


def new_fig(w_mm, h_mm):
    fig = plt.figure(figsize=(w_mm * MM, h_mm * MM)); fig._W, fig._H = w_mm, h_mm; fig._boxes = []
    return fig


def _rect(fig, x, y, w, h):
    return [x / fig._W, 1 - (y + h) / fig._H, w / fig._W, h / fig._H]


def head(fig, x, y, letter, text, sub=None):
    """Panel letter, bold title and grey subtitle at the top-left corner (x, y in mm from the top-left of the figure)."""
    if letter: fig.text((x + 0.3) / fig._W, 1 - (y + 0.2) / fig._H, letter, fontsize=12, fontweight="bold", ha="left", va="top", color=INK)
    tx = x + (5.6 if letter else 0.3)
    fig.text(tx / fig._W, 1 - (y + 0.75) / fig._H, text, fontsize=8.5, fontweight="bold", ha="left", va="top", color=INK)
    if sub: fig.text(tx / fig._W, 1 - (y + 4.6) / fig._H, sub, fontsize=6.6, ha="left", va="top", color=GREY, linespacing=1.15)


def panel(fig, box, letter, text, sub=None, l=10, r=2, t=None, b=8, **kw):
    """One axes inside box = (x, y, w, h) mm. l, r, t, b = room for labels; t defaults to the header height."""
    x, y, w, h = box
    if t is None: t = (9.5 + 2.9 * sub.count("\n") if sub else 6.5) if text else 1
    if text or letter: head(fig, x, y, letter, text, sub)
    ax = fig.add_axes(_rect(fig, x + l, y + t, w - l - r, h - t - b), **kw)
    fig._boxes.append((letter or text, box, [ax])); return ax


def panel_grid(fig, box, letter, text, sub=None, ncol=2, nrow=1, l=10, r=2, t=None, b=8, gap=3, vgap=8, sharey=True, widths=None):
    x, y, w, h = box
    if t is None: t = (9.5 if sub else 6.5) if text else 1
    if text or letter: head(fig, x, y, letter, text, sub)
    widths = np.array(widths or [1] * ncol, float); aw = (w - l - r - gap * (ncol - 1)) * widths / widths.sum()
    ah = (h - t - b - vgap * (nrow - 1)) / nrow; axes = []
    for i in range(nrow):
        for j in range(ncol):
            ax = fig.add_axes(_rect(fig, x + l + aw[:j].sum() + gap * j, y + t + i * (ah + vgap), aw[j], ah), sharey=axes[0] if (sharey and axes) else None)
            if sharey and j > 0: ax.tick_params(labelleft=False)
            axes.append(ax)
    fig._boxes.append((letter or text, box, axes)); return axes


def legend_at(fig, x, y, handles, ncol=None, loc="upper center", **kw):
    """Legend anchored at (x, y) mm from the top-left of the figure."""
    kw = {"columnspacing": 1.6, "handlelength": 1.4, **kw}
    return fig.legend(handles=handles, loc=loc, bbox_to_anchor=(x / fig._W, 1 - y / fig._H), ncol=ncol or len(handles), **kw)


def cbar_at(fig, sm, box, label, vertical=False, ticks=None):
    cax = fig.add_axes(_rect(fig, *box)); cb = fig.colorbar(sm, cax=cax, orientation="vertical" if vertical else "horizontal", ticks=ticks)
    cb.outline.set_visible(False); cb.ax.tick_params(length=2, labelsize=6.2)
    if vertical: cb.ax.set_title(label, fontsize=6.4, pad=3, loc="left")
    else: cb.set_label(label, fontsize=6.6)
    return cb


def grid(ax, axis="both"):
    ax.grid(True, axis=axis, color="#EFEFEF", linewidth=0.6)


def finish(fig, name):
    """Save PDF, TIFF, PNG; then check that nothing is clipped by the page or runs into a neighbouring panel."""
    fig.canvas.draw(); r = fig.canvas.get_renderer(); W, H = fig.bbox.width, fig.bbox.height; k = W / fig._W
    probs = []
    arts = [(t, "text '%s'" % t.get_text()[:25]) for t in fig.texts] + [(l, "legend") for l in fig.legends] + [(a, "axes #%d" % i) for i, a in enumerate(fig.axes)]
    for art, what in arts:
        bb = art.get_tightbbox(r) if not isinstance(art, matplotlib.text.Text) else art.get_window_extent(r)
        if bb is None: continue
        if bb.x0 < -1 or bb.y0 < -1 or bb.x1 > W + 1 or bb.y1 > H + 1:
            probs.append(f"OFF PAGE: {what} ({bb.x0 / k:.0f},{(H - bb.y1) / k:.0f})-({bb.x1 / k:.0f},{(H - bb.y0) / k:.0f}) mm")
    for n, (x, y, w, h), axs in fig._boxes:
        for a in axs:
            bb = a.get_tightbbox(r); e = 0.6 * k
            if bb.x0 < x * k - e or bb.x1 > (x + w) * k + e or bb.y1 > H - y * k + e or bb.y0 < H - (y + h) * k - e:
                probs.append(f"OUT OF BOX: panel {n}: axes ({bb.x0 / k:.0f},{(H - bb.y1) / k:.0f})-({bb.x1 / k:.0f},{(H - bb.y0) / k:.0f}) vs box ({x},{y})-({x + w},{y + h}) mm")
    tb = [(t, t.get_window_extent(r)) for t in fig.texts if t.get_text()] + [(l, l.get_window_extent(r)) for l in fig.legends]
    for i in range(len(tb)):
        for j in range(i + 1, len(tb)):
            if tb[i][1].overlaps(tb[j][1]): probs.append(f"OVERLAP: '{getattr(tb[i][0], 'get_text', lambda: 'legend')()[:20]}' and '{getattr(tb[j][0], 'get_text', lambda: 'legend')()[:20]}'")
        for n, _, axs in fig._boxes:
            for a in axs:
                if tb[i][1].overlaps(a.get_tightbbox(r)): probs.append(f"OVERLAP: '{getattr(tb[i][0], 'get_text', lambda: 'legend')()[:25]}' and axes of {n}")
    p = os.path.join(OUT, name)
    fig.savefig(p + ".png", dpi=300, facecolor="white"); fig.savefig(p + ".pdf", facecolor="white")
    im = Image.open(p + ".png").convert("RGB"); im.save(p + ".tiff", compression="tiff_lzw", dpi=(300, 300))
    im.resize((int(fig.get_figwidth() * 110), int(fig.get_figheight() * 110)), Image.LANCZOS).save(os.path.join(OUT, "preview_" + name + ".png"))
    plt.close(fig)
    print(name, "OK" if not probs else "\n  " + "\n  ".join(sorted(set(probs))))


def heat(ax, M, cmap=DIV, vlim=None, fmt="{:+.1f}", stars=None, fs=6.3, norm=None, na="·"):
    """Annotated heatmap from a DataFrame (rows top to bottom)."""
    V = M.values.astype(float)
    if norm is None:
        norm = Normalize(-vlim, vlim) if vlim else Normalize(np.nanmin(V), np.nanmax(V))
    ax.imshow(np.ma.masked_invalid(V), cmap=cmap, norm=norm, aspect="auto")
    ax.set_facecolor("#F4F4F4")
    for i in range(V.shape[0]):
        for j in range(V.shape[1]):
            v = V[i, j]
            if np.isnan(v):
                ax.text(j, i, na, ha="center", va="center", fontsize=fs, color="#9A9A9A"); continue
            rgba = cmap(norm(v)); lum = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
            s = fmt.format(v) + (stars.values[i, j] if stars is not None else "")
            ax.text(j, i, s, ha="center", va="center", fontsize=fs, color="white" if lum < 0.5 else INK)
    ax.set_xticks(range(V.shape[1])); ax.set_xticklabels(M.columns)
    ax.set_yticks(range(V.shape[0])); ax.set_yticklabels(M.index)
    ax.set_xticks(np.arange(-.5, V.shape[1]), minor=True); ax.set_yticks(np.arange(-.5, V.shape[0]), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.4); ax.tick_params(which="both", length=0)
    for s in ax.spines.values(): s.set_visible(False)
    return matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)


def forest(ax, y, est, lo, hi, color=INK, ms=4.5, lw=1.3, marker="o", **kw):
    ax.hlines(y, lo, hi, color=color, linewidth=lw, zorder=3)
    ax.plot(est, y, marker, color=color, markersize=ms, linestyle="none", zorder=4, **kw)


def zero(ax, v=True):
    (ax.axvline if v else ax.axhline)(0, color="#333333", linewidth=0.7, zorder=1)


def repel(ax, xs, ys, labels, fs=5.8, color=INK, dx=4, dy=4, **kw):
    """Greedy label placement: try offsets around each point, keep the first that does not overlap earlier labels."""
    fig = ax.figure; fig.canvas.draw(); r = fig.canvas.get_renderer()
    placed = []
    pts = ax.transData.transform(np.c_[xs, ys])
    cand = [(dx, dy, "left", "bottom"), (-dx, dy, "right", "bottom"), (dx, -dy, "left", "top"), (-dx, -dy, "right", "top"),
            (dx * 2.2, 0, "left", "center"), (-dx * 2.2, 0, "right", "center"), (0, dy * 2.2, "center", "bottom"), (0, -dy * 2.2, "center", "top"),
            (dx * 3, dy * 3, "left", "bottom"), (-dx * 3, dy * 3, "right", "bottom"), (dx * 3, -dy * 3, "left", "top"), (-dx * 3, -dy * 3, "right", "top")]
    axbb = ax.get_window_extent(r)
    for x, y, lab in zip(xs, ys, labels):
        best = None
        for ox, oy, ha, va in cand:
            t = ax.annotate(lab, (x, y), xytext=(ox, oy), textcoords="offset points", ha=ha, va=va, fontsize=fs, color=color, **kw)
            bb = t.get_window_extent(r).expanded(1.05, 1.15)
            ok = not any(bb.overlaps(b) for b in placed) and bb.x0 > axbb.x0 and bb.x1 < axbb.x1 and bb.y0 > axbb.y0 and bb.y1 < axbb.y1
            ok = ok and not any(bb.contains(px, py) for px, py in pts)
            if ok:
                best = t; placed.append(bb)
                if abs(ox) > dx * 2.5 or abs(oy) > dy * 2.5:
                    t.arrowprops = None
                    ax.annotate("", (x, y), xytext=(ox * 0.8, oy * 0.8), textcoords="offset points", arrowprops=dict(arrowstyle="-", lw=0.4, color="#888888"))
                break
            t.remove()
    return placed


def pfmt(p):
    if p < 1e-4:
        e = int(np.floor(np.log10(p))); m = p / 10 ** e
        return f"{m:.0f}×10$^{{{e}}}$"
    return f"{p:.3f}" if p < 0.01 else f"{p:.2f}"
