#!/usr/bin/env python3
"""
plot_paper2_figures_mpl.py

Paper 2 (artepilina C) figures, re-plotted with matplotlib following the
academic-figure-skill rules (Okabe-Ito palette per microstate, sans-serif
type, despined axes, vector EPS master + 300 dpi PNG preview).

Data-extraction helpers are imported from plot_paper2_figures.py (Veusz
version) so both renderings read exactly the same files.

Usage (veusz_env has matplotlib + numpy):
    LC_ALL=C ~/.conda/envs/veusz_env/bin/python plot_paper2_figures_mpl.py \
        --project-root ~/work_turbomole/artepillin_C_studies \
        --outdir ~/work_turbomole/artepillin_C_studies/paper_2/Figures/eps \
        [--only fig02 fig09 ...]

Widths follow the revtex4-2 (aps) column of paper_2_v2.0.1.tex: 86 mm single column,
178 mm double column. Layout coordinates are kept from the Veusz version
(defined on a 180 mm page) and rescaled horizontally to 178 mm.
"""

import argparse
import os
import re
import sys

import matplotlib as mpl

mpl.use("Agg")

# --- Academic Figure Skill typography + export baseline (COPY VERBATIM from skill refs)
mpl.rcParams.update({
    # serif per user preference (Times New Roman or metric-compatible fallback); marked in the skill
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "Liberation Serif", "Nimbus Roman", "STIX"],
    "mathtext.fontset": "stix",
    "font.size": 10,
    "axes.titlesize": 10,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "axes.spines.top": True,      # closed frame (user preference 2026-10-15)
    "axes.spines.right": True,
    "axes.linewidth": 0.6,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": "tight",
    "savefig.dpi": 300,
    "mathtext.default": "regular",
})

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
# data helpers shared with the Veusz builder (single source of data logic)
from plot_paper2_figures import (  # noqa: E402
    STATES, SNAME, COL, GRAY, HC, read_csv, sdir, whitespace_table, state_spectrum,
    loo_mean_control, sticks,
)
# Origin default sequence (user, 2026-10-15): 1 black, 2 red, 3 blue, 4 magenta.
# Microstates keep their order N, A, B, D -> black, red, blue, magenta.
COL = {"N": "#000000", "A": "#FF0000", "B": "#0000FF", "D": "#FF00FF"}

COL_W = 86.0      # single column (mm)
DBL_W = 178.0     # double column (mm)
BASE_W = 180.0    # page width the layout coordinates were written for


# ----------------------------------------------------------------------------- helpers
def new_fig(width_mm, height_mm):
    return plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4))


def add_ax(fig, l, t, w, h, width_mm, height_mm):
    """Axes from cm offsets (top-left origin, as in the Veusz layout) on the 18 cm grid, rescaled to width_mm."""
    l, t, w, h = (10 * v for v in (l, t, w, h))  # Veusz layout is in cm
    k = width_mm / BASE_W
    return fig.add_axes([l * k / width_mm, 1 - (t + h) / height_mm, w * k / width_mm, h / height_mm])


def add_ax_col(fig, t_cm, h_cm, H_mm):
    """Axes in a single-column page (86 mm wide, height H_mm): left 1.5 cm, width 6.6 cm.
    Positions are absolute (no rescaling) so the text keeps its point size."""
    Wc, Hc = COL_W / 10, H_mm / 10
    return fig.add_axes([1.5 / Wc, 1 - (t_cm + h_cm) / Hc, 6.6 / Wc, h_cm / Hc])


def style(ax, xlabel="", ylabel="", xlim=None, ylim=None, xlog=False, ylog=False):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")
    ax.tick_params(length=2.5, width=0.6, pad=1.5)


def letter(ax, s):
    # above the top-left corner: clear of the y-axis label and of the data
    ax.text(-0.02, 1.02, f"{s})", transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left")


def cat_axis(ax, positions, labels, rotation=0, fontsize=9):
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=rotation, ha="right" if rotation else "center", fontsize=fontsize)


def note(ax, s, x=0.03, y=0.04, fontsize=9, ha="left"):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=fontsize, ha=ha, va="bottom", color="#333333")


LEGEND_CANDIDATES = ["upper left", "upper right", "lower right", "lower left",
                     "center right", "center left", "upper center", "lower center"]


def _data_display_points(ax):
    """All plotted data (lines densified, markers/errorbars) in display pixels."""
    pts = []
    for ln in ax.lines:
        xy = np.asarray(ln.get_xydata(), float)
        xy = xy[np.isfinite(xy).all(axis=1)]
        if len(xy) == 0:
            continue
        d = ax.transData.transform(xy)
        if len(d) > 1:  # densify so a line crossing the legend is caught
            t = np.linspace(0, 1, 16)[:, None]
            d = np.vstack([d] + [d[i] * (1 - t) + d[i + 1] * t for i in range(len(d) - 1)])
        pts.append(d)
    for col in ax.collections:
        off = np.asarray(col.get_offsets(), float)
        if len(off):
            off = off[np.isfinite(off).all(axis=1)]
            if len(off):
                pts.append(ax.transData.transform(off))
    return np.vstack(pts) if pts else np.empty((0, 2))


def place_legend(ax, **kw):
    """Legend inside the frame, at the corner that touches no data.

    Tries each corner; if none is clear, adds headroom above the data once (linear y only).
    Falls back to the least-overlapping corner and reports it.
    """
    kw = dict(kw)
    kw.pop("loc", None)
    kw.pop("bbox_to_anchor", None)
    fig = ax.figure
    renderer = fig.canvas.get_renderer()
    best = None
    for attempt in range(4):
        for loc in LEGEND_CANDIDATES:
            leg = ax.legend(loc=loc, **kw)
            bb = leg.get_window_extent(renderer)
            axbb = ax.get_window_extent(renderer)
            inside = (bb.x0 >= axbb.x0) and (bb.x1 <= axbb.x1) and (bb.y0 >= axbb.y0) and (bb.y1 <= axbb.y1)
            P = _data_display_points(ax)
            hits = int(np.sum((P[:, 0] >= bb.x0 - 3) & (P[:, 0] <= bb.x1 + 3) &
                              (P[:, 1] >= bb.y0 - 3) & (P[:, 1] <= bb.y1 + 3))) if len(P) else 0
            score = hits + (0 if inside else 10000)
            if best is None or score < best[0]:
                best = (score, loc)
            if score == 0:
                return leg
            leg.remove()
        if attempt < 3 and ax.get_yscale() == "linear":
            lo, hi = ax.get_ylim()
            ax.set_ylim(lo, hi + 0.35 * (hi - lo))  # headroom above the data, grown until the legend fits
        elif attempt < 3 and ax.get_yscale() == "log":
            lo, hi = ax.get_ylim()
            ax.set_ylim(lo, hi * (hi / lo) ** 0.35)  # same headroom, in log space
        else:
            break
    leg = ax.legend(loc=best[1], **kw)
    lb = leg.get_window_extent(renderer); ab = ax.get_window_extent(renderer)
    print(f"  !! no clear legend corner on '{ax.get_xlabel()}' using {best[1]} (score {best[0]}); "
          f"legend {lb.width:.0f}x{lb.height:.0f}px vs axes {ab.width:.0f}x{ab.height:.0f}px")
    return leg
    return ax.legend(loc=best[1], **kw)


def save(fig, name, out):
    base = os.path.join(out, name)
    fig.savefig(base + ".eps", format="eps")
    fig.savefig(base + ".pdf")  # pdflatex companion (paper_2_v2.0.1.tex includes PDF)
    fig.savefig(base + ".png", dpi=300)
    plt.close(fig)
    print("saved", base + ".eps")


def dots(ax, x, y, color, marker="o", ms=2.5, label=None, zorder=3):
    ax.plot(x, y, linestyle="none", marker=marker, ms=ms, mfc=color, mec=color, mew=0.5,
            label=label, zorder=zorder)


def line(ax, x, y, color="black", lw=1.2, ls="-", label=None):
    ax.plot(x, y, color=color, lw=lw, ls=ls, label=label, solid_capstyle="butt")


def hz(ax, y, x0, x1, color=GRAY, ls="--", lw=0.8):
    ax.plot([x0, x1], [y, y], color=color, lw=lw, ls=ls)


def vt(ax, x, y0, y1, color=GRAY, ls="--", lw=0.8):
    ax.plot([x, x], [y0, y1], color=color, lw=lw, ls=ls)


# ----------------------------------------------------------------------------- figures
def fig01_workflow(root, out):
    """Funnel schematic; numbers are the real counts (boxes/arrows drawn in data coordinates)."""
    W, H = DBL_W, 110.0
    fig = new_fig(W, H)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    boxes = [
        ("b1", 0, 0, "CREST iMTD-GC / GFN2-xTB", "4 microstates\n337 | 322 | 37 | 21 conformers", "#d9e8f5"),
        ("b2", 1, 0, "Boltzmann filter (95%, 298 K)", "weight x CREST degeneracy\nsurvivors: 170 | 143 | 11 | 3", "#d9e8f5"),
        ("b3", 2, 0, "r2SCAN-3c re-optimisation", "327 conformers; re-ranked\n(top conformer changes in 3/4 states)", "#d9e8f5"),
        ("b4", 2, 1, "Layer A: PBE0/def2-SVP", "IR, UV-Vis; population-weighted\nscreening spectra (327)", "#e6f2e6"),
        ("b5", 1, 1, "Layer B: HSE06-D3(BJ)", "def2-TZVP; 54 labelled + extensions\n(same level as Paper 1)", "#e6f2e6"),
        ("b6", 0, 1, "Delta-learning surrogates", "B1 $\\nu$(C=O): GFN2-xTB → HSE06\nB2 S1 (E, f): sTDA → HSE06", "#e6f2e6"),
        ("b7", 0, 2, "Statistical audit", "LOO $Q^2$ + permutation test\nnegative results kept", "#fbeadb"),
        ("b8", 1, 2, "Thermodynamic cycle → pKa", "freeh + COSMO(water)\nclass-wise calibration vs references", "#fbeadb"),
        ("b9", 2, 2, "Outputs", "ensemble IR/UV-Vis per microstate\nspeciation, spectra vs pH", "#f3e0ef"),
    ]
    bw, bh = 0.27, 0.22
    xs = [0.03, 0.365, 0.70]
    ys = [0.70, 0.40, 0.10]
    pos = {}
    for bid, c, r, title, sub, colr in boxes:
        x0, y0 = xs[c], ys[r]
        pos[bid] = (x0, y0)
        ax.add_patch(FancyBboxPatch((x0, y0), bw, bh, boxstyle="round,pad=0,rounding_size=0.01",
                                    fc=colr, ec="#444444", lw=0.6))
        ax.text(x0 + bw / 2, y0 + bh * 0.66, title, ha="center", va="center", fontsize=9, fontweight="bold")
        ax.text(x0 + bw / 2, y0 + bh * 0.28, sub, ha="center", va="center", fontsize=9, linespacing=1.3)
    arrows = [("b1", "b2"), ("b2", "b3"), ("b3", "b4"), ("b4", "b5"), ("b5", "b6"), ("b6", "b7"),
              ("b7", "b8"), ("b8", "b9")]
    for a, b in arrows:
        xa, ya = pos[a]
        xb, yb = pos[b]
        if abs(ya - yb) < 1e-9:
            p0, p1 = ((xa + bw, ya + bh / 2), (xb, yb + bh / 2)) if xb > xa else \
                     ((xa, ya + bh / 2), (xb + bw, yb + bh / 2))
        else:
            p0, p1 = (xa + bw / 2, ya), (xb + bw / 2, yb + bh)
        ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=7, lw=0.9, color="#222222"))
    ax.text(0.5, 0.97, "Funnel: ensemble × protonation state → spectra and pKa (artepilin C, C$_{19}$H$_{24}$O$_3$)",
            ha="center", va="top", fontsize=9, fontweight="bold")
    save(fig, "Fig01_workflow", out)


def fig02_ensembles(root, out):
    W, H = DBL_W, 132.0  # two rows: (a) on top, (b) and (c) below, aligned to (a)
    fig = new_fig(W, H)
    # (a) cumulative population vs rank (r2SCAN-3c, degeneracy-weighted)
    ax = add_ax(fig, 1.6, 0.6, 16.0, 5.0, W, H)
    kept = {}
    for s in STATES:
        rows = read_csv(os.path.join(sdir(root, s), "layer_a", "population_table_dft.csv"))
        p = np.array(sorted([float(r["r2scan3c_population"]) for r in rows], reverse=True))
        p = p / p.sum()
        n = len(p)
        # round dots on a log-spaced subset of ranks, so the dense curve stays readable
        mark = np.unique(np.logspace(0, np.log10(n - 1), 12).astype(int))
        ax.plot(np.arange(1, n + 1), np.cumsum(p), color=COL[s], lw=1.2, marker="o", ms=2.6,
                mfc=COL[s], mec=COL[s], mew=0.4, markevery=mark, label=SNAME[s])
        kept[s] = n
    # y-range extended to 1.4 so the legend sits in the empty band above the curves (max 1.0)
    style(ax, "Conformer rank", "Cumulative population", (0.8, 200), (0, 1.4), xlog=True)
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    # one column inside the axes: two columns overflowed the right spine
    place_legend(ax, loc="upper left", ncol=1, handlelength=2.0)
    letter(ax, "a")
    # (b) GFN2 vs r2SCAN-3c relative energy (N and A)
    ax = add_ax(fig, 1.6, 7.2, 7.0, 5.0, W, H)
    allx, ally = [], []
    for s in ("N", "A"):
        rows = read_csv(os.path.join(sdir(root, s), "layer_a", "population_table_dft.csv"))
        x = [float(r["gfn2_relative_energy_kcal"]) for r in rows]
        y = [float(r["r2scan3c_relative_energy_kcal"]) for r in rows]
        allx += x
        ally += y
        dots(ax, x, y, COL[s], ms=2.2, label={"N": "Neutral", "A": "Monoanion A"}[s])  # short: narrow panel
    rx = np.argsort(np.argsort(allx))
    ry = np.argsort(np.argsort(ally))
    rho = np.corrcoef(rx, ry)[0, 1]
    style(ax, "GFN2-xTB rel. E (kcal mol$^{-1}$)", "r2SCAN-3c rel. E (kcal mol$^{-1}$)", (0, 1.1), (0, 3.0))
    place_legend(ax, loc="upper left")
    note(ax, f"$\\rho$ = {rho:.2f} (N+A)", x=1.0, y=1.02, ha="right")
    letter(ax, "b")
    # (c) conformer counts
    ax = add_ax(fig, 10.4, 7.2, 7.2, 5.0, W, H)
    crest = {"N": 337, "A": 322, "B": 37, "D": 21}
    lab = {"N": 20, "A": 20, "B": 11, "D": 3}
    xs = [1, 2, 3, 4]
    # markers joined by lines: one series per filtering stage
    for vals, color, lbl in (([crest[s] for s in STATES], "#999999", "CREST"),
                             ([kept[s] for s in STATES], "#222222", "95% filter"),
                             ([lab[s] for s in STATES], "#FF0000", "HSE06 labelled")):
        ax.plot(xs, vals, color=color, lw=0.9, marker="o", ms=3.2, mfc=color, mec=color, label=lbl)
    style(ax, "", "Number of conformers", (0.4, 4.6), (0, 400))
    cat_axis(ax, xs, ["N", "A", "B", "D"])
    place_legend(ax, loc="upper right")
    letter(ax, "c")
    save(fig, "Fig02_ensembles", out)


def fig03_functional(root, out):
    W, H = DBL_W, 70.0
    fig = new_fig(W, H)
    fc = read_csv(os.path.join(root, "functional_comparison.csv"))
    names = [r["functional"].upper().replace("CAM-B3LYP", "CAM-\nB3LYP") for r in fc]
    e = [float(r["energy_eV"]) for r in fc]
    f = [float(r["osc_len"]) for r in fc]
    xs = [1, 2, 3]
    ax = add_ax(fig, 1.6, 0.6, 4.6, 5.2, W, H)
    dots(ax, xs, e, "#000000", marker="s", ms=3.5)
    style(ax, "", "S1 vertical energy (eV)", (0.4, 3.6), (3.6, 4.4))
    cat_axis(ax, xs, names, fontsize=9)
    letter(ax, "a")
    ax = add_ax(fig, 7.6, 0.6, 4.4, 5.2, W, H)
    dots(ax, xs, f, "#FF0000", ms=3.5)
    style(ax, "", "S1 oscillator strength $f$ (length)", (0.4, 3.6), (0.4, 0.9))
    cat_axis(ax, xs, names, fontsize=9)
    note(ax, "$f$ spread = CT proxy (Paper 1)", x=1.0, y=1.02, ha="right")
    letter(ax, "b")
    # (c) sTDA vs HSE06 S1 energy, neutral, character matched
    nd = os.path.join(sdir(root, "N"), "layer_b")
    st = {r["job_dir"]: r for r in read_csv(os.path.join(nd, "s1_stda_character_matched_FIXED.csv"))}
    hs = {r["job_dir"]: r for r in read_csv(os.path.join(nd, "s1_hse06.csv"))}
    common = [j for j in hs if j in st]
    xe = [float(st[j]["s1_energy_eV_stda"]) for j in common]
    ye = [float(hs[j]["s1_energy_eV"]) for j in common]
    ax = add_ax(fig, 13.6, 0.6, 4.0, 5.2, W, H)
    dots(ax, xe, ye, "#000000", ms=2.2)
    style(ax, "sTDA S1 (eV)", "HSE06 S1 (eV)", (2.6, 3.0), (3.8, 4.1))
    note(ax, "sTDA ~1.2 eV low", x=1.0, y=1.02, ha="right")
    letter(ax, "c")
    save(fig, "Fig03_functional_and_sTDA", out)

    # companion: f scatter (log-log), used for the B2 delta-learning claim
    # wider page so the long legend fits inside the frame (lower-right corner is empty here)
    fig = new_fig(DBL_W, 80.0)
    ax = fig.add_axes([18 / DBL_W, 1 - (6 + 60) / 80, 130 / DBL_W, 60 / 80])
    xf = [float(st[j]["s1_osc_stda"]) for j in common]
    yf = [float(hs[j]["s1_osc_len"]) for j in common]
    dots(ax, xf, yf, "#000000", ms=2.2, label="neutral conformers shown; fit uses N=28")
    xx = np.logspace(-1.3, 0.18, 20)
    line(ax, xx, 10 ** (-0.2027 + 0.5794 * np.log10(xx)), color="#FF0000", lw=1.2,
         label="log$_{10}f$(HSE06) = −0.2027 + 0.5794 log$_{10}f$(sTDA)")
    line(ax, xx, xx, color=GRAY, lw=0.8, ls="--", label="$y = x$")
    style(ax, "sTDA S1 $f$", "HSE06 S1 $f$", (0.05, 1.5), (0.05, 1.5), xlog=True, ylog=True)
    place_legend(ax, loc="lower right", fontsize=9)
    save(fig, "Fig03b_sTDA_vs_HSE06_f", out)


def fig04_b1(root, out):
    W, H = DBL_W, 120.0
    fig = new_fig(W, H)
    td = {s: read_csv(os.path.join(sdir(root, s), "layer_b", "b1_training_data.csv")) for s in STATES}
    # (a) strip plot of delta per state
    ax = add_ax(fig, 1.8, 0.6, 5.0, 4.6, W, H)
    rng = np.random.default_rng(3)
    for i, s in enumerate(STATES, 1):
        d = [float(r["delta_cm1"]) for r in td[s]]
        dots(ax, i + rng.uniform(-0.15, 0.15, len(d)), d, COL[s], ms=2.2)
    style(ax, "", "HSE06 − GFN2-xTB $\\nu$(C=O) (cm$^{-1}$)", (0.4, 4.6), (-100, 60))
    cat_axis(ax, [1, 2, 3, 4], ["N", "A", "B", "D"])
    letter(ax, "a")
    # (b) monoA: delta vs Bader charge of carbonyl O
    feat = {r["job_dir"]: float(r["bader_charge_carbonyl_o"])
            for r in read_csv(os.path.join(sdir(root, "A"), "layer_b", "population_feature.csv"))}
    xa = [feat[r["job_dir"]] for r in td["A"] if r["job_dir"] in feat]
    ya = [float(r["delta_cm1"]) for r in td["A"] if r["job_dir"] in feat]
    ax = add_ax(fig, 8.2, 0.6, 4.6, 4.6, W, H)
    dots(ax, xa, ya, COL["A"], ms=2.2, label="labelled (N=20)")
    xx = np.linspace(-1.6034, -1.1553, 10)
    line(ax, xx, 4.7223 * xx - 73.9249, label="OLS, $Q^2$ = 0.277")
    style(ax, "Bader $q$(O carbonyl) (e)", "$\\Delta\\nu$(C=O) (cm$^{-1}$)", (-1.65, -1.1), (-82, -76.2))  # headroom for the legend
    place_legend(ax, loc="upper left", fontsize=9)
    letter(ax, "b")
    # (c) monoB: delta vs O...H distance
    oh = {r["job_dir"]: float(r["oh_distance_angstrom"])
          for r in read_csv(os.path.join(sdir(root, "B"), "layer_b", "oh_distance_feature.csv"))}
    xb = [oh[r["job_dir"]] for r in td["B"] if r["job_dir"] in oh]
    yb = [float(r["delta_cm1"]) for r in td["B"] if r["job_dir"] in oh]
    ax = add_ax(fig, 14.0, 0.6, 3.7, 4.6, W, H)
    dots(ax, xb, yb, COL["B"], ms=2.2, label="N=11")
    xx = np.linspace(min(xb), max(xb), 5)
    line(ax, xx, 1392.13 * xx - 3114.30, label="OLS, $Q^2$ = 0.594")
    pad = 0.1 * (max(yb) - min(yb) + 1e-6)
    style(ax, "O···H distance (Å)", "$\\Delta\\nu$(C=O) (cm$^{-1}$)", (2.2605, 2.2635),
          (min(yb) - pad, max(yb) + pad + 0.35 * (max(yb) - min(yb))))  # headroom for the legend
    place_legend(ax, loc="upper left", fontsize=9)
    letter(ax, "c")
    # (d) neutral: LOO prediction vs truth (flat offset)
    loo = read_csv(os.path.join(sdir(root, "N"), "layer_b", "b1_loo_predictions.csv"))
    t = [float(r["true_delta_cm1"]) for r in loo]
    p = [float(r["loo_predicted_delta_cm1"]) for r in loo]
    ax = add_ax(fig, 1.8, 6.8, 5.0, 4.6, W, H)
    line(ax, [27, 30], [27, 30], color=GRAY, lw=0.8, ls="--")
    dots(ax, t, p, COL["N"], ms=2.2)
    style(ax, "True $\\Delta\\nu$ (cm$^{-1}$)", "LOO predicted $\\Delta\\nu$ (cm$^{-1}$)", (27, 30), (27, 30))
    note(ax, "Neutral: no feature beats a flat offset\n28.222 ± 0.676 cm$^{-1}$", x=0.92, y=1.02, ha="right")  # inset so it does not reach the e) label
    letter(ax, "d")
    # (e) HSE06 harmonic nu(C=O) per microstate, unscaled (mean +/- sd)
    ax = add_ax(fig, 8.2, 6.8, 9.5, 4.6, W, H)
    means = [np.mean([float(r["nu_CO_hse06_cm1"]) for r in td[s]]) for s in STATES]
    sds = [np.std([float(r["nu_CO_hse06_cm1"]) for r in td[s]], ddof=1) for s in STATES]
    for i, s in enumerate(STATES, 1):
        ax.errorbar([i], [means[i - 1]], yerr=[sds[i - 1]], fmt="s", ms=3, color=COL[s],
                    ecolor=COL[s], elinewidth=0.8, capsize=2)
    style(ax, "", "HSE06 harmonic $\\nu$(C=O) (cm$^{-1}$)", (0.4, 4.6), (1550, 1850))
    cat_axis(ax, [1, 2, 3, 4], ["N", "A", "B", "D (N=3)"])
    note(ax, "unscaled harmonic; the C=O mode of A and D is the carboxylate stretch", x=1.0, y=1.02, ha="right")
    letter(ax, "e")
    save(fig, "Fig04_B1_carbonyl_shift", out)


def fig05_audit(root, out):
    W, H = COL_W, 113.0  # single column, 2 rows: (a) on top, (b) below
    fig = new_fig(W, H)
    rows = []
    for s in STATES:
        L = os.path.join(sdir(root, s), "layer_b")
        for fn, label in (("b1_permutation_test.dat", "B1"), ("b2_permutation_test_delta_energy_eV.dat", "B2 $\\Delta E$"),
                          ("b2_permutation_test_delta_log10_f.dat", "B2 $\\Delta\\log f$")):
            p = os.path.join(L, fn)
            if os.path.isfile(p):
                lines = [ln for ln in open(p).read().splitlines() if ln.strip()]
                v = lines[1].split()
                rows.append((s, label, int(v[1]), float(v[4]), float(v[6]), float(v[7]), float(v[9]), float(v[11])))
    n = len(rows)
    xs = np.arange(1, n + 1)
    ax = add_ax_col(fig, 0.8, 4.2, H)
    hz(ax, 0, 0.3, n + 0.7, color="black", ls="-", lw=0.6)
    ax.errorbar(xs, [r[4] for r in rows], yerr=[r[5] for r in rows], fmt="s", ms=2.5, color="#999999",
                ecolor="#999999", elinewidth=0.8, capsize=1.5, label="perm. null (± sd)")
    dots(ax, xs, [r[6] for r in rows], "#999999", marker="_", ms=5, label="null 95th pct.")
    dots(ax, xs, [r[3] for r in rows], "#FF0000", marker="D", ms=3, label="observed $Q^2$")
    style(ax, "", "Leave-one-out $Q^2$", (0.3, n + 0.7), (-0.5, 0.75))
    cat_axis(ax, xs, [f"{r[0]} {r[1]}" for r in rows], rotation=30, fontsize=9)
    # p-values are not drawn here (8 labels do not fit one 86 mm column); they go in the caption
    place_legend(ax, loc="upper left", fontsize=9)
    letter(ax, "a")
    # (b) LOO mean-of-the-others control on the neutral B1 labels -> rho = -1 exactly
    y = [float(r["delta_cm1"]) for r in read_csv(os.path.join(sdir(root, "N"), "layer_b", "b1_training_data.csv"))]
    ctrl = loo_mean_control(y)
    ax = add_ax_col(fig, 6.6, 4.2, H)
    dots(ax, y, ctrl, "#000000", ms=2.2)
    style(ax, "True $\\Delta\\nu$ (cm$^{-1}$)", "LOO control (cm$^{-1}$)",
          (min(y) - 0.1, max(y) + 0.1), (min(ctrl) - 0.02, max(ctrl) + 0.02))
    note(ax, "$\\rho$ = −1 (artefact)", x=1.0, y=1.02, ha="right")
    letter(ax, "b")
    save(fig, "Fig05_statistical_audit", out)


def fig06_uvvis_states(root, out):
    W, H = COL_W, 113.0  # single column, 2 rows: (a) on top, (b) below
    fig = new_fig(W, H)
    grid = np.linspace(2.0, 6.0, 800)
    ax = add_ax_col(fig, 0.8, 4.2, H)
    lam = HC / grid
    for s in STATES:
        prof, n, ntot, cov = state_spectrum(root, s, grid)
        line(ax, lam, prof, color=COL[s], lw=1.1,
             label=f"{s}: {n}/{ntot} conf., {cov * 100:.0f}% pop.")  # short: single-column legend
    style(ax, "Wavelength (nm)", "Absorption (a.u.)", (220, 480), (0, 1.25))
    place_legend(ax, loc="upper right", fontsize=9)
    letter(ax, "a")
    cons = read_csv(os.path.join(root, "b2_consolidated_dataset.csv"))
    ax = add_ax_col(fig, 6.6, 4.2, H)
    names = {"neutral": "N", "monoanion_a": "A", "monoanion_b": "B", "dianion": "D"}
    for st, s in names.items():
        for src, mk in (("hse06_ground_truth", "o"), ("surrogate_max_f", "^")):
            sel = [r for r in cons if r["protonation_state"] == st and r["source"] == src]
            if sel:
                dots(ax, [float(r["energy_eV"]) for r in sel], [float(r["osc_f"]) for r in sel], COL[s],
                     marker=mk, ms=2.5, label=f"{s} {'HSE06' if mk == 'o' else 'surrogate'}")
    style(ax, "S1 energy (eV)", "S1 oscillator strength $f$", (2.2, 4.6), (0, 1.0))
    place_legend(ax, loc="upper left", fontsize=9, ncol=2)
    letter(ax, "b")
    save(fig, "Fig06_uvvis_per_microstate", out)


def fig07_exp_and_ph(root, out):
    W, H = COL_W, 113.0  # single column, 2 rows: (a) on top, (b) below
    fig = new_fig(W, H)
    # experiment: Camuri et al., Spectrochim. Acta A 198 (2018) 71 (vector data extracted from the PDF)
    dd = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Figures", "data")
    def read_exp(name):  # '#' comment block, then one header line, then data (nan = outside the curve)
        with open(os.path.join(dd, name)) as f:
            ls = [ln.split() for ln in f if ln.strip() and not ln.startswith("#")]
        return {h: np.array([float(r[i]) for r in ls[1:]]) for i, h in enumerate(ls[0])}
    sp2 = read_exp("exp_camuri2018_fig2_spectra_vs_pH.dat")
    tit = read_exp("exp_camuri2018_fig3_lmax_vs_pH.dat")
    # (a) calculated population-weighted spectra of N and A (same 0.4 eV broadening as Fig. 6) vs the
    # first (pH 2.6) and last (pH 9.2) curves of the experimental titration; all peak-normalised
    ax = add_ax_col(fig, 0.8, 4.2, H)
    g = np.linspace(2.0, 6.0, 1600)
    lam = HC / g
    o = np.argsort(lam)
    for s_, ecol in (("N", "s00"), ("A", "s23")):
        prof = state_spectrum(root, s_, g)[0]
        pr = prof[o] / prof.max()
        line(ax, lam[o], pr, color=COL[s_], lw=1.2, label=f"{s_} calc.")
        e = sp2[ecol]
        ok = ~np.isnan(e)
        line(ax, sp2["wavelength_nm"][ok], e[ok], color=COL[s_], lw=1.0, ls=":", label=f"{s_} exp.")
    style(ax, "Wavelength (nm)", "Normalised absorbance", (250, 400), (0, 1.1))
    place_legend(ax, loc="upper right", fontsize=9)
    letter(ax, "a")
    # (b) band maximum and centroid vs pH from per-state spectra mixed by speciation; the "COSMO-shifted" curve
    # moves the N and A spectra rigidly by the vacuum->COSMO shift of their band maximum measured on the same
    # conformers (anion_band_cosmo.py); B and D have no COSMO-TDDFT and stay in vacuum
    from anion_band_cosmo import run as cosmo_shifts
    sp = whitespace_table(os.path.join(root, "speciation_v4_5refs.dat"))
    grid = np.linspace(2.0, 6.0, 1200)
    spec = {s: state_spectrum(root, s, grid)[0] for s in STATES}
    lam_c = {(r["state"], r["medium"]): r["lambda_max_nm"] for r in cosmo_shifts(root, 0.4, 3.0, 4.96)}
    dE = {s: HC / lam_c[(s, "cosmo")] - HC / lam_c[(s, "vacuum")] for s in ("N", "A")}
    spec_c = {s: (np.interp(grid - dE[s], grid, spec[s]) if s in dE else spec[s]) for s in STATES}
    w = (grid >= 3.0) & (grid <= 5.5)
    phs = np.arange(1.0, 13.01, 0.25)
    mx, mxc, cen = [], [], []
    for ph in phs:
        i = int(np.argmin(np.abs(sp["pH"] - ph)))
        prof = sum(sp[f"f_{s}"][i] * spec[s] for s in STATES)
        profc = sum(sp[f"f_{s}"][i] * spec_c[s] for s in STATES)
        mx.append(HC / grid[w][np.argmax(prof[w])])
        mxc.append(HC / grid[w][np.argmax(profc[w])])
        cen.append(HC / (np.sum(grid[w] * prof[w]) / np.sum(prof[w])))
    ax = add_ax_col(fig, 6.6, 4.2, H)
    line(ax, phs, mx, color="#FF0000", lw=1.2, ls="--", label="calc. max., vacuum")
    line(ax, phs, mxc, color="#FF0000", lw=1.2, label="calc. max., COSMO-shifted")
    line(ax, phs, cen, color="#000000", lw=1.0, label="calc. centroid, vacuum")
    ax.plot(tit["pH"], tit["lambda_max_nm"], ls="none", marker="o", ms=3.2, mfc="white", mec="#000000",
            mew=0.8, label="exp. peak")
    vt(ax, 4.65, 250, 335)
    style(ax, "pH", "Band position (nm)", (1, 13), (250, 425))  # headroom for the 4-entry legend (data end at 322 nm)
    note(ax, "exp. p$K_a$ 4.65", x=1.0, y=1.02, ha="right")
    ax.legend(loc="upper right", fontsize=9)
    letter(ax, "b")
    save(fig, "Fig07_uvvis_vs_experiment_and_pH", out)


def fig08_ir(root, out):
    W, H = COL_W, 113.0  # single column, 2 rows: (a) on top, (b) below
    fig = new_fig(W, H)
    ax = add_ax_col(fig, 0.8, 4.2, H)
    ax2 = add_ax_col(fig, 6.6, 4.2, H)
    for s in STATES:
        d = whitespace_table(os.path.join(sdir(root, s), "layer_a", "ir_population_weighted.dat"))
        k = list(d.keys())
        x, y = d[k[0]], d[k[1]]
        y = y / max(y.max(), 1e-12)
        line(ax, x, y, color=COL[s], lw=0.8, label=SNAME[s])
        m = (x >= 1550) & (x <= 1850)
        line(ax2, x[m], y[m] / max(y[m].max(), 1e-12), color=COL[s], lw=1.0)
    # headroom above the tallest peak (1.0) so the legend sits inside the frame, in two columns
    style(ax, "Wavenumber (cm$^{-1}$)", "Population-weighted IR (norm.)", (400, 4000), (0, 1.5))  # taller page: one-column legend in the headroom
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    place_legend(ax, loc="upper left", ncol=1, handlelength=1.8)
    letter(ax, "a")
    style(ax2, "Wavenumber (cm$^{-1}$)", "IR (norm.)", (1550, 1850), (0, 1.05))
    letter(ax2, "b")
    save(fig, "Fig08_ir_per_microstate", out)


def fig09_pka(root, out):
    W, H = COL_W, 113.0  # single column, 2 rows: (a) on top, (b) below
    fig = new_fig(W, H)
    refs = []
    with open(os.path.join(root, "pka_calibration_5refs.csv")) as f:
        for ln in f:
            m = re.match(r"^(\S+)\s+class=(\S+)\s+raw=([\d.-]+)\s+exp=([\d.-]+)\s+residual=([\d.-]+)", ln)
            if m:
                refs.append((m.group(1), m.group(2), float(m.group(5))))
    seen = []
    for r in refs:
        if r not in seen:
            seen.append(r)
    refs = seen
    n = len(refs)
    ax = add_ax_col(fig, 0.8, 4.2, H)
    xs = np.arange(1, n + 1)
    for i, (nm, cl, res) in enumerate(refs, 1):
        dots(ax, [i], [res], "#000000" if cl == "0to-1" else "#FF0000", ms=3.5)
    short = {"acetic_acid": "acetic", "phenol": "phenol", "4-hydroxybenzoate_phenol": "4-HBA",
             "26-dimethylphenol": "2,6-DMP", "4-coumarate_phenol": "4-coum."}
    c1 = np.mean([r[2] for r in refs if r[1] == "0to-1"])
    c2 = np.mean([r[2] for r in refs if r[1] == "-1to-2"])
    hz(ax, c1, 0.4, n + 0.6, color="#000000")
    hz(ax, c2, 0.4, n + 0.6, color="#FF0000")
    style(ax, "", "pKa residual (exp − raw)", (0.4, n + 0.6), (-5.2, -2.2))
    cat_axis(ax, xs, [short.get(r[0], r[0]) for r in refs], fontsize=9)
    note(ax, "black: 0 → −1 (N=3);  red: −1 → −2 (N=2)", x=1.0, y=1.02, ha="right")
    letter(ax, "a")
    # (b) calibrated artepilin C pKa: vacuum ensemble (production) vs water-extended.
    # The water-extended series uses the consistent definition (G_aq Boltzmann weights with
    # CREST degeneracy and G_ens, scheme 5 of Fig12c), not calib_C of pka_vacuum_vs_water.csv
    # (G_aq weights without degeneracy), which mixed definitions (2026-10-15).
    vw = {r["pair"]: r for r in read_csv(os.path.join(root, "pka_vacuum_vs_water.csv"))}
    dk = {r["pair"]: r for r in read_csv(os.path.join(root, "vac_vs_water_decomp_pka.csv"))
          if r["weights"] == "gaqg" and r["state_G"] == "ens" and r["set"] == "all"}
    pairs = ["N_to_A", "N_to_B", "A_to_D", "B_to_D"]
    xs = np.array([1, 2, 3, 4])
    ax = add_ax_col(fig, 6.6, 4.2, H)
    ax.errorbar(xs - 0.12, [float(vw[p]["calib_A"]) for p in pairs], yerr=[0.41, 0.41, 1.15, 1.15],
                fmt="s", ms=3.5, color="#000000", ecolor="#000000", elinewidth=0.8, capsize=2,
                label="vacuum-CREST ensemble (± class scatter)")
    dots(ax, xs + 0.12, [float(dk[p]["calib_pKa"]) for p in pairs], "#FF0000", ms=3,
         label="+ water-guided conformers, $G_{ens}$")
    dots(ax, [1], [4.65], "black", marker="*", ms=5, label="experiment (carboxyl)")
    style(ax, "", "Calibrated pKa", (0.4, 4.6), (3, 14))
    cat_axis(ax, xs, ["N→A", "N→B", "A→D", "B→D"], fontsize=9)
    place_legend(ax, loc="lower right", fontsize=9)
    note(ax, "−1 → −2 pairs: low confidence (N=2, sd 1.15)", x=1.0, y=1.02, ha="right")
    letter(ax, "b")
    save(fig, "Fig09_pKa_calibration", out)


def fig10_speciation(root, out):
    W, H = COL_W, 113.0  # single column, 2 rows: (a) on top, (b) below
    fig = new_fig(W, H)
    for letter_, fn, title, l in (("a", "speciation_v4_5refs.dat", "vacuum-CREST ensemble", 1.8),
                                  ("b", "speciation_v5_gens_water_ext.dat", "+ water-guided conformers ($G_{ens}$)", 10.2)):
        sp = whitespace_table(os.path.join(root, fn))
        ax = add_ax_col(fig, 0.8 if letter_ == "a" else 6.6, 4.2, H)
        for s in STATES:
            line(ax, sp["pH"], sp[f"f_{s}"], color=COL[s], lw=1.4, label=SNAME[s])
        vt(ax, 4.65, 0, 1.02)
        # headroom above 1.0 (fractions never exceed 1) so the legend sits inside the frame
        style(ax, "pH", "Mole fraction", (0, 14), (0, 1.25))
        ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        note(ax, title, x=0.92, y=1.02, ha="right", fontsize=9)  # inset so it clears the next panel label
        letter(ax, letter_)
        if letter_ == "a":
            place_legend(ax, loc="upper center", ncol=1, handlelength=1.8, fontsize=9)
    save(fig, "Fig10_speciation_vs_pH", out)


def fig11_vac_water(root, out):
    W, H = DBL_W, 140.0  # two rows: (a) on top, (b) and (c) below, aligned to (a)
    fig = new_fig(W, H)
    d = [ln.strip().split(",") for ln in open(os.path.join(root, "vacuum_vs_water_ensembles.dat")) if ln.strip()]
    hdr, data = d[0], d[1:]
    idx = {h: i for i, h in enumerate(hdr)}
    st = [r[idx["state"]] for r in data]
    lab = {"dianion": "D", "monoanion_b": "B", "monoanion_a": "A"}
    xs = np.arange(1, len(st) + 1)
    names = [lab[s] for s in st]
    ax = add_ax(fig, 1.8, 0.6, 16.0, 5.4, W, H)
    dots(ax, xs, [float(r[idx["water_pop_in_kept"]]) for r in data], "#000000", marker="s", ms=3.5)
    hz(ax, 0.95, 0.4, len(st) + 0.6)
    style(ax, "", "Water pop. retained (95%)", (0.4, len(st) + 0.6), (0.6, 1.0))
    cat_axis(ax, xs, names)
    letter(ax, "a")
    ax = add_ax(fig, 1.8, 7.6, 7.0, 5.4, W, H)
    dots(ax, xs - 0.12, [float(r[idx["n_for_95_vac"]]) for r in data], "#999999", marker="s", ms=3.5, label="vacuum")
    dots(ax, xs + 0.12, [float(r[idx["n_for_95_water"]]) for r in data], "#FF0000", ms=3.5,
         label="water (GFN2/ALPB)")
    style(ax, "", "Conformers for 95% population", (0.4, len(st) + 0.6), (1, 400), ylog=True)
    cat_axis(ax, xs, names)
    place_legend(ax, loc="lower right", fontsize=9)
    letter(ax, "b")
    ax = add_ax(fig, 10.6, 7.6, 7.2, 5.4, W, H)
    dots(ax, xs, [float(r[idx["tv_distance"]]) for r in data], "#0000FF", marker="D", ms=3.5)
    style(ax, "", "Total-variation distance", (0.4, len(st) + 0.6), (0, 0.3))
    cat_axis(ax, xs, names)
    letter(ax, "c")
    save(fig, "Fig11_vacuum_vs_water_ensembles", out)


def fig12_water_diagnosis(root, out):
    W, H = DBL_W, 140.0  # two rows: (a) on top, (b) and (c) below, aligned to (a)
    fig = new_fig(W, H)
    conf = read_csv(os.path.join(root, "vac_vs_water_decomp_conformers.csv"))
    alp = {(r["state"], r["set"], r["job"]): float(r["dG_solv_alpb_kcal"])
           for r in read_csv(os.path.join(root, "xtb_solvation_on_dft_geometries.csv"))}
    ax = add_ax(fig, 1.8, 0.6, 16.0, 5.4, W, H)
    for st in "NABD":
        rows_ = [r for r in conf if r["state"] == st and r["dG_solv_COSMO_kcal"] != ""
                 and (r["state"], r["set"], r["job"]) in alp]
        ref = [r for r in rows_ if float(r["G_aq_rel_kcal"]) == 0.0][0]
        x = [float(r["dG_solv_COSMO_kcal"]) for r in rows_]
        y = [alp[(r["state"], r["set"], r["job"])] - alp[(ref["state"], ref["set"], ref["job"])] for r in rows_]
        rr = np.corrcoef(x, y)[0, 1]
        dots(ax, x, y, COL[st], ms=3, label=f"{SNAME[st].split(' ')[0]} {st}  (r = {rr:.2f})")
    style(ax, "$\\Delta G_{solv}$ COSMO (kcal mol$^{-1}$)", "$\\Delta G_{solv}$ GFN2/ALPB (kcal mol$^{-1}$)",
          (-0.3, 3.7), (-0.5, 1.6))
    place_legend(ax, loc="upper left", fontsize=9)
    letter(ax, "a")
    # (b) monoB decomposition
    rb = [r for r in conf if r["state"] == "B" and r["dG_solv_COSMO_kcal"] != ""]
    rb.sort(key=lambda r: float(r["G_aq_rel_kcal"]))
    n = len(rb)
    xs = np.arange(1, n + 1)
    ax = add_ax(fig, 1.8, 7.6, 7.0, 5.4, W, H)
    hz(ax, 0, 0.3, n + 0.7, color="black", ls="-", lw=0.6)
    dots(ax, xs, [float(r["dE_gas_HSE06_kcal"]) for r in rb], "#999999", marker="s", ms=2.8, label="gas HSE06")
    dots(ax, xs, [float(r["dG_solv_COSMO_kcal"]) for r in rb], "#0000FF", ms=2.8, label="solvation COSMO")
    dots(ax, xs, [float(r["G_aq_rel_kcal"]) for r in rb], "black", marker="D", ms=2.8, label="total $G_{aq}$")
    new = [i for i, r in enumerate(rb, 1) if r["set"] == "new"]
    dots(ax, new, [-1.0] * len(new), "#FF0000", marker="^", ms=3.5, label="water-guided (new)")
    style(ax, "Monoanion B conformer (sorted)", "Relative to best $G_{aq}$ (kcal mol$^{-1}$)", (0.3, n + 0.7), (-1.2, 5.0))  # headroom for the legend
    place_legend(ax, loc="upper left", fontsize=9)
    letter(ax, "b")
    # (c) N->B calibrated pKa under successive definitions
    pk = read_csv(os.path.join(root, "vac_vs_water_decomp_pka.csv"))

    def pick(pair, w, g, st):
        return float([r for r in pk if r["pair"] == pair and r["weights"] == w and r["state_G"] == g and r["set"] == st][0]["calib_pKa"])

    vals = [pick("N_to_B", "pop", "mean", "old"), pick("N_to_B", "gaq", "mean", "old"),
            pick("N_to_B", "gaqg", "ens", "old"), pick("N_to_B", "gaq", "mean", "all"),
            pick("N_to_B", "gaqg", "ens", "all")]
    labs = ["prod.", "$G_{aq}$ w.", "+degen.", "+water", "+water*"]
    ax = add_ax(fig, 10.6, 7.6, 7.2, 5.4, W, H)
    line(ax, range(1, 6), vals, color="#000000", lw=1.0)
    dots(ax, range(1, 6), vals, "#000000", marker="s", ms=3.5)
    style(ax, "", "Calibrated pKa (N→B)", (0.4, 5.6), (6, 9))
    cat_axis(ax, range(1, 6), labs, rotation=30, fontsize=9)
    letter(ax, "c")
    save(fig, "Fig12_vacuum_vs_water_diagnosis", out)


def fig13_cosmo_band(root, out):
    """Vacuum vs COSMO(water) band of N and A (same HSE06/def2-TZVP conformers and weights) against the measured curve."""
    from anion_band_cosmo import cosmo_jobs
    from epsilon_abs_table import subset_spectrum, load_exp
    W, H = DBL_W, 68.0
    fig = new_fig(W, H)
    exp = load_exp(os.path.join(HERE, "..", "Figures", "data", "exp_camuri2018_fig2_spectra_vs_pH.dat"))
    g = np.linspace(HC / 400, HC / 250, 1600)   # 250-400 nm, the experimental band
    lam = HC / g
    o = np.argsort(lam)
    for k, (s_, ecol, left) in enumerate((("N", "s00", 1.8), ("A", "s23", 10.2))):
        ax = add_ax(fig, left, 1.3, 7.4, 4.6, W, H)
        jobs, pop = cosmo_jobs(root, s_)
        w = [pop[j] for j, _, _ in jobs]
        res = {}
        for medium, idx, ls_, lab in (("vacuum", 1, "--", "vacuum"), ("cosmo", 2, "-", "COSMO")):
            prof = subset_spectrum([j[idx] for j in jobs], w, g, 0.4)
            res[medium] = lam[np.argmax(prof)]
            line(ax, lam[o], (prof / prof.max())[o], color=COL[s_], lw=1.2, ls=ls_, label=f"{lab} calc.")
        e = exp[ecol]
        ok = ~np.isnan(e) & (exp["wavelength_nm"] >= 250) & (exp["wavelength_nm"] <= 400)
        line(ax, exp["wavelength_nm"][ok], e[ok] / np.nanmax(e[ok]), color=COL[s_], lw=1.0, ls=":", label="exp.")
        vt(ax, exp["wavelength_nm"][ok][np.nanargmax(e[ok])], 0, 1.0)
        style(ax, "Wavelength (nm)", "Normalised absorbance" if k == 0 else "", (250, 400), (0, 1.3))
        note(ax, f"$\\lambda_{{max}}$ exp. {exp['wavelength_nm'][ok][np.nanargmax(e[ok])]:.0f}, vac. {res['vacuum']:.0f}, "
                 f"COSMO {res['cosmo']:.0f} nm", x=1.0, y=1.02, fontsize=9, ha="right")
        ax.legend(loc="upper right", fontsize=9)  # corner above the red/black tails (checked on the PNG)
        letter(ax, "ab"[k])
    save(fig, "Fig13_cosmo_band", out)


BUILDERS = {
    "fig01": fig01_workflow, "fig02": fig02_ensembles, "fig03": fig03_functional, "fig04": fig04_b1,
    "fig05": fig05_audit, "fig06": fig06_uvvis_states, "fig07": fig07_exp_and_ph, "fig08": fig08_ir,
    "fig09": fig09_pka, "fig10": fig10_speciation, "fig11": fig11_vac_water, "fig12": fig12_water_diagnosis, "fig13": fig13_cosmo_band,
}


def main():
    ap = argparse.ArgumentParser(description="Build all Paper 2 figures (matplotlib, EPS master).")
    ap.add_argument("--project-root", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--only", nargs="*", default=None)
    a = ap.parse_args()
    root = os.path.abspath(os.path.expanduser(a.project_root))
    out = os.path.abspath(os.path.expanduser(a.outdir))
    os.makedirs(out, exist_ok=True)
    for k, fn in BUILDERS.items():
        if a.only and k not in a.only:
            continue
        try:
            fn(root, out)
        except Exception as e:  # keep building the others; report loudly
            print(f"!! {k} FAILED: {type(e).__name__}: {e}")


if __name__ == "__main__":
    main()
