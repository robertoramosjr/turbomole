#!/usr/bin/env python3
"""
plot_paper2_figures.py

Builds every Paper 2 (artepilin C, ensemble x protonation-state study) figure
as a self-contained Veusz document (.vsz, data embedded) plus .pdf/.png
exports, from the files already produced in artepillin_C_studies/.

Usage (needs the veusz_env python):
    ~/.conda/envs/veusz_env/bin/python plot_paper2_figures.py \
        --project-root ~/work_turbomole/artepillin_C_studies \
        --outdir ~/work_turbomole/artepillin_C_studies/paper_2/Figures \
        [--only fig02 fig09 ...]

Every builder only reads existing result files; nothing is recomputed except
(i) per-state UV-Vis envelopes from the stored escf.out sticks (same recipe
as predict_uvvis_vs_ph.py) and (ii) the leave-one-out "mean of the others"
control used in the statistical-audit figure.
"""

import argparse
import csv
import math
import os
import re
import shutil

import numpy as np

if not os.environ.get("DISPLAY"):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import veusz.embed as veusz  # noqa: E402

STATES = ["N", "A", "B", "D"]
SDIR = {"N": "neutral", "A": "monoanion_a", "B": "monoanion_b", "D": "dianion"}
SNAME = {"N": "Neutral (0)", "A": "Monoanion A (carboxylate)", "B": "Monoanion B (phenolate)",
         "D": "Dianion (-2)"}
# Okabe-Ito colour-blind-safe palette, one colour per microstate everywhere
COL = {"N": "#0072B2", "A": "#D55E00", "B": "#009E73", "D": "#CC79A7"}
GRAY = "#666666"
HC = 1239.84198  # eV*nm


# ----------------------------------------------------------------------------- helpers
class Fig:
    """One Veusz document (one page) with panel() graphs placed in cm."""

    def __init__(self, name, width_cm, height_cm):
        self.name, self.W, self.H = name, width_cm, height_cm
        self.g = veusz.Embedded(name, hidden=True)
        self.g.To(self.g.Add("page", name="page1"))
        self.g.Set("width", f"{width_cm}cm")
        self.g.Set("height", f"{height_cm}cm")
        self.n = 0
        self.key_done = set()

    def panel(self, letter, l, t, w, h, xlabel="", ylabel="", xmin=None, xmax=None, ymin=None,
              ymax=None, xlog=False, ylog=False, noxticks=False):
        g = self.g
        g.To("/page1")
        pname = f"p{letter}"
        g.To(g.Add("graph", name=pname))
        g.Set("leftMargin", f"{l}cm")
        g.Set("topMargin", f"{t}cm")
        g.Set("rightMargin", f"{self.W - l - w}cm")
        g.Set("bottomMargin", f"{self.H - t - h}cm")
        g.Set("x/label", xlabel)
        g.Set("y/label", ylabel)
        for ax in ("x", "y"):
            g.Set(f"{ax}/Label/size", "9pt")
            g.Set(f"{ax}/TickLabels/size", "8pt")
            g.Set(f"{ax}/MajorTicks/number", 6)
        g.Set("x/log", xlog)
        g.Set("y/log", ylog)
        for ax, lo, hi in (("x", xmin, xmax), ("y", ymin, ymax)):
            if lo is not None:
                g.Set(f"{ax}/min", float(lo))
            if hi is not None:
                g.Set(f"{ax}/max", float(hi))
        if noxticks:
            g.Set("x/TickLabels/hide", True)
            g.Set("x/MajorTicks/hide", True)
            g.Set("x/MinorTicks/hide", True)
        self.cur = pname
        self.text(letter and f"({letter})", 0.02, 0.93, size="11pt", bold=True)
        return pname

    def to(self, pname):
        self.g.To(f"/page1/{pname}")
        self.cur = pname

    def xy(self, x, y, color="black", width="1.2pt", marker="none", size="2.5pt", line=True,
           key=None, style="solid", fill=True):
        g = self.g
        self.n += 1
        nm = f"d{self.n}"
        g.SetData(f"{self.cur}_{nm}x", [float(v) for v in x])
        g.SetData(f"{self.cur}_{nm}y", [float(v) for v in y])
        g.Add("xy", name=nm)
        g.Set(f"{nm}/xData", f"{self.cur}_{nm}x")
        g.Set(f"{nm}/yData", f"{self.cur}_{nm}y")
        g.Set(f"{nm}/marker", marker)
        g.Set(f"{nm}/markerSize", size)
        g.Set(f"{nm}/MarkerFill/color", color)
        g.Set(f"{nm}/MarkerFill/hide", not fill)
        g.Set(f"{nm}/MarkerLine/color", color)
        g.Set(f"{nm}/PlotLine/hide", not line)
        g.Set(f"{nm}/PlotLine/color", color)
        g.Set(f"{nm}/PlotLine/width", width)
        g.Set(f"{nm}/PlotLine/style", style)
        if key:
            g.Set(f"{nm}/key", key)

    def err(self, x, y, yerr, color="black", marker="circle", size="3pt", key=None):
        g = self.g
        self.n += 1
        nm = f"d{self.n}"
        g.SetData(f"{self.cur}_{nm}x", [float(v) for v in x])
        g.SetData(f"{self.cur}_{nm}y", [float(v) for v in y], symerr=[float(v) for v in yerr])
        g.Add("xy", name=nm)
        g.Set(f"{nm}/xData", f"{self.cur}_{nm}x")
        g.Set(f"{nm}/yData", f"{self.cur}_{nm}y")
        g.Set(f"{nm}/marker", marker)
        g.Set(f"{nm}/markerSize", size)
        g.Set(f"{nm}/MarkerFill/color", color)
        g.Set(f"{nm}/MarkerLine/color", color)
        g.Set(f"{nm}/PlotLine/hide", True)
        g.Set(f"{nm}/ErrorBarLine/color", color)
        g.Set(f"{nm}/errorStyle", "bar")
        if key:
            g.Set(f"{nm}/key", key)

    def vline(self, x, y0, y1, color=GRAY, style="dashed"):
        self.xy([x, x], [y0, y1], color=color, width="0.8pt", style=style)

    def hline(self, y, x0, x1, color=GRAY, style="dashed"):
        self.xy([x0, x1], [y, y], color=color, width="0.8pt", style=style)

    def text(self, s, x, y, size="8pt", color="black", bold=False, align="left", angle=0):
        g = self.g
        self.n += 1
        nm = f"t{self.n}"
        g.Add("label", name=nm)
        g.Set(f"{nm}/label", s.replace("\n", "\\\\"))
        g.Set(f"{nm}/xPos", float(x))
        g.Set(f"{nm}/yPos", float(y))
        g.Set(f"{nm}/positioning", "relative")
        g.Set(f"{nm}/alignHorz", align)
        g.Set(f"{nm}/angle", float(angle))
        g.Set(f"{nm}/Text/size", size)
        g.Set(f"{nm}/Text/color", color)
        g.Set(f"{nm}/Text/bold", bold)

    def key(self, hpos="right", vpos="top"):
        g = self.g
        g.Add("key", name="key1")
        g.Set("key1/horzPosn", hpos)
        g.Set("key1/vertPosn", vpos)
        g.Set("key1/Border/hide", True)
        g.Set("key1/Text/size", "7pt")
        g.Set("key1/keyLength", "0.6cm")

    def xticks(self, positions, labels, y=-0.07, size="8pt", angle=0):
        """Category labels under a numeric x axis."""
        lo, hi = self._xrange
        self.g.Set("x/MajorTicks/manualTicks", [float(p) for p in positions])
        self.g.Set("x/TickLabels/hide", True)
        for p, lab in zip(positions, labels):
            self.text(lab, (p - lo) / (hi - lo), y, size=size, align="centre", angle=angle)

    def save(self, outdir):
        base = os.path.join(outdir, self.name)
        self.g.Save(base + ".vsz")
        self.g.Export(base + ".pdf")
        self.g.Export(base + ".png", dpi=200)
        print("saved", base)
        self.g.Close()


def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return float("nan")


def sdir(root, s):
    return os.path.join(root, f"artepillin_c_{SDIR[s]}_d3_hse06")


def whitespace_table(path):
    """Veusz-style .dat with a 'descriptor' header line -> dict of numpy columns."""
    with open(path) as f:
        head = f.readline().strip().replace("descriptor", "").replace(",", " ").split()
        rows = [[float(v) for v in ln.split()] for ln in f if ln.strip() and not ln.startswith("#")]
    a = np.array(rows)
    return {h: a[:, i] for i, h in enumerate(head)}


# ----------------------------------------------------------------------------- fig builders
def fig01_workflow(root, out):
    """Workflow schematic (boxes + arrows); numbers are the real counts."""
    F = Fig("Fig01_workflow", 18.0, 11.0)
    g = F.g
    g.To("/page1")
    boxes = [
        # (id, col, row, title, subtitle, colour)
        ("b1", 0, 0, "CREST iMTD-GC / GFN2-xTB", "4 microstates\n337 | 322 | 37 | 21 conformers", "#d9e8f5"),
        ("b2", 1, 0, "Boltzmann filter (95%, 298 K)", "weight x CREST degeneracy\nsurvivors: 170 | 143 | 11 | 3", "#d9e8f5"),
        ("b3", 2, 0, "r2SCAN-3c re-optimisation", "327 conformers; re-ranked\n(top conformer changes in 3/4 states)", "#d9e8f5"),
        ("b4", 2, 1, "Layer A: PBE0/def2-SVP IR, UV-Vis", "population-weighted\nscreening spectra (327)", "#e6f2e6"),
        ("b5", 1, 1, "Layer B: HSE06-D3(BJ)/def2-TZVP", "54 labelled + extensions\n(same level as Paper 1)", "#e6f2e6"),
        ("b6", 0, 1, "Delta-learning surrogates", "B1 nu(C=O): GFN2-xTB -> HSE06\nB2 S1 (E, f): sTDA -> HSE06", "#e6f2e6"),
        ("b7", 0, 2, "Statistical audit", "LOO Q^{2} + permutation test\nnegative results kept", "#fbeadb"),
        ("b8", 1, 2, "Thermodynamic cycle -> pKa", "freeh + COSMO(water)\nclass-wise calibration vs references", "#fbeadb"),
        ("b9", 2, 2, "Outputs", "ensemble IR/UV-Vis per microstate\nspeciation, spectra vs pH", "#f3e0ef"),
    ]
    bw, bh = 0.29, 0.20
    xs = [0.04, 0.355, 0.67]
    ys = [0.74, 0.43, 0.12]
    pos = {}
    for bid, c, r, title, sub, colr in boxes:
        x0, y0 = xs[c], ys[r]
        pos[bid] = (x0, y0)
        g.Add("rect", name=bid)
        g.Set(f"{bid}/xPos", [x0 + bw / 2])
        g.Set(f"{bid}/yPos", [y0 + bh / 2])
        g.Set(f"{bid}/width", [bw])
        g.Set(f"{bid}/height", [bh])
        g.Set(f"{bid}/Fill/color", colr)
        g.Set(f"{bid}/Border/width", "0.8pt")
        g.Set(f"{bid}/Border/color", "#444444")
        g.Add("label", name=bid + "_t")
        g.Set(f"{bid}_t/label", title)
        g.Set(f"{bid}_t/xPos", [x0 + bw / 2])
        g.Set(f"{bid}_t/yPos", [y0 + bh * 0.62])
        g.Set(f"{bid}_t/alignHorz", "centre")
        g.Set(f"{bid}_t/Text/size", "8pt")
        g.Set(f"{bid}_t/Text/bold", True)
        g.Add("label", name=bid + "_s")
        g.Set(f"{bid}_s/label", sub.replace("\n", "\\\\"))
        g.Set(f"{bid}_s/xPos", [x0 + bw / 2])
        g.Set(f"{bid}_s/yPos", [y0 + bh * 0.12])
        g.Set(f"{bid}_s/alignHorz", "centre")
        g.Set(f"{bid}_s/Text/size", "7pt")
    arrows = [("b1", "b2"), ("b2", "b3"), ("b3", "b4"), ("b4", "b5"), ("b5", "b6"), ("b6", "b7"),
              ("b7", "b8"), ("b8", "b9")]
    for i, (a, b) in enumerate(arrows):
        xa, ya = pos[a]
        xb, yb = pos[b]
        # same row -> horizontal arrow; row change -> from bottom of a to top of b
        if abs(ya - yb) < 1e-6:
            xs0, ys0, xs1, ys1 = (xa + bw, ya + bh / 2, xb, yb + bh / 2) if xb > xa else (xa, ya + bh / 2, xb + bw, yb + bh / 2)
        else:
            xs0, ys0, xs1, ys1 = xa + bw / 2, ya, xb + bw / 2, yb + bh
        nm = f"arr{i}"
        g.Add("line", name=nm)
        g.Set(f"{nm}/mode", "point-to-point")
        g.Set(f"{nm}/xPos", [xs0])
        g.Set(f"{nm}/yPos", [ys0])
        g.Set(f"{nm}/xPos2", [xs1])
        g.Set(f"{nm}/yPos2", [ys1])
        g.Set(f"{nm}/arrowright", "arrow")
        g.Set(f"{nm}/arrowSize", "7pt")
        g.Set(f"{nm}/Line/width", "1.2pt")
    g.Add("label", name="title")
    g.Set("title/label", "Funnel: ensemble x protonation state -> spectra and pKa (artepilin C, C19H24O3)")
    g.Set("title/xPos", [0.5])
    g.Set("title/yPos", [0.96])
    g.Set("title/alignHorz", "centre")
    g.Set("title/Text/size", "10pt")
    g.Set("title/Text/bold", True)
    F.save(out)


def fig02_ensembles(root, out):
    F = Fig("Fig02_ensembles", 18.0, 7.0)
    # (a) cumulative population vs rank (r2SCAN-3c, degeneracy-weighted)
    F.panel("a", 1.6, 0.6, 5.0, 5.2, "Conformer rank", "Cumulative population", 0.8, 200, 0, 1.02, xlog=True)
    kept = {}
    for s in STATES:
        rows = read_csv(os.path.join(sdir(root, s), "layer_a", "population_table_dft.csv"))
        p = np.array(sorted([float(r["r2scan3c_population"]) for r in rows], reverse=True))
        p = p / p.sum()
        F.xy(np.arange(1, len(p) + 1), np.cumsum(p), color=COL[s], width="1.4pt", key=SNAME[s])
        kept[s] = len(p)
    F.key("right", "bottom")
    # (b) GFN2 vs r2SCAN-3c relative energy
    F.panel("b", 8.2, 0.6, 4.2, 5.2, "GFN2-xTB rel. E (kcal/mol)", "r2SCAN-3c rel. E (kcal/mol)", 0, 1.1, 0, 3.0)
    for s in ("N", "A"):
        rows = read_csv(os.path.join(sdir(root, s), "layer_a", "population_table_dft.csv"))
        x = [float(r["gfn2_relative_energy_kcal"]) for r in rows]
        y = [float(r["r2scan3c_relative_energy_kcal"]) for r in rows]
        F.xy(x, y, color=COL[s], marker="circle", size="2pt", line=False, key=SNAME[s])
    F.key("left", "top")
    allx, ally = [], []
    for s_ in ("N", "A"):
        rr = read_csv(os.path.join(sdir(root, s_), "layer_a", "population_table_dft.csv"))
        allx += [float(r["gfn2_relative_energy_kcal"]) for r in rr]
        ally += [float(r["r2scan3c_relative_energy_kcal"]) for r in rr]
    rx = np.argsort(np.argsort(allx)); ry = np.argsort(np.argsort(ally))
    F.text(f"Spearman rho = {np.corrcoef(rx, ry)[0, 1]:.2f} (N + A pooled)", 0.04, 0.04, size="7pt")
    # (c) conformer counts
    F.panel("c", 14.0, 0.6, 3.6, 5.2, "", "Number of conformers", 0.4, 4.6, 0, 400)
    crest = {"N": 337, "A": 322, "B": 37, "D": 21}
    lab = {"N": 20, "A": 20, "B": 11, "D": 3}
    xs = [1, 2, 3, 4]
    F._xrange = (0.4, 4.6)
    F.xy(xs, [crest[s] for s in STATES], color="#999999", marker="square", size="3pt", line=False, key="CREST")
    F.xy(xs, [kept[s] for s in STATES], color="#222222", marker="circle", size="3pt", line=False, key="95% filter")
    F.xy(xs, [lab[s] for s in STATES], color="#D55E00", marker="diamond", size="3pt", line=False, key="HSE06 labelled")
    F.xticks(xs, ["N", "A", "B", "D"])
    F.key("right", "top")
    F.save(out)


def fig03_functional(root, out):
    F = Fig("Fig03_functional_and_sTDA", 18.0, 7.0)
    fc = read_csv(os.path.join(root, "functional_comparison.csv"))
    names = [r["functional"].upper() for r in fc]
    e = [float(r["energy_eV"]) for r in fc]
    f = [float(r["osc_len"]) for r in fc]
    F.panel("a", 1.6, 0.6, 4.6, 5.2, "", "S1 vertical energy (eV)", 0.4, 3.6, 3.6, 4.4)
    xs = [1, 2, 3]
    F._xrange = (0.4, 3.6)
    F.xy(xs, e, color="#0072B2", marker="square", size="4pt", line=False)
    F.xticks(xs, [n.replace("CAM-B3LYP", "CAM-B3LYP") for n in names], size="7pt")
    F.panel("b", 7.6, 0.6, 4.4, 5.2, "", "S1 oscillator strength f (length)", 0.4, 3.6, 0.4, 0.9)
    F._xrange = (0.4, 3.6)
    F.xy(xs, f, color="#D55E00", marker="circle", size="4pt", line=False)
    F.xticks(xs, names, size="7pt")
    F.text("f discrepancy across functionals = CT proxy (Paper 1)", 0.04, 0.05, size="7pt")
    # (c) sTDA vs HSE06 S1 energy, neutral, character matched
    nd = os.path.join(sdir(root, "N"), "layer_b")
    st = {r["job_dir"]: r for r in read_csv(os.path.join(nd, "s1_stda_character_matched_FIXED.csv"))}
    hs = {r["job_dir"]: r for r in read_csv(os.path.join(nd, "s1_hse06.csv"))}
    common = [j for j in hs if j in st]
    xe = [float(st[j]["s1_energy_eV_stda"]) for j in common]
    ye = [float(hs[j]["s1_energy_eV"]) for j in common]
    xf = [float(st[j]["s1_osc_stda"]) for j in common]
    yf = [float(hs[j]["s1_osc_len"]) for j in common]
    F.panel("c", 13.6, 0.6, 4.0, 5.2, "sTDA S1 (eV)", "HSE06 S1 (eV)", 2.6, 3.0, 3.8, 4.1)
    F.xy(xe, ye, color="#0072B2", marker="circle", size="2pt", line=False)
    F.text("Neutral, character-matched: sTDA is ~1.2 eV low", 0.04, 0.05, size="7pt")
    F.save(out)
    # companion: f scatter as its own figure (log-log), used for the B2 delta-learning claim
    G = Fig("Fig03b_sTDA_vs_HSE06_f", 9.0, 7.0)
    G.panel("a", 1.8, 0.6, 6.6, 5.4, "sTDA S1 f", "HSE06 S1 f", 0.05, 1.5, 0.05, 1.5, xlog=True, ylog=True)
    G.xy(xf, yf, color="#0072B2", marker="circle", size="2pt", line=False, key="neutral conformers shown; fit uses N=28")
    xx = np.logspace(-1.3, 0.18, 20)
    G.xy(xx, 10 ** (-0.2027 + 0.5794 * np.log10(xx)), color="#D55E00", width="1.4pt",
         key="log10 f(HSE06) = -0.2027 + 0.5794 log10 f(sTDA)")
    G.xy(xx, xx, color=GRAY, width="0.8pt", style="dashed", key="y = x")
    G.key("left", "top")
    G.save(out)


def fig04_b1(root, out):
    F = Fig("Fig04_B1_carbonyl_shift", 18.0, 12.0)
    td = {s: read_csv(os.path.join(sdir(root, s), "layer_b", "b1_training_data.csv")) for s in STATES}
    # (a) strip plot of delta per state
    F.panel("a", 1.8, 0.6, 5.0, 4.6, "", "HSE06 - GFN2-xTB  nu(C=O) (cm-1)", 0.4, 4.6, -100, 60)
    F._xrange = (0.4, 4.6)
    rng = np.random.default_rng(3)
    for i, s in enumerate(STATES, 1):
        d = [float(r["delta_cm1"]) for r in td[s]]
        F.xy(i + rng.uniform(-0.15, 0.15, len(d)), d, color=COL[s], marker="circle", size="2pt", line=False)
    F.xticks([1, 2, 3, 4], ["N", "A", "B", "D"])
    # (b) monoA: delta vs Bader charge of carbonyl O
    A = sdir(root, "A")
    feat = {r["job_dir"]: float(r["bader_charge_carbonyl_o"]) for r in read_csv(os.path.join(A, "layer_b", "population_feature.csv"))}
    xa = [feat[r["job_dir"]] for r in td["A"] if r["job_dir"] in feat]
    ya = [float(r["delta_cm1"]) for r in td["A"] if r["job_dir"] in feat]
    F.panel("b", 8.2, 0.6, 4.6, 4.6, "Bader q(O carbonyl) (e)", "Delta nu(C=O) (cm-1)", -1.65, -1.1, -82, -78)
    F.xy(xa, ya, color=COL["A"], marker="circle", size="2pt", line=False, key="labelled (N=20)")
    xx = np.linspace(-1.6034, -1.1553, 10)
    F.xy(xx, 4.7223 * xx - 73.9249, color="black", width="1.2pt", key="OLS, Q^{2} (LOO) = 0.277, p = 0.017")
    F.key("left", "top")
    # (c) monoB: delta vs O...H distance
    B = sdir(root, "B")
    oh = {r["job_dir"]: float(r["oh_distance_angstrom"]) for r in read_csv(os.path.join(B, "layer_b", "oh_distance_feature.csv"))}
    xb = [oh[r["job_dir"]] for r in td["B"] if r["job_dir"] in oh]
    yb = [float(r["delta_cm1"]) for r in td["B"] if r["job_dir"] in oh]
    F.panel("c", 14.0, 0.6, 3.7, 4.6, "O...H distance (A)", "Delta nu(C=O) (cm-1)", 2.2605, 2.2635, -0.5, 0.5)
    ymin, ymax = min(yb), max(yb)
    pad = 0.1 * (ymax - ymin + 1e-6)
    F.g.Set("y/min", float(ymin - pad))
    F.g.Set("y/max", float(ymax + pad))
    F.xy(xb, yb, color=COL["B"], marker="circle", size="2pt", line=False, key="N=11")
    xx = np.linspace(min(xb), max(xb), 5)
    F.xy(xx, 1392.13 * xx - 3114.30, color="black", width="1.2pt", key="OLS, Q^{2} (LOO) = 0.594")
    F.key("left", "top")
    # (d) neutral: LOO prediction vs truth (flat offset)
    N = sdir(root, "N")
    loo = read_csv(os.path.join(N, "layer_b", "b1_loo_predictions.csv"))
    t = [float(r["true_delta_cm1"]) for r in loo]
    p = [float(r["loo_predicted_delta_cm1"]) for r in loo]
    F.panel("d", 1.8, 6.8, 5.0, 4.6, "True delta (cm-1)", "LOO predicted delta (cm-1)", 27, 30, 27, 30)
    F.xy([27, 30], [27, 30], color=GRAY, width="0.8pt", style="dashed")
    F.xy(t, p, color=COL["N"], marker="circle", size="2pt", line=False)
    F.text("Neutral: no feature beats a flat offset\n28.222 +/- 0.676 cm-1", 0.05, 0.12, size="7pt")
    # (e) end result: HSE06 harmonic nu(C=O) (carboxylate: C-O stretch) per microstate, unscaled
    F.panel("e", 8.2, 6.8, 9.5, 4.6, "", "HSE06 harmonic nu(C=O) (cm-1)", 0.4, 4.6, 1550, 1850)
    F._xrange = (0.4, 4.6)
    means = [np.mean([float(r["nu_CO_hse06_cm1"]) for r in td[s]]) for s in STATES]
    sds = [np.std([float(r["nu_CO_hse06_cm1"]) for r in td[s]], ddof=1) for s in STATES]
    for i_, s in enumerate(STATES, 1):
        F.err([i_], [means[i_ - 1]], [sds[i_ - 1]], color=COL[s], marker="square")
    F.xticks([1, 2, 3, 4], ["N", "A", "B", "D (N=3)"], size="7pt")
    F.text("unscaled harmonic; the C=O mode of A and D is the carboxylate stretch", 0.04, 0.05, size="6pt")
    F.save(out)


def loo_mean_control(y):
    """Leave-one-out prediction = mean of the other points (Spearman is -1 by construction)."""
    y = np.asarray(y, float)
    return np.array([(y.sum() - v) / (len(y) - 1) for v in y])


def fig05_audit(root, out):
    F = Fig("Fig05_statistical_audit", 18.0, 7.5)
    rows = []
    for s in STATES:
        L = os.path.join(sdir(root, s), "layer_b")
        for fn, label in (("b1_permutation_test.dat", "B1"), ("b2_permutation_test_delta_energy_eV.dat", "B2 dE"),
                          ("b2_permutation_test_delta_log10_f.dat", "B2 dlogf")):
            p = os.path.join(L, fn)
            if os.path.isfile(p):
                lines = [ln for ln in open(p).read().splitlines() if ln.strip()]
                v = lines[1].split()
                rows.append((s, label, int(v[1]), float(v[4]), float(v[6]), float(v[7]), float(v[9]), float(v[11])))
    n = len(rows)
    F.panel("a", 1.8, 0.8, 9.8, 5.6, "", "Leave-one-out Q^{2}", 0.3, n + 0.7, -0.5, 0.75)
    F._xrange = (0.3, n + 0.7)
    xs = list(range(1, n + 1))
    F.hline(0, 0.3, n + 0.7)
    F.err(xs, [r[4] for r in rows], [r[5] for r in rows], color="#999999", marker="square", size="2pt",
          key="permutation null (mean +/- sd)")
    F.xy(xs, [r[6] for r in rows], color="#999999", marker="linehorz", size="4pt", line=False, key="null 95th pct")
    F.xy(xs, [r[3] for r in rows], color="#D55E00", marker="diamond", size="4pt", line=False, key="observed best-feature Q^{2}")
    F.xticks(xs, [f"{r[0]} {r[1]}" for r in rows], size="7pt", y=-0.1, angle=30)
    for x, r in zip(xs, rows):
        F.text(f"p={r[7]:.3f}", (x - 0.3) / (n + 0.4), 0.04, size="7pt", align="centre")
    F.key("left", "top")
    # (b) LOO mean-of-the-others control on the neutral B1 labels -> rho = -1 exactly
    y = [float(r["delta_cm1"]) for r in read_csv(os.path.join(sdir(root, "N"), "layer_b", "b1_training_data.csv"))]
    ctrl = loo_mean_control(y)
    lo, hi = min(y) - 0.1, max(y) + 0.1
    F.panel("b", 13.0, 0.8, 4.6, 5.6, "True delta (cm-1)", "LOO 'mean of others' (cm-1)", lo, hi,
            min(ctrl) - 0.02, max(ctrl) + 0.02)
    F.xy(y, ctrl, color="#0072B2", marker="circle", size="2pt", line=False)
    F.text("Spearman = -1 for ANY data set\n(algebraic artefact of LOO)", 0.04, 0.08, size="7pt")
    F.save(out)


def sticks(path):
    rx = re.compile(r"^\s*\|\s*(\d+)\s+a\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|"
                    r"\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|\s*([-\d.]+)\s*\|")
    out = []
    with open(path) as f:
        for ln in f:
            m = rx.match(ln)
            if m:
                out.append((float(m.group(3)), float(m.group(7))))
    return out


def state_spectrum(root, s, grid, sigma=0.4):
    d = sdir(root, s)
    pop = {r["job_dir"]: float(r["r2scan3c_population"]) for r in read_csv(os.path.join(d, "layer_a", "population_table_dft.csv"))}
    found = {}
    for sub in ("labeled_jobs", "confirmation_batch_jobs", "random_control_batch_jobs"):
        base = os.path.join(d, "layer_b", sub)
        if os.path.isdir(base):
            for j in os.listdir(base):
                p = os.path.join(base, j, "escf.out")
                if os.path.isfile(p) and j not in found and j in pop:
                    found[j] = p
    tot = sum(pop[j] for j in found)
    prof = np.zeros_like(grid)
    for j, p in found.items():
        for e, f in sticks(p):
            prof += pop[j] / tot * f * np.exp(-((grid - e) ** 2) / (2 * sigma ** 2))
    return prof, len(found), len(pop), tot / sum(pop.values())


def fig06_uvvis_states(root, out):
    F = Fig("Fig06_uvvis_per_microstate", 18.0, 7.0)
    grid = np.linspace(2.0, 6.0, 800)
    F.panel("a", 1.8, 0.6, 8.6, 5.4, "Wavelength (nm)", "Population-weighted absorption (a.u.)", 220, 480, 0, 1.25)
    lam = HC / grid
    for s in STATES:
        prof, n, ntot, cov = state_spectrum(root, s, grid)
        F.xy(lam, prof, color=COL[s], width="1.4pt", key=f"{SNAME[s]} ({n}/{ntot} conf., {cov * 100:.0f}% pop.)")
    F.key("right", "top")
    cons = read_csv(os.path.join(root, "b2_consolidated_dataset.csv"))
    F.panel("b", 12.0, 0.6, 5.6, 5.4, "S1 energy (eV)", "S1 oscillator strength f", 2.2, 4.6, 0, 1.0)
    names = {"neutral": "N", "monoanion_a": "A", "monoanion_b": "B", "dianion": "D"}
    for st, s in names.items():
        for src, mk in (("hse06_ground_truth", "circle"), ("surrogate_max_f", "triangle")):
            sel = [r for r in cons if r["protonation_state"] == st and r["source"] == src]
            if sel:
                F.xy([float(r["energy_eV"]) for r in sel], [float(r["osc_f"]) for r in sel], color=COL[s],
                     marker=mk, size="2pt", line=False, key=f"{s} {'HSE06' if mk == 'circle' else 'surrogate'}")
    F.key("left", "top")
    F.save(out)


def fig07_exp_and_ph(root, out):
    F = Fig("Fig07_uvvis_vs_experiment_and_pH", 18.0, 7.0)
    d = whitespace_table(os.path.join(root, "uvvis_conf1_vs_exp.dat"))
    wl = d["wavelength_nm"]
    F.panel("a", 1.8, 0.6, 7.6, 5.4, "Wavelength (nm)", "Normalised absorbance", 250, 400, 0, 1.1)
    F.xy(wl, d["calc_conf1_norm"], color=GRAY, width="1.0pt", style="dashed", key="calc., conformer 1")
    F.xy(wl, d["calc_ensemble_norm"], color="#0072B2", width="1.4pt", key="calc., ensemble")
    ok = ~np.isnan(d["exp_norm"])
    F.xy(wl[ok], d["exp_norm"][ok], color="black", width="1.4pt", key="experiment, pH 2.6 neutral (Camuri 2018)")
    F.key("right", "top")
    # (b) band maximum and centroid vs pH from per-state spectra mixed by speciation
    sp = whitespace_table(os.path.join(root, "speciation_v4_5refs.dat"))
    grid = np.linspace(2.0, 6.0, 1200)
    spec = {s: state_spectrum(root, s, grid)[0] for s in STATES}
    w = (grid >= 3.0) & (grid <= 5.5)
    phs = np.arange(1.0, 13.01, 0.25)
    mx, cen = [], []
    for ph in phs:
        i = int(np.argmin(np.abs(sp["pH"] - ph)))
        prof = sum(sp[f"f_{s}"][i] * spec[s] for s in STATES)
        mx.append(HC / grid[w][np.argmax(prof[w])])
        cen.append(HC / (np.sum(grid[w] * prof[w]) / np.sum(prof[w])))
    F.panel("b", 11.6, 0.6, 6.0, 5.4, "pH", "Predicted band position (nm)", 1, 13, 250, 330)
    F.xy(phs, mx, color="#D55E00", width="1.4pt", key="band maximum")
    F.xy(phs, cen, color="#0072B2", width="1.4pt", key="f-weighted centroid")
    F.vline(4.65, 250, 330)
    F.text("exp. pKa 4.65", 0.30, 0.05, size="7pt")
    F.key("right", "top")
    F.save(out)


def fig08_ir(root, out):
    F = Fig("Fig08_ir_per_microstate", 18.0, 7.0)
    F.panel("a", 1.8, 0.6, 9.0, 5.4, "Wavenumber (cm-1)", "Population-weighted IR (norm.)", 400, 4000, 0, 1.05)
    F2 = None
    data = {}
    for s in STATES:
        d = whitespace_table(os.path.join(sdir(root, s), "layer_a", "ir_population_weighted.dat"))
        k = list(d.keys())
        x, y = d[k[0]], d[k[1]]
        data[s] = (x, y / max(y.max(), 1e-12))
        F.xy(x, data[s][1], color=COL[s], width="1.0pt", key=SNAME[s])
    F.key("left", "top")
    F.panel("b", 12.2, 0.6, 5.4, 5.4, "Wavenumber (cm-1)", "IR (norm.)", 1550, 1850, 0, 1.05)
    for s in STATES:
        x, y = data[s]
        m = (x >= 1550) & (x <= 1850)
        F.xy(x[m], y[m] / max(y[m].max(), 1e-12), color=COL[s], width="1.2pt")
    F.save(out)


def fig09_pka(root, out):
    F = Fig("Fig09_pKa_calibration", 18.0, 7.0)
    refs, order = [], []
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
    F.panel("a", 1.8, 1.2, 6.6, 4.6, "", "pKa residual (exp - raw)", 0.4, n + 0.6, -5.2, -2.2)
    F._xrange = (0.4, n + 0.6)
    for i, (nm, cl, res) in enumerate(refs, 1):
        F.xy([i], [res], color="#0072B2" if cl == "0to-1" else "#D55E00", marker="circle", size="4pt", line=False)
    short = {"acetic_acid": "acetic", "phenol": "phenol", "4-hydroxybenzoate_phenol": "4-HBA",
             "26-dimethylphenol": "2,6-DMP", "4-coumarate_phenol": "4-coum."}
    F.xticks(list(range(1, n + 1)), [short.get(r[0], r[0]) for r in refs], size="7pt", y=-0.1, angle=0)
    c1 = np.mean([r[2] for r in refs if r[1] == "0to-1"])
    c2 = np.mean([r[2] for r in refs if r[1] == "-1to-2"])
    F.hline(c1, 0.4, n + 0.6, color="#0072B2")
    F.hline(c2, 0.4, n + 0.6, color="#D55E00")
    F.text("blue: 0 -> -1 (N=3);  orange: -1 -> -2 (N=2)", 0.04, 0.05, size="7pt")
    # (b) calibrated artepilin C pKa: vacuum ensemble vs water-extended
    vw = {r["pair"]: r for r in read_csv(os.path.join(root, "pka_vacuum_vs_water.csv"))}
    pairs = ["N_to_A", "N_to_B", "A_to_D", "B_to_D"]
    F.panel("b", 10.4, 1.2, 7.2, 4.6, "", "Calibrated pKa", 0.4, 4.6, 3, 14)
    F._xrange = (0.4, 4.6)
    xs = [1, 2, 3, 4]
    F.err([x - 0.12 for x in xs], [float(vw[p]["calib_A"]) for p in pairs], [0.41, 0.41, 1.15, 1.15], color="#0072B2", marker="square", size="4pt", key="vacuum-CREST ensemble (+/- class scatter)")
    F.xy([x + 0.12 for x in xs], [float(vw[p]["calib_C"]) for p in pairs], color="#D55E00", marker="circle", size="3pt", line=False, key="+ water-guided conformers, G(aq) weights")
    F.xy([1], [4.65], color="black", marker="star", size="4pt", line=False, key="experiment (carboxyl)")
    F.xticks(xs, ["N->A", "N->B", "A->D", "B->D"], size="7pt", y=-0.09)
    F.key("right", "bottom")
    F.text("-1 -> -2 pairs: low confidence (N=2, sd 1.15)", 0.04, 0.05, size="7pt")
    F.save(out)


def fig10_speciation(root, out):
    F = Fig("Fig10_speciation_vs_pH", 18.0, 7.0)
    for letter, fn, title, l in (("a", "speciation_v4_5refs.dat", "vacuum-CREST ensemble", 1.8),
                                 ("b", "speciation_v4_5refs_water_ext.dat", "+ water-guided conformers", 10.2)):
        sp = whitespace_table(os.path.join(root, fn))
        F.panel(letter, l, 0.6, 7.2, 5.4, "pH", "Mole fraction", 0, 14, 0, 1.02)
        for s in STATES:
            F.xy(sp["pH"], sp[f"f_{s}"], color=COL[s], width="1.6pt", key=SNAME[s])
        F.vline(4.65, 0, 1.02)
        F.text(title, 0.40, 0.12, size="7pt")
        if letter == "a":
            F.key("right", "centre")
    F.save(out)


def fig11_vac_water(root, out):
    F = Fig("Fig11_vacuum_vs_water_ensembles", 18.0, 7.0)
    rows = read_csv(os.path.join(root, "vacuum_vs_water_ensembles.dat")) if False else None
    d = [ln.strip().split(",") for ln in open(os.path.join(root, "vacuum_vs_water_ensembles.dat")) if ln.strip()]
    hdr, data = d[0], d[1:]
    idx = {h: i for i, h in enumerate(hdr)}
    st = [r[idx["state"]] for r in data]
    lab = {"dianion": "D", "monoanion_b": "B", "monoanion_a": "A"}
    xs = list(range(1, len(st) + 1))
    F.panel("a", 1.8, 0.6, 5.4, 5.4, "", "Water population retained by vacuum 95% set", 0.4, len(st) + 0.6, 0.6, 1.0)
    F._xrange = (0.4, len(st) + 0.6)
    F.xy(xs, [float(r[idx["water_pop_in_kept"]]) for r in data], color="#0072B2", marker="square", size="4pt", line=False)
    F.hline(0.95, 0.4, len(st) + 0.6)
    F.xticks(xs, [lab[s] for s in st])
    F.panel("b", 9.0, 0.6, 4.2, 5.4, "", "Conformers for 95% population", 0.4, len(st) + 0.6, 1, 400, ylog=True)
    F._xrange = (0.4, len(st) + 0.6)
    F.xy([x - 0.12 for x in xs], [float(r[idx["n_for_95_vac"]]) for r in data], color="#999999", marker="square", size="4pt", line=False, key="vacuum")
    F.xy([x + 0.12 for x in xs], [float(r[idx["n_for_95_water"]]) for r in data], color="#D55E00", marker="circle", size="4pt", line=False, key="water (GFN2/ALPB)")
    F.xticks(xs, [lab[s] for s in st])
    F.key("left", "top")
    F.panel("c", 14.6, 0.6, 3.0, 5.4, "", "Total-variation distance", 0.4, len(st) + 0.6, 0, 0.3)
    F._xrange = (0.4, len(st) + 0.6)
    F.xy(xs, [float(r[idx["tv_distance"]]) for r in data], color="#009E73", marker="diamond", size="4pt", line=False)
    F.xticks(xs, [lab[s] for s in st])
    F.save(out)


def fig12_water_diagnosis(root, out):
    """Diagnosis of the vacuum-vs-water pKa discrepancy (outputs of analyze_vacuum_vs_water_pka.py
    and xtb_solvation_on_dft_geometries.py)."""
    F = Fig("Fig12_vacuum_vs_water_diagnosis", 18.0, 7.0)
    conf = read_csv(os.path.join(root, "vac_vs_water_decomp_conformers.csv"))
    alp = {(r["state"], r["set"], r["job"]): float(r["dG_solv_alpb_kcal"])
           for r in read_csv(os.path.join(root, "xtb_solvation_on_dft_geometries.csv"))}
    F.panel("a", 1.8, 0.6, 5.2, 5.4, "dG(solv) COSMO (kcal/mol)",
            "dG(solv) GFN2/ALPB (kcal/mol)", -0.3, 3.7, -0.5, 1.6)
    for st in "NABD":
        rows_ = [r for r in conf if r["state"] == st and r["dG_solv_COSMO_kcal"] != ""
                 and (r["state"], r["set"], r["job"]) in alp]
        ref = [r for r in rows_ if float(r["G_aq_rel_kcal"]) == 0.0][0]
        x = [float(r["dG_solv_COSMO_kcal"]) for r in rows_]
        y = [alp[(r["state"], r["set"], r["job"])] - alp[(ref["state"], ref["set"], ref["job"])] for r in rows_]
        rr = np.corrcoef(x, y)[0, 1]
        F.xy(x, y, color=COL[st], marker="circle", size="3pt", line=False, key=f"{SNAME[st].split(' ')[0]} {st}  (r = {rr:.2f})")
    F.key("left", "top")
    # (b) monoB decomposition
    rb = [r for r in conf if r["state"] == "B" and r["dG_solv_COSMO_kcal"] != ""]
    rb.sort(key=lambda r: float(r["G_aq_rel_kcal"]))
    n = len(rb)
    F.panel("b", 9.0, 0.6, 4.6, 5.4, "Monoanion B conformer (sorted)", "Relative to best G_aq (kcal/mol)", 0.3, n + 0.7, -1.2, 3.6)
    xs = list(range(1, n + 1))
    F.hline(0, 0.3, n + 0.7)
    F.xy(xs, [float(r["dE_gas_HSE06_kcal"]) for r in rb], color="#999999", marker="square", size="3pt", line=False, key="gas HSE06")
    F.xy(xs, [float(r["dG_solv_COSMO_kcal"]) for r in rb], color="#009E73", marker="circle", size="3pt", line=False, key="solvation COSMO")
    F.xy(xs, [float(r["G_aq_rel_kcal"]) for r in rb], color="black", marker="diamond", size="3pt", line=False, key="total G_aq")
    new = [i for i, r in enumerate(rb, 1) if r["set"] == "new"]
    F.xy(new, [-1.0] * len(new), color="#D55E00", marker="triangle", size="4pt", line=False, key="water-guided (new)")
    F.key("left", "top")
    # (c) N->B calibrated pKa under successive definitions
    pk = read_csv(os.path.join(root, "vac_vs_water_decomp_pka.csv"))
    def pick(pair, w, g, st):
        return float([r for r in pk if r["pair"] == pair and r["weights"] == w and r["state_G"] == g and r["set"] == st][0]["calib_pKa"])
    labs = ["prod.", "G_aq w.", "+degen.\n+S_conf", "+water\nextra\n(no deg.)", "+water\nextra\n(correct)"]
    vals = [pick("N_to_B", "pop", "mean", "old"), pick("N_to_B", "gaq", "mean", "old"), pick("N_to_B", "gaqg", "ens", "old"),
            pick("N_to_B", "gaq", "mean", "all"), pick("N_to_B", "gaqg", "ens", "all")]
    F.panel("c", 14.6, 0.6, 3.0, 5.4, "", "Calibrated pKa (N->B)", 0.4, 5.6, 6, 9)
    F._xrange = (0.4, 5.6)
    F.xy(range(1, 6), vals, color="#0072B2", marker="square", size="4pt", line=True, width="1pt")
    F.xticks(range(1, 6), ["1", "2", "3", "4", "5"], size="7pt")
    F.save(out)


BUILDERS = {
    "fig01": fig01_workflow, "fig02": fig02_ensembles, "fig03": fig03_functional, "fig04": fig04_b1,
    "fig05": fig05_audit, "fig06": fig06_uvvis_states, "fig07": fig07_exp_and_ph, "fig08": fig08_ir,
    "fig09": fig09_pka, "fig10": fig10_speciation, "fig11": fig11_vac_water, "fig12": fig12_water_diagnosis,
}


def main():
    ap = argparse.ArgumentParser(description="Build all Paper 2 Veusz figures.")
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
