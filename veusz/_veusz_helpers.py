"""Shared helpers for the azobenzene Veusz plot scripts (headless-safe)."""

import os

if not os.environ.get("DISPLAY"):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import veusz.embed as veusz  # noqa: E402,F401


def new_graph(name, width="14cm", height="9cm", xlabel="", ylabel="", xmin=None, xmax=None,
              ymin=None, ymax=None):
    """Create an embedded Veusz doc with one page/graph and labelled axes."""
    g = veusz.Embedded(name, hidden=True)
    g.To(g.Add("page", name="page1"))
    g.Set("width", width)
    g.Set("height", height)
    g.To(g.Add("graph", name="graph1"))
    g.Set("x/label", xlabel)
    g.Set("y/label", ylabel)
    for axis in ("x", "y"):
        g.Set(f"{axis}/Label/size", "10pt")
        g.Set(f"{axis}/TickLabels/size", "9pt")
    for axis, lo, hi in (("x", xmin, xmax), ("y", ymin, ymax)):
        if lo is not None:
            g.Set(f"{axis}/min", float(lo))
        if hi is not None:
            g.Set(f"{axis}/max", float(hi))
    return g


def add_xy(g, name, x, y, color="black", width="1.2pt", marker="none", size="3pt",
           line=True, key=None, style="solid"):
    """Add an xy widget from numpy arrays (data embedded in the .vsz)."""
    g.SetData(f"{name}_x", [float(v) for v in x])
    g.SetData(f"{name}_y", [float(v) for v in y])
    g.Add("xy", name=name)
    g.Set(f"{name}/xData", f"{name}_x")
    g.Set(f"{name}/yData", f"{name}_y")
    g.Set(f"{name}/marker", marker)
    g.Set(f"{name}/markerSize", size)
    g.Set(f"{name}/MarkerFill/color", color)
    g.Set(f"{name}/MarkerLine/color", color)
    g.Set(f"{name}/PlotLine/hide", not line)
    g.Set(f"{name}/PlotLine/color", color)
    g.Set(f"{name}/PlotLine/width", width)
    g.Set(f"{name}/PlotLine/style", style)
    if key:
        g.Set(f"{name}/key", key)


def add_vline(g, name, x, ymin, ymax, color="gray", style="dashed", key=None):
    add_xy(g, name, [x, x], [ymin, ymax], color=color, width="1pt", style=style, key=key)


def finish(g, output, key_pos=None):
    k = g.Add("key", name="key1")  # noqa: F841
    if key_pos:
        g.Set("key1/horzPosn", key_pos[0])
        g.Set("key1/vertPosn", key_pos[1])
    g.Set("key1/Border/hide", True)
    g.Set("key1/Text/size", "7pt")
    g.Set("key1/keyLength", "0.6cm")
    g.Save(f"{output}.vsz")
    g.Export(f"{output}.pdf")
    g.Export(f"{output}.png", dpi=200)
    print(f"Saved -> {output}.vsz / .pdf / .png")
    g.Close()
