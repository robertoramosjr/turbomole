#!/usr/bin/env python3
"""
plot_bond_order_artepillinC.py

Molecule-specific variant of plot_bond_order.py for Artepillin C: bar
plot of the top-N strongest atom-pair shared electron numbers (SEN),
color-coded by chemical role (aromatic ring / carbonyl C=O / C=C double
bonds / other conjugation+prenyl bonds) with text labels on the bonds
called out in the paper figure. The generic plot_bond_order.py has no
notion of bond chemistry, so this hardcodes the atom-pair -> category
mapping for THIS molecule's atom numbering (from coord/coord.xyz) --
don't reuse it for a different molecule or a re-numbered structure
without re-checking bond_order_dataset_labels.csv.

Usage:
    python ~/scripts/veusz/plot_bond_order_artepillinC.py \
        --labels bond_order_dataset_labels.csv --top 20 \
        --output bond_order_figure
"""

import argparse
import csv
import re
import veusz.embed as veusz

# Category -> (color, legend text). Order matters for legend order.
CATEGORIES = {
    "aromatic": ("black", "Aromatic ring"),
    "carbonyl": ("blue", "Carbonyl C=O"),
    "cc_double": ("red", "C=C double bonds"),
    "other": ("grey", "Other (conjugation, prenyl)"),
}

# Aromatic ring carbons (from coord numbering): hexagon 4-6-5-10-11-9-4
AROMATIC_RING_ATOMS = {4, 5, 6, 9, 10, 11}

# Pair labels (as they appear in bond_order_dataset_labels.csv) that are
# genuine C=C double bonds: the cinnamic-acid vinyl bond + the two
# prenyl C=C bonds (diprenyl substitution).
CC_DOUBLE_PAIRS = {"16c-21c", "12c-14c", "13c-15c"}

# The carboxylic C=O (O3 is the carbonyl oxygen; O2 is the -OH oxygen,
# which has a much lower SEN and isn't in the top-N).
CARBONYL_PAIRS = {"3o-22c"}

# Pairs to annotate with a text label above the bar (subset shown in
# the reference figure -- not every bar needs a label).
LABELED_PAIRS = {"16c-21c", "13c-15c", "12c-14c", "3o-22c", "21c-22c"}


def format_label(pair_label):
    """'16c-21c' -> ('C16', '21C'... ) -> display string, e.g. 'C16=C21'."""
    a, b = pair_label.split("-")
    parts = []
    for tok in (a, b):
        m = re.match(r"(\d+)([a-z]+)", tok)
        idx, elem = m.groups()
        parts.append(f"{elem.upper()}{idx}")
    return parts


def classify(pair_label):
    if pair_label in CC_DOUBLE_PAIRS:
        return "cc_double"
    if pair_label in CARBONYL_PAIRS:
        return "carbonyl"
    a, b = pair_label.split("-")
    atoms = []
    for tok in (a, b):
        m = re.match(r"(\d+)([a-z]+)", tok)
        idx, elem = m.groups()
        atoms.append((int(idx), elem))
    if atoms[0][1] == "c" and atoms[1][1] == "c" and \
       atoms[0][0] in AROMATIC_RING_ATOMS and atoms[1][0] in AROMATIC_RING_ATOMS:
        return "aromatic"
    return "other"


def read_labels_csv(path, top_n):
    rows = []
    with open(path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({
                "pair_label": row["pair_label"],
                "sen": float(row["shared_electron_number"]),
            })
    rows.sort(key=lambda r: r["sen"], reverse=True)
    return rows[:top_n]


def main():
    parser = argparse.ArgumentParser(
        description="Plot top-N SEN bond-order pairs, color-coded by bond chemistry, for Artepillin C.")
    parser.add_argument("--labels", required=True, help="Path to bond_order_dataset_labels.csv")
    parser.add_argument("--top", type=int, default=20, help="Number of strongest pairs to display")
    parser.add_argument("--output", default="bond_order_figure", help="Output filename prefix")
    args = parser.parse_args()

    top_pairs = read_labels_csv(args.labels, args.top)
    n = len(top_pairs)

    print(f"Plotting top {n} pairs by SEN:")
    # One length-N column per category, zero-padded outside that
    # category's ranks, plotted as a single stacked bar widget sharing
    # one posn=1..N dataset. (Splitting into separate bar widgets per
    # category is broken in Veusz: a category with a single rank makes
    # its posn dataset length 1, which triggers Veusz's bar-width
    # fallback of using the *entire graph width* for that bar -- it
    # then covers every other bar.)
    by_category = {cat: [0.0] * n for cat in CATEGORIES}
    label_annotations = []
    for i, p in enumerate(top_pairs):
        cat = classify(p["pair_label"])
        by_category[cat][i] = p["sen"]
        print(f"  {i + 1}. {p['pair_label']} [{cat}]: {p['sen']:.4f}")
        if p["pair_label"] in LABELED_PAIRS:
            elem1, elem2 = format_label(p["pair_label"])
            sep = "=" if cat in ("cc_double", "carbonyl") else "-"
            label_annotations.append((i + 1, p["sen"], f"{elem1}{sep}{elem2}"))

    g = veusz.Embedded("bond_order_figure", hidden=True)

    g.To(g.Add("page", name="page1"))
    g.Set("width", "16cm")
    g.Set("height", "12cm")
    g.To(g.Add("graph", name="graph1"))

    g.Set("x/label", "Bond Pair (ranked by SEN)")
    g.Set("x/min", 0)
    g.Set("x/max", n + 1)

    g.Set("y/label", "Shared electron number (SEN)")
    g.Set("y/min", 0)
    # Explicit headroom above the tallest (labeled) bars so the rotated
    # text labels and the legend box don't collide with the bars.
    g.Set("y/max", max(p["sen"] for p in top_pairs) * 1.3)

    g.SetData("rank_all", list(range(1, n + 1)))

    active_cats = [cat for cat in CATEGORIES if any(by_category[cat])]
    length_names = []
    fills = []
    keys = []
    for cat in active_cats:
        color, legend_text = CATEGORIES[cat]
        g.SetData(f"sen_{cat}", by_category[cat])
        length_names.append(f"sen_{cat}")
        fills.append(("solid", color, False))
        keys.append(legend_text)

    g.Add("bar", name="bar1")
    g.Set("bar1/posn", "rank_all")
    g.Set("bar1/lengths", length_names)
    g.Set("bar1/mode", "stacked")
    g.Set("bar1/barfill", 0.9)
    g.Set("bar1/BarFill/fills", fills)
    g.Set("bar1/keys", keys)

    for rank, sen, text in label_annotations:
        name = f"label_{rank}"
        g.Add("label", name=name, label=text)
        g.Set(f"{name}/positioning", "axes")
        g.Set(f"{name}/xPos", rank)
        g.Set(f"{name}/yPos", sen + 0.05)
        g.Set(f"{name}/alignHorz", "left")
        g.Set(f"{name}/alignVert", "centre")
        g.Set(f"{name}/angle", 90)
        g.Set(f"{name}/Text/size", "7pt")

    # 2 columns so the long "Other (conjugation, prenyl)" entry doesn't
    # force the box wide enough to reach over the rank 1-3 labels.
    g.Add("key", name="key1")
    g.Set("key1/horzPosn", "right")
    g.Set("key1/vertPosn", "top")
    g.Set("key1/columns", 2)
    g.Set("key1/Text/size", "7pt")

    g.Save(f"{args.output}.vsz")
    g.Export(f"{args.output}.pdf")
    print(f"Saved -> {args.output}.vsz and {args.output}.pdf")
    g.Close()


if __name__ == "__main__":
    main()
