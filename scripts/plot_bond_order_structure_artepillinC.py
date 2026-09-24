#!/usr/bin/env python3
"""
plot_bond_order_structure_artepillinC.py

Combined figure for Artepillin C: (a) the SEN bar chart, (b) the
DFT-relaxed molecular structure with the same bonds color-coded
directly on the geometry, instead of only via small rotated text
labels on the bars.

Builds an actual RDKit Mol (correct elements, bond orders, aromaticity)
from the DFT-relaxed geometry + SEN-derived connectivity, so bonds
render with proper cheminformatics conventions (double bonds as
parallel lines, valence-correct aromatic ring, element symbols)
instead of hand-drawn matplotlib circles/line segments. The 2D layout
is still OUR PCA projection of the true relaxed coordinates (not
RDKit's own 2D coordinate generator), so the depicted geometry stays
faithful to the actual DFT structure -- RDKit is used here purely as
the rendering engine. Bond categories/colors match
plot_bond_order_artepillinC.py exactly. Molecule-specific (atom
numbering, ring/double-bond/carbonyl atom indices) -- re-derive for
any other molecule.

Requires rdkit + matplotlib, which the base env lacks -- installed
into the vasp_env conda env instead.

Usage:
    ~/miniconda3/envs/vasp_env/bin/python3 plot_bond_order_structure_artepillinC.py \
        --xyz coord.xyz --labels bond_order_dataset_labels.csv --top 20 \
        --output bond_order_structure_figure
"""
import argparse
import csv
import io
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point3D

CATEGORIES = {
    "aromatic": ((0, 0, 0), "Aromatic ring"),
    "carbonyl": ((0.10, 0.10, 0.95), "Carbonyl C=O"),
    "cc_double": ((0.90, 0.10, 0.10), "C=C double bonds"),
    "other": ((0.55, 0.55, 0.55), "Other (conjugation, prenyl)"),
}
AROMATIC_RING_ATOMS = {4, 5, 6, 9, 10, 11}
CC_DOUBLE_PAIRS = {"16c-21c", "12c-14c", "13c-15c"}
CARBONYL_PAIRS = {"3o-22c"}
LABELED_PAIRS = {"16c-21c", "13c-15c", "12c-14c", "3o-22c", "21c-22c"}
OTHER_HIGHLIGHTED = {
    "11c-16c", "21c-22c", "7c-12c", "4c-7c", "15c-20c",
    "15c-19c", "14c-17c", "14c-18c", "8c-13c", "5c-8c",
}
REAL_BOND_SEN_THRESHOLD = 0.5
ELEM_SYMBOL = {"c": "C", "o": "O", "h": "H"}


def format_label(pair_label):
    a, b = pair_label.split("-")
    out = []
    for tok in (a, b):
        idx, elem = re.match(r"(\d+)([a-z]+)", tok).groups()
        out.append(f"{elem.upper()}{idx}")
    return out


def classify(pair_label):
    if pair_label in CC_DOUBLE_PAIRS:
        return "cc_double"
    if pair_label in CARBONYL_PAIRS:
        return "carbonyl"
    a, b = pair_label.split("-")
    atoms = []
    for tok in (a, b):
        idx, elem = re.match(r"(\d+)([a-z]+)", tok).groups()
        atoms.append((int(idx), elem))
    if atoms[0][1] == "c" and atoms[1][1] == "c" and \
       atoms[0][0] in AROMATIC_RING_ATOMS and atoms[1][0] in AROMATIC_RING_ATOMS:
        return "aromatic"
    return "other"


def read_xyz(path):
    with open(path) as f:
        n = int(f.readline())
        f.readline()
        elements, coords = [], []
        for _ in range(n):
            parts = f.readline().split()
            elements.append(parts[0].lower())
            coords.append([float(x) for x in parts[1:4]])
    return elements, np.array(coords)


def read_pairs(path):
    rows = []
    with open(path) as f:
        for row in csv.DictReader(f):
            rows.append({
                "pair_label": row["pair_label"],
                "a": int(row["atom1_index"]), "b": int(row["atom2_index"]),
                "sen": float(row["shared_electron_number"]),
            })
    return rows


def pca_project(coords):
    centered = coords - coords.mean(axis=0)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    basis = vt[:2]
    return centered @ basis.T


def build_mol(elements, xy, real_bonds):
    rw = Chem.RWMol()
    for elem in elements:
        a = Chem.Atom(ELEM_SYMBOL[elem])
        a.SetNoImplicit(True)
        rw.AddAtom(a)

    bond_cat = {}
    for p in real_bonds:
        cat = classify(p["pair_label"])
        i, j = p["a"] - 1, p["b"] - 1
        if cat == "cc_double" or cat == "carbonyl":
            bt = Chem.BondType.DOUBLE
        elif cat == "aromatic":
            bt = Chem.BondType.AROMATIC
        else:
            bt = Chem.BondType.SINGLE
        rw.AddBond(i, j, bt)
        bond_cat[frozenset((i, j))] = (cat, p)

    for atom_idx in AROMATIC_RING_ATOMS:
        rw.GetAtomWithIdx(atom_idx - 1).SetIsAromatic(True)
    for p in real_bonds:
        if classify(p["pair_label"]) == "aromatic":
            b = rw.GetBondBetweenAtoms(p["a"] - 1, p["b"] - 1)
            b.SetIsAromatic(True)

    mol = rw.GetMol()
    Chem.SanitizeMol(mol)

    conf = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        x, y = xy[i]
        conf.SetAtomPosition(i, Point3D(float(x), float(y), 0.0))
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)

    return mol, bond_cat


def draw_bar_panel(ax, top_pairs, label_annotations):
    n = len(top_pairs)
    for cat, (color, _) in CATEGORIES.items():
        rr = [i + 1 for i, p in enumerate(top_pairs) if classify(p["pair_label"]) == cat]
        vv = [top_pairs[i - 1]["sen"] for i in rr]
        ax.bar(rr, vv, width=0.9, color=color, edgecolor="black", linewidth=0.6)

    for rank, sen, text in label_annotations:
        ax.text(rank, sen + 0.04, text, rotation=90, ha="center", va="bottom", fontsize=7)

    ax.set_xlim(0, n + 1)
    ax.set_ylim(0, max(p["sen"] for p in top_pairs) * 1.3)
    ax.set_xlabel("Bond Pair (ranked by SEN)")
    ax.set_ylabel("Shared electron number (SEN)")
    ax.set_title("(a)", loc="left", fontweight="bold")

    handles = [Line2D([0], [0], marker="s", linestyle="", markersize=10,
                       markerfacecolor=color, markeredgecolor="black")
               for color, _ in CATEGORIES.values()]
    labels = [text for _, text in CATEGORIES.values()]
    ax.legend(handles, labels, loc="upper right", fontsize=7, ncol=1, framealpha=0.9)


def render_rdkit_structure(mol, bond_cat, elements, size=(900, 900)):
    d = rdMolDraw2D.MolDraw2DCairo(*size)
    opts = d.drawOptions()
    opts.addAtomIndices = False
    opts.bondLineWidth = 2
    opts.padding = 0.12
    opts.updateAtomPalette({8: (0, 0, 0)})  # oxygen: black, not the default red -- easier to read

    highlight_bonds, highlight_bond_colors = [], {}
    for bond in mol.GetBonds():
        key = frozenset((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        if key in bond_cat:
            cat, p = bond_cat[key]
            if cat in ("cc_double", "carbonyl", "aromatic") or p["pair_label"] in OTHER_HIGHLIGHTED:
                highlight_bonds.append(bond.GetIdx())
                highlight_bond_colors[bond.GetIdx()] = CATEGORIES[cat][0]

    rdMolDraw2D.PrepareAndDrawMolecule(
        d, mol, highlightAtoms=[], highlightBonds=highlight_bonds,
        highlightBondColors=highlight_bond_colors)

    # Atom-index labels are drawn ourselves (not via RDKit's atomNote,
    # which places them right on top of the atom/bonds) so they can be
    # pushed radially outward, clear of the bonds, without moving the
    # bonds/atoms themselves.
    px = [d.GetDrawCoords(i) for i in range(mol.GetNumAtoms())]
    px = np.array([[p.x, p.y] for p in px])
    centroid = px.mean(axis=0)
    label_px = []
    for i, elem in enumerate(elements):
        if elem == "h":
            continue
        neighbor_idx = [nb.GetIdx() for nb in mol.GetAtomWithIdx(i).GetNeighbors()]
        if neighbor_idx:
            # Place the label in the widest angular gap between bonds to
            # this atom, rather than a simple "away from neighbor
            # centroid" direction -- the latter can point straight down
            # a bond line when the neighbors aren't symmetric (as
            # happens with the real, non-idealized DFT geometry), which
            # dropped labels right on top of highlighted bonds.
            vecs = px[neighbor_idx] - px[i]
            angles = sorted(np.arctan2(v[1], v[0]) for v in vecs)
            angles.append(angles[0] + 2 * np.pi)
            gaps = [(angles[k + 1] - angles[k], angles[k], angles[k + 1])
                    for k in range(len(angles) - 1)]
            _, a0, a1 = max(gaps, key=lambda g: g[0])
            bisector = (a0 + a1) / 2
            direction = np.array([np.cos(bisector), np.sin(bisector)])
        else:
            direction = px[i] - centroid
            direction /= np.linalg.norm(direction)
        label_px.append((px[i] + direction * 24, str(i + 1)))

    d.FinishDrawing()
    png_bytes = d.GetDrawingText()
    img = plt.imread(io.BytesIO(png_bytes), format="png")
    return img, label_px


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xyz", required=True)
    ap.add_argument("--labels", required=True)
    ap.add_argument("--top", type=int, default=20)
    ap.add_argument("--output", default="bond_order_structure_figure")
    args = ap.parse_args()

    elements, coords = read_xyz(args.xyz)
    xy = pca_project(coords)

    all_pairs = read_pairs(args.labels)
    real_bonds = [p for p in all_pairs if p["sen"] > REAL_BOND_SEN_THRESHOLD]
    top_pairs = sorted(all_pairs, key=lambda r: r["sen"], reverse=True)[:args.top]

    label_annotations = []
    for i, p in enumerate(top_pairs, start=1):
        if p["pair_label"] in LABELED_PAIRS:
            elem1, elem2 = format_label(p["pair_label"])
            cat = classify(p["pair_label"])
            sep = "=" if cat in ("cc_double", "carbonyl") else "-"
            label_annotations.append((i, p["sen"], f"{elem1}{sep}{elem2}"))

    mol, bond_cat = build_mol(elements, xy, real_bonds)
    img, label_px = render_rdkit_structure(mol, bond_cat, elements)

    fig, (ax_bar, ax_struct) = plt.subplots(1, 2, figsize=(13, 6),
                                             gridspec_kw={"width_ratios": [1.15, 1]})
    draw_bar_panel(ax_bar, top_pairs, label_annotations)

    ax_struct.imshow(img)
    for (x, y), text in label_px:
        ax_struct.text(x, y, text, fontsize=7.5, ha="center", va="center",
                        color="dimgrey",
                        bbox=dict(boxstyle="round,pad=0.1", facecolor="white",
                                  edgecolor="none", alpha=0.7))
    ax_struct.axis("off")
    ax_struct.set_title("(b) DFT-relaxed structure", loc="left", fontweight="bold")

    fig.tight_layout()
    fig.savefig(f"{args.output}.pdf", bbox_inches="tight")
    fig.savefig(f"{args.output}.png", dpi=200, bbox_inches="tight")
    print(f"Saved -> {args.output}.pdf and {args.output}.png")


if __name__ == "__main__":
    main()
