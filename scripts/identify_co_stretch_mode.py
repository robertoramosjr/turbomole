#!/usr/bin/env python3
"""
identify_co_stretch_mode.py

Identifies which vibrational mode in a frequency calculation is the
carbonyl (C=O) stretch, for B1's delta-learning target/proxy. Picking
"the highest-intensity band in a frequency window" is not reliable --
other modes (aromatic ring stretches, C-H bends) can be just as
IR-intense in the same 1500-1800 cm-1 region, especially once the
carboxyl deprotonates into a carboxylate (which shifts and splits into
symmetric/antisymmetric stretches). Instead, this projects each mode's
displacement onto the actual C=O bond vector: for the carbonyl carbon
and oxygen specifically, the *relative* displacement (oxygen minus
carbon) is projected onto the bond direction, and the mode with the
largest |projection| above --min-wavenumber is the C=O stretch --
this directly measures how much a mode stretches/compresses that one
bond, regardless of what else lights up nearby in frequency space.

Carbonyl atom detection reuses the same geometry-based logic as
build_protonation_state.py (O bonded to a C that also carries a second
O within carbonyl bonding distance) rather than hardcoded indices, so
this works on any conformer/protonation state of this molecule without
edits, and should generalize to other carboxylic-acid-bearing molecules.

Input format: Gaussian-98-style frequency output (g98.out from
`xtb --hess`, or an equivalent export) -- 3 modes per block, atom-by-
atom XYZ displacement vectors.

Usage:
    python ~/work_turbomole/scripts/identify_co_stretch_mode.py \
        --g98 g98.out --xyz structure.xyz [--min-wavenumber 800]
"""

import argparse
import re

import numpy as np

CARBONYL_BOND_MAX = 1.35
CO_BOND_MAX = 1.75


def read_xyz(path):
    with open(path) as f:
        lines = f.readlines()
    n = int(lines[0].split()[0])
    atoms = []
    for line in lines[2:2 + n]:
        parts = line.split()
        atoms.append((parts[0], np.array(list(map(float, parts[1:4])))))
    return atoms


def find_carbonyl_pair(atoms):
    """Return (o_index, c_index), 0-based, for the carbonyl O=C bond
    (the C=O oxygen with no attached H, whose carbon also carries the
    hydroxyl/carboxylate oxygen) -- same distance-based logic as
    build_protonation_state.py's hydroxyl detection."""
    def dist(i, j):
        return np.linalg.norm(atoms[i][1] - atoms[j][1])

    for oi, (sym_o, _) in enumerate(atoms):
        if sym_o != "O":
            continue
        # carbonyl O has no H neighbor
        if any(atoms[h][0] == "H" and dist(oi, h) < 1.15
               for h in range(len(atoms))):
            continue
        for ci, (sym_c, _) in enumerate(atoms):
            if sym_c != "C" or dist(oi, ci) >= CARBONYL_BOND_MAX:
                continue
            # confirm this C also has a second O (the hydroxyl/carboxylate one)
            other_o = [j for j, (s, _) in enumerate(atoms)
                       if s == "O" and j != oi and dist(ci, j) < CO_BOND_MAX]
            if other_o:
                return oi, ci
    return None, None


def parse_g98_modes(path, n_atoms):
    """Return list of dicts: frequency_cm1, ir_intensity, displacements
    (n_atoms x 3 array), in file order."""
    with open(path) as f:
        lines = f.readlines()

    modes = []
    i = 0
    while i < len(lines):
        if lines[i].strip().startswith("Frequencies --"):
            freqs = [float(x) for x in lines[i].split("--")[1].split()]
            inten_line = next(
                l for l in lines[i:i + 5] if l.strip().startswith("IR Inten")
            )
            intensities = [float(x) for x in inten_line.split("--")[1].split()]
            header_idx = i + 1
            while not lines[header_idx].strip().startswith("Atom"):
                header_idx += 1
            disp = [[] for _ in freqs]
            for a in range(n_atoms):
                parts = lines[header_idx + 1 + a].split()
                vals = list(map(float, parts[2:]))
                for m in range(len(freqs)):
                    disp[m].append(vals[3 * m:3 * m + 3])
            for m in range(len(freqs)):
                modes.append({
                    "frequency_cm1": freqs[m],
                    "ir_intensity": intensities[m],
                    "displacement": np.array(disp[m]),
                })
            i = header_idx + 1 + n_atoms
        else:
            i += 1
    return modes


def main():
    parser = argparse.ArgumentParser(
        description="Identify the C=O stretch mode via bond-vector "
                    "displacement projection.")
    parser.add_argument("--g98", required=True, help="Gaussian-98-style frequency output")
    parser.add_argument("--xyz", required=True, help="Geometry (.xyz) matching the g98 file")
    parser.add_argument("--min-wavenumber", type=float, default=800.0,
                         help="Ignore modes below this wavenumber [default: 800 cm-1]")
    args = parser.parse_args()

    atoms = read_xyz(args.xyz)
    o_idx, c_idx = find_carbonyl_pair(atoms)
    if o_idx is None:
        raise SystemExit("Could not identify a carbonyl C=O pair in this geometry.")

    bond_vec = atoms[o_idx][1] - atoms[c_idx][1]
    bond_dir = bond_vec / np.linalg.norm(bond_vec)

    modes = parse_g98_modes(args.g98, len(atoms))

    best_mode = None
    for m in modes:
        if m["frequency_cm1"] < args.min_wavenumber:
            continue
        rel_disp = m["displacement"][o_idx] - m["displacement"][c_idx]
        projection = abs(np.dot(rel_disp, bond_dir))
        m["co_projection"] = projection
        if best_mode is None or projection > best_mode["co_projection"]:
            best_mode = m

    if best_mode is None:
        raise SystemExit("No mode found above --min-wavenumber.")

    print(f"Carbonyl atoms (0-indexed): O{o_idx + 1}, C{c_idx + 1}, "
          f"bond length {np.linalg.norm(bond_vec):.3f} A")
    print(f"C=O stretch mode: {best_mode['frequency_cm1']:.2f} cm-1, "
          f"IR intensity {best_mode['ir_intensity']:.2f}, "
          f"bond-projection {best_mode['co_projection']:.4f}")

    top5 = sorted(modes, key=lambda m: -m.get("co_projection", 0))[:5]
    print("\nTop 5 modes by C=O bond-vector projection (sanity check -- "
          "should show one clear winner, not several close together):")
    for m in top5:
        print(f"  {m['frequency_cm1']:8.2f} cm-1  "
              f"proj={m.get('co_projection', 0):.4f}  "
              f"IR={m['ir_intensity']:.2f}")


if __name__ == "__main__":
    main()
