#!/usr/bin/env python3
"""
parse_dos.py

Parses Turbomole $pop ... dos output files (total DOS and per-element
projected DOS from separate atoms-restricted runs) and combines them
into a single Veusz-importable dataset (descriptor format).

Now also extracts the s/p/d sub-shell contributions per element
(not just the total), needed to reproduce panel (d) of Fig. 3 in
Januario & Cabral (2026): pDOS-Cs / pDOS-Cp.

Usage:
    python ~/scripts/parse_dos.py \
        --total dos_total --carbon dos_carbon \
        --hydrogen dos_hydrogen --oxygen dos_oxygen \
        --output dos_dataset.dat
"""

import argparse
import numpy as np

HARTREE_TO_EV = 27.211386245988


def read_dos_file(path):
    column_names = None
    rows = []

    with open(path, "r") as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                tokens = stripped.lstrip("#").split()
                if tokens and tokens[0].lower() == "energy":
                    column_names = tokens
                continue
            rows.append([float(x) for x in stripped.split()])

    data = np.array(rows)
    n_data_columns = data.shape[1]

    if column_names is None:
        raise ValueError(f"Could not find a header line starting with 'energy' in {path}")

    if len(column_names) < n_data_columns:
        n_missing = n_data_columns - len(column_names)
        print(f"WARNING [{path}]: header names {len(column_names)} columns "
              f"but data has {n_data_columns}. Labeling {n_missing} extra "
              f"column(s) as extra_1, extra_2, ...")
        column_names = column_names + [f"extra_{i+1}" for i in range(n_missing)]

    return column_names, data


def main():
    parser = argparse.ArgumentParser(
        description="Combine Turbomole total DOS and element-projected DOS "
                    "(including s/p sub-shell breakdown) into a Veusz dataset.")
    parser.add_argument("--total", required=True)
    parser.add_argument("--carbon", required=True)
    parser.add_argument("--hydrogen", required=True)
    parser.add_argument("--oxygen", required=True)
    parser.add_argument("--output", default="dos_dataset.dat")
    args = parser.parse_args()

    names_total, data_total = read_dos_file(args.total)
    names_c, data_c = read_dos_file(args.carbon)
    names_h, data_h = read_dos_file(args.hydrogen)
    names_o, data_o = read_dos_file(args.oxygen)

    energy_hartree = data_total[:, 0]
    for label, data in [("carbon", data_c), ("hydrogen", data_h), ("oxygen", data_o)]:
        if data.shape[0] != energy_hartree.shape[0] or \
           not np.allclose(data[:, 0], energy_hartree, atol=1e-6):
            print(f"WARNING: energy grid in '{label}' file does not match the total DOS grid.")

    energy_ev = energy_hartree * HARTREE_TO_EV

    tdos = data_total[:, names_total.index("total")]
    dos_s = data_total[:, names_total.index("s")]
    dos_p = data_total[:, names_total.index("p")]
    dos_d = data_total[:, names_total.index("d")]

    pdos_c_total = data_c[:, names_c.index("total")]
    pdos_c_s = data_c[:, names_c.index("s")]
    pdos_c_p = data_c[:, names_c.index("p")]

    pdos_h_total = data_h[:, names_h.index("total")]
    pdos_o_total = data_o[:, names_o.index("total")]

    descriptor_names = ["energy_eV", "TDOS", "DOS_s", "DOS_p", "DOS_d",
                         "pDOS_C", "pDOS_Cs", "pDOS_Cp", "pDOS_H", "pDOS_O"]
    descriptor_line = "descriptor " + ",".join(descriptor_names)

    with open(args.output, "w") as f:
        f.write(descriptor_line + "\n")
        for i in range(len(energy_ev)):
            row = [f"{energy_ev[i]:.6f}", f"{tdos[i]:.6f}", f"{dos_s[i]:.6f}",
                   f"{dos_p[i]:.6f}", f"{dos_d[i]:.6f}", f"{pdos_c_total[i]:.6f}",
                   f"{pdos_c_s[i]:.6f}", f"{pdos_c_p[i]:.6f}",
                   f"{pdos_h_total[i]:.6f}", f"{pdos_o_total[i]:.6f}"]
            f.write(" ".join(row) + "\n")

    print(f"Veusz dataset saved -> {args.output}")


if __name__ == "__main__":
    main()