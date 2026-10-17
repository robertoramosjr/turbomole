#!/usr/bin/env python3
"""Collect Bader cube-resolution test results (h_<spacing>/ACF.dat) in the cwd into a CSV."""
import argparse
import csv
from pathlib import Path

ATOMIC_NUMBER = {
    "h": 1, "he": 2, "li": 3, "be": 4, "b": 5, "c": 6, "n": 7, "o": 8, "f": 9, "ne": 10,
    "na": 11, "mg": 12, "al": 13, "si": 14, "p": 15, "s": 16, "cl": 17, "ar": 18,
    "k": 19, "ca": 20, "sc": 21, "ti": 22, "v": 23, "cr": 24, "mn": 25, "fe": 26,
    "co": 27, "ni": 28, "cu": 29, "zn": 30, "br": 35, "i": 53,
}
METALS = {"sc", "ti", "v", "cr", "mn", "fe", "co", "ni", "cu", "zn"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pattern", default="h_*", help="Glob for the per-spacing directories")
    parser.add_argument("--coord", default="coord", help="TURBOMOLE coord file name inside each directory")
    parser.add_argument("--acf", default="ACF.dat", help="Bader atomic charge file name")
    parser.add_argument("--grid-info", default="grid_info.csv", help="Per-directory grid/timing file")
    parser.add_argument("--output", default="bader_grid_convergence.csv", help="Output CSV")
    return parser.parse_args()


def read_elements(path: Path) -> list[str]:
    elements, inside = [], False
    for line in path.read_text().splitlines():
        if line.startswith("$coord"):
            inside = True
            continue
        if inside and line.startswith("$"):
            break
        if inside and line.split():
            elements.append(line.split()[3].lower())
    return elements


def read_acf(path: Path) -> list[float]:
    charges = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) >= 7 and fields[0].isdigit():
            charges.append(float(fields[4]))
    return charges


def read_grid_info(path: Path) -> dict:
    if not path.exists():
        return {}
    h, nx, ny, nz, t_proper, t_bader = path.read_text().strip().split(",")
    return {"npoints": int(nx) * int(ny) * int(nz), "proper_s": int(t_proper), "bader_s": int(t_bader)}


def summarize(directory: Path, args: argparse.Namespace) -> dict:
    elements = read_elements(directory / args.coord)
    charges = read_acf(directory / args.acf)
    if len(elements) != len(charges):
        raise ValueError(f"{directory}: {len(elements)} atoms in coord but {len(charges)} in ACF")
    net = [ATOMIC_NUMBER[e] - q for e, q in zip(elements, charges)]
    n_ref = sum(ATOMIC_NUMBER[e] for e in elements)
    row = {
        "spacing_bohr": float(directory.name.split("_", 1)[1]),
        "n_electrons_bader": round(sum(charges), 5),
        "n_electrons_reference": n_ref,
        "missing_electrons": round(n_ref - sum(charges), 5),
        "sum_net_charge_e": round(sum(net), 5),
    }
    for idx, (e, q) in enumerate(zip(elements, net), start=1):
        if e in METALS:
            row[f"net_charge_{e}{idx}_e"] = round(q, 5)
    row.update(read_grid_info(directory / args.grid_info))
    return row


def main() -> None:
    args = parse_args()
    cwd = Path.cwd()
    rows = [summarize(d, args) for d in sorted(cwd.glob(args.pattern)) if (d / args.acf).exists()]
    rows.sort(key=lambda r: -r["spacing_bohr"])
    if not rows:
        raise SystemExit(f"no {args.acf} found under {args.pattern}")
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with (cwd / args.output).open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
