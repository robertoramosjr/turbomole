#!/usr/bin/env python3
"""Compare COSMO TD-DFT spectra computed with def2-TZVPD (diffuse) vs def2-TZVP.

For each conformer in a diffuse joblist, read both `exspectrum` files and report
the bright state (max oscillator strength, length gauge), the lowest root with
f > --fmin, and the shifts (TZVPD - TZVP) in energy and f.

Usage:
    python compare_diffuse_vs_tzvp.py --joblist todo_v12/diffuse/joblist.txt \
        --tzvp-root water_tddft --output todo_v12/analysis/diffuse_vs_tzvp.csv
Joblist lines: <.../STATE/cNNNN>, STATE in N,A,B,D.
"""
import argparse
import csv
import os
import sys

STATE_DIR = {"N": "neutral", "A": "monoanion_a", "B": "monoanion_b", "D": "dianion"}


def read_exspectrum(path):
    """Return list of (root, E_eV, f_len)."""
    roots = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            roots.append((int(p[0]), float(p[3]), float(p[-1])))
    return roots


def summarize(roots, fmin):
    bright = max(roots, key=lambda r: r[2])
    low = next((r for r in roots if r[2] > fmin), bright)
    return bright, low


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--joblist", required=True)
    ap.add_argument("--tzvp-root", required=True, help="water_tddft directory (neutral/, monoanion_a/, ...)")
    ap.add_argument("--fmin", type=float, default=0.05)
    ap.add_argument("--output", required=True)
    a = ap.parse_args()
    rows = []
    for d in (l.strip() for l in open(a.joblist) if l.strip()):
        state, conf = d.rstrip("/").split("/")[-2:]
        ref = os.path.join(a.tzvp_root, STATE_DIR[state], conf, "exspectrum")
        new = os.path.join(d, "exspectrum")
        if not (os.path.isfile(ref) and os.path.isfile(new)):
            print(f"missing: {state}/{conf} ref={os.path.isfile(ref)} new={os.path.isfile(new)}", file=sys.stderr)
            continue
        (bz, lz), (bd, ld) = summarize(read_exspectrum(ref), a.fmin), summarize(read_exspectrum(new), a.fmin)
        rows.append(dict(state=state, conf=conf,
                         bright_root_tzvp=bz[0], E_bright_tzvp=bz[1], f_bright_tzvp=bz[2],
                         bright_root_tzvpd=bd[0], E_bright_tzvpd=bd[1], f_bright_tzvpd=bd[2],
                         dE_bright=bd[1] - bz[1], df_bright=bd[2] - bz[2],
                         E_low_tzvp=lz[1], E_low_tzvpd=ld[1], dE_low=ld[1] - lz[1],
                         f_low_tzvp=lz[2], f_low_tzvpd=ld[2]))
    with open(a.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()})
    for r in rows:
        print(f"{r['state']}/{r['conf']}: bright {r['E_bright_tzvp']:.3f}->{r['E_bright_tzvpd']:.3f} eV "
              f"(d={r['dE_bright']:+.3f}), f {r['f_bright_tzvp']:.3f}->{r['f_bright_tzvpd']:.3f}; "
              f"lowest f>{a.fmin}: d={r['dE_low']:+.3f} eV")


if __name__ == "__main__":
    main()
