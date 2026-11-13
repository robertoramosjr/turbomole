#!/usr/bin/env python3
"""riper_bandgap.py -- band gap and occupation diagnostics from riper's riper_bands.txt.

The 'Band gap' line in riper.out/job.last is unreliable with smearing (it reported ~0.001 Eh for an
insulator with a 1.9 eV gap). This reads the per-k band energies/occupations instead.

Usage:
  riper_bandgap.py runs/pbesol_sp/riper_bands.txt [runs/other/riper_bands.txt ...] [--nocc 192] [--per-k]

Reports VBM/CBM (by occupation > 1 for closed shells), indirect and minimum direct gap, the k-points
where they occur, and every state with fractional occupation (thr..2-thr), which flags a Fermi level
sitting inside a band. With --nocc it also reports the gap between band nocc and nocc+1 regardless
of occupations (the 'aufbau' gap).
"""
import argparse

HARTREE_EV = 27.211386245988


def read_bands(path):
    kpts, cur = [], None
    for line in open(path):
        s = line.split()
        if line.startswith("k-point"):
            cur = (tuple(float(x) for x in s[1:4]), [])
            kpts.append(cur)
        elif cur is not None and len(s) == 2:
            try:
                cur[1].append((float(s[0]), float(s[1])))
            except ValueError:
                pass
    return kpts


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("files", nargs="+")
    p.add_argument("--nocc", type=int, help="number of doubly occupied bands for the aufbau gap")
    p.add_argument("--thr", type=float, default=0.01, help="fractional-occupation threshold")
    p.add_argument("--per-k", action="store_true", help="print HOMO/LUMO per k-point")
    a = p.parse_args()
    for f in a.files:
        ks = read_bands(f)
        occ = [(e, k) for k, b in ks for e, o in b if o > 1.0]
        un = [(e, k) for k, b in ks for e, o in b if o <= 1.0]
        vbm, kv = max(occ)
        cbm, kc = min(un)
        direct = min((min(e for e, o in b if o <= 1.0) - max(e for e, o in b if o > 1.0), k) for k, b in ks)
        frac = [(k, i + 1, e, o) for k, b in ks for i, (e, o) in enumerate(b) if a.thr < o < 2 - a.thr]
        nocc_k = sorted({sum(1 for e, o in b if o > 1.0) for k, b in ks})
        print(f"== {f}  (nk={len(ks)}, nbands={len(ks[0][1])}, occupied bands per k: {nocc_k})")
        print(f"   VBM {vbm*HARTREE_EV:8.3f} eV at {kv}   CBM {cbm*HARTREE_EV:8.3f} eV at {kc}")
        print(f"   gap (occupation) {(cbm-vbm)*HARTREE_EV:.3f} eV ; min direct {direct[0]*HARTREE_EV:.3f} eV at {direct[1]}")
        if a.nocc:
            v = max(b[a.nocc - 1][0] for k, b in ks)
            c = min(b[a.nocc][0] for k, b in ks)
            print(f"   aufbau gap (band {a.nocc}->{a.nocc+1}) {(c-v)*HARTREE_EV:.3f} eV")
        print(f"   fractional occupations: {len(frac)}")
        for k, i, e, o in frac[:10]:
            print(f"     k={k} band {i} E={e*HARTREE_EV:.3f} eV occ={o:.4f}")
        if a.per_k:
            for k, b in ks:
                h = max(e for e, o in b if o > 1.0)
                l = min(e for e, o in b if o <= 1.0)
                print(f"     {k}  HOMO {h*HARTREE_EV:8.3f}  LUMO {l*HARTREE_EV:8.3f}  dir {(l-h)*HARTREE_EV:.3f}")


if __name__ == "__main__":
    main()
