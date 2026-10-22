#!/usr/bin/env python3
"""TDOS from a TURBOMOLE run via tm2molden + Multiwfn (noGUI, batch).

Acts on the current working directory (needs control, coord, basis and mos, or
alpha/beta for UHF). Writes dos_tdos.dat, dos_levels.dat and dos_gap.dat.
"""
import argparse
import os
import subprocess
from pathlib import Path

import numpy as np


def run(cmd, stdin, log, env):
    with open(log, "w") as fh:
        subprocess.run(cmd, input=stdin, text=True, stdout=fh, stderr=subprocess.STDOUT, env=env, check=False)


HARTREE_EV = 27.211386


def uhf_levels(molden):
    """Per-spin (energies_eV, occupations) from a Molden file, or None if closed shell."""
    spins = {"alpha": ([], []), "beta": ([], [])}
    cur = None
    for l in open(molden):
        t = l.split()
        if l.lstrip().startswith("Ene="):
            e = float(t[1]) * HARTREE_EV
        elif l.lstrip().startswith("Spin="):
            cur = t[1].lower()
        elif l.lstrip().startswith("Occup=") and cur:
            spins[cur][0].append(e)
            spins[cur][1].append(float(t[1]))
    if not spins["beta"][0]:
        return None
    return {k: (np.array(v[0]), np.array(v[1])) for k, v in spins.items()}


def write_uhf(a, lv):
    """Per-spin levels, HOMO/LUMO and Gaussian-broadened alpha/beta DOS (UHF)."""
    grid = np.arange(a.emin, a.emax + a.estep / 2, a.estep)
    sigma = a.fwhm / 2.3548200450309493
    homo, lumo = {}, {}
    total = np.zeros_like(grid)
    for k, (e, occ) in lv.items():
        order = np.argsort(e)
        e, occ = e[order], occ[order]
        np.savetxt(f"dos_levels_{k}.dat", e, fmt="%.5f", header=f"mo_energy_eV ascending, {k}")
        nocc = int(round(occ.sum()))
        homo[k], lumo[k] = e[nocc - 1], e[nocc]
        dos = np.exp(-0.5 * ((grid[:, None] - e[None, :]) / sigma) ** 2).sum(axis=1)
        total += dos
        np.savetxt(f"dos_{k}.dat", np.c_[grid, dos], fmt="%.6g",
                   header=f"energy_eV dos_{k}_arb (Gaussian FWHM={a.fwhm} eV, unit-height levels)")
        print(f"{k}: nocc={nocc} HOMO={homo[k]:.3f} eV LUMO={lumo[k]:.3f} eV")
    # Multiwfn's TDOS is unreliable for UHF molden: rebuild it, on its RHF scale (0.3133 per level; /2 so alpha+beta ~ one spatial set)
    np.savetxt("dos_tdos.dat", np.c_[grid, 0.3133 * total / 2], fmt="%.6g",
               header=f"energy_eV tdos_arb (Gaussian FWHM={a.fwhm} eV, (alpha+beta)/2 from Molden levels)")
    h, l = max(homo.values()), min(lumo.values())
    np.savetxt("dos_gap.dat", [[homo["alpha"], lumo["alpha"], homo["beta"], lumo["beta"], h, l, l - h]],
               fmt="%.5f", header="homo_a lumo_a homo_b lumo_b homo_eV lumo_eV gap_eV (global = max HOMO, min LUMO)")
    print(f"global HOMO={h:.3f} eV LUMO={l:.3f} eV gap={l - h:.3f} eV")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mos", default="mos", help="TURBOMOLE MO file passed to tm2molden")
    ap.add_argument("--molden", default="mol.molden")
    ap.add_argument("--tm-env", type=Path, default=Path.home() / "software/turbomole_env.sh")
    ap.add_argument("--mwfn-env", type=Path, default=Path.home() / "software/multiwfn_env.sh")
    ap.add_argument("--mwfn-exe", default="Multiwfn_noGUI")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--emin", type=float, default=-30.0, help="eV")
    ap.add_argument("--emax", type=float, default=5.0, help="eV")
    ap.add_argument("--estep", type=float, default=0.02, help="eV")
    ap.add_argument("--fwhm", type=float, default=0.3, help="eV")
    ap.add_argument("--nocc", type=int, default=0, help="occupied MOs (0 = parse from Multiwfn log)")
    a = ap.parse_args()
    cwd = Path.cwd()

    def shell(cmd, stdin, log):
        script = f"source {a.tm_env}; source {a.mwfn_env}; {cmd}"
        with open(log, "w") as fh:
            subprocess.run(["bash", "-c", script], input=stdin, text=True, stdout=fh,
                           stderr=subprocess.STDOUT, check=False)

    if not Path(a.molden).exists():
        shell("tm2molden", "\n\n\n\n", "tm2molden.log")
        src = cwd / "molden.input"
        if not src.exists():
            raise SystemExit("tm2molden failed, see tm2molden.log")
        src.rename(a.molden)
        for junk in ("gradient",):
            if (cwd / junk).stat().st_size == 0:
                (cwd / junk).unlink()
    # Multiwfn menu 10 (DOS): 8 = eV, 2 = range+step, 3 = FWHM, 0 = plot, 3 = export curve/lines
    keys = f"{a.molden}\n10\n8\n2\n{a.emin} {a.emax} {a.estep}\n3\n{a.fwhm}\n0\n3\n0\n-10\nq\n"
    shell(f"{a.mwfn_exe} -nt {a.threads}", keys, "dos_multiwfn.log")
    if not (cwd / "DOS_curve.txt").exists():
        raise SystemExit("Multiwfn produced no DOS_curve.txt, see dos_multiwfn.log")

    curve = np.loadtxt("DOS_curve.txt")
    np.savetxt("dos_tdos.dat", curve, fmt="%.6g",
               header=f"energy_eV tdos_arb (Gaussian FWHM={a.fwhm} eV, Multiwfn TDOS)")
    lines = np.loadtxt("DOS_line.txt")
    levels = np.unique(np.round(lines[:, 0], 5))
    np.savetxt("dos_levels.dat", levels, fmt="%.5f", header="mo_energy_eV (ascending)")
    if uhf_levels(a.molden) is not None:
        write_uhf(a, uhf_levels(a.molden))
        for f in ("DOS_curve.txt", "DOS_line.txt"):
            Path(f).unlink()
        return
    nocc = a.nocc
    if not nocc:
        for l in open("dos_multiwfn.log"):
            if "are occupied" in l:
                nocc = int(l.split("to")[1].split()[0])
    homo, lumo = levels[nocc - 1], levels[nocc]
    np.savetxt("dos_gap.dat", [[homo, lumo, lumo - homo]], fmt="%.5f", header="homo_eV lumo_eV gap_eV")
    for f in ("DOS_curve.txt", "DOS_line.txt"):
        Path(f).unlink()
    print(f"nocc={nocc} HOMO={homo:.3f} eV LUMO={lumo:.3f} eV gap={lumo - homo:.3f} eV")


if __name__ == "__main__":
    main()
