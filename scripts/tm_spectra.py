#!/usr/bin/env python3
"""Broaden TURBOMOLE IR (vibspectrum) and UV-vis (exspectrum) stick spectra.

Acts on the current working directory. Writes ir_sticks.dat, ir_broadened.dat,
uvvis_sticks.dat and uvvis_broadened.dat (header + units) for plotting.
"""
import argparse
import re
from pathlib import Path

import numpy as np

HC_EV_NM = 1239.841984  # eV*nm


def parse_vibspectrum(path):
    """Return arrays (wavenumber cm-1, IR intensity km/mol, IR-active flag)."""
    rows = []
    for line in path.read_text().splitlines():
        if line.startswith(("$", "#")) or not line.strip():
            continue
        parts = line.split()
        if parts[0].isdigit() and len(parts) >= 4:
            if parts[1].isalpha():  # symmetry label present
                wn, inten, ir = float(parts[2]), float(parts[3]), parts[4]
            else:
                wn, inten, ir = float(parts[1]), float(parts[2]), parts[3]
            rows.append((wn, inten, ir.upper() == "YES"))
    if not rows:
        raise SystemExit(f"No modes found in {path}")
    return np.array([r[0] for r in rows]), np.array([r[1] for r in rows])


def parse_exspectrum(path):
    """Return arrays (energy eV, wavelength nm, osc. strength length gauge, velocity gauge)."""
    ev, nm, f_len, f_vel = [], [], [], []
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        p = line.replace("D", "E").split()
        if p[0].isdigit():
            off = 0 if p[1].isalpha() else -1  # irrep label shifts columns
            ev.append(float(p[3 + off]))
            nm.append(float(p[5 + off]))
            f_vel.append(float(p[6 + off]))
            f_len.append(float(p[7 + off]))
    if not ev:
        raise SystemExit(f"No excitations found in {path}")
    return np.array(ev), np.array(nm), np.array(f_len), np.array(f_vel)


def lorentz(x, x0, w):
    return (w / 2) ** 2 / ((x - x0) ** 2 + (w / 2) ** 2)


def gauss(x, x0, w):
    sigma = w / (2 * np.sqrt(2 * np.log(2)))
    return np.exp(-((x - x0) ** 2) / (2 * sigma**2))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vibspectrum", type=Path, default=Path("vibspectrum"))
    ap.add_argument("--exspectrum", type=Path, default=Path("exspectrum"))
    ap.add_argument("--ir-fwhm", type=float, default=20.0, help="IR FWHM (cm-1)")
    ap.add_argument("--ir-min", type=float, default=0.0, help="IR axis min (cm-1)")
    ap.add_argument("--ir-max", type=float, default=4000.0, help="IR axis max (cm-1)")
    ap.add_argument("--ir-step", type=float, default=1.0, help="IR grid step (cm-1)")
    ap.add_argument("--ir-scale", type=float, default=1.0, help="frequency scaling factor")
    ap.add_argument("--uv-fwhm", type=float, default=0.3, help="UV-vis FWHM (eV)")
    ap.add_argument("--uv-min-nm", type=float, default=200.0)
    ap.add_argument("--uv-max-nm", type=float, default=800.0)
    ap.add_argument("--uv-step-nm", type=float, default=0.5)
    ap.add_argument("--uv-gauge", choices=["length", "velocity"], default="length")
    ap.add_argument("--uv-shape", choices=["gauss", "lorentz"], default="gauss")
    ap.add_argument("--outdir", type=Path, default=Path.cwd())
    a = ap.parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)

    # --- IR ---
    if a.vibspectrum.exists():
        wn, inten = parse_vibspectrum(a.vibspectrum)
        keep = wn > 1.0  # drop translations/rotations (0.00)
        wn, inten = wn[keep] * a.ir_scale, inten[keep]
        np.savetxt(a.outdir / "ir_sticks.dat", np.c_[wn, inten],
                   header="wavenumber_cm-1 ir_intensity_km_mol", fmt="%.4f")
        x = np.arange(a.ir_min, a.ir_max + a.ir_step, a.ir_step)
        y = sum(i * lorentz(x, w, a.ir_fwhm) for w, i in zip(wn, inten))
        np.savetxt(a.outdir / "ir_broadened.dat", np.c_[x, y],
                   header=f"wavenumber_cm-1 ir_intensity_km_mol (Lorentz FWHM={a.ir_fwhm} cm-1, scale={a.ir_scale})",
                   fmt="%.6g")
        print(f"IR: {len(wn)} modes, strongest {wn[inten.argmax()]:.1f} cm-1 ({inten.max():.1f} km/mol)")

    # --- UV-vis ---
    if a.exspectrum.exists():
        ev, nm, f_len, f_vel = parse_exspectrum(a.exspectrum)
        f = f_len if a.uv_gauge == "length" else f_vel
        np.savetxt(a.outdir / "uvvis_sticks.dat", np.c_[nm, ev, f],
                   header=f"wavelength_nm energy_eV oscillator_strength_{a.uv_gauge}", fmt="%.6g")
        wl = np.arange(a.uv_min_nm, a.uv_max_nm + a.uv_step_nm, a.uv_step_nm)
        e = HC_EV_NM / wl
        shape = gauss if a.uv_shape == "gauss" else lorentz
        y = sum(fi * shape(e, ei, a.uv_fwhm) for ei, fi in zip(ev, f))
        np.savetxt(a.outdir / "uvvis_broadened.dat", np.c_[wl, e, y],
                   header=f"wavelength_nm energy_eV intensity_arb ({a.uv_shape} FWHM={a.uv_fwhm} eV, f_{a.uv_gauge})",
                   fmt="%.6g")
        ok = (nm >= a.uv_min_nm) & (nm <= a.uv_max_nm)
        print(f"UV-vis: {len(ev)} states ({ok.sum()} in window); lowest-energy bright (f>0.05): "
              f"{[(round(n,1), round(fi,3)) for n, fi in zip(nm, f) if fi > 0.05][-3:]}")


if __name__ == "__main__":
    main()
