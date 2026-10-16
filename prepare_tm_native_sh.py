#!/usr/bin/env python3
"""Monta uma trajetoria de surface hopping NATIVO do TURBOMOLE (frog/jobex -md) a partir de um IC do Newton-X.

Uso: prepare_tm_native_sh.py <TRAJ_NX> <destino> [--gap EV] [--nsteps N] [--dt-fs 0.5] [--seed S]
  <TRAJ_NX> = pasta com geom (bohr, col: sym Z x y z massa), veloc (a.u.) e JOB_AD/{basis,auxbasis,mos,control}.
Gera: coord, control (CAM-B3LYP/SVP TDA, $surface_hopping, $nacme), mdmaster, mdlog.1, basis/auxbasis/mos.
Obs.: veloc do NX e' em t=0; frog espera t-dt/2 (erro O(dt/2) na velocidade inicial, aceito no teste).
"""
import argparse, re, shutil
from pathlib import Path

FS = 41.341373  # a.u. por fs

ap = argparse.ArgumentParser()
ap.add_argument("traj"); ap.add_argument("dest")
ap.add_argument("--gap", type=float, default=None, help="$gap_threshold em eV (omitido = default do TM)")
ap.add_argument("--nsteps", type=int, default=80)
ap.add_argument("--dt-fs", type=float, default=0.5)
ap.add_argument("--seed", type=int, default=179)
ap.add_argument("--scf-profile", default="iter300", choices=["iter300", "baseline", "damp_shift", "strong", "fermi"],
                help="SCF profile. NOTE: the post-hop SCF failure was caused by `$denconv` in the control (any value, even 1d-7, stalls ridft near the crossing); it is no longer written. Next-geometry test (2026-10-05): baseline (100 iterations) FAILS, iter300 converges in 114, damp_shift in 226; iter300 is the default.")
a = ap.parse_args()

src, dst = Path(a.traj), Path(a.dest)
dst.mkdir(parents=True, exist_ok=True)
geom = [l.split() for l in (src / "geom").read_text().splitlines() if l.strip()]
vel = [l.split() for l in (src / "veloc").read_text().splitlines() if l.strip()]
assert len(geom) == len(vel)
dt = a.dt_fs * FS

for f in ("basis", "auxbasis", "mos"):
    shutil.copy(src / "JOB_AD" / f, dst / f)

coord = "".join(f"{float(g[2]):20.12f}{float(g[3]):20.12f}{float(g[4]):20.12f}     {g[0].lower()}\n" for g in geom)
(dst / "coord").write_text("$coord\n" + coord + "$end\n")

ctrl = (src / "JOB_AD" / "control").read_text()
# remove blocos que o NX/egrad antigos deixaram (linha numerica orfa, $grad, $scfinstab/$soes/$rpacor/$denconv antigos)
out = []
for l in ctrl.splitlines():
    if re.match(r"^\s*-?\d+\.\d+\s+-?\d+\.\d+", l):  # linha numerica orfa
        continue
    if re.match(r"^\$(grad|end|scfinstab|soes|rpacor|denconv|dft|gridsize|scfiterlimit|scfdamp|scforbitalshift|scfconv|fermi)\b", l):
        continue
    if re.match(r"^\s+(functional|gridsize)\b", l) or re.match(r"^ a\s+3\s*$", l):
        continue
    out.append(l)
ctrl = "\n".join(out) + "\n"
SCF = {"iter300": ("300", "start=0.300  step=0.050  min=0.100", "automatic=.1", False),
       "baseline": ("100", "start=0.300  step=0.050  min=0.100", "automatic=.1", False),
       "damp_shift": ("300", "start=1.000  step=0.050  min=0.200", "automatic=0.4", False),
       "strong": ("500", "start=1.500  step=0.100  min=0.300", "automatic=0.8", False),
       "fermi": ("300", "start=0.300  step=0.050  min=0.100", "automatic=.1", True)}[a.scf_profile]
ctrl += f"$scfiterlimit {SCF[0]}\n$scfdamp   {SCF[1]}\n$scforbitalshift  {SCF[2]}\n$scfconv   7\n"
if SCF[3]:
    ctrl += "$fermi tmstrt=300.00 tmend=50.00 tmfac=0.9 hlcrit=5.d-2 stop=1.d-3\n"
ctrl += """$dft
   functional cam-b3lyp
   gridsize   m4
   weight derivatives
$scfinstab ciss
$soes
 a 1
$exopt all
$rpacor 4000
$rpaconv 8
$nacme
$grad file=gradient
$exgradient1 file=grad1
$nacvector01 file=nacv01
$excitation
$currentstate 1
$active
 cycle= 1
   1
$surface_hopping
$energy file=energy
$end
"""
(dst / "control").write_text(ctrl)

gap = f"$gap_threshold {a.gap}\n" if a.gap is not None else ""
(dst / "mdmaster").write_text(f"""# A I M D : master file (native TM SH test)
$nsteps      {a.nsteps}
$natoms        {len(geom)}
$surface_hopping
{gap}$log file=mdlog.1
$current file=mdlog.1
$turbomole file=control
$seed                        {a.seed}
$title
$log_history
  maxlog      256
       0                    mdlog.1
$ke_control
  length     50
  response      1
$barrier
$constraints
  algorithm   none
  tolerance  0.1000000E-01
  corpercyc  0.2500000
$md_status
  free
  from t= {-dt:.10f}       until t=  0.00000000000
$md_action
$sh_seed                              {a.seed}
$au_seed                              {a.seed}
$end
""")

lines = ["# AIMD log file", "$log", "$current", "t=   0.000000000"]
lines += [f"{float(g[2]):20.12f}{float(g[3]):20.12f}{float(g[4]):20.12f}     {g[0].lower()}" for g in geom]
lines.append(f"  {dt:.10f}       0.00000000000")
lines += [f"{g[0].lower():<2}{float(v[0]):22.12E}{float(v[1]):22.12E}{float(v[2]):22.12E}" for g, v in zip(geom, vel)]
lines.append("$end")
(dst / "mdlog.1").write_text("\n".join(lines) + "\n")
print(f"ok {dst}  natoms={len(geom)} dt={dt:.4f} a.u. nsteps={a.nsteps} gap={a.gap}")
