#!/usr/bin/env python3
"""
prepare_tm_s0_continuation_v2.py

Builds a pure GROUND-STATE molecular-dynamics continuation (TURBOMOLE frog/jobex -md, no surface hopping, no
TDDFT) from the LAST frame of a native TM surface-hopping run, to follow a trajectory after its S1->S0 hop and
classify the photoproduct. Positions/velocities/time come from the last frame of `mdlog.1` of the source
directory; the electronic setup (basis, auxbasis, control) is copied and stripped of every excited-state /
SH keyword; the SCF settings are replaced by a robust profile (the plain ones failed right after the hop,
close to the S1/S0 crossing).

v2 (2026-10-16): also drops `$actual step` (stale marker from a phase-1 ridft/dscf that died; it made the S0
continuation of tests t3/t9 die in 2 s with "longjmp ... ridft_omp terminated").

Usage:
    python prepare_tm_s0_continuation_v2.py <sh_dir> <out_dir> [--nsteps 300] [--scf-profile strong] [--natoms 24]
Then in <out_dir>:  jobex -ri -md -c <nsteps+5>
"""
import argparse
import os
import re
import shutil

PROFILES = {
    "iter300":    ("300", "start=0.300  step=0.050  min=0.100", "automatic=.1", "7", False),
    "baseline":   ("100", "start=0.300  step=0.050  min=0.100", "automatic=.1", "7", False),
    "damp_shift": ("300", "start=1.000  step=0.050  min=0.200", "automatic=0.4", "7", False),
    "strong":     ("500", "start=1.500  step=0.100  min=0.300", "automatic=0.8", "7", False),
    "fermi":      ("300", "start=0.300  step=0.050  min=0.100", "automatic=.1", "7", True),
}
DROP = (r"^\$(scfiterlimit|scfdamp|scforbitalshift|scfconv|scfinstab|soes|exopt|nacme|exgradient1|nacvector01|"
        r"excitation|currentstate|active|surface_hopping|sh_coeffs|population_el|prob_fssh|rpacor|rpaconv|denconv|"
        r"actual|grad|end|energy|fermi|gap_threshold|nac|subenergy|gsenergy|esenergy|couplingvector|esgrad|egradmon|last|xctype|dipole|charge)\b")


def last_frame(mdlog, nat):
    lines = open(mdlog).read().split("\n")
    idx = [i for i, l in enumerate(lines) if l.startswith("t=")]
    for i in reversed(idx):  # last COMPLETE frame (coords, dt line, velocities)
        blk = lines[i:i + 1 + nat + 1 + nat]
        body = blk[1:1 + nat] + blk[2 + nat:]  # coordinates and velocities (3 numbers); line nat+1 is dt
        if len(blk) == 1 + 2 * nat + 1 and all(len(l.split()) >= 3 for l in body):
            return blk
    raise SystemExit("no complete frame in " + mdlog)


def clean_control(path, profile):
    it, damp, shift, conv, fermi = PROFILES[profile]
    keep, skip = [], False
    for l in open(path).read().split("\n"):
        if l.startswith("$"):
            skip = bool(re.match(DROP, l))
        if skip:
            continue
        keep.append(l)
    txt = "\n".join(keep).rstrip() + "\n"
    txt += f"$scfiterlimit {it}\n$scfdamp   {damp}\n$scforbitalshift  {shift}\n$scfconv   {conv}\n"
    if fermi:
        txt += "$fermi tmstrt=300.00 tmend=50.00 tmfac=0.9 hlcrit=5.d-2 stop=1.d-3\n"
    return txt + "$grad file=gradient\n$energy file=energy\n$end\n"


def main():
    ap = argparse.ArgumentParser(description="Ground-state MD continuation from a native TM SH run.")
    ap.add_argument("sh_dir")
    ap.add_argument("out_dir")
    ap.add_argument("--nsteps", type=int, default=300)
    ap.add_argument("--scf-profile", default="iter300", choices=list(PROFILES))
    ap.add_argument("--natoms", type=int, default=24)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    for f in ("basis", "auxbasis", "mos"):
        shutil.copy(os.path.join(a.sh_dir, f), a.out_dir)
    frame = last_frame(os.path.join(a.sh_dir, "mdlog.1"), a.natoms)
    with open(os.path.join(a.out_dir, "mdlog.1"), "w") as f:
        f.write("# AIMD log file (S0 continuation)\n$log\n$current\n" + "\n".join(frame) + "\n$end\n")
    with open(os.path.join(a.out_dir, "coord"), "w") as f:
        f.write("$coord\n" + "\n".join(frame[1:1 + a.natoms]) + "\n$end\n")
    with open(os.path.join(a.out_dir, "control"), "w") as f:
        f.write(clean_control(os.path.join(a.sh_dir, "control"), a.scf_profile))
    dt = float(frame[1 + a.natoms].split()[0])
    with open(os.path.join(a.out_dir, "mdmaster"), "w") as f:
        f.write(f"""# A I M D : master file (ground-state continuation after S1->S0 hop)
$nsteps      {a.nsteps}
$natoms        {a.natoms}
$log file=mdlog.1
$current file=mdlog.1
$turbomole file=control
$seed                        179
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
  from t= {-dt:.10f}       until t=  {float(frame[0][2:]):.10f}
$md_action
$end
""")
    print(f"ok {a.out_dir}: t0={float(frame[0][2:]) / 41.341373:.2f} fs, nsteps={a.nsteps}, profile={a.scf_profile}")


if __name__ == "__main__":
    main()
