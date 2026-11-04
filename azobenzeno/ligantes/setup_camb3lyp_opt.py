#!/usr/bin/env python3
"""Monta uma pasta CAM-B3LYP-D3(BJ)/def2-TZVP + $senex pronta para jobex -> ridft -> NumForce.

Baseado em scripts/setup_hse06_tzvp_template.py (mesmo define automatizado), com: functional
cam-b3lyp, gridsize 5, sem bloco escf, $scfconv 10, $denconv 1d-10, $scfiterlimit 1500,
$forceconv 7, $ricore 4000 e $statpt apertado do passo 1a (thrmaxgrad 2e-5: 1e-5 é piso de
ruído da grade 5). Requer TURBODIR no ambiente (define no PATH) e pexpect.
Uso: setup_camb3lyp_opt.py --coord anchored/trans_sh/coord --charge 0 --output-dir lig_trans_sh
"""
import argparse
import os
import re
import shutil
import time

import pexpect

STATPT = ("$statpt\n   threchange 1.0d-9\n   thrmaxgrad 2.0d-5\n   thrrmsgrad 5.0d-6\n"
          "   thrmaxdispl 1.0d-4\n   thrrmsdispl 1.0d-4\n")
EXTRA = ("$senex\n$rij\n$disp3 bj\n$scfconv 10\n$denconv 1d-10\n$scfiterlimit 1500\n"
         "$ricore 4000\n$forceconv 7\n" + STATPT)


def run_define(d, charge, grid, log):
    child = pexpect.spawn("define", cwd=d, timeout=20, encoding="utf-8")
    with open(log, "w") as lf:
        child.logfile = lf

        def send(s, w=0.8):
            child.sendline(s)
            time.sleep(w)
        for s in ("", "", "a coord", "*", "no"):
            send(s)
        send("b all def2-TZVP", 2.0)
        for s in ("*",):
            send(s, 1.0)
        for s in ("eht", "", str(charge), "", "dft", "func cam-b3lyp", f"grid {grid}", "on", ""):
            send(s)
        send("*", 2.0)
        for _ in range(5):
            i = child.expect(["GENERAL MENU", pexpect.TIMEOUT], timeout=3)
            if i == 0:
                break
            child.sendline("")
        send("*", 2.0)
        try:
            child.expect(pexpect.EOF, timeout=25)
        finally:
            child.close()


ap = argparse.ArgumentParser()
ap.add_argument("--coord", required=True)
ap.add_argument("--charge", type=int, default=0)
ap.add_argument("--grid", default="5")
ap.add_argument("--output-dir", required=True)
a = ap.parse_args()
if os.path.exists(a.output_dir):
    raise SystemExit(f"{a.output_dir} já existe")
os.makedirs(a.output_dir)
shutil.copy(a.coord, os.path.join(a.output_dir, "coord"))
run_define(a.output_dir, a.charge, a.grid, os.path.join(a.output_dir, "define_session.log"))
cp = os.path.join(a.output_dir, "control")
t = open(cp).read()
if "cam-b3lyp" not in t:
    raise SystemExit("control sem cam-b3lyp: ver define_session.log")
t = re.sub(r"^\$(scfconv|denconv)\s.*\n", "", t, flags=re.M)
t = t.replace("$end", EXTRA + "$end", 1)
open(cp, "w").write(t)
print("ok", a.output_dir)
