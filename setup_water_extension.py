#!/usr/bin/env python3
"""
setup_water_extension.py

Stage A of the water-driven ensemble extension
(protonation_and_conformer_investigation.md Sec. 9.39/9.40): takes the
CREST conformers that were dropped by the 95% VACUUM GFN2 cut but are
needed to reach 95% of the GFN2/ALPB(water) population
(analyze_alpb_rerank.py), and builds r2SCAN-3c jobex directories for
them EXACTLY like the original Layer A (same template, same start
geometry = the vacuum CREST structure, same 10-cycle cap in
job_r2scan3c_optimize.sh), so their r2SCAN-3c energies are directly
comparable with the existing population_table_dft.csv.

Writes <state_dir>/layer_a/water_ext/{conformer_jobs/job_NNNN,
ensemble.xyz, map.csv (job_dir -> crest_index), joblist.txt}.

Usage:
    python setup_water_extension.py <state_dir> <crest_index> [<crest_index> ...]
"""

import os
import subprocess
import sys


def main():
    state_dir = os.path.abspath(sys.argv[1])
    idx = [int(x) for x in sys.argv[2:]]
    with open(os.path.join(state_dir, "crest_conformers.xyz")) as f:
        lines = f.readlines()
    nat = int(lines[0].split()[0])
    out = os.path.join(state_dir, "layer_a", "water_ext")
    os.makedirs(out)
    with open(os.path.join(out, "ensemble.xyz"), "w") as f:
        for i in idx:
            f.write("".join(lines[(i - 1) * (nat + 2): i * (nat + 2)]))
    here = os.path.dirname(os.path.abspath(__file__))
    fan = os.path.join(here, "scripts", "fanout_conformer_jobs.py")
    if not os.path.isfile(fan):
        fan = os.path.join(here, "fanout_conformer_jobs.py")
    subprocess.run(["python3", fan, "--template", os.path.join(state_dir, "layer_a", "template"),
                    "--ensemble", os.path.join(out, "ensemble.xyz"),
                    "--output-dir", os.path.join(out, "conformer_jobs")], check=True)
    with open(os.path.join(out, "map.csv"), "w") as f:
        f.write("job_dir,crest_index\n")
        for k, i in enumerate(idx, start=1):
            f.write(f"job_{k:04d},{i}\n")
    with open(os.path.join(out, "joblist.txt"), "w") as f:
        for k in range(1, len(idx) + 1):
            f.write(os.path.join(out, "conformer_jobs", f"job_{k:04d}") + "\n")
    print("done:", out)


if __name__ == "__main__":
    main()
