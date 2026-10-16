#!/usr/bin/env python3
"""Report how close capped TURBOMOLE jobex optimizations are to convergence.

For every job directory (one per conformer) it reads the `energy` file and the
`CONVERGENCE INFORMATION` block of `job.last`, and reports: status marker
(GEO_OPT_CONVERGED / GEO_OPT_FAILED / none), cycles run, |dE| over the last 1
and last 5 cycles, total drift from cycle 1, MAX gradient and how many of the
5 statpt criteria are met. Use it to decide whether a run that hit the cycle
cap (`jobex -c N`) must be continued before its energies feed Boltzmann weights.

Usage:
    python3 reopt_convergence_report.py reopt/neutral reopt/dianion
    python3 reopt_convergence_report.py reopt/* --csv reopt_status.csv --tol-kcal 0.05
"""
import argparse
import csv
import glob
import os
import re
import statistics

HARTREE_TO_KCAL = 627.5095
CRITERIA = ("Energy change", "RMS of displacement", "RMS of gradient",
            "MAX displacement", "MAX gradient")


def read_energies(path):
    """Return the list of total energies (Eh) from a TURBOMOLE `energy` file."""
    energies = []
    with open(path) as fh:
        for line in fh:
            if re.match(r"^\s*\d+\s+-?\d", line):
                energies.append(float(line.split()[1]))
    return energies


def read_criteria(path):
    """Return (max_gradient, n_criteria_met) from the last CONVERGENCE INFORMATION block."""
    if not os.path.isfile(path):
        return float("nan"), 0
    block = open(path, errors="replace").read().split("CONVERGENCE INFORMATION")[-1]
    met = sum(bool(re.search(rf"^\s+{c}\s+yes", block, re.M)) for c in CRITERIA)
    m = re.search(r"MAX gradient\s+(?:yes|no)\s+(\S+)", block)
    return (float(m.group(1)) if m else float("nan")), met


def job_status(job_dir):
    for marker in ("GEO_OPT_CONVERGED", "GEO_OPT_FAILED", "GEO_OPT_RUNNING"):
        if os.path.isfile(os.path.join(job_dir, marker)):
            return marker.replace("GEO_OPT_", "").lower()
    return "none"


def analyse(job_dir):
    energies = read_energies(os.path.join(job_dir, "energy"))
    if len(energies) < 2:
        return None
    gmax, met = read_criteria(os.path.join(job_dir, "job.last"))
    k = min(5, len(energies) - 1)
    return {
        "job_dir": job_dir,
        "status": job_status(job_dir),
        "cycles": len(energies),
        "dE_last_kcal": abs(energies[-1] - energies[-2]) * HARTREE_TO_KCAL,
        "dE_last5_kcal": abs(energies[-1] - energies[-1 - k]) * HARTREE_TO_KCAL,
        "drift_kcal": (energies[-1] - energies[0]) * HARTREE_TO_KCAL,
        "max_grad": gmax,
        "criteria_met": met,
        "energy_Eh": energies[-1],
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("groups", nargs="+", help="directories holding one job dir per conformer")
    ap.add_argument("--pattern", default="c*", help="glob for job dirs inside each group (default: c*)")
    ap.add_argument("--tol-kcal", type=float, default=0.05,
                    help="|dE| over the last 5 cycles regarded as 'flat' (default 0.05 kcal/mol)")
    ap.add_argument("--csv", help="optional per-job CSV output")
    args = ap.parse_args()

    all_rows = []
    for group in args.groups:
        rows = [r for d in sorted(glob.glob(os.path.join(group, args.pattern)))
                if os.path.isdir(d) and (r := analyse(d))]
        if not rows:
            print(f"{group}: no job dirs with an energy file")
            continue
        all_rows += rows
        count = lambda f: sum(1 for r in rows if f(r))
        print(f"{group}: n={len(rows)} converged={count(lambda r: r['status'] == 'converged')} "
              f"failed(cap)={count(lambda r: r['status'] == 'failed')} | "
              f"flat dE5<{args.tol_kcal}: {count(lambda r: r['dE_last5_kcal'] < args.tol_kcal)} | "
              f"median dE5={statistics.median(r['dE_last5_kcal'] for r in rows):.3f} "
              f"max dE5={max(r['dE_last5_kcal'] for r in rows):.3f} kcal/mol | "
              f"median drift={statistics.median(r['drift_kcal'] for r in rows):.2f} kcal/mol")
    if args.csv and all_rows:
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(all_rows[0]))
            w.writeheader()
            w.writerows(all_rows)
        print(f"wrote {args.csv}")


if __name__ == "__main__":
    main()
