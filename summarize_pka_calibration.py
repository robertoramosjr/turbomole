#!/usr/bin/env python3
"""
summarize_pka_calibration.py

Multi-reference version of calibrate_pka_pair.py
(protonation_and_conformer_investigation.md Sec. 9.34/9.35): instead of
one residual per pair, groups the reference compounds by CLASS (the
kind of deprotonation step, e.g. "0to-1" / "-1to-2"), reports each
class's mean offset +/- sample standard deviation, and applies the
class mean to each artepilin C target pair of that class, with the
reference SCATTER as its uncertainty. A class with a single reference
gets no scatter (reported as N=1, uncertainty unknown -- never faked).

Also propagates the spread of the experimental pKa itself when a
reference has a literature range (exp_pka_low/exp_pka_high): the class
mean is evaluated at the low, nominal and high experimental values and
the extremes are reported as a separate "experimental-range" term, so
the two sources of imprecision (chemical-class scatter vs. reference
data uncertainty) are not mixed.

--refs CSV columns:   name,class,acid_dir,base_dir,exp_pka,exp_pka_low,exp_pka_high
--targets CSV columns: name,class,acid_thermo,base_thermo,sconf_shift_kcal
(sconf_shift_kcal = ensemble conformational free energy difference
base-acid, Sec. 9.34; applied only with --apply-sconf.)

Usage:
    python ~/work_turbomole/scripts/summarize_pka_calibration.py \
        --refs pka_refs.csv --targets pka_targets.csv \
        --output pka_calibration_summary.dat
"""

import argparse
import csv
import math
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compute_pka_thermodynamic_cycle import compute_cycle, R_KCAL_MOL_K, T_KELVIN  # noqa: E402
from calibrate_pka_pair import read_reference_species, read_state_thermo  # noqa: E402


def sd(values):
    return statistics.stdev(values) if len(values) >= 2 else float("nan")


def fmt(x, nd=3):
    return "nan" if x != x else f"{x:.{nd}f}"


def main():
    parser = argparse.ArgumentParser(
        description="Class-averaged pKa calibration from several reference "
                    "compounds, applied to the artepilin C pairs.")
    parser.add_argument("--refs", required=True)
    parser.add_argument("--targets", required=True)
    parser.add_argument("--apply-sconf", action="store_true",
                         help="add each target's ensemble conformational "
                              "free-energy shift (sensitivity, off by default)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    refs = []
    with open(args.refs) as f:
        for r in csv.DictReader(f):
            a = read_reference_species(r["acid_dir"])
            b = read_reference_species(r["base_dir"])
            raw = compute_cycle(a[0], a[1], a[2], b[0], b[1], b[2])[-1]
            exp = float(r["exp_pka"])
            lo = float(r["exp_pka_low"] or exp)
            hi = float(r["exp_pka_high"] or exp)
            refs.append({"name": r["name"], "cls": r["class"], "raw": raw,
                          "res": exp - raw, "res_lo": lo - raw, "res_hi": hi - raw,
                          "exp": exp})

    classes = {}
    for r in refs:
        classes.setdefault(r["cls"], []).append(r)

    lines = ["# per-reference residuals (pKa_exp - pKa_raw)"]
    print("=== references ===")
    for r in refs:
        s = (f"{r['name']:28s} class={r['cls']:7s} raw={r['raw']:.3f} exp={r['exp']:.3f} "
             f"residual={r['res']:+.3f} pKa ({r['res'] * R_KCAL_MOL_K * T_KELVIN * math.log(10):+.3f} kcal/mol)")
        print(s); lines.append(s)

    stats = {}
    print("\n=== class offsets ===")
    for cls, rs in classes.items():
        vals = [r["res"] for r in rs]
        mean = statistics.mean(vals)
        mean_lo = statistics.mean([r["res_lo"] for r in rs])
        mean_hi = statistics.mean([r["res_hi"] for r in rs])
        stats[cls] = (mean, sd(vals), len(rs), min(mean_lo, mean_hi), max(mean_lo, mean_hi))
        s = (f"class {cls}: N={len(rs)} mean offset = {mean:+.3f}, scatter (sample sd) = "
             f"{fmt(sd(vals))}, offset range from experimental-pKa spread = "
             f"[{min(mean_lo, mean_hi):+.3f}, {max(mean_lo, mean_hi):+.3f}]")
        print(s); lines.append(s)

    lines.append("# targets")
    print("\n=== artepilin C targets ===")
    with open(args.targets) as f:
        for t in csv.DictReader(f):
            ga, sa, _ = read_state_thermo(t["acid_thermo"])
            gb, sb, _ = read_state_thermo(t["base_thermo"])
            raw = compute_cycle(ga, ga + sa, 0.0, gb, gb + sb, 0.0)[-1]
            if args.apply_sconf and t.get("sconf_shift_kcal"):
                raw += float(t["sconf_shift_kcal"]) / (R_KCAL_MOL_K * T_KELVIN * math.log(10))
            mean, scatter, n, lo, hi = stats[t["class"]]
            s = (f"{t['name']:14s} class={t['class']:7s} raw={raw:.3f} calibrated={raw + mean:.3f} "
                 f"+/- {fmt(scatter, 2)} (class scatter, N={n}); experimental-range term "
                 f"[{raw + lo:.3f}, {raw + hi:.3f}]")
            print(s); lines.append(s)

    with open(args.output, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
