#!/usr/bin/env python3
"""
aggregate_sh.py

Aggregates per-trajectory CSVs from summarize_sh_traj.py into the fraction of trajectories that switched isomer,
with a Wilson 95 % interval. The denominator is ONLY trajectories with status hop_ok; no_hop / died_S1 / died_S0
are counted separately. This is a switching fraction of the model, NOT a quantum yield.

Usage:
    python3 aggregate_sh.py --from cis|trans <summary.csv> [<summary.csv> ...] [--out table.md]
"""
import argparse
import csv
import math


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (c - h, c + h)


def aggregate(rows, frm):
    seen = set()
    for r in rows:
        t = r.get("traj")
        if t is not None:
            if t in seen:
                raise ValueError(f"duplicate traj: {t}")
            seen.add(t)
    rows = [r for r in rows if r["status"] != "not_run"]
    ok = [r for r in rows if r["status"] == "hop_ok"]
    n_cis = sum(r["product"] == "cis" for r in ok)
    n_trans = sum(r["product"] == "trans" for r in ok)
    k = n_trans if frm == "cis" else n_cis
    out = {"n_total": len(rows), "n_hop_ok": len(ok),
           "n_no_hop": sum(r["status"] == "no_hop" for r in rows),
           "n_died_S1": sum(r["status"] == "died_S1" for r in rows),
           "n_died_S0": sum(r["status"] == "died_S0" for r in rows),
           "n_cis": n_cis, "n_trans": n_trans,
           "frac_switch": None, "ci_low": None, "ci_high": None}
    if ok:
        out["frac_switch"] = k / len(ok)
        out["ci_low"], out["ci_high"] = wilson(k, len(ok))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csvs", nargs="+")
    ap.add_argument("--from", dest="frm", required=True, choices=["cis", "trans"])
    ap.add_argument("--out")
    a = ap.parse_args()
    lines = ["| set | n | hop_ok | no_hop | died_S1 | died_S0 | cis | trans | switch frac | Wilson 95 % |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for path in a.csvs:
        rows = [{**r, "product": r["product"] or None} for r in csv.DictReader(open(path))]
        g = aggregate(rows, a.frm)
        f = "-" if g["frac_switch"] is None else f"{g['frac_switch']:.2f}"
        ci = "-" if g["ci_low"] is None else f"{g['ci_low']:.2f}–{g['ci_high']:.2f}"
        lines.append(f"| {path} | {g['n_total']} | {g['n_hop_ok']} | {g['n_no_hop']} | {g['n_died_S1']} | "
                     f"{g['n_died_S0']} | {g['n_cis']} | {g['n_trans']} | {f} | {ci} |")
    txt = "\n".join(lines) + "\n"
    print(txt)
    if a.out:
        open(a.out, "w").write(txt)


if __name__ == "__main__":
    main()
