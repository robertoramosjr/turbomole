#!/usr/bin/env python3
"""
seed_sh_log_groups.py

Pre-seeds the four SH log groups ($population_el, $prob_fssh, $sh_coeffs, $egradmon) at the top of a TURBOMOLE
`control` with a chosen `maxlog`. The TM creates them itself on the first MD step with maxlog=256, and the roll-over
of these logs at step 256->257 truncates the control (P00a); the `maxlog` in mdmaster's $log_history does NOT
govern them. If the groups already exist only their `maxlog=` value is changed. Idempotent.

Usage: python3 seed_sh_log_groups.py <control> --maxlog N
"""
import argparse
import re

GROUPS = ("population_el", "prob_fssh", "sh_coeffs", "egradmon")


def block(g, n):
    return f"${g} log=true\n logbase={g}log\n hist=none\n current=  1\n maxlog={n:4d}\n"


def seed(text, n):
    for g in reversed(GROUPS):  # prepend => final order = GROUPS order, as the TM writes it
        m = re.search(rf"^\${g} log=true[^\n]*\n(?:(?! *\$)[^\n]*\n)*?( maxlog=[^\n]*\n)", text, re.M)
        if m:
            text = text[:m.start(1)] + f" maxlog={n:4d}\n" + text[m.end(1):]
        else:
            text = block(g, n) + text
    return text


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("control")
    ap.add_argument("--maxlog", type=int, required=True)
    a = ap.parse_args()
    txt = open(a.control).read()
    with open(a.control, "w") as f:
        f.write(seed(txt, a.maxlog))


if __name__ == "__main__":
    main()
