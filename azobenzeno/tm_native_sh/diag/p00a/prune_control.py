#!/usr/bin/env python3
"""Poda $nac e $excitation do control mantendo só os últimos K ciclos (numeração original preservada). Uso: prune_control.py in out [K=2]"""
import re, sys
src, dst = sys.argv[1], sys.argv[2]
K = int(sys.argv[3]) if len(sys.argv) > 3 else 2
lines = open(src).read().split("\n")
out, i = [], 0
while i < len(lines):
    l = lines[i]
    if l.startswith("$excitation") or l.startswith("$nac "):
        j = i + 1
        while j < len(lines) and not lines[j].startswith("$"):
            j += 1
        body = lines[i + 1:j]
        starts = [k for k, b in enumerate(body) if re.match(r"\s*cycle\s*=", b)]
        keep = starts[-K:] if len(starts) > K else starts
        out.append(l)
        if keep:
            out.extend(body[keep[0]:])
        i = j
        continue
    out.append(l); i += 1
open(dst, "w").write("\n".join(out))
