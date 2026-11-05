#!/usr/bin/env python3
"""Build the CAM-B3LYP range-separation diagnostic cases (aoforce vs finite-difference curvature).

Each case copies a converged CAM-B3LYP trans-azobenzeno minimum, cleans the control (no $grad/$hessian/
$vibrational/$energy/...), and changes ONE ingredient: omega (w in $xctype RSH-GGA), RI-J off, or $senex off.
Writes <outdir>/<case>/{control,coord,basis,auxbasis,mos,variant.txt} and <outdir>/cases.txt.

Usage (from the study folder):
  setup_rsh_diagnostic.py --base azo_trans_t4_grid6_aoforce --outdir rsh_diag --gridsize 5 \
      --case nosenex:nosenex --case norij:norij --case w010:w=0.10 --case w050:w=0.50 --case w020:w=0.20
"""
import argparse
import os
import re
import shutil

DROP = ('$grad', '$hessian', '$vibrational', '$optinfo', '$hessapprox', '$subenergy', '$energy', '$last',
        '$actual', '$orbital_max_rnorm', '$dipgrad', '$scfinstab', '$soes', '$nvibro', '$dipole', '$charge')


def clean_control(text):
    out, skip = [], False
    for line in text.splitlines():
        if line.startswith('$'):
            key = line.split()[0]
            skip = key in DROP and key != '$end'
        if not skip:
            out.append(line)
    return '\n'.join(out) + '\n'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--base', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--gridsize', default='5')
    ap.add_argument('--case', action='append', required=True, help='name:nosenex | name:norij | name:w=<omega>')
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    names = []
    for spec in a.case:
        name, var = spec.split(':', 1)
        d = os.path.join(a.outdir, name)
        os.makedirs(d, exist_ok=True)
        for f in ('coord', 'basis', 'auxbasis', 'mos'):
            shutil.copy(os.path.join(a.base, f), d)
        c = clean_control(open(os.path.join(a.base, 'control')).read())
        c = re.sub(r'(gridsize\s+)\S+', rf'\g<1>{a.gridsize}', c)
        if var == 'nosenex':
            c = re.sub(r'^\$(e?senex).*\n', '', c, flags=re.M)
        elif var == 'norij':
            c = re.sub(r'^\$(rij|jbas|ricore).*\n', '', c, flags=re.M)
        elif var.startswith('w='):
            c, n = re.subn(r'(\$xctype RSH-GGA w=\s*)\S+', rf'\g<1>{float(var[2:]):.5f}', c)
            assert n == 1, 'xctype line not found'
        else:
            raise SystemExit(f'unknown variant {var}')
        open(os.path.join(d, 'control'), 'w').write(c)
        open(os.path.join(d, 'variant.txt'), 'w').write(var + '\n')
        names.append(name)
    open(os.path.join(a.outdir, 'cases.txt'), 'w').write('\n'.join(names) + '\n')
    print(len(names), 'cases in', a.outdir)


if __name__ == '__main__':
    main()
