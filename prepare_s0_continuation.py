#!/usr/bin/env python3
"""Prepara continuacao em S0 (NX nstat=1/nstatdyn=1, DFT fundamental) de trajetorias
Newton-X mortas por EDIFF/EDRIFT (aproximacao S1/S0).
Uso: prepare_s0_continuation.py <dir 5-DYN/xxx_run> [--tmax 150] [--only 1,2,3]
Cria <run>/../S0_cont/TRAJECTORIES/TRAJn com geom (ultimo frame ESCRITO, bohr), veloc (a.u.),
JOB_AD fundamental (sem $scfinstab/$soes/$rpacor/$denconv; mos do ultimo QM da trajetoria)
e user_config.nml (nstat=1, nstatdyn=1). Imprime a lista de N preparados."""
import sys, os, re, shutil, argparse, glob
A2B = 1/0.529177210903
ap = argparse.ArgumentParser(); ap.add_argument('run'); ap.add_argument('--tmax', type=float, default=150.0)
ap.add_argument('--only', default=''); a = ap.parse_args()
run = os.path.abspath(a.run); out = os.path.join(os.path.dirname(run), 'S0_cont', 'TRAJECTORIES')
only = {int(x) for x in a.only.split(',') if x}
done = []
for d in sorted(glob.glob(run + '/TRAJECTORIES/TRAJ*'), key=lambda s: int(re.findall(r'\d+$', s)[0])):
    n = int(re.findall(r'\d+$', d)[0])
    if only and n not in only: continue
    try: log = open(d + '/md.log').read()
    except OSError: continue
    if 'KILLED' not in log or 'ABNORMAL TERMINATION' not in log: continue   # so mortas
    try:
        gx = open(d + '/RESULTS/geometries.xyz').read().split('\n')
        vl = open(d + '/RESULTS/velocities.dat').read().split('\n')
    except OSError: continue
    nat = int(gx[0]); blk = nat + 2
    nfr = len([l for l in gx if l.strip()]) // blk
    xyz = gx[(nfr-1)*blk+2:(nfr-1)*blk+2+nat]; tfs = float(gx[(nfr-1)*blk+1].split()[0])
    vb = vl[(nfr-1)*(nat+2)+2:(nfr-1)*(nat+2)+2+nat]
    assert float(vl[(nfr-1)*(nat+2)+1].split()[0]) == tfs, (n, 'tempo geom != veloc')
    g0 = open(d + '/geom').read().split('\n')[:nat]
    t = os.path.join(out, 'TRAJ%d' % n)
    if os.path.exists(t): shutil.rmtree(t)
    os.makedirs(t)
    with open(t + '/geom', 'w') as f:
        for g, l in zip(g0, xyz):
            p = g.split(); x = [float(v)*A2B for v in l.split()[1:4]]
            assert p[0].lower() == l.split()[0].lower()
            f.write(' %s %6.1f %14.8f %14.8f %14.8f %14.8f\n' % (p[0], float(p[1]), x[0], x[1], x[2], float(p[5])))
    open(t + '/veloc', 'w').write('\n'.join(vb) + '\n')
    shutil.copytree(d + '/JOB_AD', t + '/JOB_AD')
    shutil.copy(d + '/TEMP/turbomole.job/mos', t + '/JOB_AD/mos')
    c = open(t + '/JOB_AD/control').read()
    c = re.sub(r'\$scfinstab.*\n', '', c); c = re.sub(r'\$soes\n.*\n', '', c)
    c = re.sub(r'\$rpacor.*\n', '', c); c = re.sub(r'\$denconv.*\n', '', c)
    c = c.replace('$scfiterlimit       30', '$scfiterlimit       100')
    open(t + '/JOB_AD/control', 'w').write(c)
    u = open(d + '/user_config.nml').read()
    for k, v in (('nstat', '1'), ('nstatdyn', '1'), ('tmax', '%.1f' % a.tmax), ('init_time', '0.0'), ('init_step', '0'), ('thres', '0.0')):
        u = re.sub(r'(?m)^ %s = .*$' % k, ' %s = %s' % (k, v), u)
    open(t + '/user_config.nml', 'w').write(u)
    open(t + '/START_INFO', 'w').write('S1 morta em t=%.1f fs (traj %d); start S0 do ultimo frame escrito\n' % (tfs, n))
    done.append(n)
print(','.join(map(str, done)))
