#!/bin/bash
#SBATCH -J bader_grid
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --exclusive
#SBATCH -o slurm_bader_grid_%j.out

# Bader cube-resolution convergence test: rebuild td.cub from converged MOs
# (ridft -proper, no SCF) on a finer uniform grid and rerun bader.
# Box = default TURBOMOLE cube box (taken from the reference proper.out).
#
# Usage: sbatch -p <partition> -t <time> job_bader_grid_test.sh <ref_dir> <spacing_bohr> [<spacing_bohr> ...]
#   <ref_dir> must hold control, coord, basis, auxbasis, mos (or alpha/beta)

set -u
ulimit -s unlimited
REF=$(readlink -f "${1:?ref_dir}"); shift
[ $# -ge 1 ] || { echo "need at least one spacing"; exit 1; }
BOX=${BOX:-"11.31327 44.97012 12.39716 43.33180 10.74064 41.50817"}

export PARA_ARCH=SMP PARNODES=${SLURM_CPUS_PER_TASK:-48} OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-48}
source "$HOME/software/turbomole_env.sh"
export TURBOTMPDIR="${TMPDIR:-/tmp}/tm_bader_${SLURM_JOB_ID:-local}_$$"; mkdir -p "$TURBOTMPDIR"
trap 'rm -rf "$TURBOTMPDIR"' EXIT
BADER_BIN=$(command -v bader || echo "$HOME/local_bin/bader")

read -r x1 x2 y1 y2 z1 z2 <<< "$BOX"
for h in "$@"; do
    W="$(dirname "$REF")/h_${h}"
    mkdir -p "$W"; cp -p "$REF"/* "$W"/; cd "$W"
    nx=$(awk -v a="$x1" -v b="$x2" -v h="$h" 'BEGIN{print int((b-a)/h)+1}')
    ny=$(awk -v a="$y1" -v b="$y2" -v h="$h" 'BEGIN{print int((b-a)/h)+1}')
    nz=$(awk -v a="$z1" -v b="$z2" -v h="$h" 'BEGIN{print int((b-a)/h)+1}')
    actual -r > /dev/null 2>&1
    sed -i '/^\$pointval/d; /^\$restart/d' control
    sed -i "/^\$end/i \$pointval dens fmt=cub\n   grid1 vector 1 0 0 range $x1,$x2 points $nx\n   grid2 vector 0 1 0 range $y1,$y2 points $ny\n   grid3 vector 0 0 1 range $z1,$z2 points $nz\n   origin 0 0 0" control
    echo "=== h=$h bohr: $nx x $ny x $nz = $((nx*ny*nz)) points ($(date)) ==="
    t0=$(date +%s); ridft -proper > proper.out 2>&1; t1=$(date +%s)
    grep -q "ended normally" proper.out || { echo "ridft -proper FAILED for h=$h"; tail -5 proper.out; continue; }
    grep "calculation at a grid" proper.out
    cub=$(ls -t *.cub | head -1); echo "cube: $cub ($(du -h "$cub" | cut -f1)), proper ${t1}-${t0} = $((t1-t0)) s"
    "$BADER_BIN" "$cub" > bader_run.out 2>&1; t2=$(date +%s)
    grep "NUMBER OF ELECTRONS" bader_run.out; echo "bader $((t2-t1)) s"
    echo "$h,$nx,$ny,$nz,$((t1-t0)),$((t2-t1))" > grid_info.csv
    rm -f "$cub"   # cubes are GBs; ACF.dat is what we keep
done
