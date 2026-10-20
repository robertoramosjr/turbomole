#!/bin/bash
#SBATCH -t 02:00:00
#SBATCH -N 1 --ntasks-per-node=8
#SBATCH -n 8
#SBATCH --partition=short
#SBATCH --mem=16G
#SBATCH -o slurm_%A_%a.out
#SBATCH -e slurm_%A_%a.error
#SBATCH -J w4_sp

# W4: HSE06/def2-TZVP + COSMO(water) single point (ridft) on one water-reoptimised geometry.
# Dirs prepared by prepare_w4_sp.py.  Marker: "ridft : all done" in ridft.out; validate by TEXT, not exit code.
# Usage: sbatch --array=1-N job_w4_sp.sh <joblist_file>   (absolute dir paths, one per line)
set -u
ulimit -s unlimited
CONDA_BASE=/opt/gridunesp/dist/miniconda/24.4.0/libmamba
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate paper2_funnel
JOBLIST=${1:?"Usage: sbatch --array=1-N job_w4_sp.sh <joblist_file>"}
WORK_DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")
[ -n "$WORK_DIR" ] && [ -d "$WORK_DIR" ] || { echo "ERROR: no valid dir on line $SLURM_ARRAY_TASK_ID of $JOBLIST"; exit 1; }
cd "$WORK_DIR"
export PARA_ARCH=SMP PARNODES=$SLURM_NTASKS OMP_NUM_THREADS=$SLURM_NTASKS OMP_STACKSIZE=2G
ridft > ridft.out 2>&1
grep -q "ridft : all done" ridft.out && ! grep -qi "ended abnormally\|did not converge" ridft.out && echo W4_DONE > W4_DONE
echo "w4 exit for $WORK_DIR done=$([ -f W4_DONE ] && echo 1 || echo 0)"
