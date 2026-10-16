#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH -N 1 --ntasks-per-node=8
#SBATCH -n 8
#SBATCH --partition=short
#SBATCH --mem=16G
#SBATCH -o slurm_%A_%a.out
#SBATCH -e slurm_%A_%a.error
#SBATCH -J w_reopt

# Re-optimises ONE tier-2 survivor (line SLURM_ARRAY_TASK_ID of <survivors file>) with r2SCAN-3c/def2-mTZVP + COSMO(water):
# jobex -ri -c <cycles>. Resumable: an existing GEO_OPT_CONVERGED is skipped; a partial or capped (GEO_OPT_FAILED) run continues from its last geometry.
# Validate by TEXT (GEO_OPT_CONVERGED / energy), never by exit code.
#
# Usage: sbatch --array=1-N%T job_water_reopt.sh <state> <survivors.txt> [cycles]   (cwd = artepillin_C_studies)

set -u
source /opt/gridunesp/dist/miniconda/24.4.0/libmamba/etc/profile.d/conda.sh; conda activate paper2_funnel
ulimit -s unlimited
export TURBODIR=$HOME/software/TURBOMOLE; source $TURBODIR/Config_turbo_env
export PARA_ARCH=SMP PARNODES=$SLURM_NTASKS OMP_NUM_THREADS=$SLURM_NTASKS
export PATH=$TURBODIR/scripts:$TURBODIR/bin/$(sysname):$PATH

ST=${1:?state}; LIST=${2:?survivors file}; CYC=${3:-40}
R=/home/rramos/work_turbomole/artepillin_C_studies
IDX=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$R/$LIST")
[ -n "$IDX" ] || { echo "ERROR: no index on line $SLURM_ARRAY_TASK_ID of $LIST"; exit 1; }
OUT=$R/reopt/$ST/c$(printf %04d "$IDX")
[ -f "$OUT/GEO_OPT_CONVERGED" ] && { echo "already converged $OUT"; exit 0; }
if [ ! -f "$OUT/control" ]; then
  python3 /home/rramos/work_turbomole/scripts/water_reopt_conformer.py --xyz $R/water_ensemble/$ST/crest_conformers.xyz \
      --template $R/tier2/$ST/template --index "$IDX" --outdir "$OUT" || exit 1
fi
cd "$OUT"
# A previous run that hit the cycle cap leaves GEO_OPT_FAILED/not.converged; clear them and keep its log, then continue from coord.
if [ -f GEO_OPT_FAILED ]; then
  rm -f GEO_OPT_FAILED not.converged
  cat jobex.out >> jobex.prev.out 2>/dev/null
fi
jobex -ri -c "$CYC" > jobex.out 2>&1
echo "jobex exit $? for $OUT conv=$([ -f GEO_OPT_CONVERGED ] && echo 1 || echo 0)"
