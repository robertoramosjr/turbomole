#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH -N 1 --ntasks-per-node=8
#SBATCH -n 8
#SBATCH --partition=short
#SBATCH --mem=16G
#SBATCH -o slurm_%A_%a.out
#SBATCH -e slurm_%A_%a.error
#SBATCH -J solv_relax

# Relaxes ONE conformer in implicit water: copies an r2SCAN-3c-optimised conformer directory
# (control/basis/mos/coord), adds COSMO(eps=78.4) with cosmoprep (setup_cosmo_job.py, retried
# because cosmoprep occasionally "ends normally" without writing $cosmo), then runs
# jobex -ri (r2SCAN-3c/def2-mTZVP + COSMO) with a cycle cap. Answers whether the conformer
# ranking found with COSMO single points on gas-phase geometries survives relaxation in solvent.
# Verify results by TEXT (energy file / GEO_OPT_CONVERGED), never by exit code alone.
#
# Usage: sbatch --array=1-N job_solvent_relax.sh <joblist: "source_dir output_dir" per line> [max_cycles]

set -u
ulimit -s unlimited
CONDA_BASE=/opt/gridunesp/dist/miniconda/24.4.0/libmamba
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate paper2_funnel

JOBLIST=${1:?"Usage: sbatch --array=1-N job_solvent_relax.sh <joblist> [max_cycles]"}
CYC=${2:-25}
read SRC OUT < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")
[ -d "$SRC" ] || { echo "ERROR: bad source $SRC"; exit 1; }

if [ ! -f "$OUT/control" ] || ! grep -q '^\$cosmo$' "$OUT/control"; then
  for attempt in 1 2 3; do
    rm -rf "$OUT"
    python3 /home/rramos/work_turbomole/scripts/setup_cosmo_job.py --source-dir "$SRC" --output-dir "$OUT" \
        --epsilon 78.4 > "${OUT}.setup.log" 2>&1
    if grep -q '^\$cosmo$' "$OUT/control" 2>/dev/null && ! grep -q "amat file=\*" "$OUT/control"; then break; fi
  done
fi
grep -q '^\$cosmo$' "$OUT/control" || { echo "FAIL: no \$cosmo in $OUT"; exit 1; }

cd "$OUT"
export PARA_ARCH=SMP PARNODES=$SLURM_NTASKS
jobex -ri -c "$CYC" > jobex.out 2>&1
echo "jobex exit $? for $OUT"
