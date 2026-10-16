#!/bin/bash
#SBATCH -t 00:30:00
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH -n 1
#SBATCH --partition=short
#SBATCH --mem=2G
#SBATCH -o slurm_%A_%a.out
#SBATCH -e slurm_%A_%a.error
#SBATCH -J cosmo_prep

# COSMO(water) job-dir preparation for ONE labeled conformer, as a
# SLURM array (cosmoprep is light on CPU but slow on wall time because
# setup_cosmo_job.py's pexpect loop waits out a timeout per prompt;
# 1 core each, all conformers in parallel). Retries up to 3x: cosmoprep
# was seen to report "ended normally" WITHOUT writing $cosmo in ~1/40
# runs (Sec. 9.33) -- the $cosmo group is verified after every attempt.
# Output dir = source dir with /labeled_jobs/ -> /cosmo_jobs/.
#
# Usage: sbatch --array=1-N job_cosmo_prep.sh <joblist of labeled job dirs, absolute>

set -u
CONDA_BASE=/opt/gridunesp/dist/miniconda/24.4.0/libmamba
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate paper2_funnel

JOBLIST=${1:?"Usage: sbatch --array=1-N job_cosmo_prep.sh <joblist>"}
SRC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")
if [ -z "$SRC" ] || [ ! -d "$SRC" ]; then
    echo "ERROR: no valid directory on line $SLURM_ARRAY_TASK_ID of $JOBLIST"; exit 1
fi
OUT="${SRC/\/labeled_jobs\//\/cosmo_jobs\/}"
mkdir -p "$(dirname "$OUT")"

for attempt in 1 2 3; do
    rm -rf "$OUT"
    python3 /home/rramos/work_turbomole/scripts/setup_cosmo_job.py \
        --source-dir "$SRC" --output-dir "$OUT" --epsilon 78.4 > "${OUT}.setup.log" 2>&1
    if [ -f "$OUT/control" ] && grep -q '^\$cosmo$' "$OUT/control" && ! grep -q "amat file=\*" "$OUT/control"; then
        echo "COSMO prep OK (attempt $attempt): $OUT"; exit 0
    fi
    echo "attempt $attempt failed for $OUT"
done
echo "FAIL: COSMO prep after 3 attempts: $OUT"; exit 1
