#!/bin/bash
#SBATCH -t 04:00:00
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH -n 1
#SBATCH --partition=short
#SBATCH --mem=4G
#SBATCH -o slurm_%A_%a.out
#SBATCH -e slurm_%A_%a.error
#SBATCH -J xtb_alpb

# Re-optimizes every CREST conformer of one protonation state with
# GFN2-xTB + ALPB(water) (the CREST ensembles were built in VACUUM --
# protonation_and_conformer_investigation.md Sec. 9.38), one array task
# per chunk (conformer i goes to task ((i-1) % NCHUNKS)+1). The total
# GFN2/ALPB energy (includes the solvation free energy) is what
# analyze_alpb_rerank.py compares with the vacuum GFN2 ranking.
# Runs on 1 core/task; light (seconds-minutes per conformer).
#
# Usage: sbatch --array=1-<NCHUNKS> job_xtb_alpb_rerank.sh <state_dir> <NCHUNKS>

set -u
STATE_DIR=${1:?"Usage: sbatch --array=1-N job_xtb_alpb_rerank.sh <state_dir> <NCHUNKS>"}
NCHUNKS=${2:?"need NCHUNKS"}
TASK=$SLURM_ARRAY_TASK_ID
CHARGE=$(cat "$STATE_DIR/CHARGE.txt")
export OMP_NUM_THREADS=1
export OMP_STACKSIZE=1G
ulimit -s unlimited

OUT="$STATE_DIR/alpb_rerank"
mkdir -p "$OUT"
NAT=$(head -1 "$STATE_DIR/crest_conformers.xyz" | tr -d ' ')
NCONF=$(( $(wc -l < "$STATE_DIR/crest_conformers.xyz") / (NAT + 2) ))

for (( i=1; i<=NCONF; i++ )); do
    if [ $(( (i-1) % NCHUNKS + 1 )) -ne "$TASK" ]; then continue; fi
    D=$(printf "%s/conf_%03d" "$OUT" "$i")
    if [ -f "$D/xtbopt.xyz" ] && grep -q "TOTAL ENERGY" "$D/xtb.out" 2>/dev/null; then continue; fi
    mkdir -p "$D"
    tail -n +$(( (i-1)*(NAT+2) + 1 )) "$STATE_DIR/crest_conformers.xyz" | head -n $(( NAT + 2 )) > "$D/start.xyz"
    ( cd "$D" && xtb start.xyz --gfn 2 --alpb water --chrg "$CHARGE" --opt tight > xtb.out 2>&1 )
    echo "conf $i done: $(grep 'TOTAL ENERGY' "$D/xtb.out" | tail -1)"
done
