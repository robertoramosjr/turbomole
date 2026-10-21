#!/bin/bash
#SBATCH -J tm_sh_s0
#SBATCH -t 16:00:00
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH -o slurm_%A_%a.out

# One native-TURBOMOLE surface-hopping trajectory, end to end:
#   phase 1  S1 dynamics (jobex -ri -ex -md, frog $surface_hopping + $gap_threshold) until the forced/FSSH
#            S1->S0 hop is detected in event.yaml (a watcher then touches `stop` 2 steps after the hop);
#            the run may also end earlier by itself (e.g. SCF failure near the crossing) -- both paths go on.
#   phase 2  ground-state MD continuation from the LAST mdlog frame (prepare_tm_s0_continuation.py), cheap
#            (~5 s/step), NSTEPS_S0 steps; then classify the product by CNNC (classify_tm_md.py).
# v3 (2026-10-20): SHLOG seeding (+SHLOG_RESET watcher) and S0 ladder (K mos modo perfil); PHASE=s0only skips phase 1.
# v2 (2026-10-16): sem salto => fase 2 pulada; result.txt vem de summarize_sh_traj.py (hop_ok|no_hop|died_S1|died_S0).
# Verify results by TEXT (result.txt, mdlog.1), never by exit code.
#
# Usage: sbatch --array=1-N job_tm_native_sh_s0.sh <tasks.txt: "<traj_dir> <TM version 78|801>" per line> [NSTEPS_SH=400] [NSTEPS_S0=300] [SCF_PROFILE=iter300]   env opcional: SEGMENT=<passos> SHLOG=<maxlog, 0=off> PHASE=<full|s0only>
#        tasks line: "<traj_dir> <TM 78|801> [<JOB_AD/mos path, needed when a ladder entry uses mos=init>]"
set -u
ulimit -s unlimited
TASKS=${1:?"tasks file"}; NSH=${2:-400}; NS0=${3:-300}; PROF=${4:-iter300}
SHLOG=${SHLOG:-1000}; PHASE=${PHASE:-full}
LADDER=( "0 sh omp78 iter300" "2 sh omp78 iter300" "0 sh omp78 damp_shift" )   # K mos mode profile (ledger 2026-10-20, T3)
SEGMENT=${SEGMENT:-}   # fase 1 em blocos (usado pela Tarefa 2; vazio = bloco único)
read DIR TM MOSINIT < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASKS")
if [ "$TM" = "801" ]; then export TURBODIR=$HOME/software/TURBOMOLE801; else export TURBODIR=$HOME/software/TURBOMOLE; fi
export PARA_ARCH=SMP PARNODES=$SLURM_CPUS_ON_NODE OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE MKL_NUM_THREADS=$SLURM_CPUS_ON_NODE
source $TURBODIR/Config_turbo_env
export PATH=$TURBODIR/scripts:$TURBODIR/bin/$(sysname):$PATH
SC=/home/rramos/work_turbomole/scripts
cd "$DIR" || exit 1
echo "TM=$TM host=$(hostname) $(date)" > run_info.txt

if [ "$PHASE" = full ]; then
  [ "$SHLOG" != 0 ] && python3 $SC/seed_sh_log_groups.py control --maxlog $SHLOG >> run_info.txt 2>&1
  sed -i "s/^\$nsteps .*/\$nsteps      $NSH/" mdmaster
  ( while sleep 30; do
      if [ -f event.yaml ] && grep -q 'event: hop' event.yaml; then
        hs=$(grep -m1 -A1 'event: hop' event.yaml | grep -o 'steps: *[0-9]*' | grep -o '[0-9]*')
        nf=$(grep -c '^t=' mdlog.1)
        if [ "$nf" -ge $((hs + 3)) ]; then touch stop; break; fi
      fi
      if [ "$SHLOG" != 0 ] && [ "$(grep -c 'maxlog= *'$SHLOG'\b' control)" -lt 4 ] && [ "$(grep -c '^t=' mdlog.1)" -ge 2 ]; then
        touch stop; echo SHLOG_RESET >> run_info.txt; break
      fi
      pgrep -u $USER -f "jobex -ri -ex -md" > /dev/null || break
    done ) &
  WATCH=$!
  jobex -ri -ex -md -c $((NSH + 5)) > jobex.out 2>&1
  kill $WATCH 2>/dev/null; rm -f stop
  echo "phase1 done $(date) frames=$(grep -c '^t=' mdlog.1)" >> run_info.txt
fi

if [ "$PHASE" = s0only ] || grep -q 'event: hop' event.yaml 2>/dev/null; then
  for ENTRY in "${LADDER[@]}"; do
    read K MOS MODE PROFILE <<< "$ENTRY"
    if [ -d s0 ]; then n=1; while [ -d s0_fail$n ]; do n=$((n+1)); done; mv s0 s0_fail$n; fi
    MOSARG=(); [ "$MOS" = init ] && MOSARG=(--mos-from "$MOSINIT")
    python3 $SC/prepare_tm_s0_continuation_v3.py . s0 --nsteps $NS0 --scf-profile $PROFILE --frame-back $K "${MOSARG[@]}" >> run_info.txt 2>&1
    (
      case "$MODE" in
        serial78) unset PARA_ARCH; export PARNODES=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 ;;
        omp801) export TURBODIR=$HOME/software/TURBOMOLE801; source $TURBODIR/Config_turbo_env; export PATH=$TURBODIR/scripts:$TURBODIR/bin/$(sysname):$PATH ;;
      esac
      cd s0 && jobex -ri -md -c $((NS0 + 5)) > jobex.out 2>&1
    )
    ST=$(python3 $SC/summarize_sh_traj.py . 2>&1 | grep -o 'status=[a-z_0-9A-Z]*' | head -1)
    echo "s0 attempt (K=$K mos=$MOS $MODE $PROFILE): $ST frames=$(grep -c '^t=' s0/mdlog.1 2>/dev/null)" >> run_info.txt
    [ "$ST" = "status=hop_ok" ] && break
  done
  echo "phase2 done $(date)" >> run_info.txt
else
  echo "phase2 skipped (no hop)" >> run_info.txt
fi
python3 $SC/summarize_sh_traj.py . > result.txt 2>&1
cat result.txt >> run_info.txt
