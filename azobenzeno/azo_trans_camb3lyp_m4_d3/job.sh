#!/bin/bash
#SBATCH --job-name=azo_trans_camb3lyp_m4_d3
#SBATCH --time=02:00:00
#SBATCH --output=job.out
#SBATCH --error=job.err
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4

# Funcional-benchmark de estados excitados (rota B1, etapa 2): single-point
# ridft + escf na geometria HSE06 ja convergida. CAM-B3LYP usa senex (troca
# seminumerica) em vez de RI-K -- Turbomole recusa RI-K pra hibridos
# range-separated ("range-separated hybrids not possible with RI-K!").
# Pasta gerada por prepare_functional_benchmark.py --enable-senex.

set -e
ulimit -s unlimited

export TURBODIR=$HOME/software/TURBOMOLE
export PATH=$TURBODIR/scripts:$PATH
export PARA_ARCH=SMP
export PARNODES=$SLURM_CPUS_PER_TASK
source $TURBODIR/Config_turbo_env

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd $SLURM_SUBMIT_DIR

# 1. SCF na geometria fixa (RI-J + senex, D3(BJ), grid m4)
ridft > ridft.out

# 2. TD-DFT -- 10 estados singleto, pra capturar n->pi* (escuro) e pi->pi* (brilhante)
escf > escf.out
