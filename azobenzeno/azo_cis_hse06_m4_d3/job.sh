#!/bin/bash
#SBATCH --job-name=azo_cis
#SBATCH --time=24:00:00
#SBATCH --output=job.out
#SBATCH --error=job.err
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4

# Pipeline completo (Estagio 1 do MANUAL.md em ~/work_turbomole/scripts),
# para uma pasta recem-saida do define, com $optimize presente no control
# (por isso comeca com jobex, como no job.pbs do alizarina).
# Ajuste --partition/--cpus-per-task/--mem conforme o tamanho da molecula
# (short = ate 1 dia, medium = ate 7 dias, long = ate 30 dias).

set -e
ulimit -s unlimited

# --- Ambiente Turbomole (mesmo setup do ~/.bashrc, explicito para nao
#     depender do bash sourcear o rc em job nao-interativo) ---
export TURBODIR=$HOME/software/TURBOMOLE
export PATH=$TURBODIR/scripts:$PATH
export PARA_ARCH=SMP
export PARNODES=$SLURM_CPUS_PER_TASK
source $TURBODIR/Config_turbo_env

# Evita oversubscription de threads em bibliotecas linkadas (MKL/OpenMP)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK  
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd $SLURM_SUBMIT_DIR

# 1. Otimizacao de geometria (HSE06/senex, criterios do $statpt)
#    Atualiza coord e mos ao final, na geometria relaxada
jobex > jobex.out

# 2. Cubo de densidade para Bader (na geometria ja otimizada)
ridft -proper > proper.out

# 3. Bader charges (fora do Turbomole)
BADER_BIN=$(command -v bader || echo ~/local_bin/bader)
$BADER_BIN td.cub > bader_run.out

# 4. Reconverge densidade apertada antes do TD-DFT
ridft > ridft_tight.out

# 5. TD-DFT (excitacoes -> JDOS + optico)
escf > escf.out

# 6. Hessiana / IR (geometria ja relaxada -- sem frequencia imaginaria esperada)
aoforce > aoforce.out

# 7. Populacao (Mulliken + SEN/bond-order)
dscf -proper > pop.out
