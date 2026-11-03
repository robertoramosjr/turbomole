#!/bin/bash
#SBATCH -t 02:00:00
#SBATCH -N 1 --ntasks-per-node=16
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH -J scf_chain
#SBATCH -o slurm_scf_chain_%j.out
ulimit -s unlimited
export TURBODIR=$HOME/software/TURBOMOLE; source $TURBODIR/Config_turbo_env
export PARA_ARCH=SMP PARNODES=2 OMP_NUM_THREADS=2
export PATH=$TURBODIR/scripts:$TURBODIR/bin/$(sysname):$PATH
cd $SLURM_SUBMIT_DIR
IC=../newtonx_cis_camb3lyp/5-DYN/_test_tda/TRAJECTORIES
for v in d9 d8 d7 d9_strong; do
python3 scf_chain_test.py --traj t0_ic1_gap0.1_tm78 --ic-mos $IC/TRAJ1/JOB_AD/mos --start-frame 50 --variants $v --outdir scf_chain2/t0_$v &
python3 scf_chain_test.py --traj t1_ic2_gap0.1_tm78 --ic-mos $IC/TRAJ2/JOB_AD/mos --start-frame 64 --variants $v --outdir scf_chain2/t1_$v &
done
wait
