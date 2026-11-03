#!/bin/bash
# usage: make_batch.sh <series_name> <NX TRAJECTORIES dir> <gap_eV> <profile>
S=$1; NX=$2; GAP=$3; PROF=${4:-iter300}
OUT=batch_${S}_gap${GAP}; mkdir -p $OUT; : > $OUT/tasks.txt
for t in $(ls $NX | grep '^TRAJ' | sort -V); do
  python3 ../../scripts/prepare_tm_native_sh.py $NX/$t $OUT/$t --gap $GAP --nsteps 400 --scf-profile $PROF > /dev/null && echo "$PWD/$OUT/$t 78" >> $OUT/tasks.txt
done
wc -l $OUT/tasks.txt
