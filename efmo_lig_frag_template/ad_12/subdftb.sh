#!/bin/bash
sed -i 's/ngroup=4/ngroup=8/gI' *inp
sed -i 's/ngrfmo(1)=4,4/ngrfmo(1)=8,8/gI' *inp
xgms=/lustre/hdd/LAS/mgordon-lab/tsatta/gamess-w2/gms
for i in *inp; do
 log=${i%.inp}.log
 $xgms -p 64 -logn 8 -ppn 32 -w 10:00:00 $i -l $log
done
# -logn 8
