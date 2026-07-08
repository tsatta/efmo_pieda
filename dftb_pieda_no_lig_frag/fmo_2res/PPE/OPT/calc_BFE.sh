#!/bin/bash
gfortran -o boltzman_HO boltzman_HO.F90 

rm *binding

#protein=4r1y
#basis=6311+Gd
#EFMO=efp_noD
# for lig in ac_67
#do
# ligand=${lig}
#./boltzman_HO < ${protein}_${ligand}_PPE_${basis}_${EFMO} > ${protein}_${ligand}_${basis}_${EFMO}_binding
#done

 #for lig in $f1; do
 for lig in *30HO ; do
 echo $lig
 ./boltzman_HO < $lig > ${lig}_binding
 done

