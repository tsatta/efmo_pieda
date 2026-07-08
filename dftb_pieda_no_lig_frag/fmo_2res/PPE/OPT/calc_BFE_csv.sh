#!/bin/bash
# obtain free energy of ligands from qm-vm2 using vm2 procedure
 lig_FEpath=/work/LAS/mgordon-lab/tsatta/PLQM/ligand_qmvm2_results
 lig_FE=$lig_FEpath/hiv_ligand_qmvm2_freeEnergy_vm2opt.txt
# extract list of ligands
 ligand=`sed 's/,/ /gI' $lig_FE | awk '{print $1}'`

# print binding FE using feprocess and VM2 for ligands for PPE qm-vm2
 output=binding_FE_fullprotein.csv
 echo 'lig, qmvm2_feprocess, qmvm2_vm2' > $output

# loop over ligands
 for lig in $ligand
 do
  f1=2i0d_${lig}_2res_OPT_30HO_binding
  #f1=${lig}_newopttol_OPT_30HO.out
#------------ binding FE from feprocess ------------
  E_host=`grep 'E_host' $f1 | awk '{print $2}'`
  E_lig=`grep 'E_guest' $f1 | awk '{print $2}'`
  E_complx=`grep 'E_complx' $f1 | awk '{print $2}'`
  bindFE=`echo $E_complx - $E_host - $E_lig | bc -l`

#----- binding FE using qmvm2_vm2 for lig  ------------
# lig free energy is stored in the 2nd column of lig_FE
#
  E_lig_vm2=`grep "$lig"  $lig_FE | awk '{print $2}'`
  bindFE_VM2=`echo $E_complx - $E_host - $E_lig_vm2 | bc -l`
  echo $lig, $bindFE, $bindFE_VM2 >> $output

 done
