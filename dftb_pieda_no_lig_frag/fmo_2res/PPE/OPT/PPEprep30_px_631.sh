#!/bin/bash

# grep "TOTAL FREE ENERGY IN SOLVENT =" unit will be in kcal/mol
#efmo=fmo
#MM_path=/work/LAS/mgordon-lab/tsatta/PLQM/2024_tnks2_cutout/cutout_real6_live4/feprocess
#PPE_path=/work/LAS/mgordon-lab/tsatta/PLQM/2024_tnks2_cutout/full_plqm_real6_live4/OPT
#hope=/work/LAS/mgordon-lab/tsatta/PLQM/2024_tnks2_cutout/cutout_real6_live4/hope_tnks2.txt
 ligdir=`awk '{print $1}' lig_list`
protein=2i0d
 for ligand in $ligdir
do
PPE_file=${protein}_${ligand}_PPE
##xyz_lig=$MM_path/lig_${ligand}_TS_info.txt
##xyz_pro=$MM_path/protein_TS_info.txt
##xyz_com=${MM_path}/complex_${ligand}_TS_info.txt
#xyz_lig=$MM_path/ligands/$ligand/${ligand}_*vm2_conformers.xyz
#xyz_pro=$MM_path/protein/2i0d_p4a_tleap/*_vm2_conformers.xyz
#xyz_com=$MM_path/complexes/$ligand/*_vm2_conformers.xyz

##Protein
#protein_path=${PPE_path}/protein
#nrank=` ls $protein_path/*log | wc -l`
#HO=`grep 'protein' $hope | awk '{print $2}'`
#echo 'Protein' $nrank $HO > $PPE_file
#
# i=1
#    while [ $i -le $nrank ]
#    do
#     grep "Rank=${i} " $xyz_pro >> $PPE_file
#    ((i++))
#    done
#
##nrank=1
#echo 'DFTB OPT' >> $PPE_file
#
#i=1
#    while [ $i -le $nrank ]
#    do
#    protein_log=$protein_path/2i0d_FMO2DFTB4H_OPT_rank${i}.log
#    num="$(grep 'The best FMO energy is' $protein_log | awk '{print $6}')"
#    #echo $num  >> $PPE_file
#    echo $num*627.509 | bc -l  >> $PPE_file
#    ((i++))
#    done
##-------------------
##ligand
##-------------------
# lig_path=$PPE_path/ligands/${ligand}
# echo 'lig_path', $lig_path
#nrank=` ls $lig_path/*log | wc -l`
#HO=`grep "lig ${ligand}" $hope | awk '{print $3}'`
#echo 'ligand' $nrank $HO >> $PPE_file
#
# i=1
#    while [ $i -le $nrank ]
#    do
#     grep "Rank=${i} " $xyz_lig >> $PPE_file
#    ((i++))
#    done
#
#echo 'DFTB OPT' >> $PPE_file
#
#i=1
#    while [ $i -le $nrank ]
#    do
#    lig_log=$lig_path/${ligand}_DFTB4H_OPT_rank${i}.log
#    num="$(grep 'TOTAL FREE ENERGY IN SOLVENT =' $lig_log | awk '{print $7}')"
#    echo $num  >> $PPE_file
#    #echo $num*627.509 | bc -l  >> $PPE_file
#    ((i++))
#    done
#-------------------
#complex
#-------------------
 #com_path=$PPE_path/complexes/${ligand}
 com_path=/work/LAS/mgordon-lab/meganthegreat/hiv_2res/OPT/complexes/${ligand}
nrank=` ls $com_path/*log | wc -l`
#HO=`grep "complex ${ligand}" $hope | awk '{print $3}'`
#echo 'complex' $nrank $HO >> $PPE_file

# i=1
#    while [ $i -le $nrank ]
#    do
#     grep "Rank=${i} " $xyz_com >> $PPE_file
#    ((i++))
#    done

echo 'DFTB OPT' >> $PPE_file

i=1
    while [ $i -le $nrank ]
    do
    com_log=$com_path/2i0d_${ligand}_FMO2DFTB4H_OPT_rank${i}.log
    #num="$(grep 'EunD+so(2)=' $com_log  | awk '{print $7}')"
    num="$(grep 'The best FMO energy is' $com_log |tail -n 1 | awk '{print $6}')"
    #echo $num  >> $PPE_file
    echo $num*627.509 | bc -l  >> $PPE_file
    ((i++))
    done

 sed -i 's/=/ /gI' $PPE_file
 echo 'done' $PPE_file
done
