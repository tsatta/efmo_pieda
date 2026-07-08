#!/bin/bash

# grep "TOTAL FREE ENERGY IN SOLVENT =" unit will be in kcal/mol
efmo=fmo
MM_path=/mnt/lustre/koa/lab/tsqc_group/tsatta/hiv_2025/MM_xyz/TS_info
PPE_path=/mnt/lustre/koa/lab/tsqc_group/loreab/fmo_2res
hope=/mnt/lustre/koa/lab/tsqc_group/tsatta/hiv_2025/MM_xyz/hope_hiv.txt
 ligdir=`awk '{print $1}' lig_list`
protein=2i0d
 for ligand in $ligdir
do
PPE_file=${protein}_${ligand}_PPE
xyz_lig=$MM_path/lig_${ligand}_TS_info.txt
xyz_pro=$MM_path/protein_TS_info.txt
xyz_com=${MM_path}/complex_${ligand}_TS_info.txt

#Protein
protein_path=${PPE_path}/protein
nrank=` ls $protein_path/*log | wc -l`
HO=`grep 'protein' $hope | awk '{print $2}'`
echo 'Protein' $nrank $HO > $PPE_file

 i=1
    while [ $i -le $nrank ]
    do
     grep "Rank=${i} " $xyz_pro >> $PPE_file
    ((i++))
    done

#nrank=1
echo 'DFTB SPE' >> $PPE_file

i=1
    while [ $i -le $nrank ]
    do
    protein_log=$protein_path/2i0d_FMO2DFTB4H_SPE_rank${i}.log
    num="$(grep 'The best FMO energy is' $protein_log | awk '{print $6}')"
    #echo $num  >> $PPE_file
    echo $num*627.509 | bc -l  >> $PPE_file
    ((i++))
    done
#-------------------
#ligand
#-------------------
 lig_path=$PPE_path/ligands/${ligand}
 echo 'lig_path', $lig_path
nrank=` ls $lig_path/*log | wc -l`
HO=`grep "lig ${ligand}" $hope | awk '{print $3}'`
echo 'ligand' $nrank $HO >> $PPE_file

 i=1
    while [ $i -le $nrank ]
    do
     grep "Rank=${i} " $xyz_lig >> $PPE_file
    ((i++))
    done

echo 'DFTB SPE' >> $PPE_file

i=1
    while [ $i -le $nrank ]
    do
    lig_log=$lig_path/${ligand}_DFTB4H_SPE_rank${i}.log
    num="$(grep 'TOTAL FREE ENERGY IN SOLVENT =' $lig_log | awk '{print $7}')"
    echo $num  >> $PPE_file
    #echo $num*627.509 | bc -l  >> $PPE_file
    ((i++))
    done
#-------------------
#complex
#-------------------
 com_path=$PPE_path/complexes_spe/${ligand}
nrank=` ls $com_path/*log | wc -l`
HO=`grep "complex ${ligand}" $hope | awk '{print $3}'`
echo 'complex' $nrank $HO >> $PPE_file

 i=1
    while [ $i -le $nrank ]
    do
     grep "Rank=${i} " $xyz_com >> $PPE_file
    ((i++))
    done

echo 'DFTB SPE' >> $PPE_file

i=1
    while [ $i -le $nrank ]
    do
    com_log=$com_path/2i0d_${ligand}_FMO2DFTB4H_SPE_rank${i}.log
    #num="$(grep 'EunD+so(2)=' $com_log  | awk '{print $7}')"
    num="$(grep 'The best FMO energy is' $com_log | awk '{print $6}')"
    #echo $num  >> $PPE_file
    echo $num*627.509 | bc -l  >> $PPE_file
    ((i++))
    done

 sed -i 's/=/ /gI' $PPE_file
 echo 'done' $PPE_file
done
