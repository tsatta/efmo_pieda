#!/bin/bash
data=experimental_data.csv
ligdir=`cat $data | awk -F',' '{print $1}'`
protein=2i0d
workdir=`pwd`
out=binding.txt
echo lig,exp,1res_cal > $out
for lig in $ligdir; do
FE=`grep 'binding' -i $lig/*binding | awk '{print $NF}'`
exp=`grep "$lig" $data | awk -F',' '{print $NF}'`
echo $lig,$exp,$FE >> $out
done

