#!/bin/bash
ligdir=`ls *binding`
for lig in $ligdir; do
awk '/Complex/,/E_/' $lig | sed '$d' | tail -n +2 > temp
out=${lig}_complex_prob.txt
echo rankID, probability'%' >> $out
awk '{print $2", "$3}' temp >> $out
done
rm temp
