#!/bin/bash
bz2=`ls *bz2`
for i in $bz2; do
tar xvf $i
done
mkdir TARFILE
mv *bz2 TARFILE
