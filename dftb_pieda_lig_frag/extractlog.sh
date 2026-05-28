#!/bin/bash
ligdir=`ls -d a* k*`
protein=2i0d
workdir=`pwd`
for lig in $ligdir; do
   cd $lig
   sed -n -ie '/Two-body FMO properties/,/Done/p' *log
 rm *loge
 cd $workdir
done

