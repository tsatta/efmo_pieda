#!/bin/bash
dir=`ls -d a* k*`
for i in $dir ; do
 cd $i
 sed -n -ie '/Two-body FMO properties/,/Done with FMO/p' *log
 rm *loge
 cd ../
 tar -cvjf ${i}.tar.bz2 $i
done

