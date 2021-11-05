#!/bin/sh
#
set -x
here=`pwd`
cd ..
base=`pwd`
cd ../../exec
execDir=`pwd`
cd $base/work
workDir=`pwd`

cd $here

mv $workDir/ocndas.ts.mpi.x $execDir/godas_dw_mom3das

