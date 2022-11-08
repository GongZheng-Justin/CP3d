#!/bin/bash
cmpStr="gcc"
appStr="mpirun"
p_row=4
p_col=2
if [[ -n $1 ]]; then
  cmpStr=$1
fi
if [[ -n $2 ]]; then
  appStr=$2
fi
if [[ -n $3 ]]; then
  p_row=$3
fi
if [[ -n $4 ]]; then
  p_col=$4
fi
let nproc=p_row*p_col

cd ../../../
chmod a+x mymake.sh
./mymake.sh -exe-dem -cmp-"$cmpStr"_MPI -CompileThirdParty-0 -deleteCompileFile-1 -DEM_DEFS_Add
chmod a+x ./dem
$appStr -n $nproc ./dem ./Input/CFDACM/TestScaling/CFDACMStrong.dem $p_row $p_col
