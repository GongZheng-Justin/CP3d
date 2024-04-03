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

cd ../../
chmod a+x mymake.sh
./mymake.sh -exe-channelLPT -cmp-gcc_MPI -CompileThirdParty-0 -deleteCompileFile-1 -CFD_DEFS_Add -LPT_DEFS_Add -CFDOrder-2 -Coupling-2

for loop in 4 3 2 1
do
  $appStr -n $nproc ./channelLPT ./Input/CFDLPT_TwoWay/Channel4th_LPT.Lee"$loop" ./Input/CFDLPT_TwoWay/LPT_Channel4th.Lee"$loop" \
  $p_row $p_col
done
cd ./Input/CFDLPT_TwoWay/
