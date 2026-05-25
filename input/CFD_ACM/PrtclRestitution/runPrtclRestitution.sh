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
mkdir -p ./ACM/Restart/ 2> /dev/null

./mymake.sh -exe-cfd_ACM -cmp-"$cmpStr"_MPI -CompileThirdParty-0 -deleteCompileFile-1 -CFD_DEFS_Add -ACM_DEFS_Add -CFDACM_DEFS_Add-DSeveralSphereInfo

for loop in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20
do
  cp ./input/CFD_ACM/PrtclRestitution/SpheresCoord.Case"$loop" ./ACM/Restart/SpheresCoord.bin
  $appStr -n $nproc ./cfd_ACM ./input/CFD_ACM/PrtclRestitution/PrtclRestitutionCase"$loop".cfd ./input/CFD_ACM/PrtclRestitution/PrtclRestitutionCase"$loop".acm $p_row $p_col
done