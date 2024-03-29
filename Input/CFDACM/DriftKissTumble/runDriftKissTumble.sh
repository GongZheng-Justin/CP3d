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

./mymake.sh -exe-channelACM -cmp-"$cmpStr"_MPI -CompileThirdParty-0 -deleteCompileFile-1 -CFD_DEFS_Add -ACM_DEFS_Add-DDriftKissTumbleBreugem \
 -CFDACM_DEFS_Add-DSeveralSphereInfo

cp ./Input/CFDACM/DriftKissTumble/SpheresCoord.DKT ./ACM/Restart/SpheresCoord.dat
chmod a+x ./channelACM
$appStr -n $nproc ./channelACM ./Input/CFDACM/DriftKissTumble/DriftKissTumble.cfd ./Input/CFDACM/DriftKissTumble/DriftKissTumble.acm $p_row $p_col
cd ./Input/CFDACM/DriftKissTumble/
