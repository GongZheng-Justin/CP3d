#!/bin/bash
nproc=8
appName="mpirun"

cd ../../../
chmod a+x mymake.sh
./mymake.sh -exe-channelACM -cmp-gcc_MPI -CompileThirdParty-0 -deleteCompileFile-1 -CFD_DEFS_Add -ACM_DEFS_Add -CFDACM_DEFS_Add-DSeveralSphereInfo

mkdir -p ./ACM/Restart/ 2> /dev/null
cp ./Input/CFDACM/DriftKissTumble/SpheresCoord.DKT ./ACM/Restart/SpheresCoord.dat
chmod a+x ./channelACM
$appName -n $nproc ./channelACM ./Input/CFDACM/DriftKissTumble/DriftKissTumble.cfd ./Input/CFDACM/DriftKissTumble/DriftKissTumble.acm
cd ./Input/CFDACM/DriftKissTumble/
