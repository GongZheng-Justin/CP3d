#!/bin/bash
nproc=8
appName="mpirun"

cd ../../../
chmod a+x mymake.sh
./mymake.sh -exe-channelACM -cmp-gcc_MPI -CompileThirdParty-0 -deleteCompileFile-1 -CFD_DEFS_Add -ACM_DEFS_Add -CFDACM_DEFS_Add-DIBMDistributeLinear

mkdir -p ./ACM/Restart/ 2> /dev/null
cp ./Input/CFDACM/TestRestart/SpheresCoord.Ji ./ACM/Restart/SpheresCoord.dat

$appName -n $nproc ./channelACM ./Input/CFDACM/TestRestart/SphereJi.cfd01 ./Input/CFDACM/TestRestart/SphereJi.acm01
rm -rf CFD/Restart/RestartForJiSedPOF020000000000
rm -rf ACM/Restart/RestartForJiSedPOF_P0000000000
mv CFD/Restart/RestartForJiSedPOF010000000010 CFD/Restart/RestartForJiSedPOF020000000000
mv ACM/Restart/RestartForJiSedPOF_P0000000010 ACM/Restart/RestartForJiSedPOF_P0000000000
$appName -n $nproc ./channelACM ./Input/CFDACM/TestRestart/SphereJi.cfd02 ./Input/CFDACM/TestRestart/SphereJi.acm02
cp ./CFD/Results/JiSedPOF02.log JiSedPOF02_backup.txt
