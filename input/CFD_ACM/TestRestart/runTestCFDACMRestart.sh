#!/bin/bash
nproc=8
appName="mpirun"

cd ../../../
chmod a+x mymake.sh
./mymake.sh -exe-cfd_ACM -cmp-gcc_MPI -CompileThirdParty-0 -deleteCompileFile-1 -CFD_DEFS_Add -ACM_DEFS_Add -CFDACM_DEFS_Add-DIBMDistributeLinear

mkdir -p ./ACM/Restart/ 2> /dev/null
cp ./input/CFD_ACM/TestRestart/SpheresCoord.Ji ./ACM/Restart/SpheresCoord.bin

$appName -n $nproc ./cfd_ACM ./input/CFD_ACM/TestRestart/SphereJi.cfd01 ./input/CFD_ACM/TestRestart/SphereJi.acm01
rm -rf CFD/Restart/RestartForJiSedPOF020000000000
rm -rf ACM/Restart/RestartForJiSedPOF_P0000000000
mv CFD/Restart/RestartForJiSedPOF010000000010 CFD/Restart/RestartForJiSedPOF020000000000
mv ACM/Restart/RestartForJiSedPOF_P0000000010 ACM/Restart/RestartForJiSedPOF_P0000000000
$appName -n $nproc ./cfd_ACM ./input/CFD_ACM/TestRestart/SphereJi.cfd02 ./input/CFD_ACM/TestRestart/SphereJi.acm02
cp ./CFD/Results/JiSedPOF02.log JiSedPOF02_backup.txt