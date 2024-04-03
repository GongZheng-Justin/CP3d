#!/bin/bash
nproc=8
appName="mpirun"

cd ../../../
chmod a+x mymake.sh
./mymake.sh -exe-channel2nd -cmp-gcc_MPI -CompileThirdParty-0 -deleteCompileFile-1 -CFD_DEFS_Add-DFFTW_1D

$appName -n $nproc ./channel2nd ./Input/CFD_2nd/TestRestart/TurbCha0180_2nd.prm01
rm -rf CFD/Restart/RestartForCha180_2nd020000000000
mv CFD/Restart/RestartForCha180_2nd010000000100 CFD/Restart/RestartForCha180_2nd020000000000
$appName -n $nproc ./channel2nd ./Input/CFD_2nd/TestRestart/TurbCha0180_2nd.prm02
cp ./CFD/Results/Cha180_2nd02.log Cha180_2nd02_backup.txt
