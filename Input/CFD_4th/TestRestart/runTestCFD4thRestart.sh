#!/bin/bash
nproc=8
appName="mpirun"

cd ../../../
chmod a+x mymake.sh
./mymake.sh -exe-channel4th -cmp-gcc_MPI -CompileThirdParty-0 -deleteCompileFile-1 -CFD_DEFS_Add-DFFTW_1D -IsSolveScalar-1

$appName -n $nproc ./channel4th ./Input/CFD_4th/TestRestart/SF_Ri018.prm01
rm -rf CFD/Restart/RestartForRi018_020000000000
mv CFD/Restart/RestartForRi018_010000000100 CFD/Restart/RestartForRi018_020000000000
$appName -n $nproc ./channel4th ./Input/CFD_4th/TestRestart/SF_Ri018.prm02
cp ./CFD/Results/Ri018_02.log Ri018_02_backup.txt
