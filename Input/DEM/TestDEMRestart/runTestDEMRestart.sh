#!/bin/bash
nproc=8
appName="mpirun"

cd ../../../
chmod a+x mymake.sh
./mymake.sh -exe-dem -cmp-gcc_MPI -CompileThirdParty-0 -deleteCompileFile-1 -DEM_DEFS_Add-DTestDEMRestart

$appName -n $nproc ./dem ./Input/DEM/TestDEMRestart/TestDEMRestart.prm01
rm -rf DEM/Restart/RestartForSettling020000000000
mv DEM/Restart/RestartForSettling010000000100 DEM/Restart/RestartForSettling020000000000
$appName -n $nproc ./dem ./Input/DEM/TestDEMRestart/TestDEMRestart.prm02

cp ./DEM/Results/Settling02.log Settling02_backup.log
