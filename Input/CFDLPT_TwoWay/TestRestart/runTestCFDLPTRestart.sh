#!/bin/bash
nproc=8
appName="mpirun"

cd ../../../
chmod a+x mymake.sh
./mymake.sh -exe-channelLPT -cmp-gcc_MPI -CompileThirdParty-0 -deleteCompileFile-1 -CFD_DEFS_Add-DFFTW_1D -LPT_DEFS_Add -CFDOrder-2 -Coupling-2

$appName -n $nproc ./channelLPT ./Input/CFDLPT_TwoWay/TestRestart/Channel4th_LPT.prm01 ./Input/CFDLPT_TwoWay/TestRestart/LPT_Channel4th.prm01
rm -rf CFD/Restart/RestartForTowWayLee020000000000
rm -rf LPT/Restart/RestartForTowWayLee_P0000000000
mv CFD/Restart/RestartForTowWayLee010000000100 CFD/Restart/RestartForTowWayLee020000000000
mv LPT/Restart/RestartForTowWayLee_P0000000100 LPT/Restart/RestartForTowWayLee_P0000000000
$appName -n $nproc ./channelLPT ./Input/CFDLPT_TwoWay/TestRestart/Channel4th_LPT.prm02 ./Input/CFDLPT_TwoWay/TestRestart/LPT_Channel4th.prm02
cp ./CFD/Results/TowWayLee02.log TowWayLee02_backup.txt
