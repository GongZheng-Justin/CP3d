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
mkdir -p ./DEM/Restart/ 2> /dev/null

./mymake.sh -exe-dem -cmp-"$cmpStr"_MPI -CompileThirdParty-0 -deleteCompileFile-1 -DEM_DEFS_Add
cp ./input/CFD_DEM/bedload/FixedSpheresCoord.bedload ./DEM/Restart/FixedSpheresCoord.bin
$appStr -n $nproc ./dem ./input/CFD_DEM/bedload/Bedload.dem
octave ./input/CFD_DEM/bedload/DivideBed.m
mv SpheresCoord.bin ./DEM/Restart/SpheresCoord.bin
./mymake.sh -exe-cfd_DEM -cmp-"$cmpStr"_MPI -CompileThirdParty-0 -deleteCompileFile-1 -CFD_DEFS_Add -DEM_DEFS_Add -CFDDEM_DEFS_Add

$appStr -n $nproc ./cfd_DEM ./input/CFD_DEM/bedload/cfd_DEM_bedload.01 ./input/CFD_DEM/bedload/DEMChannel_bedload.01

mv ./CFD/Restart/RestartForBedloadMove01_F0000005000 ./CFD/Restart/RestartForBedloadMove02_F0000000000
$appStr -n $nproc ./cfd_DEM ./input/CFD_DEM/bedload/cfd_DEM_bedload.02 ./input/CFD_DEM/bedload/DEMChannel_bedload.02

mv ./CFD/Restart/RestartForBedloadMove02_F0000004000 ./CFD/Restart/RestartForBedloadMove03_F0000000000
mv ./DEM/Restart/RestartForBedloadMove02_0000004000 ./DEM/Restart/RestartForBedloadMove03_0000000000
mv ./DEM/Restart/FixEdSpheresRestart0000004000 ./DEM/Restart/FixEdSpheresRestart0000000000
$appStr -n $nproc ./cfd_DEM ./input/CFD_DEM/bedload/cfd_DEM_bedload.03 ./input/CFD_DEM/bedload/DEMChannel_bedload.03