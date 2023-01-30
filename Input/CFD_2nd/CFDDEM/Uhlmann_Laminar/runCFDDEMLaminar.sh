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
cp ./Input/CFDDEM/Uhlmann_Laminar/FixedSpheresCoord.UhlmannLaminar ./DEM/Restart/FixedSpheresCoord.dat
$appStr -n $nproc ./dem ./Input/CFDDEM/Uhlmann_Laminar/Uhlmann_Laminar.dem
octave ./Input/CFDDEM/Uhlmann_Laminar/DivideBed.m
mv SpheresCoord.dat ./DEM/Restart/SpheresCoord.dat
./mymake.sh -exe-channelDEM -cmp-"$cmpStr"_MPI -CompileThirdParty-0 -deleteCompileFile-1 -CFD_DEFS_Add -DEM_DEFS_Add -CFDDEM_DEFS_Add

for loop in 01 02 03 04 05 06 07 08 09 10
do
  $appStr -n $nproc ./channelDEM ./Input/CFDDEM/Uhlmann_Laminar/ChannelDEM_Laminar.case"$loop" ./Input/CFDDEM/Uhlmann_Laminar/DEMChannel_Laminar.case"$loop"
done
