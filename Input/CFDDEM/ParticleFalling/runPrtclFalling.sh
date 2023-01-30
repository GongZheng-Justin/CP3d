#!/bin/bash
cmpStr="gcc"
appStr="mpirun"
p_row=1
p_col=1
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

cp ./Input/CFDDEM/ParticleFalling/SpheresCoord.falling ./DEM/Restart/SpheresCoord.dat
./mymake.sh -exe-channelDEM -cmp-"$cmpStr"_MPI -CompileThirdParty-0 -deleteCompileFile-1 -CFD_DEFS_Add -DEM_DEFS_Add -CFDDEM_DEFS_Add-DSeveralSphereInfo

for loop in 01 02 04 05 06 08
do
  $appStr -n $nproc ./channelDEM ./Input/CFDDEM/ParticleFalling/ChannelDEM_falling.case"$loop" ./Input/CFDDEM/ParticleFalling/DEMChannel_falling.case"$loop"
done
