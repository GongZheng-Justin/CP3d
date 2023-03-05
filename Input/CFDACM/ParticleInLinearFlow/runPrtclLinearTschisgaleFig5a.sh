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
mkdir -p ./ACM/Restart/ 2> /dev/null

# Fig.5a
./mymake.sh -exe-channelACM -cmp-"$cmpStr"_MPI -CompileThirdParty-0 -deleteCompileFile-1 -CFD_DEFS_Add -ACM_DEFS_Add-DRotateOnly -CFDACM_DEFS_Add-DSeveralSphereInfo

cp -rf ./Input/CFDACM/ParticleInLinearFlow/SpheresCoord.TschisgaleFig5 ./ACM/Restart/SpheresCoord.dat
$appStr -n $nproc ./channelACM ./Input/CFDACM/ParticleInLinearFlow/PrtclLinearTschisgaleFig5a.cfdT1 ./Input/CFDACM/ParticleInLinearFlow/PrtclLinearTschisgaleFig5a.acmT1 $p_row $p_col
$appStr -n $nproc ./channelACM ./Input/CFDACM/ParticleInLinearFlow/PrtclLinearTschisgaleFig5a.cfdT2 ./Input/CFDACM/ParticleInLinearFlow/PrtclLinearTschisgaleFig5a.acmT2 $p_row $p_col


for loop in 1 2 3
do
  cp -rf ./CFD/Restart/RestartForTschiFig05aT0000015000 ./CFD/Restart/RestartForTschiFig05a"$loop"0000000000
  $appStr -n $nproc ./channelACM ./Input/CFDACM/ParticleInLinearFlow/PrtclLinearTschisgaleFig5a.cfd"$loop" ./Input/CFDACM/ParticleInLinearFlow/PrtclLinearTschisgaleFig5a.acm"$loop" $p_row $p_col
done
