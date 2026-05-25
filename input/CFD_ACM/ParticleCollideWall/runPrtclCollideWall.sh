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

./mymake.sh -exe-cfd_ACM -cmp-"$cmpStr"_MPI -CompileThirdParty-0 -deleteCompileFile-1 -CFD_DEFS_Add -ACM_DEFS_Add -CFDACM_DEFS_Add-DSeveralSphereInfo

# Experiment
for loop in 01 02 03 04 05 06 07
do 
  cp ./input/CFD_ACM/ParticleCollideWall/SpheresCoord.Case"$loop" ./ACM/Restart/SpheresCoord.bin
  $appStr -n $nproc ./cfd_ACM ./input/CFD_ACM/ParticleCollideWall/PrtclCollideCase"$loop".cfd ./input/CFD_ACM/ParticleCollideWall/PrtclCollide.acm $p_row $p_col
done

# No lubrication test
for loop in 08 09
do 
  cp ./input/CFD_ACM/ParticleCollideWall/SpheresCoord.Case"$loop" ./ACM/Restart/SpheresCoord.bin
  $appStr -n $nproc ./cfd_ACM ./input/CFD_ACM/ParticleCollideWall/PrtclCollideCase"$loop".cfd ./input/CFD_ACM/ParticleCollideWall/PrtclCollide.acm0 $p_row $p_col
done