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

./mymake.sh -exe-channelACM -cmp-"$cmpStr"_MPI -CompileThirdParty-0 -deleteCompileFile-1 -CFD_DEFS_Add -ACM_DEFS_Add -CFDACM_DEFS_Add

# Fixed Case
cp -rf ./Input/CFDACM/JiPoF2013/FixedSpheresCoord.dat ./ACM/Restart/FixedSpheresCoord.dat
$appStr -n $nproc ./channelACM ./Input/CFDACM/JiPoF2013/FixedJi.cfd ./Input/CFDACM/JiPoF2013/FixedJi.acm $p_row $p_col

# Move Case
cp -rf ./Input/CFDACM/JiPoF2013/SpheresCoord.fixed ./ACM/Restart/FixedSpheresCoord.dat
cp -rf ./Input/CFDACM/JiPoF2013/SpheresCoord.move ./ACM/Restart/SpheresCoord.dat
mv ./CFD/Restart/RestartForFixedJiPOF0000500000 ./CFD/Restart/RestartForMoveJiPOF0000000000
$appStr -n $nproc ./channelACM ./Input/CFDACM/JiPoF2013/MoveJi.cfd ./Input/CFDACM/JiPoF2013/MoveJi.acm $p_row $p_col