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

# Compiling exes and doing other preparations
./mymake.sh -exe-cfd_ACM -cmp-"$cmpStr"_MPI -CompileThirdParty-0 -deleteCompileFile-1 -CFD_DEFS_Add -ACM_DEFS_Add-DChanBraunJFM2011 -CFDACM_DEFS_Add-DChanBraunJFM2011
chmod a+x ./cfd_ACM
cd ./tool/interpolateField/
chmod a+x ./makeInterp.sh
./makeInterp.sh -cmp-"$cmpStr"_MPI -deleteCompileFile-1 IsSolveScalar-0 -POST_DEFS_Add
chmod a+x ./interpolateField
cd ../../

cp ./input/CFD_ACM/ChanBraunJFM2011/SpheresCoord.ChanD50      ./ACM/Restart/SpheresCoord.bin
cp ./input/CFD_ACM/ChanBraunJFM2011/FixedSpheresCoord.ChanD50 ./ACM/Restart/FixedSpheresCoord.bin
chmod a+x ./cfd_ACM
$appStr -n $nproc ./cfd_ACM ./input/CFD_ACM/ChanBraunJFM2011/ChanBraun_D50.cfd ./input/CFD_ACM/ChanBraunJFM2011/ChanBraun_D50.acm $p_row $p_col