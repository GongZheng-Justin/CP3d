#!/bin/bash
npoc1=1
npoc2=4
npoc3=8
appName="mpirun"

cd ../../../
chmod a+x mymake.sh
mkdir -p ./ACM/Restart/ 2> /dev/null
./mymake.sh -exe-cfd_ACM -cmp-gcc_MPI -CompileThirdParty-0 -deleteCompileFile-1 -CFD_DEFS_Add -ACM_DEFS_Add -CFDACM_DEFS_Add-DSeveralSphereInfo

for loop in 01 02 03 04 05 06 07 08 09 10
do
  cp -rf ./input/CFD_ACM/SingleStokes/SpheresCoord.Case"$loop" ./ACM/Restart/SpheresCoord.bin
  $appName -n $npoc1 ./cfd_ACM ./input/CFD_ACM/SingleStokes/SingleStokesCase"$loop".cfd ./input/CFD_ACM/SingleStokes/SingleStokesCase"$loop".acm
done

for loop in 11 12
do
  cp -rf ./input/CFD_ACM/SingleStokes/SpheresCoord.Case"$loop" ./ACM/Restart/SpheresCoord.bin
  $appName -n $npoc2 ./cfd_ACM ./input/CFD_ACM/SingleStokes/SingleStokesCase"$loop".cfd ./input/CFD_ACM/SingleStokes/SingleStokesCase"$loop".acm
done

for loop in 13
do
  cp -rf ./input/CFD_ACM/SingleStokes/SpheresCoord.Case"$loop" ./ACM/Restart/SpheresCoord.bin
  $appName -n $npoc3 ./cfd_ACM ./input/CFD_ACM/SingleStokes/SingleStokesCase"$loop".cfd ./input/CFD_ACM/SingleStokes/SingleStokesCase"$loop".acm
done

cd ./input/CFD_ACM/SingleStokes