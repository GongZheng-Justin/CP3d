#!/bin/bash
nproc=8
appName="mpirun"

cd ../../../
chmod a+x mymake.sh
./mymake.sh -exe-channelACM -cmp-gcc_MPI -CompileThirdParty-0 -deleteCompileFile-1 -CFD_DEFS_Add -ACM_DEFS_Add -CFDACM_DEFS_Add-DTest_IBM_MPI

mkdir -p ./ACM/Restart/ 2> /dev/null
cp ./Input/CFDACM/TestIBM_MPI/SpheresCoord.test ./ACM/Restart/SpheresCoord.dat

$appName -n $nproc ./channelACM ./Input/CFDACM/TestIBM_MPI/SphereTest.cfd ./Input/CFDACM/TestIBM_MPI/SphereTest.acm
