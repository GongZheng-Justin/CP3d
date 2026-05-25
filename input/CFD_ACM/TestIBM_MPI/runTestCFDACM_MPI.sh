#!/bin/bash
nproc=8
appName="mpirun"

cd ../../../
chmod a+x mymake.sh
./mymake.sh -exe-cfd_ACM -cmp-gcc_MPI -CompileThirdParty-0 -deleteCompileFile-1 -CFD_DEFS_Add -ACM_DEFS_Add -CFDACM_DEFS_Add-DTest_IBM_MPI

mkdir -p ./ACM/Restart/ 2> /dev/null
cp ./input/CFD_ACM/TestIBM_MPI/SpheresCoord.test ./ACM/Restart/SpheresCoord.bin

$appName -n $nproc ./cfd_ACM ./input/CFD_ACM/TestIBM_MPI/SphereTest.cfd ./input/CFD_ACM/TestIBM_MPI/SphereTest.acm