#!/bin/bash
npoc1=1
npoc2=4
npoc3=8
appName="mpirun"

cd ../../../
chmod a+x mymake.sh
mkdir -p ./ACM/Restart/ 2> /dev/null
./mymake.sh -exe-channelACM -cmp-gcc_MPI -CompileThirdParty-0 -deleteCompileFile-1 -CFD_DEFS_Add -ACM_DEFS_Add -CFDACM_DEFS_Add-DSeveralSphereInfo

cp -rf ./Input/CFDACM/SingleStokes/SpheresCoord.Case01 ./ACM/Restart/SpheresCoord.dat
$appName -n $npoc1 ./channelACM ./Input/CFDACM/SingleStokes/SingleStokesCase01.cfd ./Input/CFDACM/SingleStokes/SingleStokesCase01.acm

cp -rf ./Input/CFDACM/SingleStokes/SpheresCoord.Case02 ./ACM/Restart/SpheresCoord.dat
$appName -n $npoc1 ./channelACM ./Input/CFDACM/SingleStokes/SingleStokesCase02.cfd ./Input/CFDACM/SingleStokes/SingleStokesCase02.acm

cp -rf ./Input/CFDACM/SingleStokes/SpheresCoord.Case03 ./ACM/Restart/SpheresCoord.dat
$appName -n $npoc1 ./channelACM ./Input/CFDACM/SingleStokes/SingleStokesCase03.cfd ./Input/CFDACM/SingleStokes/SingleStokesCase03.acm

cp -rf ./Input/CFDACM/SingleStokes/SpheresCoord.Case04 ./ACM/Restart/SpheresCoord.dat
$appName -n $npoc1 ./channelACM ./Input/CFDACM/SingleStokes/SingleStokesCase04.cfd ./Input/CFDACM/SingleStokes/SingleStokesCase04.acm

cp -rf ./Input/CFDACM/SingleStokes/SpheresCoord.Case05 ./ACM/Restart/SpheresCoord.dat
$appName -n $npoc1 ./channelACM ./Input/CFDACM/SingleStokes/SingleStokesCase05.cfd ./Input/CFDACM/SingleStokes/SingleStokesCase05.acm

cp -rf ./Input/CFDACM/SingleStokes/SpheresCoord.Case06 ./ACM/Restart/SpheresCoord.dat
$appName -n $npoc1 ./channelACM ./Input/CFDACM/SingleStokes/SingleStokesCase06.cfd ./Input/CFDACM/SingleStokes/SingleStokesCase06.acm

cp -rf ./Input/CFDACM/SingleStokes/SpheresCoord.Case07 ./ACM/Restart/SpheresCoord.dat
$appName -n $npoc1 ./channelACM ./Input/CFDACM/SingleStokes/SingleStokesCase07.cfd ./Input/CFDACM/SingleStokes/SingleStokesCase07.acm

cp -rf ./Input/CFDACM/SingleStokes/SpheresCoord.Case08 ./ACM/Restart/SpheresCoord.dat
$appName -n $npoc1 ./channelACM ./Input/CFDACM/SingleStokes/SingleStokesCase08.cfd ./Input/CFDACM/SingleStokes/SingleStokesCase08.acm

cp -rf ./Input/CFDACM/SingleStokes/SpheresCoord.Case09 ./ACM/Restart/SpheresCoord.dat
$appName -n $npoc1 ./channelACM ./Input/CFDACM/SingleStokes/SingleStokesCase09.cfd ./Input/CFDACM/SingleStokes/SingleStokesCase09.acm

cp -rf ./Input/CFDACM/SingleStokes/SpheresCoord.Case10 ./ACM/Restart/SpheresCoord.dat
$appName -n $npoc1 ./channelACM ./Input/CFDACM/SingleStokes/SingleStokesCase10.cfd ./Input/CFDACM/SingleStokes/SingleStokesCase10.acm

cp -rf ./Input/CFDACM/SingleStokes/SpheresCoord.Case11 ./ACM/Restart/SpheresCoord.dat
$appName -n $npoc2 ./channelACM ./Input/CFDACM/SingleStokes/SingleStokesCase11.cfd ./Input/CFDACM/SingleStokes/SingleStokesCase11.acm

cp -rf ./Input/CFDACM/SingleStokes/SpheresCoord.Case12 ./ACM/Restart/SpheresCoord.dat
$appName -n $npoc2 ./channelACM ./Input/CFDACM/SingleStokes/SingleStokesCase12.cfd ./Input/CFDACM/SingleStokes/SingleStokesCase12.acm

cp -rf ./Input/CFDACM/SingleStokes/SpheresCoord.Case13 ./ACM/Restart/SpheresCoord.dat
$appName -n $npoc3 ./channelACM ./Input/CFDACM/SingleStokes/SingleStokesCase13.cfd ./Input/CFDACM/SingleStokes/SingleStokesCase13.acm

cd ./Input/CFDACM/SingleStokes
