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

# Re=0.1
cp ./input/CFD_ACM/TwoParticleInShearFlow/SpheresCoord.Case1 ./ACM/Restart/SpheresCoord.bin
./mymake.sh -exe-cfd_ACM -cmp-"$cmpStr"_MPI -CompileThirdParty-0 -deleteCompileFile-1 -CFD_DEFS_Add -ACM_DEFS_Add-DRotateOnly -CFDACM_DEFS_Add-DSeveralSphereInfo
$appStr -n $nproc ./cfd_ACM ./input/CFD_ACM/TwoParticleInShearFlow/TwoPrtclShearCase1.cfd0 ./input/CFD_ACM/TwoParticleInShearFlow/TwoPrtclShear.acm0 $p_row $p_col

./mymake.sh -exe-cfd_ACM -cmp-"$cmpStr"_MPI -CompileThirdParty-0 -deleteCompileFile-1 -CFD_DEFS_Add -ACM_DEFS_Add -CFDACM_DEFS_Add-DSeveralSphereInfo
$appStr -n $nproc ./cfd_ACM ./input/CFD_ACM/TwoParticleInShearFlow/TwoPrtclShearCase1.cfd1 ./input/CFD_ACM/TwoParticleInShearFlow/TwoPrtclShear.acm1 $p_row $p_col


# Re=0.2
cp ./input/CFD_ACM/TwoParticleInShearFlow/SpheresCoord.Case2 ./ACM/Restart/SpheresCoord.bin
./mymake.sh -exe-cfd_ACM -cmp-"$cmpStr"_MPI -CompileThirdParty-0 -deleteCompileFile-1 -CFD_DEFS_Add -ACM_DEFS_Add-DRotateOnly -CFDACM_DEFS_Add-DSeveralSphereInfo
$appStr -n $nproc ./cfd_ACM ./input/CFD_ACM/TwoParticleInShearFlow/TwoPrtclShearCase2.cfd0 ./input/CFD_ACM/TwoParticleInShearFlow/TwoPrtclShear.acm0 $p_row $p_col

./mymake.sh -exe-cfd_ACM -cmp-"$cmpStr"_MPI -CompileThirdParty-0 -deleteCompileFile-1 -CFD_DEFS_Add -ACM_DEFS_Add -CFDACM_DEFS_Add-DSeveralSphereInfo
$appStr -n $nproc ./cfd_ACM ./input/CFD_ACM/TwoParticleInShearFlow/TwoPrtclShearCase2.cfd1 ./input/CFD_ACM/TwoParticleInShearFlow/TwoPrtclShear.acm1 $p_row $p_col