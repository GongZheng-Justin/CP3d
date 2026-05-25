#!/bin/bash
#=======================================================================
# mymake.sh example --Zheng Gong, 2020-03-04(yy/mm/dd)
#=======================================================================

# The below line is needed to be modified if necessary.
SRC="./src"
CompilingLog="CompilationLog.txt"

#-----------------------------------------------------------------------
# Normally no need to change anything below.
PathCurrent="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
CompilingLog=$PathCurrent/$CompilingLog
TimeString=$(date  "+%Y-%m-%d %H:%M:%S")
rm -rf $CompilingLog; touch $CompilingLog
echo                                                                  | tee -a $CompilingLog
echo "!========================*- CP3d -*========================!"   | tee -a $CompilingLog
echo "!                                                          !"   | tee -a $CompilingLog
echo "!          CP3d:    CFD-Particle 3d                        !"   | tee -a $CompilingLog
echo "!          Version: 1.1                                    !"   | tee -a $CompilingLog
echo "!          Author:  Zheng Gong                             !"   | tee -a $CompilingLog
echo "!          E-mail:  gongzheng_justin@outlook.com           !"   | tee -a $CompilingLog
echo "!                                                          !"   | tee -a $CompilingLog
echo "!====================*- Fortran 95/03 -*===================!"   | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog
echo "  Source  Path:   "$SRC                                         | tee -a $CompilingLog
echo "  Current Path:   "$PathCurrent                                 | tee -a $CompilingLog
echo "  Compiling Time: "$TimeString                                  | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog

# Set exe name
echo "  Which EXE do you want to compile? "                           | tee -a $CompilingLog
echo "     1: dem       "                                             | tee -a $CompilingLog
echo "     2: cfd_2nd"                                                | tee -a $CompilingLog
echo "     3: cfd_4th"                                                | tee -a $CompilingLog
echo "     4: cfd_LPT"                                                | tee -a $CompilingLog
echo "     5: cfd_ATP"                                                | tee -a $CompilingLog
echo "     6: cfd_Diffuse_IBM"                                        | tee -a $CompilingLog
echo "     7: cfd_DEM"                                                | tee -a $CompilingLog
echo "     8: cfd_ACM"                                                | tee -a $CompilingLog
echo "     9: cfd_Porous_IBM"                                         | tee -a $CompilingLog
if [[ -n $1 ]]; then
  strTemp=$1
  EXE=${strTemp:5}
else
  read -p "  Please type a EXE index(1-9): " id_exe
  echo    "  Please type a EXE index(1-9): "$id_exe >> $CompilingLog
  if [ "$id_exe" == 1 ]; then
    EXE="dem"
  elif [ "$id_exe" == 2 ]; then
    EXE="cfd_2nd"
  elif [ "$id_exe" == 3 ]; then
    EXE="cfd_4th"
  elif [ "$id_exe" == 4 ]; then
    EXE="cfd_LPT" 
  elif [ "$id_exe" == 5 ]; then
    EXE="cfd_ATP"
  elif [ "$id_exe" == 6 ]; then
    EXE="cfd_Diffuse_IBM"
  elif [ "$id_exe" == 7 ]; then
    EXE="cfd_DEM"
  elif [ "$id_exe" == 8 ]; then
    EXE="cfd_ACM"
  elif [ "$id_exe" == 9 ]; then
    EXE="cfd_Porous_IBM"
  else
    echo "  Sorry, EXE type cannot be recognized."                    | tee -a $CompilingLog
    echo "  Compiling filed"                                          | tee -a $CompilingLog
    exit 1
  fi
fi
echo "  "$EXE"  will be compiled"                                     | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog

# Set compiler
echo "  Which compiler do you use? "                                  | tee -a $CompilingLog
echo "     1: Intel MPI (mpiifort)"                                   | tee -a $CompilingLog
echo "     2: gcc MPI   (mpif90). Default"                            | tee -a $CompilingLog
if [[ -n $2 ]]; then
  strTemp=$2
  CMP=${strTemp:5}
else
  read -p "  Please type a compiler index (1 or 2): " id_cmp
  echo    "  Please type a compiler index (1 or 2): "$id_cmp >> $CompilingLog
  if [ "$id_cmp" == 1 ]; then
    CMP="intel_MPI"
  else
    CMP="gcc_MPI"
  fi
fi
echo "  "$CMP"  will be used"                                         | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog

# Compile ThirdParty or not
echo "  Do you want to recompile ThirdParty? "                        | tee -a $CompilingLog
echo "     0: No need, use the old ThirdParty compilation. Default"   | tee -a $CompilingLog
echo "     1: Yes, recompile ThirdParty (Recommended for first use)"  | tee -a $CompilingLog
if [[ -n $3 ]]; then
  strTemp=$3
  Third_flag=${strTemp:0-1}
  echo    "  Please type a choice (0 or 1): "$Third_flag              | tee -a $CompilingLog
else
  read -p "  Please type a choice (0 or 1): " Third_flag
  echo    "  Please type a choice (0 or 1): "$Third_flag >> $CompilingLog
fi
if [ "$Third_flag" == 1 ]; then
  echo                                                                | tee -a $CompilingLog
  cd $SRC/ThirdParty
  echo "Compiling ThirdParty begins."                                 | tee -a $CompilingLog
  chmod a+x ./install_thirdParty.sh
  ./install_thirdParty.sh
  echo "Compiling ThirdParty done !"                                  | tee -a $CompilingLog
  echo                                                                | tee -a $CompilingLog
  cd ../..
else
  echo "  Choose to use the old ThirdParty compilation"               | tee -a $CompilingLog
fi
echo                                                                  | tee -a $CompilingLog

# Delete temporary compiling files or not
echo "  Do you want to delete temporary compiling files? "            | tee -a $CompilingLog
echo "     0: No, save them. "                                        | tee -a $CompilingLog
echo "     1: Yes,delete them. Default"                               | tee -a $CompilingLog
if [[ -n $4 ]]; then
  strTemp=$4
  DeleteFlag=${strTemp:0-1}
  echo    "  Please type a choice (0 or 1): "$DeleteFlag              | tee -a $CompilingLog
else
  read -p "  Please type a choice (0 or 1): " DeleteFlag
  echo    "  Please type a choice (0 or 1): "$DeleteFlag >> $CompilingLog
fi
if [ "$DeleteFlag" != 0 ]; then
  echo "  Choose to DELETE temporary compiling files"                 | tee -a $CompilingLog
else
  echo "  Choose to SAVE temporary compiling files"                   | tee -a $CompilingLog
fi
echo                                                                  | tee -a $CompilingLog

# Compile begins
echo  "!==================*- Compiling begins -*=================!"   | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog
rm -fr $EXE
cd $SRC
if [ "$DeleteFlag" != 0 ]; then
  make -f $EXE".make" clean >&/dev/null
fi
if [ "$EXE" == "dem" ]; then
  DEM_DEFS_Add=""
  if [[ -n $5 ]]; then
    strTemp=$5
    DEM_DEFS_Add=${strTemp:13}
  fi  
  make -f $EXE".make" CMP=$CMP exeName=$EXE DEM_DEFS_Add=$DEM_DEFS_Add 2>&1 | tee -a $CompilingLog
elif [ "$EXE" == "cfd_2nd" ]; then
  CFD_DEFS_Add=""
  if [[ -n $5 ]]; then
    strTemp=$5
    CFD_DEFS_Add=${strTemp:13}
  fi  
  make -f $EXE".make" CMP=$CMP exeName=$EXE CFD_DEFS_Add=$CFD_DEFS_Add 2>&1 | tee -a $CompilingLog
elif [ "$EXE" == "cfd_4th" ]; then
  CFD_DEFS_Add=""
  IsSolveScalar=""
  if [[ -n $5 ]]; then
    strTemp=$5
    CFD_DEFS_Add=${strTemp:13}
  fi  
  if [[ -n $6 ]]; then
    strTemp=$6
    IsSolveScalar=${strTemp:15}
  fi
  if [ "$IsSolveScalar" == "" ]; then
    make -f $EXE".make" CMP=$CMP exeName=$EXE CFD_DEFS_Add=$CFD_DEFS_Add 2>&1 | tee -a $CompilingLog  
  else
    make -f $EXE".make" CMP=$CMP exeName=$EXE CFD_DEFS_Add=$CFD_DEFS_Add IsSolveScalar=$IsSolveScalar 2>&1 | tee -a $CompilingLog  
  fi
elif [ "$EXE" == "cfd_LPT" ]; then
  CFD_DEFS_Add=""
  LPT_DEFS_Add=""
  CFDOrder="2"
  Coupling="2"
  if [[ -n $5 ]]; then
    strTemp=$5
    CFD_DEFS_Add=${strTemp:13}
  fi  
  if [[ -n $6 ]]; then
    strTemp=$6
    LPT_DEFS_Add=${strTemp:13}
  fi 
  if [[ -n $7 ]]; then
    strTemp=$7
    CFDOrder=${strTemp:10}
  fi
  if [[ -n $8 ]]; then
    strTemp=$8
    Coupling=${strTemp:10}
  fi
  echo 
  make -f $EXE".make" CMP=$CMP exeName=$EXE CFD_DEFS_Add=$CFD_DEFS_Add LPT_DEFS_Add=$LPT_DEFS_Add \
  CFDOrder=$CFDOrder Coupling=$Coupling 2>&1 | tee -a $CompilingLog

elif [ "$EXE" == "cfd_ATP" ]; then
  CFD_DEFS_Add=""
  LPT_DEFS_Add=""
  CFDOrder="2"
  if [[ -n $5 ]]; then
    strTemp=$5
    CFD_DEFS_Add=${strTemp:13}
  fi  
  if [[ -n $6 ]]; then
    strTemp=$6
    LPT_DEFS_Add=${strTemp:13}
  fi 
  if [[ -n $7 ]]; then
    strTemp=$7
    CFDOrder=${strTemp:10}
  fi
  echo 
  make -f $EXE".make" CMP=$CMP exeName=$EXE CFD_DEFS_Add=$CFD_DEFS_Add ATP_DEFS_Add=$LPT_DEFS_Add \
  CFDOrder=$CFDOrder 2>&1 | tee -a $CompilingLog

elif [ "$EXE" == "cfd_Porous_IBM" ]; then
  CFD_DEFS_Add=""
  IBM_DEFS_Add=""
  CFDIBM_DEFS_Add=""
  if [[ -n $5 ]]; then
    strTemp=$5
    CFD_DEFS_Add=${strTemp:13}
  fi
  if [[ -n $6 ]]; then
    strTemp=$6
    IBM_DEFS_Add=${strTemp:13}
  fi
  if [[ -n $7 ]]; then
    strTemp=$7
    CFDIBM_DEFS_Add=${strTemp:16}
  fi
  make -f $EXE".make" CMP=$CMP exeName=$EXE CFD_DEFS_Add=$CFD_DEFS_Add IBM_DEFS_Add=$IBM_DEFS_Add \
  2>&1                                                                | tee -a $CompilingLog
      
elif [ "$EXE" == "cfd_DEM" ]; then
  CFD_DEFS_Add=""
  DEM_DEFS_Add=""
  CFDDEM_DEFS_Add=""
  if [[ -n $5 ]]; then
    strTemp=$5
    CFD_DEFS_Add=${strTemp:13}
  fi
  if [[ -n $6 ]]; then
    strTemp=$6
    DEM_DEFS_Add=${strTemp:13}
  fi
  if [[ -n $7 ]]; then
    strTemp=$7
    CFDDEM_DEFS_Add=${strTemp:16}
  fi
  make -f $EXE".make" CMP=$CMP exeName=$EXE CFD_DEFS_Add=$CFD_DEFS_Add DEM_DEFS_Add=$DEM_DEFS_Add \
  CFDDEM_DEFS_Add=$CFDDEM_DEFS_Add 2>&1                               | tee -a $CompilingLog

elif [ "$EXE" == "cfd_ACM" ]; then
  CFD_DEFS_Add=""
  ACM_DEFS_Add=""
  CFDACM_DEFS_Add=""
  if [[ -n $5 ]]; then
    strTemp=$5
    CFD_DEFS_Add=${strTemp:13}
  fi
  if [[ -n $6 ]]; then
    strTemp=$6
    ACM_DEFS_Add=${strTemp:13}
  fi
  if [[ -n $7 ]]; then
    strTemp=$7
    CFDACM_DEFS_Add=${strTemp:16}
  fi
  make -f $EXE".make" CMP=$CMP exeName=$EXE CFD_DEFS_Add=$CFD_DEFS_Add ACM_DEFS_Add=$ACM_DEFS_Add \
  CFDACM_DEFS_Add=$CFDACM_DEFS_Add 2>&1                               | tee -a $CompilingLog

elif [ "$EXE" == "cfd_Diffuse_IBM" ]; then
  CFD_DEFS_Add=""
  IBM_DEFS_Add=""
  CFDIBM_DEFS_Add=""
  if [[ -n $5 ]]; then
    strTemp=$5
    CFD_DEFS_Add=${strTemp:13}
  fi
  if [[ -n $6 ]]; then
    strTemp=$6
    IBM_DEFS_Add=${strTemp:13}
  fi
  if [[ -n $7 ]]; then
    strTemp=$7
    CFDIBM_DEFS_Add=${strTemp:16}
  fi
  make -f $EXE".make" CMP=$CMP exeName=$EXE CFD_DEFS_Add=$CFD_DEFS_Add IBM_DEFS_Add=$IBM_DEFS_Add \
  2>&1                                                                | tee -a $CompilingLog
  
else
  echo  $EXE" wrong, please check !!!"                                | tee -a $CompilingLog
fi
echo                                                                  | tee -a $CompilingLog
mv $EXE $PathCurrent                              
if [ $? -ne 0 ]; then
  if [ "$DeleteFlag" != 0 ]; then
    make -f $EXE".make" clean >&/dev/null
  fi
  echo  $EXE" CANNOT be compiled correctly, please check !!!"         | tee -a $CompilingLog
else
  if [ "$DeleteFlag" != 0 ]; then
    make -f $EXE".make" clean >&/dev/null
  fi
  echo  $EXE" has been compiled normally. Enjoy it !!!"               | tee -a $CompilingLog
  cd ..
  chmod a+x ./$EXE
fi
echo                                                                  | tee -a $CompilingLog
echo  "!===================*- Compiling ends -*==================!"   | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog
