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
echo "  Source  Path:   "$SRC                                         | tee -a $CompilingLog
echo "  Current Path:   "$PathCurrent                                 | tee -a $CompilingLog
echo "  Compiling Time: "$TimeString                                  | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog

# Set exe name
EXE="interpolateField"

# Set compiler
echo "  Which compiler do you use? "                                  | tee -a $CompilingLog
echo "     1: Intel MPI (mpiifort)"                                   | tee -a $CompilingLog
echo "     2: gcc MPI   (mpif90). Default"                            | tee -a $CompilingLog
if [[ -n $1 ]]; then
  strTemp=$1
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

# Delete temporary compiling files or not
echo "  Do you want to delete temporary compiling files? "            | tee -a $CompilingLog
echo "     0: No, save them. "                                        | tee -a $CompilingLog
echo "     1: Yes,delete them. Default"                               | tee -a $CompilingLog
if [[ -n $2 ]]; then
  strTemp=$2
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

# With OR W/O  Scalar
echo "  Which EXE do you want to compile? "                           | tee -a $CompilingLog
echo "     0: interpolateField, W/O  Scalar. Default"                 | tee -a $CompilingLog
echo "     1: interpolateField, With Scalar"                          | tee -a $CompilingLog
if [[ -n $3 ]]; then
  strTemp=$3
  IsSolveScalar=${strTemp:15}
  echo    "  Please type a choice (0 or 1): "$IsSolveScalar           | tee -a $CompilingLog
else
  read -p "  please type a EXE index(0 or 1): " id_exe
  echo    "  please type a EXE index(0 or 1): "$id_exe >> $CompilingLog
  if [ "$id_exe" == 1 ]; then
    IsSolveScalar="1"
  else
    IsSolveScalar="0"
  fi
fi
if [ "$IsSolveScalar" == 1 ]; then
  EXEStr="interpolateField, With Scalar"
else
  EXEStr="interpolateField, W/O  Scalar"
fi
echo "  "$EXEStr"  will be compiled"                                  | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog
                                                                
echo  "!==================*- Compiling begins -*=================!"   | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog
rm -fr $EXE
cd $SRC
if [ "$DeleteFlag" != 0 ]; then
  make clean >&/dev/null 
fi
POST_DEFS_Add=""
if [[ -n $4 ]]; then
  strTemp=$4
  POST_DEFS_Add=${strTemp:14}
fi
make CMP=$CMP exeName=$EXE IsSolveScalar=$IsSolveScalar POST_DEFS_Add=$POST_DEFS_Add 2>&1 | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog
mv $EXE $PathCurrent                              
if [ $? -ne 0 ]; then
  if [ "$DeleteFlag" != 0 ]; then
    make clean >&/dev/null
  fi
  echo  $EXE" CANNOT be compiled correctly, please check !!!"         | tee -a $CompilingLog
else
  if [ "$DeleteFlag" != 0 ]; then
    make clean >&/dev/null
  fi
  echo  $EXE" has been compiled normally. Enjoy !!!"                  | tee -a $CompilingLog
  cd ..
  chmod a+x ./$EXE
fi
echo                                                                  | tee -a $CompilingLog
echo  "!===================*- Compiling ends -*==================!"   | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog
