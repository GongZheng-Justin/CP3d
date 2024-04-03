#!/bin/bash
#=======================================================================
# mymake.sh example --Zheng Gong, 2020-03-04(yy/mm/dd)
#=======================================================================

# The below line is needed to be modified if necessary.
SRC="./src_post/"

#-----------------------------------------------------------------------
# Normally no need to change anything below
PathCurrent="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

echo "  Which EXE do you want to compile? " 
echo "     1: readPrtclDump" 
echo "     2: chooseSaltationFile"
echo "     3: clcCntctFlagSeq"
read -p "  please type a EXE index(1, 2, 3, 4 or 5): " id_exe

if [ $id_exe -eq 1 ]; then
  CMP="gcc_MPI"
  EXE="readPrtclDump"
elif [ $id_exe -eq 2 ]; then
  CMP="gcc_MPI"
  EXE="chooseSaltationFile"
elif [ $id_exe -eq 3 ]; then
  CMP="gcc_MPI"
  EXE="clcCntctFlagSeq"
else
  echo "  Sorry, EXE type cannot be recognized." 
  echo "  Compiling filed" 
  exit 2
fi
echo "  "$EXE"  will be compiled"                                
echo                                                                  
echo  "!==================*- Compiling begins -*=================!"  
echo
rm -fr $EXE
cd $SRC
make -f "make_"$EXE clean >&/dev/null 
make -f "make_"$EXE CMP=$CMP   2>&1                       
echo                         
mv $EXE $PathCurrent                              
if [ $? -ne 0 ]; then
  make -f "make_"$EXE clean >&/dev/null
  echo  $EXE" CANNOT be compiled correctly, please check !!!"  
else
  make -f "make_"$EXE clean >&/dev/null
  echo  $EXE" has been compiled normally. Enjoy !!!" 
  cd ..
  chmod a+x ./$EXE
fi
echo                                                                 
echo  "!===================*- Compiling ends -*==================!"
echo
