#!/bin/bash
#=======================================================================
# mymake.sh example --Zheng Gong, 2020-03-04(yy/mm/dd)
#=======================================================================

# The below line is needed to be modified if necessary.
SRC="./src_PPD/"

#-----------------------------------------------------------------------
# Normally no need to change anything below
PathCurrent="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

echo "  Which compiler do you use? "                                
echo "     1: Intel MPI    ( mpiifort )"                               
echo "     2: gcc MPI      ( mpif90   )"  
echo "     3: Intel serial ( ifort    )"
echo "     4: gcc serial   ( gfortran )"            
read -p "  please type a compiler index (1, 2, 3 or 4): " id_cmp
#read -p "  please type a compiler index (1 or 2): " id_cmp

if [ $id_cmp -eq 1 ]; then
  CMP="intel_MPI"
elif [ $id_cmp -eq 2 ]; then
  CMP="gcc_MPI"
elif [ $id_cmp -eq 3 ]; then
  CMP="intel_serial"
elif [ $id_cmp -eq 4 ]; then
  CMP="gcc_serial"
else
  echo "  Sorry, compiler type cannot be recognized." 
  echo "  Compiling filed"
  exit 1
fi
echo "  "$CMP"  will be used"

echo "  Which EXE do you want to compile? " 
echo "     1: Post Particle Dump"
read -p "  please type a EXE index(1): " id_exe

if [ $id_exe -eq 1 ]; then
  EXE="PPD"
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
