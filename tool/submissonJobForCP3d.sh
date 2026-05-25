#!/bin/bash
#=======================================================================
# Job submission example for E5 HPC system of Tsinghua University
# Zheng Gong, 2020-11-12(yy/mm/dd)
#=======================================================================

#-----------------------------------------------------------------------
# Normally no need to change anything below

jobForCP3dName="jobForCP3d"
rm -rf $jobForCP3dName; touch $jobForCP3dName

echo "  Which EXE do you want to run? "
echo "     1: dem       "
echo "     2: channel3d "
echo "     3: channelLPT"
echo "     4: channelDEM"
echo "     5: channelACM"
echo "     6: channel3d_4th"
read -p "  please type a EXE index(1, 2, 3, 4 5 or 6): " id_exe

if [ $id_exe -eq 1 ]; then
  EXE="dem"
elif [ $id_exe -eq 2 ]; then
  EXE="channel3d"
elif [ $id_exe -eq 3 ]; then
  EXE="channelLPT"
elif [ $id_exe -eq 4 ]; then
  EXE="channelDEM"
elif [ $id_exe -eq 5 ]; then
  EXE="channelACM"
elif [ $id_exe -eq 6 ]; then
  EXE="channel3d_4th"
else
  echo "  Sorry, EXE type cannot be recognized."
  echo "  Please check"
  exit 1
fi

if [ ! -e $EXE ]; then
  echo "  Sorry, "$EXE " EXE not exist."
  echo "  Please check"
  exit 2
fi

echo
read -p "  How many number of nodes do you need: " nnode
if [ $nnode -lt 1 ]; then
  echo "  Sorry, node number should be bigger than 0  "
  echo "  Please check"
  exit 3
elif [ $nnode -gt 10 ]; then
  echo "  Sorry, node number should be smaller than 5"
  echo "  Please check"
  exit 4
fi

echo
let mproc=nnode*28
echo
echo "  Note: maximum processor number =  $mproc"
read -p "  How many number of processors do you need: " nproc
if [ $nproc -ge 1 -a $nproc -le $mproc ]; then
  echo "  nproc = $nproc OK! "
else
  echo
  echo "  Sorry, processor number should be bigger than 0 and smaller than $mproc"
  echo "  Please check"
  exit 6
fi

echo
read -p "  please type a job name: " jobName

echo "#!/bin/bash"                       >> $jobForCP3dName
echo "  "                                >> $jobForCP3dName
echo "#SBATCH -J "$jobName               >> $jobForCP3dName
echo "#SBATCH -p cnall"                  >> $jobForCP3dName
echo "#SBATCH -N "$nnode                 >> $jobForCP3dName
echo "#SBATCH -o "$jobName"Out.txt"      >> $jobForCP3dName
echo "#SBATCH -e "$jobName"Err.txt"      >> $jobForCP3dName
echo "#SBATCH --no-requeue"              >> $jobForCP3dName
echo "#SBATCH --ntasks-per-node=28"      >> $jobForCP3dName
echo "##SBATCH -w c06b01n[08-15,17-18]"  >> $jobForCP3dName
echo "##SBATCH -x c01b01n01"             >> $jobForCP3dName
echo "  "                                >> $jobForCP3dName

#echo "module load compiles/intel/2019/u4/config"  >> $jobForCP3dName
echo "module load compiles/intel/2018/u1/config"  >> $jobForCP3dName
echo "mpiexec.hydra -n $nproc ./$EXE"  >> $jobForCP3dName
chmod a+x ./$jobForCP3dName

echo
read -p "  submit the job now?(1=Yes, others=No): " submitFlag
if [ $submitFlag -eq 1 ]; then
  echo "  $jobForCP3dName will be submitted immediately! "
  sbatch $jobForCP3dName
else
    echo "  Not submit $jobForCP3dName now"
    echo "  Bye !!!  "  
fi
