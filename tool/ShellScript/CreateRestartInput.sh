#!/bin/bash
RunName="Cha180_2ndttt"
RestartDir="./CFD/Restart"
ResultsDir="./CFD/Results"
BaseFile='TurbCha0180_2nd.standard'
if [[ -n $1 ]]; then
  RunName=$1
fi
if [[ -n $2 ]]; then
  RestartDir=$2
fi
if [[ -n $3 ]]; then
  ResultsDir=$3
fi
if [[ -n $4 ]]; then
  BaseFile=$4
fi

# Put the file name in array "file_list" 
index=-1
for file_a in ${RestartDir}/*
do
  let index=index+1
  file_list[index]=`basename $file_a`
done

# Find the maximum nRestart
nRestart=0
for index in ${!file_list[@]};
do
  filename=${file_list[index]}
  if [[ $filename =~ "RestartFor"$RunName ]] ;then
    filename=${filename:0-10}
    iRestart=`echo $filename|awk '{print int($0)}'`
    if [ "$iRestart" -gt "$nRestart" ];then
      nRestart=$iRestart
    fi
  fi
done

# Create the new input file
NewFile=$BaseFile"_"$nRestart
rm -rf $NewFile
touch $NewFile
IFS=$(echo -en "\n") # save the blank space
while read -r StrRead
do 
  StrLine=$StrRead
  if [[ $StrRead =~ "RestartFlag" ]] ;then
    StrLine="  RestartFlag = T ! Restart or not"  
  fi
  if [[ $StrRead =~ "ifirst" ]] ;then
    let ifirst=nRestart+1
    StrLine="  ifirst= "$ifirst" ! Restart or not"  
  fi
  echo $StrLine >> $NewFile
done < $BaseFile
IFS=$SAVEIFS

echo $NewFile

# Copy the log file
NewLog=$ResultsDir/$RunName.log$nRestart
if [ ! -f $NewLog ]; then
  cp $ResultsDir/$RunName.log $NewLog
fi
