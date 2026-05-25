#!/bin/bash

nrow=4
ncol=2
RunName="Cha180_2nd"
RestartDir="./CFD/Restart"
ResultsDir="./CFD/Results"
BaseFile="./Input/CFD_2nd/TurbCha0180_2nd.standard"
ExeName="./channel2nd"

echo " "
echo "  Creating the restart input begins !"
echo " "
chmod a+x ./CreateRestartInput.sh
NewFile=`./CreateRestartInput.sh $RunName $RestartDir $ResultsDir $BaseFile`
echo "  >>New restart input file is: "
echo "  "$NewFile
echo " "
echo "  Creating the restart input ends !"
echo " "
let nproc=nrow*ncol
mpirun -n $nproc $ExeName $NewFile $nrow $ncol

