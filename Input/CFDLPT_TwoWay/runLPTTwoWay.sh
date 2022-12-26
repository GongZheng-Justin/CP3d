#!/bin/bash

cd ../../
mpirun -n 8 ./channelLPT ./Input/LPTTwoWay/Channel4th_LPT.Lee4 ./Input/LPTTwoWay/LPT_Channel4th.Lee4
mpirun -n 8 ./channelLPT ./Input/LPTTwoWay/Channel4th_LPT.Lee3 ./Input/LPTTwoWay/LPT_Channel4th.Lee3
mpirun -n 8 ./channelLPT ./Input/LPTTwoWay/Channel4th_LPT.Lee2 ./Input/LPTTwoWay/LPT_Channel4th.Lee2
mpirun -n 8 ./channelLPT ./Input/LPTTwoWay/Channel4th_LPT.Lee1 ./Input/LPTTwoWay/LPT_Channel4th.Lee1
