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

#01. Drifting-Kissing-Tumbling
cd ./DriftKissTumble/
chmod a+x run*.sh
./runDriftKissTumble.sh $cmpStr $appStr $p_row $p_col # 5 5
cd ../

#02. Drag coefficient of flow past isolated sphere
cd ./FlowPastSingleSphere/
chmod a+x run*.sh
./runFlowPastSphere.sh $cmpStr $appStr $p_row $p_col # 11 10
cd ../

#03. Rebounded angle of Oblique Collision
cd ./ObliqueCollision/
chmod a+x run*.sh
./runObliqueCollide.sh $cmpStr $appStr $p_row $p_col # 11 10
cd ../

#04. Trajectory of the sphere when normally colliding the wall
cd ./ParticleCollideWall/
chmod a+x run*.sh
./runPrtclCollideWall.sh $cmpStr $appStr $p_row $p_col # 11 10
cd ../

#05. Trajectory of the sphere when falling
cd ./ParticleFalling/
chmod a+x run*.sh
./runParticleFalling.sh $cmpStr $appStr $p_row $p_col  # 6 6
cd ../

#06. Particle dynamics in linear flow
cd ./ParticleInLinearFlow/
chmod a+x run*.sh
./runPrtclLinearTschisgaleFig5a.sh $cmpStr $appStr $p_row $p_col # 11 10
./runPrtclLinearTschisgaleFig5b.sh $cmpStr $appStr $p_row $p_col # 11 10
cd ../

#07. Normal restitution of particle in fluid
cd ./PrtclRestitution/
chmod a+x run*.sh
./runPrtclRestitution.sh $cmpStr $appStr $p_row $p_col # 11 10
cd ../

#08. Trajectory of two particles in Couttee shear flow at low Reynolds number
cd ./TwoParticleInShearFlow/
chmod a+x run*.sh
./runTwoPrtclInShearFlow.sh $cmpStr $appStr $p_row $p_col # 9 9
cd ../
