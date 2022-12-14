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

#01. Drag coefficient of flow past isolated sphere
cd ./FlowPastSingleSphere/
chmod a+x run*.sh
./runFlowPastSphere.sh $cmpStr $appStr $p_row $p_col
cd ../

#02. Rebounded angle of Oblique Collision
cd ./ObliqueCollision/
chmod a+x run*.sh
./runObliqueCollide.sh $cmpStr $appStr $p_row $p_col
cd ../

#03. Trajectory of the sphere when normally colliding the wall
cd ./ParticleCollideWall/
chmod a+x run*.sh
./runPrtclCollideWall.sh $cmpStr $appStr $p_row $p_col
cd ../

#04. Trajectory of the sphere when falling
cd ./ParticleFalling/
chmod a+x run*.sh
./runParticleFalling.sh $cmpStr $appStr $p_row $p_col
cd ../

#05. Particle dynamics in linear flow
cd ./ParticleInLinearFlow/
chmod a+x run*.sh
./runPrtclLinearTschisgaleFig5a.sh $cmpStr $appStr $p_row $p_col
./runPrtclLinearTschisgaleFig5b.sh $cmpStr $appStr $p_row $p_col
cd ../

#06. Normal restitution of particle in fluid
cd ./PrtclRestitution/
chmod a+x run*.sh
./runPrtclRestitution.sh $cmpStr $appStr $p_row $p_col
cd ../

#07. Trajectory of two particles in Couttee shear flow at low Reynolds number
cd ./TwoPartcileInShearFlow/
chmod a+x run*.sh
./runTwoPrtclInShearFlow.sh $cmpStr $appStr $p_row $p_col
cd ../
