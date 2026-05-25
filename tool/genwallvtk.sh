#!/bin/bash
GeomToVTKSrc="./GeomToVTK.f90"
CMP="gfortran" #ifort

$CMP $GeomToVTKSrc  -o wallvtk 
chmod 777 ./wallvtk
./wallvtk  
#rm -rf wallvtk
