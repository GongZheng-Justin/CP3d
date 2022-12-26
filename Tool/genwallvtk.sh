#!/bin/bash
GeomToVTKSrc="./GeomToVTK.f90"
CMP="gfortran" #ifort

$CMP $GeomToVTKSrc  -o wallvtk 
chmod a+x ./wallvtk
./wallvtk  
rm -rf wallvtk
