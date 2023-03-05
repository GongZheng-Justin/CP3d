#!/bin/bash
#=======================================================================
# install_thirdParty.sh.sh example --Zheng Gong, 2021-08-27(yy/mm/dd)
#=======================================================================

blas_flag=0
fftw_flag=1
hdf5_flag=0
hypre_flag=0
lapack_flag=0
blasVersion="blas-3.10.0"
fftwVersion="fftw-3.3.9"
hdf5Version="hdf5-1.12.1"
hypreVersion="hypre-2.22.0"
lapackVersion="lapack-3.10.0"

#======================================================================#
# Normally no need to change anything below

#-------------------------------- blas --------------------------------#
if [ $blas_flag -eq 1 ]; then
echo
echo "Compiling "$blasVersion" begins (Wait several minutes).....!"
blasDir="blas"
rm -rf $blasDir $blasVersion
mkdir $blasDir
tar -xvf $blasVersion".tar.gz" >&/dev/null
cd $blasVersion
make >&/dev/null
cp -rf blas*.a ../$blasDir/libblas.a
cd ..
rm -rf $blasVersion
echo "Compiling "$blasVersion" is done !"
fi

#-------------------------------- fftw --------------------------------#
if [ $fftw_flag -eq 1 ]; then
echo
echo "Compiling "$fftwVersion" begins (Wait several minutes).....!"
fftwDir="fftw"
rm -rf $fftwDir $fftwVersion
mkdir $fftwDir
tar -xvf $fftwVersion".tar.gz" >&/dev/null
cd $fftwVersion
mkdir $fftwDir
./configure --with-our-malloc16 --enable-threads --with-combined-threads --enable-sse2 --prefix=$(pwd)/$fftwDir >&/dev/null
make install >&/dev/null
cp -rf $fftwDir/include/* $fftwDir/lib/lib*.a ../$fftwDir
cd ..
rm -rf $fftwVersion
echo "Compiling "$fftwVersion" is done !"
fi

#-------------------------------- hdf5 --------------------------------#
if [ $hdf5_flag -eq 1 ]; then
echo
echo "Compiling "$hdf5Version" begins (Wait several minutes).....!"
hdf5Dir="hdf5"
rm -rf $hdf5Dir $hdf5Version
mkdir  $hdf5Dir
tar -xvf $hdf5Version".tar.gz" >&/dev/null
cd $hdf5Version
mkdir $hdf5Dir
./configure --enable-fortran --enable-parallel --prefix=$(pwd)/$hdf5Dir >&/dev/null
make install >&/dev/null
cp -rf $hdf5Dir/include/* $hdf5Dir/lib/* ../$hdf5Dir
cd .. 
rm -rf $hdf5Version
echo "Compiling "$hdf5Version" is done !"
fi

#-------------------------------- hypre -------------------------------#
if [ $hypre_flag -eq 1 ]; then
echo
echo "Compiling "$hypreVersion" begins (Wait several minutes).....!"
hypreDir="hypre"
rm -rf $hypreDir $hypreVersion
mkdir  $hypreDir
tar -xvf $hypreVersion".tar.gz" >&/dev/null
cd $hypreVersion/src
mkdir $hypreDir
./configure --prefix=$(pwd)/$hypreDir >&/dev/null
make install >&/dev/null
cp -rf $hypreDir/include/* $hypreDir/lib/* ../../$hypreDir
cd ../../
rm -rf $hypreVersion
echo "Compiling "$hypreVersion" is done !"
fi

#------------------------------- lapack -------------------------------#
if [ $lapack_flag -eq 1 ]; then
echo
echo "Compiling "$lapackVersion" begins (Wait several minutes).....!"
lapackDir="lapack"
rm -rf $lapackDir $lapackVersion
mkdir $lapackDir
tar -xvf $lapackVersion".tar.gz" >&/dev/null
cd $lapackVersion
cp make.inc.example make.inc
make >&/dev/null
cp -rf lib*.a ../$lapackDir
cd ..
rm -rf $lapackVersion
echo "Compiling "$lapackVersion" is done !"
echo
fi
