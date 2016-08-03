#!/bin/bash

################################################################
## * This script builds available configurations of QMCPACK   ##
##   on EOS, Oak Ridge National Lab.                          ##
##                                                            ##
## * Execute this script in trunk/                            ##
##   ./config/build_eos.sh                                    ##
##                                                            ##
## Last modified: Aug 3, 2016                                 ##
################################################################


# Load required modules (assuming default settings have not been modified)
source $MODULESHOME/init/bash
if (echo $LOADEDMODULES | grep -q pgi)
then
module unload PrgEnv-pgi
fi
module load PrgEnv-intel
module load cray-hdf5
module load fftw
module load boost
module load subversion
module load cmake


# Set environment variables
export FFTW_HOME=$FFTW_DIR/..


# Set cmake variables, shared for cpu builds
CMAKE_FLAGS="-DCMAKE_C_COMPILER=cc \ 
             -DCMAKE_CXX_COMPILER=CC"


# Configure and build cpu real
echo ""
echo ""
echo "building qmcpack for cpu real"
mkdir build_cpu_real
cd build_cpu_real
cmake $CMAKE_FLAGS .. 
make -j 32
cd ..
ln -s ./build_cpu_real/bin/qmcpack ./qmcpack_cpu_real


# Configure and build cpu complex
echo ""
echo ""
echo "building qmcpack for cpu complex"
mkdir build_cpu_comp
cd build_cpu_comp
cmake -DQMC_COMPLEX=1 $CMAKE_FLAGS .. 
make -j 32
cd ..
ln -s ./build_cpu_comp/bin/qmcpack ./qmcpack_cpu_comp




