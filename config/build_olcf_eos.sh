#!/bin/bash

################################################################
## * This script builds available configurations of QMCPACK   ##
##   on EOS, at Oak Ridge Leadership Computing Facility,      ##
##   Oak Ridge National Lab.                                  ##
##                                                            ##
## * Execute this script in trunk/                            ##
##   ./config/build_olcf_eos.sh                               ##
##                                                            ##
## Settings should be consistent with test scripts e.g.       ##
##  ./tests/test_automation/nightly_ornl_olcf_eos.job         ##
##                                                            ##
## Last modified: Jan 6, 2017                                 ##
################################################################


# Load required modules (assuming default settings have not been modified)
source $MODULESHOME/init/bash
if (echo $LOADEDMODULES | grep -q pgi)
then
module unload PrgEnv-pgi
fi
if (echo $LOADEDMODULES | grep -q gnu)
then
module unload PrgEnv-gnu
fi
module load PrgEnv-intel
module load gcc
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
mkdir build_eos_cpu_real
cd build_eos_cpu_real
cmake $CMAKE_FLAGS .. 
make -j 32
cd ..
ln -s ./build_eos_cpu_real/bin/qmcpack ./qmcpack_eos_cpu_real


# Configure and build cpu complex
echo ""
echo ""
echo "building qmcpack for cpu complex"
mkdir build_eos_cpu_comp
cd build_eos_cpu_comp
cmake -DQMC_COMPLEX=1 $CMAKE_FLAGS .. 
make -j 32
cd ..
ln -s ./build_eos_cpu_comp/bin/qmcpack ./qmcpack_eos_cpu_comp




