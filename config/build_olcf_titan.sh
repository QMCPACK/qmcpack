#!/bin/bash

################################################################
## * This script builds available configurations of QMCPACK   ##
##   on Titan, Oak Ridge National Lab.                        ##
##                                                            ##
## * Execute this script in trunk/                            ##
##   ./config/build_titan.sh                                  ##
##                                                            ##
## (c) Scientific Computing Group, NCCS, ORNL.                ##
## Last modified: Sep 25, 2017                                ##
################################################################


# Load required modules (assuming default settings have not been modified)
source $MODULESHOME/init/bash
if (echo $LOADEDMODULES | grep -q pgi)
then
module unload PrgEnv-pgi
fi
if (echo $LOADEDMODULES | grep -q cuda)
then
module unload cudatoolkit
fi
if (echo $LOADEDMODULES | grep -q hdf5)
then
module unload cray-hdf5
fi
module load PrgEnv-gnu
module load cray-hdf5-parallel
module load fftw
module load boost
module load subversion
module load cmake3/3.6.1

# always dynamic linking
export CRAYPE_LINK_TYPE=dynamic

# Set environment variables
export FFTW_HOME=$FFTW_DIR/..

export CC=cc
export CXX=CC
XT_FLAGS="-DHAVE_AMDLIBM=1"

# Set cmake variables, shared for cpu builds
CMAKE_FLAGS="-D QMC_INCLUDE=/sw/xk7/amdlibm/include \
             -D QMC_EXTRA_LIBS=/sw/xk7/amdlibm/lib/static/libamdlibm.a"

# Configure and build cpu real
echo ""
echo ""
echo "building qmcpack for cpu real"
mkdir build_cpu_real
cd build_cpu_real
cmake -D CMAKE_C_FLAGS="$XT_FLAGS" \
      -D CMAKE_CXX_FLAGS="$XT_FLAGS" \
      $CMAKE_FLAGS ..
make -j 32
cd ..
ln -s ./build_cpu_real/bin/qmcpack ./qmcpack_cpu_real


# Configure and build cpu complex
echo ""
echo ""
echo "building qmcpack for cpu complex"
mkdir build_cpu_comp
cd build_cpu_comp
cmake -D CMAKE_C_FLAGS="$XT_FLAGS" \
      -D CMAKE_CXX_FLAGS="$XT_FLAGS" \
      -D QMC_COMPLEX=1 $CMAKE_FLAGS ..
make -j 32
cd ..
ln -s ./build_cpu_comp/bin/qmcpack ./qmcpack_cpu_comp

# load cuda toolkit for gpu build
module load cudatoolkit

# Configure and build gpu real
echo ""
echo ""
echo "building qmcpack for gpu real"
mkdir build_gpu_real
cd build_gpu_real
cmake -D QMC_CUDA=1 ..
cmake -D QMC_CUDA=1 ..
make -j 32
cd ..
ln -s ./build_gpu_real/bin/qmcpack ./qmcpack_gpu_real

# Configure and build gpu complex
echo ""
echo ""
echo "building qmcpack for gpu complex"
mkdir build_gpu_comp
cd build_gpu_comp
cmake -D QMC_CUDA=1 -D QMC_COMPLEX=1 ..
cmake -D QMC_CUDA=1 -D QMC_COMPLEX=1 ..
make -j 32
cd ..
ln -s ./build_gpu_comp/bin/qmcpack ./qmcpack_gpu_comp
