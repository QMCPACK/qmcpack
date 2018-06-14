#!/bin/bash

################################################################
## * This script builds available configurations of QMCPACK   ##
##   on Titan, Oak Ridge National Lab.                        ##
##                                                            ##
## * Execute this script in trunk/                            ##
##   ./config/build_olcf_titan.sh                             ##
##                                                            ##
## (c) Scientific Computing Group, NCCS, ORNL.                ##
## Last modified: Dec 7, 2017                                 ##
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



# Configure and build cpu real AoS
echo ""
echo ""
echo "building qmcpack for cpu AoS real for titan"
mkdir -p build_titan_cpu_real
cd build_titan_cpu_real
cmake -D CMAKE_C_FLAGS="$XT_FLAGS" \
      -D CMAKE_CXX_FLAGS="$XT_FLAGS" \
      $CMAKE_FLAGS ..
make -j 32
cd ..
ln -sf ./build_titan_cpu_real/bin/qmcpack ./qmcpack_titan_cpu_real


# Configure and build cpu complex AoS
echo ""
echo ""
echo "building qmcpack for cpu AoS complex for titan"
mkdir -p build_titan_cpu_comp
cd build_titan_cpu_comp
cmake -D CMAKE_C_FLAGS="$XT_FLAGS" \
      -D CMAKE_CXX_FLAGS="$XT_FLAGS" \
      -D QMC_COMPLEX=1 $CMAKE_FLAGS ..
make -j 32
cd ..
ln -sf ./build_titan_cpu_comp/bin/qmcpack ./qmcpack_titan_cpu_comp

# Configure and build cpu real SoA
echo ""
echo ""
echo "building qmcpack for cpu SoA real for titan"
mkdir -p build_titan_cpu_real_SoA
cd build_titan_cpu_real_SoA
cmake -D CMAKE_C_FLAGS="$XT_FLAGS" \
      -D CMAKE_CXX_FLAGS="$XT_FLAGS" \
      -D ENABLE_SOA=1 \
      $CMAKE_FLAGS ..
make -j 32
cd ..
ln -sf ./build_titan_cpu_real_SoA/bin/qmcpack ./qmcpack_titan_cpu_real_SoA


# Configure and build cpu complex SoA
echo ""
echo ""
echo "building qmcpack for cpu SoA complex for titan"
mkdir -p build_titan_cpu_comp_SoA
cd build_titan_cpu_comp_SoA
cmake -D CMAKE_C_FLAGS="$XT_FLAGS" \
      -D CMAKE_CXX_FLAGS="$XT_FLAGS" \
      -D ENABLE_SOA=1 \
      -D QMC_COMPLEX=1 $CMAKE_FLAGS ..
make -j 32
cd ..
ln -sf ./build_titan_cpu_comp_SoA/bin/qmcpack ./qmcpack_titan_cpu_comp_SoA


# load cuda toolkit for gpu build
module load cudatoolkit

# Configure and build gpu real for titan
echo ""
echo ""
echo "building qmcpack for gpu real for titan"
mkdir -p build_titan_gpu_real
cd build_titan_gpu_real
cmake -D QMC_CUDA=1 ..
cmake -D QMC_CUDA=1 ..
make -j 32
cd ..
ln -sf ./build_titan_gpu_real/bin/qmcpack ./qmcpack_titan_gpu_real

# Configure and build gpu complex
echo ""
echo ""
echo "building qmcpack for gpu complex for titan"
mkdir -p build_titan_gpu_comp
cd build_titan_gpu_comp
cmake -D QMC_CUDA=1 -D QMC_COMPLEX=1 ..
cmake -D QMC_CUDA=1 -D QMC_COMPLEX=1 ..
make -j 32
cd ..
ln -sf ./build_titan_gpu_comp/bin/qmcpack ./qmcpack_titan_gpu_comp
