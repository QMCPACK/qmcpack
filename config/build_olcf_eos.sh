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
## Last modified: Sep 25, 2017                                ##
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
if (echo $LOADEDMODULES | grep -q hdf5)
then
module unload cray-hdf5
fi
if (echo $LOADEDMODULES | grep -q libsci)
then
module unload cray-libsci
fi
module load PrgEnv-intel
module unload gcc
module load gcc
module load cray-hdf5-parallel
module load fftw
module load boost
module load subversion
module load cmake

# always dynamic linking
export CRAYPE_LINK_TYPE=dynamic

# Set environment variables
export FFTW_HOME=$FFTW_DIR/..


# Set cmake variables, shared for cpu builds
CMAKE_FLAGS="-DCMAKE_C_COMPILER=cc \ 
             -DCMAKE_CXX_COMPILER=CC"


# Configure and build cpu real AoS
echo ""
echo ""
echo "building qmcpack for cpu AoS real for eos"
mkdir -p build_cpu_AoS_real_eos
cd build_cpu_AoS_real_eos
cmake $CMAKE_FLAGS .. 
make -j 32
cd ..
ln -sf ./build_cpu_AoS_real_eos/bin/qmcpack ./qmcpack_cpu_AoS_real_eos


# Configure and build cpu complex AoS
echo ""
echo ""
echo "building qmcpack for cpu AoS complex for eos"
mkdir -p build_cpu_AoS_comp_eos
cd build_cpu_AoS_comp_eos
cmake -DQMC_COMPLEX=1 $CMAKE_FLAGS .. 
make -j 32
cd ..
ln -sf ./build_cpu_AoS_comp_eos/bin/qmcpack ./qmcpack_cpu_AoS_comp_eos


# Configure and build cpu real SoA
echo ""
echo ""
echo "building qmcpack for cpu SoA real for eos"
mkdir -p build_cpu_SoA_real_eos
cd build_cpu_SoA_real_eos
cmake -DENABLE_SOA=1 $CMAKE_FLAGS .. 
make -j 32
cd ..
ln -sf ./build_cpu_SoA_real_eos/bin/qmcpack ./qmcpack_cpu_SoA_real_eos

# Configure and build cpu complex SoA
echo ""
echo ""
echo "building qmcpack for cpu SoA complex for eos"
mkdir -p build_cpu_SoA_comp_eos
cd build_cpu_SoA_comp_eos
cmake -DQMC_COMPLEX=1 -DENABLE_SOA=1 $CMAKE_FLAGS .. 
make -j 32
cd ..
ln -sf ./build_cpu_SoA_comp_eos/bin/qmcpack ./qmcpack_cpu_SoA_comp_eos







