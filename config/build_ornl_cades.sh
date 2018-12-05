#!/bin/bash

################################################################
## * This script builds available configurations of QMCPACK   ##
##   on CADES SHPC Condos , at Oak Ridge National Lab.        ##
##                                                            ##
## * Execute this script in trunk/                            ##
##   ./config/build_ornl_cades.sh                             ##
##                                                            ##
## Last verified: Dec 5, 2018                                 ##
################################################################

source $MODULESHOME/init/bash
module purge
module load PE-intel
module load xalt/0.7.6
module load hdf5-parallel/1.8.17
module load boost/1.64.0
module load cmake/3.12.0
module swap intel intel/18.0.0
module swap openmpi openmpi/1.10.3
module load gcc/6.3.0
module list

# Set environment variables
export FFTW_HOME=$FFTW_DIR/..

CMAKE_FLAGS="-DCMAKE_C_COMPILER=mpicc \
             -DCMAKE_CXX_COMPILER=mpicxx \
             -DCMAKE_C_FLAGS=-xCOMMON-AVX512 \
             -DCMAKE_CXX_FLAGS=-xCOMMON-AVX512"

# Configure and build cpu real AoS
echo ""
echo ""
echo "building qmcpack for cpu AoS real for CADES SHPC Condo"
mkdir -p build_cades_cpu_real
cd build_cades_cpu_real
cmake $CMAKE_FLAGS ..
make -j 16
cd ..
ln -sf ./build_cades_cpu_real/bin/qmcpack ./qmcpack_cades_cpu_real

# Configure and build cpu complex AoS
echo ""
echo ""
echo "building qmcpack for cpu AoS complex for CADES SHPC Condo"
mkdir -p build_cades_cpu_comp
cd build_cades_cpu_comp
cmake -DQMC_COMPLEX=1 $CMAKE_FLAGS ..
make -j 16
cd ..
ln -sf ./build_cades_cpu_comp/bin/qmcpack ./qmcpack_cades_cpu_comp

# Configure and build cpu real SoA
echo ""
echo ""
echo "building qmcpack for cpu SoA real for CADES SHPC Condo"
mkdir -p build_cades_cpu_real_SoA
cd build_cades_cpu_real_SoA
cmake -DENABLE_SOA=1 $CMAKE_FLAGS ..
make -j 16
cd ..
ln -sf ./build_cades_cpu_real_SoA/bin/qmcpack ./qmcpack_cades_cpu_real_SoA

# Configure and build cpu complex SoA
echo ""
echo ""
echo "building qmcpack for cpu SoA complex for CADES SHPC Condo"
mkdir -p build_cades_cpu_comp_SoA
cd build_cades_cpu_comp_SoA
cmake -DQMC_COMPLEX=1 -DENABLE_SOA=1 $CMAKE_FLAGS ..
make -j 16
cd ..
ln -sf ./build_cades_cpu_comp_SoA/bin/qmcpack ./qmcpack_cades_cpu_comp_SoA
