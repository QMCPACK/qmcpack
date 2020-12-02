#!/bin/bash

################################################################
## * This script builds available configurations of QMCPACK   ##
##   on CADES SHPC Condos , at Oak Ridge National Lab.        ##
##                                                            ##
## * Execute this script in trunk/                            ##
##   ./config/build_ornl_cades.sh                             ##
##                                                            ##
## Last verified: Nov 12, 2020                                ##
################################################################

# module files resulting from module imports below:
# Currently Loaded Modulefiles:
#   1) python/3.6.3           3) openmpi/3.1.5          5) gcc/6.3.0              7) hdf5_parallel/1.10.3   9) cmake/3.12.0          11) libxml2/2.9.9
#   2) intel/18.0.0           4) PE-intel/3.0           6) intel/19.0.3           8) fftw/3.3.5            10) boost/1.67.0

source $MODULESHOME/init/bash
module purge
module load python
module load PE-intel/3.0
module load gcc/6.3.0
module load intel/19.0.3
module load hdf5_parallel/1.10.3
module load fftw/3.3.5
module load cmake/3.12.0
module load boost/1.67.0
module load libxml2/2.9.9
module list

# Set environment variables
export FFTW_HOME=$FFTW_DIR/..
export BOOST_ROOT=$BOOST_DIR

CMAKE_FLAGS="-DCMAKE_C_COMPILER=mpicc \
             -DCMAKE_CXX_COMPILER=mpicxx \
             -DCMAKE_C_FLAGS=-xCOMMON-AVX512 \
             -DCMAKE_CXX_FLAGS=-xCOMMON-AVX512"

# Configure and build cpu real SoA. Build targets skylake nodes.
echo ""
echo ""
echo "building QMCPACK for cpu SoA real for CADES SHPC Condo -- Using AVX512 for Skylake nodes"
mkdir -p build_cades_cpu_real_skylake
cd build_cades_cpu_real_skylake
cmake -DENABLE_SOA=1 $CMAKE_FLAGS ..
make -j 16
cd ..
ln -sf ./build_cades_cpu_real_skylake/bin/qmcpack ./qmcpack_cades_cpu_real_skylake

# Configure and build cpu complex SoA. Build targets skylake nodes.
echo ""
echo ""
echo "building QMCPACK for cpu SoA complex for CADES SHPC Condo -- Using AVX512 for Skylake nodes"
mkdir -p build_cades_cpu_comp_skylake
cd build_cades_cpu_comp_skylake
cmake -DENABLE_SOA=1 -DQMC_COMPLEX=1 $CMAKE_FLAGS ..
make -j 16
cd ..
ln -sf ./build_cades_cpu_comp_skylake/bin/qmcpack ./qmcpack_cades_cpu_comp_skylake

CMAKE_FLAGS="-DCMAKE_C_COMPILER=mpicc \
             -DCMAKE_CXX_COMPILER=mpicxx"

# Configure and build cpu real SoA 
echo ""
echo ""
echo "building QMCPACK for cpu SoA real for CADES SHPC Condo"
mkdir -p build_cades_cpu_real
cd build_cades_cpu_real
cmake -DENABLE_SOA=1 $CMAKE_FLAGS ..
make -j 16
cd ..
ln -sf ./build_cades_cpu_real/bin/qmcpack ./qmcpack_cades_cpu_real

# Configure and build cpu complex SoA
echo ""
echo ""
echo "building QMCPACK for cpu SoA complex for CADES SHPC Condo"
mkdir -p build_cades_cpu_comp
cd build_cades_cpu_comp
cmake -DQMC_COMPLEX=1 -DENABLE_SOA=1 $CMAKE_FLAGS ..
make -j 16
cd ..
ln -sf ./build_cades_cpu_comp/bin/qmcpack ./qmcpack_cades_cpu_comp

