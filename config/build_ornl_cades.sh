#!/bin/bash

################################################################
## * This script builds available configurations of QMCPACK   ##
##   on CADES SHPC Condos , at Oak Ridge National Lab.        ##
##                                                            ##
## * Execute this script in trunk/                            ##
##   ./config/build_ornl_cades.sh                             ##
##                                                            ##
## Last verified: Jan 27, 2022                                ##
################################################################

# module files resulting from module imports below:
# Currently Loaded Modulefiles:
#   1) python/3.6.3           3) openmpi/3.1.5          5) gcc/7.2.0              7) fftw/3.3.5             9) boost/1.70.0
#   2) intel/19.0.3           4) PE-intel/3.0           6) hdf5-parallel/1.8.21   8) cmake/3.18.4          10) libxml2/2.9.9

source $MODULESHOME/init/bash
module purge
module load python
module load PE-intel/3.0
module swap intel intel/19.0.3
module load gcc/7.2.0
module load hdf5-parallel/1.8.21
module load fftw/3.3.5
module load cmake/3.18.4
module load boost/1.70.0
module load libxml2/2.9.9
module list

# Set environment variables
export FFTW_HOME=$FFTW_DIR/..
export BOOST_ROOT=$BOOST_DIR

CMAKE_FLAGS="-DCMAKE_C_COMPILER=mpicc \
             -DCMAKE_CXX_COMPILER=mpicxx \
             -DCMAKE_C_FLAGS=-xCOMMON-AVX512 \
             -DCMAKE_CXX_FLAGS=-xCOMMON-AVX512"

# Configure and build cpu real. Build targets skylake nodes.
echo ""
echo ""
echo "building QMCPACK for cpu real for CADES SHPC Condo -- Using AVX512 for Skylake nodes"
mkdir -p build_cades_cpu_real_skylake
cd build_cades_cpu_real_skylake
cmake $CMAKE_FLAGS ..
make -j 16
cd ..
ln -sf ./build_cades_cpu_real_skylake/bin/qmcpack ./qmcpack_cades_cpu_real_skylake

# Configure and build cpu complex. Build targets skylake nodes.
echo ""
echo ""
echo "building QMCPACK for cpu complex for CADES SHPC Condo -- Using AVX512 for Skylake nodes"
mkdir -p build_cades_cpu_comp_skylake
cd build_cades_cpu_comp_skylake
cmake -DQMC_COMPLEX=1 $CMAKE_FLAGS ..
make -j 16
cd ..
ln -sf ./build_cades_cpu_comp_skylake/bin/qmcpack_complex ./qmcpack_cades_cpu_comp_skylake

CMAKE_FLAGS="-DCMAKE_C_COMPILER=mpicc \
             -DCMAKE_CXX_COMPILER=mpicxx"

# Configure and build cpu real 
echo ""
echo ""
echo "building QMCPACK for cpu real for CADES SHPC Condo"
mkdir -p build_cades_cpu_real
cd build_cades_cpu_real
cmake $CMAKE_FLAGS ..
make -j 16
cd ..
ln -sf ./build_cades_cpu_real/bin/qmcpack ./qmcpack_cades_cpu_real

# Configure and build cpu complex
echo ""
echo ""
echo "building QMCPACK for cpu complex for CADES SHPC Condo"
mkdir -p build_cades_cpu_comp
cd build_cades_cpu_comp
cmake -DQMC_COMPLEX=1 $CMAKE_FLAGS ..
make -j 16
cd ..
ln -sf ./build_cades_cpu_comp/bin/qmcpack_complex ./qmcpack_cades_cpu_comp
