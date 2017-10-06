#!/bin/bash

################################################################
## * This script builds available configurations of QMCPACK   ##
##   on OIC5, at Oak Ridge National Lab.                      ##
##                                                            ##
## * Execute this script in trunk/                            ##
##   ./config/build_ornl_oic.sh                               ##
##                                                            ##
## Last modified: Oct 6, 2017                                 ##
################################################################


module () 
{ 
    eval `/opt/modules/3.1.6/bin/modulecmd bash $*`
}
module purge
module load gcc/4.9.3
module load mpi/openmpi-1.4.5-gcc4
module load hdf5/1.8.8-gcc4-parallel

export LIBXML2_HOME=/home/j1k/share/oic5_gcc4/libxml2-2.7.6/build

CMAKE_FLAGS="-DCMAKE_C_COMPILER=mpicc \ 
             -DCMAKE_CXX_COMPILER=mpicxx"

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


