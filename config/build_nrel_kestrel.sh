#!/bin/bash

#####################################################################################
## * This script builds available configurations of QMCPACK                        ##
##   on kestrel at NREL                                                            ##
##                                                                                 ##
## * Execute this script in qmcpack git top level directory                        ##
##   ./config/build_nrel_kestrel.sh                                                ##
##                                                                                 ##
## Last verified: July 23, 2024                                                    ##
#####################################################################################

# module files resulting from module imports below:
# Currently Loaded Modules:
#  1) intel-oneapi-compilers/2023.2.0      2) intel-oneapi-mpi/2021.11.0-intel
#  3) gcc/13.1.0                           4) intel-oneapi-tbb/2021.10.0-intel
#  5) intel-oneapi-mkl/2023.2.0-intel      6) boost/1.84.0-intel-oneapi-mpi-intel
#  7) fftw/3.3.10-intel-oneapi-mpi-intel   8) hdf5/1.14.3-intel-oneapi-mpi-intel
#  9) curl/8.6.0                          10) cmake/3.27.9

source $MODULESHOME/init/bash
module purge
module load intel-oneapi-compilers/2023.2.0
module load intel-oneapi-mpi/2021.11.0-intel
module load gcc/13.1.0
module load intel-oneapi-mkl/2023.2.0-intel
module load boost/1.84.0-intel-oneapi-mpi-intel
module load fftw/3.3.10-intel-oneapi-mpi-intel
module load hdf5/1.14.3-intel-oneapi-mpi-intel
module load cmake/3.27.9
 
module list

CMAKE_FLAGS="-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx"

# Configure and build cpu real.
echo ""
echo ""
echo "building QMCPACK for cpu real for Kestrel"
mkdir -p build_cpu_real
cd build_cpu_real
cmake $CMAKE_FLAGS ..
make -j $(nproc)
cd ..
ln -sf ./build_cpu_real/bin/qmcpack ./qmcpack_cpu_real

# Configure and build cpu complex.
echo ""
echo ""
echo "building QMCPACK for cpu complex for Kestrel"
mkdir -p build_cpu_complex
cd build_cpu_complex
cmake -DQMC_COMPLEX=1 $CMAKE_FLAGS ..
make -j $(nproc)
cd ..
ln -sf ./build_cpu_complex/bin/qmcpack_complex ./qmcpack_cpu_complex

