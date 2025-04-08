#!/bin/bash

#####################################################################################
## * This script builds available configurations of QMCPACK                        ##
##   on BESMS at ORNL                                                              ##
##                                                                                 ##
## * Execute this script in qmcpack git top level directory                        ##
##   ./config/build_ornl_cades_besms.sh                                            ##
##                                                                                 ##
## Last verified: April 4, 2025                                                    ##
#####################################################################################

# module files resulting from module imports below:
# Currently Loaded Modules:
#   1) DefApps/default   2) intel/2024.1.0   3) openmpi/4.1.6   4) boost/1.85.0   5) fftw/3.3.10-mpi-omp   
#   6) hdf5/1.14.3-mpi   7) cmake/3.24.4   8) libxml2/2.10.3
source $MODULESHOME/init/bash
source /sw/cades-besms/oneapi/2024.1.0/setvars.sh 
module purge
module load DefApps/default
module load intel/2024.1.0
module load openmpi/4.1.6
module load boost/1.85.0
module load fftw/3.3.10-mpi-omp
# module load openblas/0.3.23
module load hdf5/1.14.3-mpi
module load cmake/3.24.4
module load libxml2/2.10.3
# MKLROOT is pointing to a wrong path, for now correct by hand. 
export MKLROOT=/sw/cades-besms/oneapi/2024.1.0/mkl/2024.1/


module list

CMAKE_FLAGS="-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpiCC \
             -DENABLE_TIMERS=1 -DBUILD_AFQMC=0 "             
# If QMC_DATA is available, consider adding e.g.
#             -DQMC_DATA=/gpfs/wolf2/cades/mat269/world-shared/pk7/QMC_DATA_FULL -DQMC_PERFORMANCE_NIO_MAX_ATOMS=128 -DQMC_PERFORMANCE_C_MOLECULE_MAX_ATOMS=16 -DQMC_PERFORMANCE_C_GRAPHITE_MAX_ATOMS=16

# Configure and build cpu real.
echo ""
echo ""
echo "building QMCPACK for cpu real for BESMS"
mkdir -p build_besms_cpu_real
cd build_besms_cpu_real
cmake $CMAKE_FLAGS ..
nice make -j 32
cd ..
ln -sf ./build_besms_cpu_real/bin/qmcpack ./qmcpack_besms_cpu_real

# Configure and build cpu complex.
echo ""
echo ""
echo "building QMCPACK for cpu complex for BESMS"
mkdir -p build_besms_cpu_complex
cd build_besms_cpu_complex
cmake -DQMC_COMPLEX=1 $CMAKE_FLAGS ..
nice make -j 32
cd ..
ln -sf ./build_besms_cpu_complex/bin/qmcpack_complex ./qmcpack_besms_cpu_complex

