#!/bin/bash

#####################################################################################
## * This script builds available configurations of QMCPACK                        ##
##   on baseline at ORNL                                                           ##
##   https://docs.cades.olcf.ornl.gov/baseline_user_guide/baseline_user_guide.html ##
##                                                                                 ##
## * Execute this script in qmcpack git top level directory                        ##
##   ./config/build_ornl_cades_baseline.sh                                         ##
##                                                                                 ##
## Last verified: May 23, 2024                                                     ##
#####################################################################################

# module files resulting from module imports below:
# Currently Loaded Modules:
#  1) DefApps   2) gcc/12.2.0   3) openmpi/4.0.4   4) boost/1.83.0   5) fftw/3.3.10
#  6) openblas/0.3.23   7) hdf5/1.14.3   8) cmake/3.26.3

source $MODULESHOME/init/bash
module purge
module load DefApps
module load gcc/12.2.0
module load openmpi/4.0.4
module load boost/1.83.0
module load fftw/3.3.10
module load openblas/0.3.23
module load hdf5/1.14.3
module load cmake/3.26.3
 
module list

CMAKE_FLAGS="-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpiCC \
             -DENABLE_TIMERS=1 -DBUILD_AFQMC=0 \
             -DMPIEXEC_EXECUTABLE=/usr/bin/srun -DMPIEXEC_NUMPROC_FLAG='-n' -DMPIEXEC_PREFLAGS='-c;16;--distribution=block:cyclic'"
# If QMC_DATA is available, consider adding e.g.
#             -DQMC_DATA=/gpfs/wolf2/cades/mat269/world-shared/pk7/QMC_DATA_FULL -DQMC_PERFORMANCE_NIO_MAX_ATOMS=128 -DQMC_PERFORMANCE_C_MOLECULE_MAX_ATOMS=16 -DQMC_PERFORMANCE_C_GRAPHITE_MAX_ATOMS=16

# Configure and build cpu real.
echo ""
echo ""
echo "building QMCPACK for cpu real for Baseline"
mkdir -p build_baseline_cpu_real
cd build_baseline_cpu_real
cmake $CMAKE_FLAGS ..
nice make -j 32
cd ..
ln -sf ./build_baseline_cpu_real/bin/qmcpack ./qmcpack_baseline_cpu_real

# Configure and build cpu complex.
echo ""
echo ""
echo "building QMCPACK for cpu complex for Baseline"
mkdir -p build_baseline_cpu_complex
cd build_baseline_cpu_complex
cmake -DQMC_COMPLEX=1 $CMAKE_FLAGS ..
nice make -j 32
cd ..
ln -sf ./build_baseline_cpu_complex/bin/qmcpack_complex ./qmcpack_baseline_cpu_complex

