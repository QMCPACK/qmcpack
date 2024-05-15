#!/bin/bash

################################################################
## * This script builds available configurations of QMCPACK   ##
##   on Improv, at Argonne National Lab.                      ##
##                                                            ##
## * Execute this script in trunk/                            ##
##   ./config/build_lcrc_improv.sh                            ##
##                                                            ##
## Last verified: May 15, 2024                                ##
################################################################

# module files resulting from module imports below:
# Currently Loaded Modulefiles:
#  1) xalt/3.0.1     (S)   4) gcc/13.2.0                 7) amdblis/4.2-gcc-13.2.0
#  2) cmake/3.27.4         5) openmpi/5.0.2-gcc-13.2.0   8) aocl-utils/4.2-gcc-13.2.0
#  3) .binutils/2.41 (H)   6) python/3.11.6-gcc-13.2.0   9) amdlibflame/4.2-gcc-13.2.0


source $MODULESHOME/init/bash
module purge
module load cmake/3.27.4
module load gcc/13.2.0
module load openmpi/5.0.2-gcc-13.2.0
module load amdlibflame/4.2-gcc-13.2.0
module list


CMAKE_FLAGS="-DENABLE_PPCONVERT=0 \
             -DLAPACK_LIBRARIES=/gpfs/fs1/soft/improv/software/spack-built/linux-rhel8-zen3/gcc-13.2.0/openblas-0.3.26-blohgyt/lib/libopenblas.so \
	     -DFFTW_LIBRARIES=/gpfs/fs1/soft/improv/software/spack-built/linux-rhel8-zen3/gcc-13.2.0/fftw-3.3.10-x5237xr/lib/libfftw3.so \
	     -DFFTW_INCLUDE_DIR=/gpfs/fs1/soft/improv/software/spack-built/linux-rhel8-zen3/gcc-13.2.0/fftw-3.3.10-x5237xr/include \
	     -DLIBXML2_INCLUDE_DIR=/gpfs/fs1/soft/improv/software/spack-built/linux-rhel8-zen3/gcc-13.2.0/libxml2-2.10.3-xkoaaap/include/libxml2 \
	     -DLIBXML2_LIBRARY=/gpfs/fs1/soft/improv/software/spack-built/linux-rhel8-zen3/gcc-13.2.0/libxml2-2.10.3-xkoaaap/lib/libxml2.so \
	     -DBOOST_ROOT=/gpfs/fs1/soft/improv/software/custom-built/boost/1.84.0 \
	     -DHDF5_ROOT=/gpfs/fs1/soft/improv/software/spack-built/linux-rhel8-zen3/gcc-13.2.0/hdf5-1.14.3-6qo7t6e \
             -DCMAKE_C_COMPILER=mpicc \
             -DCMAKE_CXX_COMPILER=mpicxx"


# Configure and build cpu real.
echo ""
echo ""
echo "building QMCPACK for cpu real for Improv"
mkdir -p build_improv_cpu_real
cd build_improv_cpu_real
cmake $CMAKE_FLAGS ..
make -j 32
cd ..
ln -sf ./build_improv_cpu_real/bin/qmcpack ./qmcpack_improv_cpu_real

# Configure and build cpu complex.
echo ""
echo ""
echo "building QMCPACK for cpu complex for Improv"
mkdir -p build_improv_cpu_complex
cd build_improv_cpu_complex
cmake -DQMC_COMPLEX=1 $CMAKE_FLAGS ..
make -j 32
cd ..
ln -sf ./build_improv_cpu_complex/bin/qmcpack_complex ./qmcpack_improv_cpu_complex
