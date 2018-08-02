#!/bin/bash

################################################################
##   This script builds available configurations of QMCPACK   ##
##   on Blue Waters, University of Illinois Urbana-Champaign  ##
##   highly similar to OLCF Titan build script                ##
##                                                            ##
## Last modified: Apr 18, 2018                                ##
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
# bwpy conflicts with boost and cmake?                                          
if (echo $LOADEDMODULES | grep -q bwpy)                                         
then                                                                            
module unload bwpy                                                              
fi 

# use parallel HDF5 to speed up I/O of large jobs
module load cray-hdf5-parallel

# FFT library is important for orbital splining
module load fftw

# miscellaneous
module load boost
module load libxml2
module load cmake
#module load bwpy # numpy, h5py libraries are used in ctest


# always use dynamic linking
export CRAYPE_LINK_TYPE=dynamic

# use AMD optimized math libraries (performance critical!)
XT_FLAGS="-DHAVE_AMDLIBM=1"

# Set cmake variables, shared for cpu builds; not used in gpu builds
CMAKE_FLAGS="-D CMAKE_C_FLAGS=$XT_FLAGS \
  -D CMAKE_CXX_FLAGS=$XT_FLAGS \
  -D QMC_INCLUDE=/u/staff/rmokos/libs/amdlibm/include \
  -D QMC_EXTRA_LIBS=/u/staff/rmokos/libs/amdlibm/lib/static/libamdlibm.a
"

# Set environment variables
export FFTW_HOME=$FFTW_DIR/..

export CC=cc
export CXX=CC

################################################################
## CPU Binaries                                               ##
################################################################

# Configure and build cpu real
suffix=_cpu_real
echo ""
echo ""
echo "building qmcpack for cpu real"
mkdir build$suffix
cd build$suffix
cmake $CMAKE_FLAGS ..
make -j 32
cd ..
ln -s ./build$suffix/bin/qmcpack ./qmcpack$suffix

# Configure and build cpu complex
suffix=_cpu_comp
echo ""
echo ""
echo "building qmcpack for cpu complex"
mkdir build$suffix
cd build$suffix
cmake $CMAKE_FLAGS -D QMC_COMPLEX=1 ..
make -j 32
cd ..
ln -s ./build$suffix/bin/qmcpack ./qmcpack$suffix

# Configure and build cpu Structure-of-Array (SoA) optimized binaries
suffix=_cpu_real_SoA
echo ""
echo ""
echo "building qmcpack for cpu real"
mkdir build$suffix
cd build$suffix
cmake $CMAKE_FLAGS -D ENABLE_SOA=1 ..
make -j 32
cd ..
ln -s ./build$suffix/bin/qmcpack ./qmcpack$suffix

suffix=_cpu_comp_SoA
echo ""
echo ""
echo "building qmcpack for cpu complex"
mkdir build$suffix
cd build$suffix
cmake $CMAKE_FLAGS -D QMC_COMPLEX=1 -D ENABLE_SOA=1 ..
make -j 32
cd ..
ln -s ./build$suffix/bin/qmcpack ./qmcpack$suffix

################################################################
## GPU Binaries                                               ##
################################################################

module load cudatoolkit

# Configure and build gpu real
suffix=_gpu_real
echo ""
echo ""
echo "building qmcpack for gpu real"
mkdir build$suffix
cd build$suffix
cmake -D QMC_CUDA=1 -DCUDA_HOST_COMPILER=$(which CC) ..
cmake -D QMC_CUDA=1 -DCUDA_HOST_COMPILER=$(which CC) ..
make -j 32
cd ..
ln -s ./build$suffix/bin/qmcpack ./qmcpack$suffix

# Configure and build gpu complex
suffix=_gpu_comp
echo ""
echo ""
echo "building qmcpack for gpu real"
mkdir build$suffix
cd build$suffix
cmake -D QMC_COMPLEX=1 -D QMC_CUDA=1 -DCUDA_HOST_COMPILER=$(which CC) ..
cmake -D QMC_COMPLEX=1 -D QMC_CUDA=1 -DCUDA_HOST_COMPILER=$(which CC) ..
make -j 32
cd ..
ln -s ./build$suffix/bin/qmcpack ./qmcpack$suffix
