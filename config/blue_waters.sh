#!/bin/bash

################################################################
##   This script builds available configurations of QMCPACK   ##
##   on Blue Waters, University of Illinois Urbana-Champaign  ##
##   highly similar to OLCF Titan build script                ##
##                                                            ##
## Last modified: 19 May 2020                                 ##
################################################################

# Load required modules (assuming default settings have not been modified)
source $MODULESHOME/init/bash
if (echo $LOADEDMODULES | grep -q pgi)
then
module unload PrgEnv-pgi
fi
if (echo $LOADEDMODULES | grep -q cray)
then
module unload PrgEnv-cray
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

# for C++14 support (gcc/6.3.0)
module load PrgEnv-gnu
# disable user statistics reporting (otherwise crash)
module unload darshan
# use parallel HDF5 to speed up I/O of large jobs
module load cray-hdf5-parallel

## FFT library is important for orbital splining
#module load fftw
# 2020-01-30: fftw breaks mpi? set FFTW_HOME instead of load module
## 2019-10-25: bwpy conflicts with cmake and boost
#module load boost/1.63.0
# 2020-01-30: python3 is required, set BOOST_ROOT instead of load module

# miscellaneous
#module load bwpy # numpy, h5py libraries are used in ctest
module load bwpy
module load libxml2
module load cmake/3.9.4

# always use dynamic linking
export CRAYPE_LINK_TYPE=dynamic

# use AMD optimized math libraries (performance critical!)
XT_FLAGS="-DHAVE_AMDLIBM=1"

# Set cmake variables, shared for cpu builds; not used in gpu builds
AMD_LIB_HOME=/projects/sciteam/bbak/soft/amdlibm-3-0-2
CMAKE_FLAGS="-D CMAKE_C_FLAGS=$XT_FLAGS \
  -D CMAKE_CXX_FLAGS=$XT_FLAGS \
  -D QMC_INCLUDE=$AMD_LIB_HOME/include \
  -D QMC_EXTRA_LIBS=$AMD_LIB_HOME/lib/static/libamdlibm.a
"

# Set environment variables
export FFTW_HOME=/opt/cray/fftw/3.3.4.10/interlagos
export BOOST_ROOT=/sw/xe/boost/1.63.0/sles11.3_gnu5.3.0

export CC=cc
export CXX=CC

################################################################
## CPU Binaries                                               ##
################################################################

target=qmcpack

# Configure and build cpu real
suffix=_cpu_real
echo ""
echo ""
echo "building qmcpack for cpu real"
mkdir build$suffix
cd build$suffix
cmake $CMAKE_FLAGS ..
make -j 32 $target
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
make -j 32 $target
cd ..
ln -s ./build$suffix/bin/qmcpack_complex ./qmcpack$suffix

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
make -j 32 $target
cd ..
ln -s ./build$suffix/bin/qmcpack ./qmcpack$suffix

# Configure and build gpu complex
suffix=_gpu_comp
echo ""
echo ""
echo "building qmcpack for gpu complex"
mkdir build$suffix
cd build$suffix
cmake -D QMC_COMPLEX=1 -D QMC_CUDA=1 -DCUDA_HOST_COMPILER=$(which CC) ..
cmake -D QMC_COMPLEX=1 -D QMC_CUDA=1 -DCUDA_HOST_COMPILER=$(which CC) ..
make -j 32 $target
cd ..
ln -s ./build$suffix/bin/qmcpack_complex ./qmcpack$suffix
