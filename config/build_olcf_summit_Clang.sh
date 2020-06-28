#!/bin/bash

echo "----------------------- WARNING ------------------------------------"
echo "This is **not** production ready and intended for development only!!"
echo "Use config/build_olcf_summit.sh for production on Summit."
echo "----------------------- WARNING ------------------------------------"

echo "Purging current module set"
module purge
echo "Loading QMCPACK dependency modules for summit"
module load gcc/8.1.1
module load spectrum-mpi
module load cmake
module load git
module load cuda
module load essl
module load netlib-lapack
module load hdf5
module load python/3.6.6-anaconda3-5.3.0
# private module until OLCF provides a new llvm build
module load llvm/master-latest

#the XL built fftw is buggy, use the gcc version
#module load fftw
export FFTW_HOME=/autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/gcc-6.4.0/fftw-3.3.8-5gcj2ic4el7acu3rqnfnh735jz2ez7j5
export BOOST_ROOT=/autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/gcc-6.4.0/boost-1.66.0-l3sghp3ggjzwi4vtvyb5yzsjm36npgrk

TYPE=Release
Compiler=Clang

CURRENT_FOLDER=`pwd`

for name in offload_real_MP offload_real offload_cplx offload_cplx_MP
do

CMAKE_FLAGS="-D CMAKE_BUILD_TYPE=$TYPE -D ENABLE_CUDA=1 -D CUDA_ARCH=sm_70 -D ENABLE_MASS=1 -D MASS_ROOT=/sw/summit/xl/16.1.1-5/xlmass/9.1.1 -D MPIEXEC_EXECUTABLE=/bin/sh -D MPIEXEC_NUMPROC_FLAG=$CURRENT_FOLDER/tests/scripts/jsrunhelper.sh"
if [[ $name == *"cplx"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -D QMC_COMPLEX=1"
fi

if [[ $name == *"_MP"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -D QMC_MIXED_PRECISION=1"
fi

if [[ $name == *"offload"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -D ENABLE_OFFLOAD=ON -D CUDA_HOST_COMPILER=`which gcc` -D USE_OBJECT_TARGET=ON"
fi

folder=build_summit_${Compiler}_${name}
echo "**********************************"
echo "$folder"
echo "$CMAKE_FLAGS"
echo "**********************************"
mkdir $folder
cd $folder
if [ ! -f CMakeCache.txt ] ; then
cmake $CMAKE_FLAGS -D CMAKE_C_COMPILER=mpicc -D CMAKE_CXX_COMPILER=mpicxx -D ENABLE_TIMERS=1 ..
cmake ..
fi
make -j24
cd ..

echo
done
