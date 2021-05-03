#!/bin/bash

echo "----------------------- WARNING ------------------------------------"
echo "This is **not** production ready and intended for development only!!"
echo "Use config/build_olcf_summit.sh for production on Summit."
echo "----------------------- WARNING ------------------------------------"

echo "Purging current module set"
module purge
echo "Loading QMCPACK dependency modules for summit"
module load xl
module load spectrum-mpi
module load cmake
module load git
module load cuda
module load essl
module load netlib-lapack
module load hdf5
module load python/3.6.6-anaconda3-5.3.0

#the XL built fftw is buggy, use the gcc version
#module load fftw
export FFTW_HOME=/autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/gcc-6.4.0/fftw-3.3.8-5gcj2ic4el7acu3rqnfnh735jz2ez7j5
export BOOST_ROOT=/autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/gcc-6.4.0/boost-1.66.0-l3sghp3ggjzwi4vtvyb5yzsjm36npgrk

TYPE=Release
Compiler=XL

source_folder=..

for name in offload_real_MP offload_real # offload_cplx offload_cplx_MP
do

CMAKE_FLAGS="-D CMAKE_BUILD_TYPE=$TYPE -D ENABLE_CUDA=1 -D CUDA_ARCH=sm_70 -D CUDA_HOST_COMPILER=/sw/summit/gcc/6.4.0/bin/gcc -D ENABLE_MASS=1 -D MASS_ROOT=/sw/summit/xl/16.1.1-5/xlmass/9.1.1 -D MPIEXEC_EXECUTABLE=`which jsrun` -D MPIEXEC_NUMPROC_FLAG='-n' -D MPIEXEC_PREFLAGS='-c;16;-g;1;-b;packed:16;--smpiargs=off'"

if [[ $name == *"cplx"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -D QMC_COMPLEX=1"
fi

if [[ $name == *"_MP"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -D QMC_MIXED_PRECISION=1"
fi

if [[ $name == *"offload"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -D ENABLE_OFFLOAD=1"
fi

folder=build_summit_${Compiler}_${name}
echo "**********************************"
echo "$folder"
echo "$CMAKE_FLAGS"
echo "**********************************"
mkdir $folder
cd $folder
if [ ! -f CMakeCache.txt ] ; then
cmake $CMAKE_FLAGS -D CMAKE_C_COMPILER=mpixlc -D CMAKE_CXX_COMPILER=mpixlC -D CMAKE_C_FLAGS=-qarch=pwr9 \
  -D CMAKE_CXX_FLAGS="-qarch=pwr9 -qxflag=disable__cplusplusOverride -isystem /sw/summit/gcc/6.4.0/include/c++/6.4.0/powerpc64le-none-linux-gnu -qgcc_cpp_stdinc=/sw/summit/gcc/6.4.0/include/c++/6.4.0" \
  -D CMAKE_CXX_STANDARD_LIBRARIES=/sw/summit/gcc/6.4.0/lib64/libstdc++.a \
  -D BLAS_essl_LIBRARY=$OLCF_ESSL_ROOT/lib64/libessl.so $source_folder
cmake .
fi
make -j24
cd ..

echo
done
