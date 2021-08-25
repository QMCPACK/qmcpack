#!/bin/bash

echo "----------------------- WARNING ------------------------------------"
echo "This is **not** production ready and intended for development only!!"
echo "Use config/build_olcf_summit.sh for production on Summit."
echo "----------------------- WARNING ------------------------------------"

echo "Purging current module set"
module purge
echo "Loading QMCPACK dependency modules for summit"
module load gcc/9.3.0
module load spectrum-mpi
module load cmake
module load git
# default cuda/11.0.3 causes frequent hanging with the offload code in production runs.
module load cuda/10.1.243
module load essl
module load netlib-lapack
module load hdf5
module load fftw
module load boost/1.76.0
module load python/3.8-anaconda3
# private module until OLCF provides a new llvm build
if [[ ! -d /ccs/proj/mat151/opt/modules ]] ; then
  echo "Required module folder /ccs/proj/mat151/opt/modules not found!"
  exit 1
fi
module use /ccs/proj/mat151/opt/modules
# requires cuda/10.1.243 to have safe production runs.
module load llvm/main-20210811-cuda10.1

TYPE=Release
Compiler=Clang

source_folder=..

for name in offload_cuda_real_MP offload_cuda_real offload_cuda_cplx_MP offload_cuda_cplx \
            cpu_real_MP cpu_real cpu_cplx_MP cpu_cplx
do

CMAKE_FLAGS="-D CMAKE_BUILD_TYPE=$TYPE -D QMC_MATH_VENDOR=IBM_MASS -D MASS_ROOT=/sw/summit/xl/16.1.1-10/xlmass/9.1.1 -D MPIEXEC_EXECUTABLE=`which jsrun` -D MPIEXEC_NUMPROC_FLAG='-n' -D MPIEXEC_PREFLAGS='-c;16;-g;1;-b;packed:16;--smpiargs=off'"

if [[ $name == *"cplx"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -D QMC_COMPLEX=1"
fi

if [[ $name == *"_MP"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -D QMC_MIXED_PRECISION=1"
fi

if [[ $name == *"offload"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -D ENABLE_OFFLOAD=ON -D USE_OBJECT_TARGET=ON -DOFFLOAD_ARCH=sm_70"
fi

if [[ $name == *"cuda"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -D ENABLE_CUDA=1 -D CUDA_ARCH=sm_70 -D CUDA_HOST_COMPILER=/usr/bin/gcc -D CUDA_NVCC_FLAGS='-Xcompiler;-mno-float128'"
fi

folder=build_summit_${Compiler}_${name}
echo "**********************************"
echo "$folder"
echo "$CMAKE_FLAGS"
echo "**********************************"
mkdir $folder
cd $folder
if [ ! -f CMakeCache.txt ] ; then
cmake $CMAKE_FLAGS -D CMAKE_C_COMPILER=mpicc -D CMAKE_CXX_COMPILER=mpicxx $source_folder
cmake .
fi
make -j16
cd ..

echo
done
