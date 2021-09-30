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

source_folder=~/opt/qmcpack

for name in offload_cuda_real_MP offload_cuda_real offload_cuda_cplx_MP offload_cuda_cplx \
            cpu_real_MP cpu_real cpu_cplx_MP cpu_cplx
do

CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=$TYPE -DQMC_MATH_VENDOR=IBM_MASS -DMASS_ROOT=/sw/summit/xl/16.1.1-10/xlmass/9.1.1 -DMPIEXEC_EXECUTABLE=`which jsrun` -DMPIEXEC_NUMPROC_FLAG='-n' -DMPIEXEC_PREFLAGS='-c;16;-g;1;-b;packed:16;--smpiargs=off' -DCMAKE_CXX_STANDARD_LIBRARIES=/sw/summit/gcc/9.3.0-2/lib64/libstdc++.a"

if [[ $name == *"cplx"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DQMC_COMPLEX=ON"
fi

if [[ $name == *"_MP"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DQMC_MIXED_PRECISION=ON"
fi

if [[ $name == *"offload"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DENABLE_OFFLOAD=ON -DUSE_OBJECT_TARGET=ON -DOFFLOAD_ARCH=sm_70"
fi

if [[ $name == *"cuda"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=70 -DCMAKE_CUDA_HOST_COMPILER=/usr/bin/g++"
  CUDA_FLAGS="-Xcompiler -mno-float128"
fi

folder=build_summit_${Compiler}_${name}
echo "**********************************"
echo "$folder"
echo "$CMAKE_FLAGS"
echo "**********************************"
mkdir $folder
cd $folder
if [ ! -f CMakeCache.txt ] ; then
cmake $CMAKE_FLAGS -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_CUDA_FLAGS="$CUDA_FLAGS" $source_folder
cmake .
fi
make -j16
cd ..

echo
done
