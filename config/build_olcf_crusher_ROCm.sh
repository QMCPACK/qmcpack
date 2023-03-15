#!/bin/bash

# Build script for test and development system crusher at OLCF
# See https://github.com/QMCPACK/qmcpack/pull/4123 for more details on the module file if needed

echo "Loading QMCPACK dependency modules for crusher"
module unload PrgEnv-gnu PrgEnv-cray PrgEnv-amd PrgEnv-gnu-amd PrgEnv-cray-amd
module unload amd amd-mixed gcc gcc-mixed cce cce-mixed
module load PrgEnv-amd amd/5.4.0 gcc-mixed/11.2.0
module unload cray-libsci
module load cmake/3.22.2
module load cray-fftw
module load openblas/0.3.17-omp
module load cray-hdf5-parallel
module load boost/1.78.0

module list >& module_list.txt

TYPE=Release
Compiler=rocm540

if [[ $# -eq 0 ]]; then
  source_folder=`pwd`
else
  source_folder=$1
fi

if [[ -f $source_folder/CMakeLists.txt ]]; then
  echo Using QMCPACK source directory $source_folder
else
  echo "Source directory $source_folder doesn't contain CMakeLists.txt. Pass QMCPACK source directory as the first argument."
  exit
fi

for name in offload_cuda2hip_real_MP offload_cuda2hip_real offload_cuda2hip_cplx_MP offload_cuda2hip_cplx \
            cpu_real_MP cpu_real cpu_cplx_MP cpu_cplx
do

CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=$TYPE -DMPIEXEC_EXECUTABLE=`which srun`"

if [[ $name == *"cplx"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DQMC_COMPLEX=ON"
fi

if [[ $name == *"_MP"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DQMC_MIXED_PRECISION=ON"
fi

if [[ $name == *"offload"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DENABLE_OFFLOAD=ON"
fi

if [[ $name == *"cuda2hip"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DENABLE_CUDA=ON -DQMC_CUDA2HIP=ON -DCMAKE_HIP_ARCHITECTURES=gfx90a"
fi

folder=build_crusher_${Compiler}_${name}
echo "**********************************"
echo "$folder"
echo "$CMAKE_FLAGS"
echo "**********************************"
mkdir $folder
cd $folder
if [ ! -f CMakeCache.txt ] ; then
cmake $CMAKE_FLAGS -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment \
      -DCMAKE_C_FLAGS=--gcc-toolchain=/opt/cray/pe/gcc/11.2.0/snos -DCMAKE_CXX_FLAGS=--gcc-toolchain=/opt/cray/pe/gcc/11.2.0/snos \
      $source_folder
fi
make -j16
cd ..

echo
done
