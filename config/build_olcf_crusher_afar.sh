#!/bin/bash

# Build script for test and development system crusher at OLCF
# See https://github.com/QMCPACK/qmcpack/pull/4123 for more details on the module file if needed

echo "Loading QMCPACK dependency modules for crusher"
module load cmake/3.22.2
module load cray-fftw
module load openblas/0.3.17-omp
module load boost/1.77.0-cxx17
# private module until OLCF provides MPI compiler wrappers for afar compilers.
if [[ ! -d /ccs/proj/mat189/modules/crusher ]] ; then
  echo "Required module folder /ccs/proj/mat189/modules/crusher not found!"
  exit 1
fi
module use /ccs/proj/mat189/modules/crusher
module load mpiwrappers/cray-mpich-afar
module load cray-hdf5-parallel

TYPE=Release
Compiler=afar

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
  CMAKE_FLAGS="$CMAKE_FLAGS -DENABLE_CUDA=ON -DQMC_CUDA2HIP=ON -DHIP_ARCH=gfx90a"
fi

folder=build_crusher_${Compiler}_${name}
echo "**********************************"
echo "$folder"
echo "$CMAKE_FLAGS"
echo "**********************************"
mkdir $folder
cd $folder
if [ ! -f CMakeCache.txt ] ; then
cmake $CMAKE_FLAGS -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx $source_folder
fi
make -j16
cd ..

echo
done
