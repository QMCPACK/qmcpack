#!/bin/bash
# This recipe is intended for ALCF Polaris https://www.alcf.anl.gov/polaris
# It builds all the varaints of QMCPACK in the current directory
# last revision: Aug 17th 2022
#
# How to invoke this script?
# build_alcf_polaris_Clang.sh # build all the variants assuming the current directory is the source directory.
# build_alcf_polaris_Clang.sh <source_dir> # build all the variants with a given source directory <source_dir>
# build_alcf_polaris_Clang.sh <source_dir> <install_dir> # build all the variants with a given source directory <source_dir> and install to <install_dir>

module load mpiwrappers/cray-mpich-llvm llvm/release-15.0.0
module load cudatoolkit-standalone/11.2.2
module load cray-fftw/3.3.8.13
module load cray-hdf5-parallel/1.12.1.3
module load cmake/3.23.2

export BOOST_ROOT=/soft/applications/qmcpack/boost_1_79_0
export CMAKE_PREFIX_PATH=/soft/libraries/openblas/0.3.20-omp:$CMAKE_PREFIX_PATH


echo "**********************************"
echo '$ clang -v'
clang -v
echo "**********************************"

TYPE=Release
Machine=polaris
Compiler=Clang

if [[ $# -eq 0 ]]; then
  source_folder=`pwd`
elif [[ $# -eq 1 ]]; then
  source_folder=$1
else
  source_folder=$1
  install_folder=$2
fi

if [[ -f $source_folder/CMakeLists.txt ]]; then
  echo Using QMCPACK source directory $source_folder
else
  echo "Source directory $source_folder doesn't contain CMakeLists.txt. Pass QMCPACK source directory as the first argument."
  exit
fi

for name in offload_cuda_real_MP offload_cuda_real offload_cuda_cplx_MP offload_cuda_cplx \
            cpu_real_MP cpu_real cpu_cplx_MP cpu_cplx
do

CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=$TYPE"

if [[ $name == *"cplx"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DQMC_COMPLEX=ON"
fi

if [[ $name == *"_MP"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DQMC_MIXED_PRECISION=ON"
fi

if [[ $name == *"offload"* || $name == *"cuda"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DQMC_GPU_ARCHS=sm_80"
fi

if [[ $name == *"offload"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DENABLE_OFFLOAD=ON"
fi

if [[ $name == *"cuda"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DENABLE_CUDA=ON"
fi

folder=build_${Machine}_${Compiler}_${name}

if [[ -v install_folder ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DCMAKE_INSTALL_PREFIX=$install_folder/$folder"
fi

echo "**********************************"
echo "$folder"
echo "$CMAKE_FLAGS"
echo "**********************************"

mkdir $folder
cd $folder

if [ ! -f CMakeCache.txt ] ; then
cmake $CMAKE_FLAGS -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx $source_folder
fi

if [[ -v install_folder ]]; then
  make -j16 install && chmod -R -w $install_folder/$folder
else
  make -j16
fi

cd ..

echo
done
