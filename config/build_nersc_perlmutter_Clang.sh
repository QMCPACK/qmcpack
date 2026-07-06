#!/bin/bash
# This recipe is intended for NERSC Perlmutter https://docs.nersc.gov/systems/perlmutter
# It builds all the varaints of QMCPACK in the current directory
# last revision: July 6th 2026
#
# How to invoke this script?
# build_nersc_perlmutter_Clang.sh # build all the variants assuming the current directory is the source directory.
# build_nersc_perlmutter_Clang.sh <source_dir> # build all the variants with a given source directory <source_dir>
# build_nersc_perlmutter_Clang.sh <source_dir> <install_dir> # build all the variants with a given source directory <source_dir> and install to <install_dir>

module load llvm/21.1.4
export MPICH_CC=`which clang`
export MPICH_CXX=`which clang++`

module load PrgEnv-gnu
module load cray-libsci
CRAY_LIBSCI_LIB=$CRAY_LIBSCI_PREFIX_DIR/lib/libsci_gnu_mp.so
module unload PrgEnv-gnu
module load craype cray-mpich
module load cray-fftw
module load cray-hdf5-parallel
module load cmake
export BOOST_ROOT=/global/common/software/m2113/opt/boost/1_82_0

echo "**********************************"
echo '$ clang -v'
clang -v
echo "**********************************"

TYPE=Release
Machine=perlmutter
Compiler=Clang21

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

for name in gpu_real_MP gpu_real gpu_cplx_MP gpu_cplx \
            cpu_real_MP cpu_real cpu_cplx_MP cpu_cplx
do

CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=$TYPE -DBLAS_LIBRARIES=$CRAY_LIBSCI_LIB"

if [[ $name == *"cplx"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DQMC_COMPLEX=ON"
fi

if [[ $name == *"_MP"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DQMC_MIXED_PRECISION=ON"
fi

if [[ $name == *"gpu"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DQMC_GPU_ARCHS=sm_80"
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
