#!/bin/bash

# Build script for Frontier
# It builds all the varaints of QMCPACK in the current directory
# last revision: Apr 22nd 2025

echo "Loading QMCPACK dependency modules for frontier"
for module_name in PrgEnv-gnu PrgEnv-cray PrgEnv-amd PrgEnv-gnu-amd PrgEnv-cray-amd \
                   amd amd-mixed gcc gcc-mixed gcc-native cce cce-mixed rocm
do
  if module is-loaded $module_name ; then module unload $module_name; fi
done

module load PrgEnv-amd amd/6.3.1
module unload darshan-runtime
unset HIP_PATH # it messed up clang as a HIP compiler.
module unload cray-libsci
module load cray-fftw cray-hdf5-parallel
module load Core/25.03 cmake openblas/0.3.28-omp

# edit this line if you are not a member of mat151
export BOOST_ROOT=/ccs/proj/mat151/opt/boost/1_81_0

module list >& module_list.txt

TYPE=Release
Compiler=rocm631

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

CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=$TYPE -DMPIEXEC_EXECUTABLE=`which srun`"

if [[ $name == *"cplx"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DQMC_COMPLEX=ON"
fi

if [[ $name == *"_MP"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DQMC_MIXED_PRECISION=ON"
fi

if [[ $name == *"gpu"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DQMC_GPU_ARCHS=gfx90a"
fi

folder=build_frontier_${Compiler}_${name}

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
cmake $CMAKE_FLAGS -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC \
      -DCMAKE_CXX_FLAGS="-add-runpath" \
      $source_folder
fi

if [[ -v install_folder ]]; then
  make -j16 install && chmod -R -w $install_folder/$folder
else
  make -j16
fi

cd ..

echo
done
