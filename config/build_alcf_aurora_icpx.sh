#!/bin/bash
# This recipe is intended for ALCF Aurora https://www.alcf.anl.gov/support-center/aurora-sunspot
# last revision: Apr 25th 2025
#
# How to invoke this script?
# build_alcf_aurora_icpx.sh # build all the variants assuming the current directory is the source directory.
# build_alcf_aurora_icpx.sh <source_dir> # build all the variants with a given source directory <source_dir>
# build_alcf_aurora_icpx.sh <source_dir> <install_dir> # build all the variants with a given source directory <source_dir> and install to <install_dir>

for module_name in oneapi/release oneapi/eng-compiler
do
  if module is-loaded $module_name ; then module unload $module_name; fi
done

module load oneapi/release/2025.0.5
module load cmake hdf5 boost
module list >& module_list.txt

# unset the following to desensitize CMake to modules/environment variables.
unset CPATH
unset LIBRARY_PATH
unset C_INCLUDE_PATH
unset CPLUS_INCLUDE_PATH

echo "**********************************"
echo '$ icpx -v'
icpx -v
echo "**********************************"

TYPE=Release
Machine=aurora
Compiler=oneapi2025.0.5

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

CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=$TYPE -DMPIEXEC_PREFLAGS='--cpu-bind;depth;-d;8'"
unset CMAKE_CXX_FLAGS

if [[ $name == *"cplx"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DQMC_COMPLEX=ON"
fi

if [[ $name == *"_MP"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DQMC_MIXED_PRECISION=ON"
fi

if [[ $name == *"gpu"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DQMC_GPU_ARCHS=intel_gpu_pvc"
  CMAKE_CXX_FLAGS="-mllvm -vpo-paropt-atomic-free-reduction-slm=true"
fi

folder=build_${Machine}_${Compiler}_${name}

if [[ -v install_folder ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DCMAKE_INSTALL_PREFIX=$install_folder/$folder"
fi

echo "**********************************"
echo "folder $folder"
echo "CMAKE_FLAGS: $CMAKE_FLAGS"
echo "CMAKE_CXX_FLAGS: $CMAKE_CXX_FLAGS"
echo "**********************************"

mkdir $folder
cd $folder

if [ ! -f CMakeCache.txt ] ; then
cmake $CMAKE_FLAGS -DCMAKE_CXX_FLAGS="$CMAKE_CXX_FLAGS" \
      -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx $source_folder
fi

if [[ -v install_folder ]]; then
  make -j16 install && chmod -R -w $install_folder/$folder
else
  make -j16
fi

cd ..

echo
done
