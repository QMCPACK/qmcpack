#!/bin/bash
# This recipe is intended for ALCF Sunspot https://www.alcf.anl.gov/support-center/aurora-sunspot
# last revision: Jan 8th 2023
#
# How to invoke this script?
# build_alcf_sunspot_icpx.sh # build all the variants assuming the current directory is the source directory.
# build_alcf_sunspot_icpx.sh <source_dir> # build all the variants with a given source directory <source_dir>
# build_alcf_sunspot_icpx.sh <source_dir> <install_dir> # build all the variants with a given source directory <source_dir> and install to <install_dir>

module load spack libxml2 cmake
module load cray-hdf5
module load oneapi/eng-compiler/2023.05.15.007

module list >& module_list.txt

# edit this line for your own boost header files.
export BOOST_ROOT=/home/yeluo/opt/boost_1_80_0

echo "**********************************"
echo '$ icpx -v'
icpx -v
echo "**********************************"

TYPE=Release
Machine=sunspot
Compiler=icpx20230613

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

for name in offload_sycl_real_MP offload_sycl_real offload_sycl_cplx_MP offload_sycl_cplx \
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

if [[ $name == *"offload"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DENABLE_OFFLOAD=ON"
  CMAKE_CXX_FLAGS="-mllvm -vpo-paropt-atomic-free-reduction-slm=true"
fi

if [[ $name == *"sycl"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DENABLE_SYCL=ON"
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
