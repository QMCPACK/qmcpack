#!/usr/bin/env bash                                                                                                                                                                                                

# This recipe is intended for NERSC Perlmutter https://docs.nersc.gov/systems/perlmutter/
# It builds the CPU varaints of QMCPACK in the current directory
# Last revision: 2023/01/13
#
# How to invoke this script?
# build_nersc_perlmutter.sh # build all the variants assuming the current directory is the source directory
# build_nersc_perlmutter.sh <source_dir> # build all the variants with a given source directory <source_dir>
# build_nersc_perlmutter.sh <source_dir> <install_dir> # build all the variants with a given source directory <source_dir> and install to <install_dir>

# Load the required modules
module load cpu
module load PrgEnv-gnu
module load cray-hdf5-parallel
module load cray-fftw
module load cmake
module list

# Make LibXml2 available to QMCPACK
export CMAKE_PREFIX_PATH=/global/common/software/spackecp/perlmutter/e4s-22.05/78535/spack/opt/spack/cray-sles15-zen3/gcc-11.2.0/libxml2-2.9.13-u2ai4xjq2lmljvej4p3ly7qd6hfbrz7h:$CMAKE_PREFIX_PATH

# On the date above, the loaded modules were:
#Currently Loaded Modules:
#  1) craype-x86-milan     4) xpmem/2.5.2-2.4_3.18__gd0f7936.shasta   7) cpe/22.11      10) cpu/1.0           13) cray-mpich/8.1.22      16) cray-hdf5-parallel/1.12.2.1
#  2) libfabric/1.15.2.0   5) gcc/11.2.0                              8) xalt/2.10.2    11) craype/2.7.19     14) cray-libsci/22.11.1.2  17) cray-fftw/3.3.10.2
#  3) craype-network-ofi   6) perftools-base/22.09.0                  9) darshan/3.4.0  12) cray-dsmml/0.2.2  15) PrgEnv-gnu/8.3.3       18) cmake/3.24.3

TYPE=Release
Machine=perlmutter
Compiler=GNU

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

for name in cpu_real_MP cpu_real cpu_cplx_MP cpu_cplx
do

CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=$TYPE"

if [[ $name == *"cplx"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DQMC_COMPLEX=ON"
fi

if [[ $name == *"_MP"* ]]; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DQMC_MIXED_PRECISION=ON"
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
cmake $CMAKE_FLAGS -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment $source_folder
fi

if [[ -v install_folder ]]; then
  make -j16 install && chmod -R -w $install_folder/$folder
else
  make -j16
fi

cd ..

echo
done
