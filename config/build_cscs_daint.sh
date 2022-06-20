#!/bin/bash

# module purge
#echo "Purging current module set"

#BUILD_MODULES=config/load_cscs_daint_modules.sh
#echo "Sourcing file: $BUILD_MODULES to build QMCPACK"
#. $BUILD_MODULES

echo "Loading QMCPACK dependency modules for cscs piz-daint"
echo "https://user.cscs.ch/access/running/piz_daint/"
echo
module swap PrgEnv-cray PrgEnv-intel
module load daint-gpu
module load cudatoolkit
module load EasyBuild-custom/cscs
module load cray-hdf5-parallel
module load CMake
module load cray-python
module load Boost
# install libxml2 for CrayIntel
#eb libxml2-2.9.7-CrayIntel-20.08.eb -r
#module load libxml2/2.9.7-CrayIntel-20.08
module load libxml2
module unload cray-libsci
module unload cray-libsci_acc
# make sure there is a recent gcc compiler in the path
#module load gcc/8.3.0

module list

echo "Either source $BUILD_MODULES or load these same modules to run QMCPACK"

declare -A builds=( \
["cpu"]="        -DQMC_COMPLEX=0 -DQMC_CUDA=0" \
["complex_cpu"]="-DQMC_COMPLEX=1 -DQMC_CUDA=0" \
["legacy_gpu"]="        -DQMC_COMPLEX=0 -DQMC_CUDA=1 -DCMAKE_CUDA_ARCHITECTURES=60 -DENABLE_PHDF5=On -DCUDA_PROPAGATE_HOST_FLAGS=Off -DCUDA_HOST_COMPILER=`which gcc`" \
["complex_legacy_gpu"]="-DQMC_COMPLEX=1 -DQMC_CUDA=1 -DCMAKE_CUDA_ARCHITECTURES=60 -DENABLE_PHDF5=On -DCUDA_PROPAGATE_HOST_FLAGS=Off -DCUDA_HOST_COMPILER=`which gcc`" \
)

mkdir bin

for build in "${!builds[@]}"
do
    echo "building: $build with ${builds[$build]}"
    rm bin/qmcpack_${build}
    rm -rf build_${build}
    mkdir build_${build}
    cd build_${build}
    cmake \
          -DBUILD_LMYENGINE_INTERFACE=0 \
          -DQMC_MPI=On -DQMC_OMP=On \
          -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment \
          ${builds[$build]} \
          ..
    make -j 20
    if [ $? -eq 0 ]; then
      build_dir=$(pwd)
      ln -sf ${build_dir}/bin/qmcpack ${build_dir}/../bin/qmcpack_${build}
    fi
    cd ..
done
