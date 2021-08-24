#!/bin/bash

BUILD_MODULES=config/load_olcf_summit_modules.sh

module purge
echo "Purging current module set"
echo "Sourcing file: $BUILD_MODULES to build QMCPACK"

. $BUILD_MODULES

echo "Either source $BUILD_MODULES or load these same modules to run QMCPACK"

declare -A builds=( ["cpu"]=" -DQMC_MATH_VENDOR=IBM_MASS -DMASS_ROOT=/sw/summit/xl/16.1.1-10/xlmass/9.1.1" \
                    ["complex_cpu"]="-DQMC_COMPLEX=1  -DQMC_MATH_VENDOR=IBM_MASS -DMASS_ROOT=/sw/summit/xl/16.1.1-10/xlmass/9.1.1" \
                    ["legacy_gpu"]="-DQMC_CUDA=1 -DCUDA_ARCH=sm_70 " \
		    ["complex_legacy_gpu"]="-DQMC_CUDA=1 -DQMC_COMPLEX=1 -DCUDA_ARCH=sm_70 " )

mkdir bin

for build in "${!builds[@]}"
do
    echo "building: $build with ${builds[$build]}"
    rm bin/qmcpack_${build}
    mkdir build_summit_${build}
    cd build_summit_${build}
    cmake -DCMAKE_C_COMPILER="mpicc" \
          -DCMAKE_CXX_COMPILER="mpicxx" \
          -DBUILD_LMYENGINE_INTERFACE=0 \
          ${builds[$build]} \
          ..
    make -j 20
    if [ $? -eq 0 ]; then
      build_dir=$(pwd)
      if [ -e ${build_dir}/bin/qmcpack_complex ]; then
        ln -sf ${build_dir}/bin/qmcpack_complex ${build_dir}/../bin/qmcpack_${build}
      else
        ln -sf ${build_dir}/bin/qmcpack ${build_dir}/../bin/qmcpack_${build}
      fi
    fi
    cd ..
done
