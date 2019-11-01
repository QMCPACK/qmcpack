#!/bin/bash
module purge
echo "Purging current module set"

. config/load_summit_modules.sh

declare -A builds=( ["cpu"]=" " \
                    ["complex_cpu"]="-DQMC_COMPLEX=1 " \
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
      ln -sf ${build_dir}/bin/qmcpack ${build_dir}/../bin/qmcpack_${build}
    fi
    cd ..
done
