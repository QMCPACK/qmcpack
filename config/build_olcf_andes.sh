#!/bin/bash

BUILD_MODULES=config/load_olcf_andes_modules.sh

module purge
echo "Purging current module set"
echo "Sourcing file: $BUILD_MODULES to build QMCPACK"

. $BUILD_MODULES

echo "Either source $BUILD_MODULES or load these same modules to run QMCPACK"

export BLAS_LIBS="-L$OLCF_OPENBLAS_ROOT/lib -lopenblas"
export LAPACK_LIBS="$BLAS_LIBS $OLCF_NETLIB_LAPACK_ROOT/lib64/liblapack.a"

declare -A builds=( ["cpu"]="-DBUILD_PPCONVERT=1" \
                    ["complex_cpu"]="-DQMC_COMPLEX=1" \
		  )

mkdir bin_andes

for build in "${!builds[@]}"
do
    echo "building: $build with ${builds[$build]}"
    rm bin_andes/qmcpack_${build}
    mkdir build_andes_${build}
    cd build_andes_${build}
    cmake -DCMAKE_C_COMPILER="mpicc" \
          -DCMAKE_CXX_COMPILER="mpicxx" \
          -DBUILD_LMYENGINE_INTERFACE=0 \
          ${builds[$build]} \
          ..
    make -j 20
    if [ $? -eq 0 ]; then
      build_dir=$(pwd)
      if [ -e ${build_dir}/bin/qmcpack_complex ]; then
        ln -sf ${build_dir}/bin/qmcpack_complex ${build_dir}/../bin_andes/qmcpack_${build}
      else
        ln -sf ${build_dir}/bin/qmcpack ${build_dir}/../bin_andes/qmcpack_${build}
      fi
    fi
    cd ..
done

