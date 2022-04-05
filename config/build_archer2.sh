#!/bin/bash

BUILD_MODULES=config/load_archer2.sh

#module purge
#echo "Purging current module set"
echo "Sourcing file: $BUILD_MODULES to build QMCPACK"

. $BUILD_MODULES

module list

echo "Either source $BUILD_MODULES or load these same modules to run QMCPACK"

#export BLAS_LIBS="-L$OLCF_OPENBLAS_ROOT/lib -lopenblas"
#export LAPACK_LIBS="$BLAS_LIBS $OLCF_NETLIB_LAPACK_ROOT/lib64/liblapack.a"

declare -A builds=( ["cpu"]="-DBUILD_PPCONVERT=1" \
                    ["complex_cpu"]="-DQMC_COMPLEX=1" \
#		    ["legacy_gpu"]="-DQMC_CUDA=1 " \
#		    ["complex_legacy_gpu"]="-DQMC_CUDA=1 -DQMC_COMPLEX=1 " \
		  )

mkdir bin

for build in "${!builds[@]}"
do
    echo "building: $build with ${builds[$build]}"
    rm bin/qmcpack_${build}
    mkdir build_${build}
    cd build_${build}
    cmake -DCMAKE_C_COMPILER="cc" \
          -DCMAKE_CXX_COMPILER="CC" \
	  -DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment \
          -D LIBXML2_LIBRARY="$LIBXML2_ROOT/lib/libxml2.so" \
          -D LIBXML2_INCLUDE_DIR="$LIBXML2_ROOT//include/libxml2" \
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

