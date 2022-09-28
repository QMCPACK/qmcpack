#!/bin/bash

echo "ARCHER2: Information on hardware and software"
echo "https://www.archer2.ac.uk/about/hardware.html"
echo "and documentation:"
echo "https://docs.archer2.ac.uk"
echo

echo "Loading QMCPACK dependency modules for archer2"
echo
module restore
module load PrgEnv-gnu
module load cray-hdf5-parallel
module load cray-fftw
export FFTW_ROOT=$FFTW_DIR/..
module load libxml2
module load cmake
module load boost
module load cray-python
echo
echo "Loaded moduli:"
module list

echo
echo "In the running scipt (but not in compilation) also load the following two modules:"
echo "   module load craype-network-ucx"
echo "   module load cray-mpich-ucx"
echo "which improves a lot the scaling efficiency. "
echo
echo


declare -A builds=( ["cpu"]="-DBUILD_PPCONVERT=1" \
                    ["complex_cpu"]="-DQMC_COMPLEX=1" \
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
          -D LibXml2_ROOT=$LIBXML2_ROOT \
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

