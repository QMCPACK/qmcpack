#!/bin/bash
module purge
echo "Purging current module set"

echo "Loading QMCPACK dependency modules for summit"
module load gcc
module load spectrum-mpi
module load essl
module load netlib-lapack
module load hdf5
module load fftw
export FFTW_ROOT=$OLCF_FFTW_ROOT
module load cmake
module load boost
module load cuda
module load python/2.7.15-anaconda2-5.3.0

declare -A builds=( ["legacy_gpu"]="-DQMC_CUDA=1 " \
		    ["complex_legacy_gpu"]="-DQMC_CUDA=1 -DQMC_COMPLEX=1 ")

mkdir bin

for build in "${!builds[@]}"
do
    echo "building: $build with ${builds[$build]}"
    mkdir build_summit_${build}
    cd build_summit_${build}
    cmake -DCMAKE_C_COMPILER="mpicc" \
      -DCMAKE_CXX_COMPILER="mpicxx" \
      -DBUILD_LMYENGINE_INTERFACE=0 \
      ${builds[$build]} \
      -DCUDA_ARCH="sm_70" \
      ..
    make -j 20
    cd ..
    ln -sf ./build_summit_${build}/bin/qmcpack ./bin/qmcpack_${build}
done
