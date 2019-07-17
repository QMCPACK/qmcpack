#!/bin/bash
module purge
echo "Purging current module set"

echo "Loading QMCPACK dependency modules for summit"
module load gcc
module load spectrum-mpi/10.3.0.1-20190611
module load essl
module load netlib-lapack
module load hdf5
module load fftw
export FFTW_ROOT=$OLCF_FFTW_ROOT
module load cmake
module load boost
module load cuda
module load python/2.7.15-anaconda2-5.3.0
mkdir build_summit_legacy_cuda
cd build_summit_legacy_cuda
cmake -DCMAKE_C_COMPILER="mpicc" \
      -DCMAKE_CXX_COMPILER="mpicxx" \
      -DBUILD_LMYENGINE_INTERFACE=0 \
      -DQMC_CUDA=1 \
      -DCUDA_ARCH="sm_70" \
..
make -j 20
cd ..
ln -sf ./build_summit_legacy_cuda/bin/qmcpack ./qmcpack_summit_legacy_cuda

mkdir build_summit_complex_legacy_cuda
cd build_summit_complex_legacy_cuda
cmake -DCMAKE_C_COMPILER="mpicc" \
      -DCMAKE_CXX_COMPILER="mpicxx" \
      -DBUILD_LMYENGINE_INTERFACE=0 \
      -DQMC_CUDA=1 \
      -DQMC_COMPLEX=1 \
      -DCUDA_ARCH="sm_70" \
..
make -j 20
cd ..
ln -sf ./build_summit_complex_legacy_cuda/bin/qmcpack ./qmcpack_summit_complex_legacy_cuda

mkdir build_summit_enable_cuda
cd build_summit_enable_cuda
cmake -DCMAKE_C_COMPILER="mpicc" \
      -DCMAKE_CXX_COMPILER="mpicxx" \
      -DCUDA_TOOLKIT_ROOT_DIR=${OLCF_CUDA_ROOT} \
      -DENABLE_CUDA=1 \
      -DCUDA_ARCH="sm_70" \
..
make -j 20
cd ..
ln -sf ./build_summit_enable_cuda/bin/qmcpack ./qmcpack_summit_enable_cuda

mkdir build_summit_complex_enable_cuda
cd build_summit_complex_enable_cuda
cmake -DCMAKE_C_COMPILER="mpicc" \
      -DCMAKE_CXX_COMPILER="mpicxx" \
      -DCUDA_TOOLKIT_ROOT_DIR=${OLCF_CUDA_ROOT} \
      -DENABLE_CUDA=1 \
      -DQMC_COMPLEX=1 \
      -DCUDA_ARCH="sm_70" \
..
make -j 20
cd ..
ln -sf ./build_summit_complex_enable_cuda/bin/qmcpack ./qmcpack_summit_complex_enable_cuda


