#/bin/bash

SRC="${BASH_SOURCE[0]}"
SRC_DIR="$( cd -P "$( dirname "$SRC" )" >/dev/null 2>&1 && pwd )"
. ${SRC_DIR}/start_spack_env.sh

QMC_IMMUTABLE_FLAGS="-DQMC_CUDA=1"

spack load boost@$boost_vnew%gcc@$gcc_vnew
spack load gcc@$gcc_vcuda
spack load hdf5@$hdf5_vnew%gcc@$gcc_vcuda~mpi
spack load cmake@$cmake_vnew%gcc@$gcc_vnew
spack load -r openmpi@$ompi_vnew%gcc@$gcc_vcuda
spack load fftw@$fftw_vnew%gcc@$gcc_vnew
