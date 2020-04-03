#!/bin/bash

SRC="${BASH_SOURCE[0]}"
SRC_DIR="$( cd -P "$( dirname "$SRC" )" >/dev/null 2>&1 && pwd )"
. ${SRC_DIR}/start_spack_env.sh

spack load gcc@$gcc_vnew%gcc@$gcc_vold
spack load boost@$boost_vnew%gcc@$gcc_vnew
spack load hdf5@$hdf5_vnew%gcc@$gcc_vnew~mpi
spack load openssl@1.1.1d%gcc@9.2.0+systemcerts
spack load cmake@$cmake_vnew%gcc@$gcc_vnew
spack load hwloc%gcc@$gcc_vnew
spack load libiconv%gcc@$gcc_vnew
spack load openmpi@$ompi_vnew%gcc@$gcc_vnew
spack load libxml2@$libxml2_vnew%gcc@$gcc_vnew
spack load fftw@$fftw_vnew%gcc@$gcc_vnew
spack load openblas%gcc@$gcc_vnew
spack load netlib-lapack%gcc@$gcc_vnew

#if you've installed more than one python for the new compiler this will fail
spack load python%gcc@$gcc_vnew

QMC_IMMUTABLE_FLAGS="-DBUILD_AFQMC=1"

echo ${PATH}
