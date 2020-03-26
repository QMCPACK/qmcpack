#/bin/bash

. $SPACK_ROOT/share/spack/setup-env.sh

# This is to deal with the fact that we don't know the relationship between the working directory
# and checkout directory in a relative way.
# this doesn't deal with symlinks
SRC="${BASH_SOURCE[0]}"
SRC_DIR="$( cd -P "$( dirname "$SRC" )" >/dev/null 2>&1 && pwd )"

. ${SRC_DIR}/spack_supported_package_versions.sh

QMC_IMMUTABLE_FLAGS="-DENABLE_CUDA=1"

spack load boost@$boost_vnew%gcc@$gcc_vnew
spack load gcc@$gcc_vcuda
spack load hdf5@$hdf5_vnew%gcc@$gcc_vcuda~mpi
spack load cmake@$cmake_vnew%gcc@$gcc_vnew
spack load -r openmpi@$ompi_vnew%gcc@$gcc_vcuda
spack load fftw@$fftw_vnew%gcc@$gcc_vnew
