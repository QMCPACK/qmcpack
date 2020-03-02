#/bin/bash

. $SPACK_ROOT/share/spack/setup-env.sh

. ../tests/test_automation/spack_supported_package_versions.sh

QMC_IMMUTABLE_FLAGS="-DENABLE_CUDA=1"

spack load boost@$boost_vnew%gcc@$gcc_vnew
spack load gcc@$gcc_vcuda
spack load hdf5@$hdf5_vnew%gcc@$gcc_vcuda~mpi
spack load cmake@$cmake_vnew%gcc@$gcc_vnew
spack load openmpi@$ompi_vnew%gcc@$gcc_vcuda
spack load libxml2@$libxml2_vnew%gcc@$gcc_vnew
spack load fftw@$fftw_vnew%gcc@$gcc_vnew
