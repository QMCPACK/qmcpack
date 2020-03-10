#!/bin/bash --login

. $SPACK_ROOT/share/spack/setup-env.sh

. ../tests/test_automation/spack_supported_package_versions.sh

QMC_IMMUTABLE_FLAGS="-DBUILD_AFQMC=1"

spack load gcc@$gcc_vnew
spack load boost@$boost_vnew%gcc@$gcc_vnew
spack load hdf5@$hdf5_vnew%gcc@$gcc_vnew~mpi
spack load cmake@$cmake_vnew%gcc@$gcc_vnew
spack load openmpi@$ompi_vnew%gcc@$gcc_vnew
spack load libxml2@$libxml2_vnew%gcc@$gcc_vnew
spack load fftw@$fftw_vnew%gcc@$gcc_vnew
