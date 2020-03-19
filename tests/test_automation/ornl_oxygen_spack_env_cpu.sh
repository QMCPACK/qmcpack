#/bin/bash

. $SPACK_ROOT/share/spack/setup-env.sh

. ../tests/test_automation/spack_supported_package_versions.sh

spack load llvm@$llvm_vnew
spack load boost@$boost_vnew%gcc@$gcc_vnew
spack load gcc@$gcc_vnew
spack load hdf5@$hdf5_vnew%gcc@$gcc_vnew~mpi
spack load cmake@$cmake_vnew%gcc@$gcc_vnew
spack load -r openmpi@$ompi_vnew%clang@$llvm_vnew
spack load fftw@$fftw_vnew%gcc@$gcc_vnew
GCC_TOOLCHAIN_PATH=$(spack find -p gcc@$gcc_vnew | awk -v gcc_version="gcc@$gcc_vnew" -e '{ if ($1 == gcc_version) { print $2 }}')
QMC_IMMUTABLE_FLAGS="-DBUILD_AFQMC=1 -DCMAKE_CXX_FLAGS=--gcc-toolchain=${GCC_TOOLCHAIN_PATH};-stdlib=libstdc++ -DCMAKE_EXE_LINKER_FLAGS=-L${GCC_TOOLCHAIN_PATH}/lib;-Wl,-rpath,${GCC_TOOLCHAIN_PATH}/lib"
