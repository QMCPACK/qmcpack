#/bin/bash

. $SPACK_ROOT/share/spack/setup-env.sh

# This is to deal with the fact that we don't know the relationship between the working directory
# and checkout directory in a relative way.
# This is from https://stackoverflow.com/questions/59895/how-to-get-the-source-directory-of-a-bash-script-from-within-the-script-itself
# licensing issue?
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"

. ${DIR}/tests/test_automation/spack_supported_package_versions.sh

#temporary change from clang$llvm_vnew

spack load gcc@$gcc_vnew
spack load boost@$boost_vnew%gcc@$gcc_vnew
spack load gcc@$gcc_vnew
spack load hdf5@$hdf5_vnew%gcc@$gcc_vnew~mpi
spack load cmake@$cmake_vnew%gcc@$gcc_vnew
spack load -r openmpi@$ompi_vnew%gcc@$gcc_vnew
spack load fftw@$fftw_vnew%gcc@$gcc_vnew

# GCC_TOOLCHAIN_PATH=$(spack find -p gcc@$gcc_vnew | awk -v gcc_version="gcc@$gcc_vnew" -e '{ if ($1 == gcc_version) { print $2 }}')
# CXXFLAGS="--gcc-toolchain=${GCC_TOOLCHAIN_PATH} -stdlib=libstdc++"
# CFLAGS="--gcc-toolchain=${GCC_TOOLCHAIN_PATH}"
# LDFLAGS="-L${GCC_TOOLCHAIN_PATH}/lib64 -Wl,-rpath,${GCC_TOOLCHAIN_PATH}/lib64"

QMC_IMMUTABLE_FLAGS="-DBUILD_AFQMC=1"
