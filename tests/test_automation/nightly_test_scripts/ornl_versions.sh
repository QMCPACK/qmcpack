#
# Versions should be consistent with setup script
#

# GCC
# Dates at https://gcc.gnu.org/releases.html
gcc_vnew=12.2.0 # Released 2022-08-19
gcc_vold=10.4.0 # Released 2022-06-28

#gcc_vcuda=10.2.0  # Released 2020-07-23 https://docs.nvidia.com/hpc-sdk/hpc-sdk-release-notes/index.html
#gcc_vcuda=12.1.0
#gcc_vintel=10.2.0 # Released 2020-07-23 Compiler for C++ library used by Intel compiler
gcc_vintel=10.4.0
#gcc_vnvhpc=10.2.0 # Released 2020-07-23 Use makelocalrc to configure NVHPC with this compiler 
gcc_vnvhpc=12.2.0

# LLVM 
# Dates at https://releases.llvm.org/
llvm_vnew=15.0.2 # Released 2022-10-04

# HDF5
# Dates at https://portal.hdfgroup.org/display/support/Downloads
#hdf5_vnew=1.13.0 # Released 2021-12-01 # odd versions are development versions
hdf5_vnew=1.12.2 # Released 2022-04-27
hdf5_vold=${hdf5_vnew}

# CMake 
# Dates at https://cmake.org/files/
cmake_vnew=3.24.2 # Released 2022-09-13
cmake_vold=3.18.4 # Released 2020-19-06

# OpenMPI
# Dates at https://www.open-mpi.org/software/ompi/v4.1/
ompi_vnew=4.1.4 # Released 2022-05-26
ompi_vold=3.1.6 # Released 2020-03-18

# Libxml2
libxml2_v=2.9.13 # Released 2022-02-20 See http://xmlsoft.org/sources/

# FFTW
# Dates at http://www.fftw.org/release-notes.html
fftw_vnew=3.3.10 # Released 2021-09-15
fftw_vold=3.3.8 # Released 2018-05-28

# BOOST
# Dates at https://www.boost.org/users/history/
boost_vnew=1.79.0 # Released 2022-04-13
boost_vold=1.74.0 # Released 2020-08-14

# Python
# Use a single version to reduce dependencies. Ideally the spack prefered version.
python_version=3.9.13




