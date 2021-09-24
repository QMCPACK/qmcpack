#
# Versions should be consistent with setup script
#

# GCC
# Dates at https://gcc.gnu.org/releases.html
gcc_vnew=11.2.0 # Released 2021-07-28
#gcc_vnew=10.3.0 # Released 2021-04-08
gcc_vold=9.1.0  # Released 2019-05-03

gcc_vcuda=9.1.0  # Released 2019-05-03 https://docs.nvidia.com/hpc-sdk/hpc-sdk-release-notes/index.html
gcc_vintel=9.1.0 # Released 2019-05-03 Compiler for C++ library used by Intel compiler
gcc_vpgi=9.1.0   # Released 2019-05-03 Use makelocalrc to configure PGI with this compiler 

# LLVM 
# Dates at http://releases.llvm.org/
llvm_vnew=12.0.1 # Released 2021-07-08

# HDF5
# Dates at https://portal.hdfgroup.org/display/support/Downloads
hdf5_vnew=1.12.1 # Released 2020-07-01
hdf5_vold=1.8.19  # Released 2017-06-16

# CMake 
# Dates at https://cmake.org/files/
cmake_vnew=3.21.1 # Released 2021-07-27
cmake_vold=3.15.3 # Released 2019-09-04

# OpenMPI
# Dates at https://www.open-mpi.org/software/ompi/v4.1/
ompi_vnew=4.1.1 # Released 2021-04-24
ompi_vold=3.1.6 # Released 2020-03-18

# Libxml2
#libxml2_v=2.9.12 # Released 2021-05-13 See http://xmlsoft.org/sources/
libxml2_v=2.9.10 # Released 2019-10-30 See http://xmlsoft.org/sources/

# FFTW
# Dates at http://www.fftw.org/release-notes.html
fftw_vnew=3.3.9 # Released 2020-12-13
fftw_vold=3.3.8 # Released 2018-05-28

# BOOST
# Dates at https://www.boost.org/users/history/
#boost_vnew=1.77.0 # Released 2021-08-11
boost_vnew=1.76.0 # Released 2021-04-16
boost_vold=1.68.0 # Released 2018-08-09

# Python
# Use a single version to reduce dependencies. Ideally the spack prefered version.
python_version=3.8.11

