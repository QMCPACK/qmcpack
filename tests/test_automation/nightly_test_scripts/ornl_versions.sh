#
# Versions should be consistent with setup script
#

# AMD Zen2 (Rome) well supported from gcc 9.2+, llvm 10+

# GCC
# Zen2 optimziations are only in gcc 9.1+, with improved scheduling in 9.2+
# Dates at https://gcc.gnu.org/releases.html
gcc_vnew=10.2.0 # 2020-07-23
gcc_vold=7.3.0 # 2018-01-25

#gcc_vcuda=8.2.1 #  For CUDA 10.2 compatibility  https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html
gcc_vcuda=8.3.0 #  Use 8.3.0 instead since available in spack
gcc_vintel=8.3.0 # Compiler for C++ library used by Intel compiler
gcc_vpgi=8.3.0 # Use makelocalrc to configure PGI with this compiler

# LLVM 
# Zen2 scheduling optimizations are only in LLVM 10+
# Dates at http://releases.llvm.org/
llvm_vnew=10.0.1 # 2020-08-06
llvm_vold=6.0.1 # 2018-07-05
llvm_vcuda=8.0.0 # For CUDA 10.2 compatibility https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html

# HDF5
# Dates at https://portal.hdfgroup.org/display/support/Downloads
#hdf5_vnew=1.12.0 # Released 2020-02-28 . Not supported by py-h5py 2.10.0, latest in spack on 20200604
hdf5_vnew=1.10.5 # Releeased 2019-02-28
hdf5_vold=1.8.19 # Released 2017-06-16

# CMake 
# Dates at https://cmake.org/files/
cmake_vnew=3.18.2 # Released 2020-08-20
cmake_vold=3.12.1 # Released 2018-08-09

# OpenMPI
# Dates at https://www.open-mpi.org/software/ompi/v4.0/
ompi_vnew=4.0.4 # Released 2020-06-10
ompi_vold=3.1.2 # Released 2018-08-22

# Libxml2
libxml2_vnew=2.9.10 # Released 2019-10-30 See http://xmlsoft.org/sources/
libxml2_vold=2.9.1 # Released 2013-04-19

# FFTW
# Dates at http://www.fftw.org/release-notes.html
fftw_vnew=3.3.8 # Released 2018-05-28
fftw_vold=3.3.4 # Released 2014-03-16

# BOOST
# Dates at https://www.boost.org/users/history/
boost_vnew=1.74.0 # Released 2020-08-14
boost_vold=1.68.0 # Released 2018-08-09


# Python
# Use a single version to reduce dependencies
# Can not use spack python default 3.7.x as of 20200831 due to hidden 3.8 requirement in LLVM master
python_version=3.8.5

