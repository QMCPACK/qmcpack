#
# Versions should be consistent with setup script
#

# GCC
# Dates at https://gcc.gnu.org/releases.html
gcc_vnew=10.2.0 # Released 2020-07-23
#gcc_vold=8.2.0  # Released 2018-07-26
gcc_vold=8.3.0  # Released 2019-02-22 (8.2.0 results in too many install issues 20201202)

gcc_vcuda=9.3.0  # Released 2020-03-12 9.x For CUDA 11.1 compatibility https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html
gcc_vintel=8.3.0 # Released 2019-02-22 Compiler for C++ library used by Intel compiler
gcc_vpgi=8.3.0   # Released 2019-02-22 Use makelocalrc to configure PGI with this compiler

# LLVM 
# Dates at http://releases.llvm.org/
#llvm_vnew=11.0.0 # Released 2020-10-12
llvm_vnew=10.0.1 # Released 2020-08-06
#llvm_vold=7.0.1  # Released 2018-12-21 (disabled due to SPACK install issues 20201202)
#llvm_vcuda=9.0.0 # Released 2019-09-19 9.0.0 For CUDA 11.1 compatibility https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html
llvm_vcuda=9.0.1 # Released 2019-12-20 9.0.1 (because 9.0.0 will not install 20201203)

# HDF5
# Dates at https://portal.hdfgroup.org/display/support/Downloads
#hdf5_vnew=1.12.0 # Released 2020-02-28 . Not supported by py-h5py 2.10.0, needs v3+
hdf5_vnew=1.10.7  # Released 2020-09-15
hdf5_vold=1.8.19  # Released 2017-06-16

# CMake 
# Dates at https://cmake.org/files/
cmake_vnew=3.19.2 # Released 2020-12-16
cmake_vold=3.13.2 # Released 2018-12-13

# OpenMPI
# Dates at https://www.open-mpi.org/software/ompi/v4.1/
ompi_vnew=4.1.0 # Released 2020-12-18
ompi_vold=3.1.6 # Released 2020-03-18

# Libxml2
libxml2_v=2.9.10 # Released 2019-10-30 See http://xmlsoft.org/sources/

# FFTW
# Dates at http://www.fftw.org/release-notes.html
fftw_vnew=3.3.8 # Released 2018-05-28
fftw_vold=3.3.4 # Released 2014-03-16

# BOOST
# Dates at https://www.boost.org/users/history/
boost_vnew=1.75.0 # Released 2020-12-11
boost_vold=1.68.0 # Released 2018-08-09

# Python
# Use a single version to reduce dependencies
python_version=3.8.6

