#!/bin/bash

# GCC
# Dates at https://gcc.gnu.org/releases.html
gcc_vnew=9.2.0 # 2019-08-12
gcc_vold=7.3.0 # 2018-01-25

#For Intel:
gcc_vintel=7.4.0 # 2018-12-06

#PGI 19.4
# makelocalrc configured with 8.3.0 currently
gcc_vpgi=8.3.0 # 2019-02-22

# For CUDA toolkit compatibility
gcc_vcuda=8.3.0 #  2019-02-22

# LLVM 
# Dates at http://releases.llvm.org/
llvm_vnew=9.0.0 # 2019-09-19
llvm_vold=5.0.1 # 2017-12-21
# for CUDA 10.1 update 2
llvm_vcuda=8.0.0 # 2019-03-

# HDF5
hdf5_vnew=1.10.5 # Releeased 2019-02-28
hdf5_vold=1.8.19 # Released 2017-06-16

# CMake 
# Dates at https://cmake.org/files/
cmake_vnew=3.16.2 # Released 2019-12-19
cmake_vold=3.10.2 # Released 2018-01-18

# OpenMPI
# Dates at https://www.open-mpi.org/software/ompi/v4.0/
ompi_vnew=4.0.2 # Released 2019-10-07
ompi_vold=2.1.2 # Released 2017-09-20

libxml2_vnew=2.9.9 # Released 2019-01-03 See http://xmlsoft.org/sources/
libxml2_vold=2.9.1 # Released 2013-04-19

# FFTW
# Dates at http://www.fftw.org/release-notes.html
fftw_vnew=3.3.8 # Released 2018-05-28
fftw_vold=3.3.4 # Released 2014-03-16

# BOOST
# Dates at https://www.boost.org/users/history/
boost_vnew=1.70.0 # Released 2019-04-12
boost_vold=1.65.1 # Released 2016-05-13
