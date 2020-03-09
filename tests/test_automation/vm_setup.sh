#!/bin/bash

#
# Installs compilers & libraries for QMCPACK testing via SPACK 
#
#          *** DESTROYS EXISTING INSTALLATION ***
#
#

rm -r -f $HOME/apps/spack $HOME/.spack

mkdir $HOME/.spack

cat >$HOME/.spack/config.yaml<<EOF
config:

  build_stage:
    - /scratch/$USER/spack_build_stage
    - /home/$USER/apps/spack/var/spack/stage

  build_jobs: 12

EOF
# Use flash /scratch for builds
rm -r -f /scratch/$USER/spack_build_stage
mkdir /scratch/$USER/spack_build_stage

cd $HOME/apps
git clone https://github.com/spack/spack.git
cd spack
# For reproducibility, use a specific version of Spack
# Use tagged releases https://github.com/spack/spack/releases
# git checkout v0.13.3
git checkout b9dc263801ab8b9ce46e83adec8002c299fe2e44
#Author: Justin S <3630356+codeandkey@users.noreply.github.com>
#Date:   Fri Jan 3 15:52:59 2020 -0600
#
#    py-intervaltree: new package at 3.0.2 (#14277)

module() { eval `/usr/bin/modulecmd bash $*`; }

cd bin
./spack bootstrap

export SPACK_ROOT=$HOME/apps/spack
. $SPACK_ROOT/share/spack/setup-env.sh

echo --- Spack list
spack find
echo --- Spack compilers
spack compilers
echo --- Modules list
module list
echo --- End listings

#
# Versions should be consistent with setup script
#

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

echo --- START env `date`
echo --- gcc@${gcc_vold}
spack install gcc@${gcc_vold}
spack load gcc@${gcc_vold}
spack compiler find
echo --- Spack compilers
spack compilers
spack install --no-checksum libxml2@${libxml2_vold}%gcc@${gcc_vold}
spack install cmake@${cmake_vold}%gcc@${gcc_vold}
spack install boost@${boost_vold}%gcc@${gcc_vold}
spack install openmpi@${ompi_vold}%gcc@${gcc_vold}
#spack install hdf5@${hdf5_vold}%gcc@${gcc_vold}^openmpi@${ompi_vold}%gcc@${gcc_vold}
spack install hdf5@${hdf5_vold}~mpi %gcc@${gcc_vold}
spack install fftw@${fftw_vold}%gcc@${gcc_vold}
spack unload gcc@${gcc_vold}
echo --- gcc@${gcc_vnew}
spack install gcc@${gcc_vnew}
spack load gcc@${gcc_vnew}
spack compiler add
spack install --no-checksum libxml2@${libxml2_vnew}%gcc@${gcc_vnew}
#spack install --no-checksum cmake@${cmake_vnew}%gcc@${gcc_vnew}
# BUILD LATEST CMAKE WITH OLD GCC DUE TO BUILD FAILURE WITH GCC_VNEW
spack install --no-checksum cmake@${cmake_vnew}%gcc@${gcc_vold}
spack install boost@${boost_vnew}%gcc@${gcc_vnew}
spack install openmpi@${ompi_vnew}%gcc@${gcc_vnew}^libxml2@${libxml2_vnew}%gcc@${gcc_vnew}
#spack install hdf5@${hdf5_vnew}%gcc@${gcc_vnew}^openmpi@${ompi_vnew}%gcc@${gcc_vnew}
spack install hdf5@${hdf5_vnew}~mpi %gcc@${gcc_vnew}
spack install fftw@${fftw_vnew}%gcc@${gcc_vnew}
spack unload gcc@${gcc_vnew}
echo --- llvm@${llvm_vnew}
spack install llvm@${llvm_vnew}
spack load llvm@${llvm_vnew}
spack compiler add
spack unload llvm@${llvm_vnew}
echo --- llvm@${llvm_vold}
spack install llvm@${llvm_vold}%gcc@${gcc_vold}
spack load llvm@${llvm_vold}%gcc@$gcc_vold
spack compiler add
spack unload llvm@${llvm_vold}%gcc@$gcc_vold
echo --- llvm@${llvm_vcuda}
spack install llvm@${llvm_vcuda}
spack load llvm@${llvm_vcuda}
spack compiler add
spack unload llvm@${llvm_vcuda}
echo --- gcc@${gcc_vintel}
spack install gcc@${gcc_vintel}
spack load gcc@${gcc_vintel}
spack compiler add
spack unload gcc@${gcc_vintel}
echo --- gcc@${gcc_vpgi}
spack install gcc@${gcc_vpgi}
spack load gcc@${gcc_vpgi}
spack compiler add
spack unload gcc@${gcc_vpgi}
echo --- gcc@${gcc_vcuda}
spack install gcc@${gcc_vcuda}
spack load gcc@${gcc_vcuda}
spack compiler add
spack unload gcc@${gcc_vcuda}
echo --- Convenience
spack install git
echo --- Python setup for NEXUS
spack install py-numpy@1.18.0%gcc@9.2.0
spack install py-h5py^py-numpy@1.18.0%gcc@9.2.0
spack install py-pandas@0.25.1^py-numpy@1.18.0%gcc@9.2.0
spack install py-scipy@1.4.1^py-numpy@1.18.0%gcc@9.2.0
spack activate py-numpy@1.18.0%gcc@9.2.0
spack activate py-h5py^py-numpy@1.18.0%gcc@9.2.0
spack activate py-pandas@0.25.1^py-numpy@1.18.0%gcc@9.2.0
spack activate py-scipy@1.4.1^py-numpy@1.18.0%gcc@9.2.0
# Specify py-numpy@1.16.4 to avoid python3 dependencies in later versions
#spack install py-numpy@1.16.4
#spack install py-h5py^py-numpy@1.16.4
#spack install py-pandas@0.24.2^py-numpy@1.16.4
#spack activate py-numpy@1.16.4
#spack activate py-h5py^py-numpy@1.16.4
#spack activate py-pandas@0.24.2^py-numpy@1.16.4

echo --- PGI setup reminder
echo "To configure the PGI compilers with one of the newly installed C++ libraries:"
echo "spack load gcc@8.2.0 # For example"
echo "cd /opt/pgi/linux86-64/19.4/bin"
echo "sudo ./makelocalrc -x /opt/pgi/linux86-64/19.4/ -gcc `which gcc` -gpp `which g++` -g77 `which gfortran`"
echo "gcc_vpgi is set to" $gcc_vpgi
echo --- FINISH `date`
