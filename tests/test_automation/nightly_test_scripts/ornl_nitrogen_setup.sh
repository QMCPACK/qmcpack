#!/bin/bash

#
# Installs compilers & libraries for QMCPACK testing via SPACK 
#
#          *** DESTROYS EXISTING INSTALLATION ***
#
#

echo --- START setup script `date`

rm -r -f $HOME/apps/spack $HOME/.spack

mkdir $HOME/.spack

cat >$HOME/.spack/config.yaml<<EOF
config:

  build_stage:
    - /scratch/$USER/spack_build_stage
    - /home/$USER/apps/spack/var/spack/stage

  build_jobs: 128

EOF
# Use flash /scratch for builds
rm -r -f /scratch/$USER/spack_build_stage
mkdir /scratch/$USER/spack_build_stage

cd $HOME/apps
git clone https://github.com/spack/spack.git

cd $HOME/apps/spack
# For reproducibility, use a specific version of Spack
# Use tagged releases https://github.com/spack/spack/releases
# git checkout v0.13.3
#git checkout b9dc263801ab8b9ce46e83adec8002c299fe2e44
#Author: Justin S <3630356+codeandkey@users.noreply.github.com>
#Date:   Fri Jan 3 15:52:59 2020 -0600
#
#    py-intervaltree: new package at 3.0.2 (#14277)
#git checkout v0.13.4
git log -1

module() { eval `/usr/bin/modulecmd bash $*`; }

cd bin
./spack bootstrap

export SPACK_ROOT=$HOME/apps/spack
. $SPACK_ROOT/share/spack/setup-env.sh

echo --- Spack list
spack find
echo --- Spack compilers
spack compilers
echo --- Spack compiler add
spack compiler find
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
gcc_vnew=9.3.0 # 2020-03-12

#Zen2 optimziations are only in gcc 9.1+, with improved scheduling in 9.2+
#For now, only use newer compilers

# For CUDA toolkit compatibility
gcc_vcuda=8.3.0 #  2019-02-22

# LLVM 
# Dates at http://releases.llvm.org/
llvm_vnew=10.0.0 # 2020-03-24
#Zen2 scheduling optimizations are only in LLVM 10+

# HDF5
hdf5_vnew=1.10.5 # Releeased 2019-02-28
hdf5_vold=1.8.19 # Released 2017-06-16

# CMake 
# Dates at https://cmake.org/files/
cmake_vnew=3.16.5 # Released 2020-03-04
cmake_vold=3.10.2 # Released 2018-01-18

# OpenMPI
# Dates at https://www.open-mpi.org/software/ompi/v4.0/
ompi_vnew=4.0.3 # Released 2019-10-07
#ompi_vold=2.1.2 # Released 2017-09-20

libxml2_vnew=2.9.9 # Released 2019-01-03 See http://xmlsoft.org/sources/
libxml2_vold=2.9.1 # Released 2013-04-19

# FFTW
# Dates at http://www.fftw.org/release-notes.html
fftw_vnew=3.3.8 # Released 2018-05-28
#fftw_vold=3.3.4 # Released 2014-03-16

# BOOST
# Dates at https://www.boost.org/users/history/
boost_vnew=1.72.0 # Released 2019-12-11
boost_vold=1.67.0 # Released 2018-04-14

echo --- START env `date`
echo --- gcc@${gcc_vnew}
spack install gcc@${gcc_vnew}
echo --- load gcc@${gcc_vnew}
spack load gcc@${gcc_vnew}
module list
spack compiler find
echo --- Convenience
spack install git
spack install python # Needed by libflame, good safety measure
spack load python
echo --- gcc@${gcc_vnew} consumers
spack install libxml2@${libxml2_vnew}%gcc@${gcc_vnew}
spack install cmake@${cmake_vnew}%gcc@${gcc_vnew}
spack install boost@${boost_vnew}%gcc@${gcc_vnew}
spack install openmpi@${ompi_vnew}%gcc@${gcc_vnew} ^libxml2@${libxml2_vnew}%gcc@${gcc_vnew}
#spack HDF5 package requires fortran and hl (high level) support to be specifically enabled for use with QE
#spack install hdf5@${hdf5_vnew}%gcc@${gcc_vnew} +fortran +hl ^openmpi@${ompi_vnew}%gcc@${gcc_vnew} fabrics=auto
spack install hdf5@${hdf5_vnew}%gcc@${gcc_vnew} +fortran +hl -mpi #Avoid MPI otherwise nompi build can break
#spack install fftw@${fftw_vnew}%gcc@${gcc_vnew} ^openmpi@${ompi_vnew}%gcc@${gcc_vnew} fabrics=auto
spack install fftw@${fftw_vnew}%gcc@${gcc_vnew} -mpi
spack unload gcc@${gcc_vnew}
echo --- gcc@${gcc_vcuda}
spack install gcc@${gcc_vcuda}
spack load gcc@${gcc_vcuda}
spack compiler find
spack unload gcc@${gcc_vcuda}
echo --- llvm@${llvm_vnew}
spack install llvm@${llvm_vnew}
spack load llvm@${llvm_vnew}
spack compiler add
spack unload llvm@${llvm_vnew}
echo --- Zen2 BLAS+LAPACK
# Additional development required here due to only partial AMD Rome coverage and
# spack support for optimized BLAS and LAPACK

spack install blis%gcc@${gcc_vnew}
spack install amdblis%gcc@${gcc_vnew}
# Spack has no amdlibflame package for LAPACK
# libflame builds fail 20200326 with SHA-1 collisions
#spack install libflame%gcc@${gcc_vnew} ^blis%gcc@${gcc_vnew} # Needs env python to work to build. Python is loaded above.
spack install netlib-lapack%gcc@${gcc_vnew} # Netlib failback

#spack install openblas%gcc@${gcc_vnew} # Will crash in Performance tests

echo --- Convenience
spack install git
#echo --- Python setup for NEXUS `date`
#spack install py-numpy^blis%gcc@${gcc_vnew} # will pull in libflame (problems 2020-03-27)
spack install py-numpy # Will pull in OpenBLAS
spack install py-scipy
spack install py-mpi4py^openmpi@${ompi_vnew}%gcc@${gcc_vnew} 
#spack install py-h5py^openmpi@${ompi_vnew}%gcc@${gcc_vnew}
spack install py-h5py -mpi ^hdf5@${hdf5_vnew}%gcc@${gcc_vnew} +fortran +hl -mpi #Avoid MPI otherwise nompi build can break
spack install py-pandas
spack activate py-numpy
spack activate py-scipy
spack activate py-h5py
spack activate py-pandas


echo --- FINISH setup script `date`
