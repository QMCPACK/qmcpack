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
git checkout 2e6b44328a54293cebc5e871523c974b8be5a1fb
# This is a few commits after gcc@9.2.0 was included
#commit 2e6b44328a54293cebc5e871523c974b8be5a1fb
#Author: Adam J. Stewart <ajstewart426@gmail.com>
#Date:   Tue Aug 13 21:09:25 2019 -0500
#
#    zstd: fix build linking error (#12413)

cd bin
./spack bootstrap

module() { eval `/usr/bin/modulecmd bash $*`; }

export SPACK_ROOT=$HOME/apps/spack
. $SPACK_ROOT/share/spack/setup-env.sh

echo --- Spack list
spack find
echo --- Modules list
module list
echo --- End listings

#
# Versions should be consistent with setup script
#
boost_vnew=1.70.0 # Released 2019-04-12
boost_vold=1.61.0 # Released 2016-05-13

gcc_vnew=8.3.0 # Released 2019-02-22 See https://www.gnu.org/software/gcc/releases.html
gcc_vold=5.5.0 # Released 2017-10-27
gcc_vintel=7.4.0 # Released 2018-12-06

hdf5_vnew=1.10.5 # Releeased 2019-02-28
hdf5_vold=1.8.19 # Released 2017-06-16

cmake_vnew=3.14.4 # Released 2019-05-15
cmake_vold=3.8.2 # Released 2017-05-31

ompi_vnew=4.0.1 # Released 2019-03-26
ompi_vold=2.1.1 # Released 2017-05-10

libxml2_vnew=2.9.9 # Released 2019-01-03 See http://xmlsoft.org/sources/
libxml2_vold=2.9.1  #Released 2013-04-19

fftw_vnew=3.3.8 # Released 2018-05-28
fftw_vold=3.3.4 # Released 2014-03-16

llvm_vnew=7.0.1 # Released 2018-12-21
llvm_vold=4.0.1 # Released 2017-07-04
llvm_vcuda=6.0.1 # Released 2018-07-05

echo --- START env `date`
echo --- gcc@${gcc_vnew}
spack install gcc@${gcc_vnew}
spack load gcc@${gcc_vnew}
spack compiler add
spack install --no-checksum libxml2@${libxml2_vnew}%gcc@${gcc_vnew}
spack install cmake@${cmake_vnew}%gcc@${gcc_vnew}
spack install boost@${boost_vnew}%gcc@${gcc_vnew}
spack install openmpi@${ompi_vnew}%gcc@${gcc_vnew}^libxml2@${libxml2_vnew}%gcc@${gcc_vnew}
#spack install hdf5@${hdf5_vnew}%gcc@${gcc_vnew}^openmpi@${ompi_vnew}%gcc@${gcc_vnew}
spack install hdf5@${hdf5_vnew}~mpi %gcc@${gcc_vnew}
spack install fftw@${fftw_vnew}%gcc@${gcc_vnew}
spack unload gcc@${gcc_vnew}
echo --- gcc@${gcc_vold}
spack install gcc@${gcc_vold}
spack load gcc@${gcc_vold}
spack compiler add
spack install --no-checksum libxml2@${libxml2_vold}%gcc@${gcc_vold}
spack install cmake@${cmake_vold}%gcc@${gcc_vold}
spack install boost@${boost_vold}%gcc@${gcc_vold}
spack install openmpi@${ompi_vold}%gcc@${gcc_vold}
#spack install hdf5@${hdf5_vold}%gcc@${gcc_vold}^openmpi@${ompi_vold}%gcc@${gcc_vold}
spack install hdf5@${hdf5_vold}~mpi %gcc@${gcc_vold}
spack install fftw@${fftw_vold}%gcc@${gcc_vold}
spack unload gcc@${gcc_vold}
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
echo --- Convenience
spack install git
echo --- Python setup for NEXUS
# Specify py-numpy@1.16.4 to avoid python3 dependencies in later versions
spack install py-numpy@1.16.4
spack install py-h5py^py-numpy@1.16.4
spack install py-pandas@0.24.2^py-numpy@1.16.4
spack activate py-numpy@1.16.4
spack activate py-h5py^py-numpy@1.16.4
spack activate py-pandas@0.24.2^py-numpy@1.16.4

echo --- PGI setup reminder
echo "To configure the PGI compilers with one of the newly installed C++ libraries:"
echo "spack load gcc@8.2.0 # For example"
echo "cd /opt/pgi/linux86-64/19.4/bin"
echo "sudo ./makelocalrc -x /opt/pgi/linux86-64/19.4/ -gcc `which gcc` -gpp `which g++` -g77 `which gfortran`"
echo --- FINISH `date`
