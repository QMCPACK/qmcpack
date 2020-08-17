#!/bin/bash

#
# Installs compilers & libraries for QMCPACK testing via SPACK 
#
#          *** DESTROYS EXISTING INSTALLATION ***
#
#

echo --- START setup script `date`

if [ -e `dirname "$0"`/ornl_versions.sh ]; then
    source `dirname "$0"`/ornl_versions.sh
else
    echo Did not find version numbers script ornl_versions.sh
    exit 1
fi

plat=`lscpu|grep Vendor|sed 's/.*ID:[ ]*//'`
case "$plat" in
    GenuineIntel )
	ourplatform=Intel
	;;
    AuthenticAMD )
	ourplatform=AMD
	;;
    * )
	# ARM support should be trivial, but is not yet done
	echo Unknown platform
	exit 1
	;;
esac

echo --- Installing for $ourplatform architecture
ourhostname=`hostname|sed 's/\..*//g'`
echo --- Host is $ourhostname

if [ -e $HOME/apps/spack ]; then
    rm -r -f $HOME/apps/spack
fi
if [ -e $HOME/.spack ]; then
    rm -r -f $HOME/.spack
fi
mkdir $HOME/.spack

# Setup build multiplicity and preferred directories for spack
# Choose the fastest filesytem. Don't abuse shared nodes.
case "$ourhostname" in
    nitrogen )
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
	;;
    sulfur )
	cat >$HOME/.spack/config.yaml<<EOF
config:

  build_stage:
    - /scratch/$USER/spack_build_stage
    - /home/$USER/apps/spack/var/spack/stage

  build_jobs: 96

EOF
# Use flash /scratch for builds
	rm -r -f /scratch/$USER/spack_build_stage
	mkdir /scratch/$USER/spack_build_stage
	# Workaround linux-rhel8-cascadelake problem on sulfur for GCC in spack circa 202006
#	cat >$HOME/.spack/packages.yaml<<EOF
#packages:
#  all:
#    target: [skylake_avx512]
#EOF
	cat >$HOME/.spack/packages.yaml<<EOF
packages:
  all:
    target: [x86_64]
EOF
	;;
    *)
	echo "*** WARNING: Unknown host in initial ourhostname case statement. No custom onfiguration"
	;;
esac


if [ ! -e $HOME/apps ]; then
mkdir $HOME/apps
fi

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

export SPACK_ROOT=$HOME/apps/spack
export PATH=$SPACK_ROOT/bin:$PATH
. $SPACK_ROOT/share/spack/setup-env.sh

# IMPORTANT: Install+Use a GCC toolset on Red Hat systems to use
# recent compilers with better architecture support.
# e.g. yum install gcc-toolset-9
if [ -e /opt/rh/gcc-toolset-9/root/bin/gcc ]; then
    export PATH=/opt/rh/gcc-toolset-9/root/bin/:$PATH
fi


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


echo --- START env `date`
echo --- gcc@${gcc_vnew}
spack install gcc@${gcc_vnew}
echo --- load gcc@${gcc_vnew}
spack load gcc@${gcc_vnew}
module list
spack compiler find
echo --- Convenience
spack install git%gcc@${gcc_vnew}
spack install python%gcc@${gcc_vnew} # Needed by libflame, good safety measure
spack load python%gcc@${gcc_vnew}
echo --- gcc@${gcc_vnew} consumers
spack install libxml2@${libxml2_vnew}%gcc@${gcc_vnew}
spack install cmake@${cmake_vnew}%gcc@${gcc_vnew}
spack install boost@${boost_vnew}%gcc@${gcc_vnew}
spack install openmpi@${ompi_vnew}%gcc@${gcc_vnew} ^libxml2@${libxml2_vnew}%gcc@${gcc_vnew}
#spack HDF5 package requires fortran and hl (high level) support to be specifically enabled for use with QE
spack install hdf5@${hdf5_vnew}%gcc@${gcc_vnew} +fortran +hl -mpi #Avoid MPI otherwise nompi build can break
spack install fftw@${fftw_vnew}%gcc@${gcc_vnew} -mpi #Avoid MPI for simplicity
spack unload gcc@${gcc_vnew}
echo --- gcc@${gcc_vold}
spack install gcc@${gcc_vold}
echo --- load gcc@${gcc_vold}
spack load gcc@${gcc_vold}
module list
spack compiler find
echo --- gcc@${gcc_vold} consumers
spack install --no-checksum libxml2@${libxml2_vold}%gcc@${gcc_vold}
spack install cmake@${cmake_vold}%gcc@${gcc_vold}
spack install boost@${boost_vold}%gcc@${gcc_vold}
spack install openmpi@${ompi_vold}%gcc@${gcc_vold} ^libxml2@${libxml2_vold}%gcc@${gcc_vold}
#spack HDF5 package requires fortran and hl (high level) support to be specifically enabled for use with QE
spack install hdf5@${hdf5_vold}%gcc@${gcc_vold} +fortran +hl -mpi #Avoid MPI otherwise nompi build can break
spack install fftw@${fftw_vold}%gcc@${gcc_vold} -mpi #Avoid MPI for simplicity
spack install python%gcc@${gcc_vold}
echo --- gcc@${gcc_vcuda}
spack install gcc@${gcc_vcuda}
spack load gcc@${gcc_vcuda}
spack compiler find
spack unload gcc@${gcc_vcuda}
if [ "$ourplatform" == "Intel" ]; then
echo --- gcc@${gcc_vintel}
spack install gcc@${gcc_vintel}
spack load gcc@${gcc_vintel}
spack compiler find
spack unload gcc@${gcc_vintel}
fi
echo --- gcc@${gcc_vpgi}
spack install gcc@${gcc_vpgi}
spack load gcc@${gcc_vpgi}
spack compiler find
spack unload gcc@${gcc_vpgi}
echo --- llvm@${llvm_vnew}
spack install llvm@${llvm_vnew}
spack load llvm@${llvm_vnew}
spack compiler find
spack unload llvm@${llvm_vnew}
echo --- llvm@${llvm_vold}
spack install llvm@${llvm_vold}%gcc@${gcc_vold}
spack load llvm@${llvm_vold}%gcc@$gcc_vold
spack compiler find
spack unload llvm@${llvm_vold}%gcc@$gcc_vold
echo --- llvm@${llvm_vcuda}
spack install llvm@${llvm_vcuda}
spack load llvm@${llvm_vcuda}
spack compiler find
spack unload llvm@${llvm_vcuda}
echo --- BLAS+LAPACK
# Additional development required here due to only partial AMD Rome coverage and
# spack support for optimized BLAS and LAPACK, problems building libflame etc.
spack install blis%gcc@${gcc_vnew}
if [ "$ourplatform" == "AMD" ]; then
    spack install amdblis%gcc@${gcc_vnew}
fi

# Spack has no amdlibflame package for LAPACK
# libflame builds fail 20200326 with SHA-1 collisions
#spack install libflame%gcc@${gcc_vnew} ^blis%gcc@${gcc_vnew} # Needs env python to work to build. Python is loaded above.
spack install netlib-lapack%gcc@${gcc_vnew} # Netlib failback
#spack install openblas%gcc@${gcc_vnew} # Historically crashed in Performance tests with large thread counts on AMD

echo --- Python setup for NEXUS `date`
echo --- New python modules
#spack install py-numpy^blis%gcc@${gcc_vnew} # will pull in libflame (problems 2020-03-27)
spack install py-numpy%gcc@${gcc_vnew} # Will pull in OpenBLAS
spack install py-scipy%gcc@${gcc_vnew} 
spack install py-mpi4py%gcc@${gcc_vnew} ^openmpi@${ompi_vnew}%gcc@${gcc_vnew} ^libxml2@${libxml2_vnew}%gcc@${gcc_vnew}
#spack install py-mpi4py^openmpi@${ompi_vnew}%gcc@${gcc_vnew} 
#spack install py-h5py^openmpi@${ompi_vnew}%gcc@${gcc_vnew}
spack install py-h5py%gcc@${gcc_vnew}  -mpi ^hdf5@${hdf5_vnew}%gcc@${gcc_vnew} +fortran +hl -mpi #Avoid MPI otherwise nompi build can break
spack install py-pandas%gcc@${gcc_vnew} 
spack install py-lxml%gcc@${gcc_vnew} 
spack activate py-numpy%gcc@${gcc_vnew} 
spack activate py-scipy%gcc@${gcc_vnew} 
spack activate py-h5py%gcc@${gcc_vnew} 
spack activate py-pandas%gcc@${gcc_vnew} 
spack activate py-lxml%gcc@${gcc_vnew} 
spack unload python
echo --- Old python modules
spack load python%gcc@${gcc_vold} 
spack install py-numpy%gcc@${gcc_vold} # Will pull in OpenBLAS
spack install py-scipy%gcc@${gcc_vold} 
spack install py-mpi4py%gcc@${gcc_vold} ^openmpi@${ompi_vold}%gcc@${gcc_vold} ^libxml2@${libxml2_vold}%gcc@${gcc_vold}
#spack install py-mpi4py^openmpi@${ompi_vold}%gcc@${gcc_vold} 
#spack install py-h5py^openmpi@${ompi_vold}%gcc@${gcc_vold}
spack install py-h5py%gcc@${gcc_vold}  -mpi ^hdf5@${hdf5_vold}%gcc@${gcc_vold} +fortran +hl -mpi #Avoid MPI otherwise nompi build can break
spack install py-pandas%gcc@${gcc_vold} 
spack install py-lxml%gcc@${gcc_vold} 
spack activate py-numpy%gcc@${gcc_vold} 
spack activate py-scipy%gcc@${gcc_vold} 
spack activate py-h5py%gcc@${gcc_vold} 
spack activate py-pandas%gcc@${gcc_vold} 
spack activate py-lxml%gcc@${gcc_vold} 
spack unload python
echo --- PGI setup reminder
echo "To configure the PGI compilers with one of the newly installed C++ libraries:"
echo "spack load gcc@${gcc_vpgi} # For example"
echo "cd /opt/pgi/linux86-64/19.10/bin"
echo "sudo ./makelocalrc -x /opt/pgi/linux86-64/19.10/ -gcc `which gcc` -gpp `which g++` -g77 `which gfortran`"
echo "gcc_vpgi is set to" ${gcc_vpgi}
echo --- FINISH setup script `date`
