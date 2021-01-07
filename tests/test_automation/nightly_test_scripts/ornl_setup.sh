#!/bin/bash

#
# Installs compilers & libraries for QMCPACK testing via SPACK 
#
#          *** DESTROYS EXISTING INSTALLATION ***
#
#

echo --- START setup script `date`

# Bug avoidance 20200902
if [ -e /usr/lib/aomp/bin/clang ]; then
    echo AOMP Clang install detected. This will break llvm install.
    echo Suggest temporarily: sudo chmod og-rx /usr/lib/aomp
    exit 1
fi


if [ -e `dirname "$0"`/ornl_versions.sh ]; then
    echo --- Contents of ornl_versions.sh
    cat `dirname "$0"`/ornl_versions.sh
    echo --- End of contents of ornl_versions.sh
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

	#Use system installed SSL. See https://spack.readthedocs.io/en/latest/getting_started.html#openssl
	cat >>$HOME/.spack/packages.yaml<<EOF
packages:
    openssl:
        externals:
        - spec: openssl@1.1.1g
          prefix: /usr
          buildable: False
EOF
#	cat >>$HOME/.spack/packages.yaml<<EOF
#packages:
#EOF
# Experiment to See if clang builds and other problems clear 20200827
	cat >>$HOME/.spack/packages.yaml<<EOF
    all:
        target: [x86_64]
EOF
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
	#Use system installed SSL. See https://spack.readthedocs.io/en/latest/getting_started.html#openssl
	#Use system RHEL gcc8 compiler for cmake to solve bootstrap problems
	cat >>$HOME/.spack/packages.yaml<<EOF
packages:
    openssl:
        externals:
        - spec: openssl@1.1.1g
          prefix: /usr
          buildable: False
EOF
# Workaround cmake bootstrap error with old versions/old spack compilers. e.g.
# spack install cmake@3.13.2%gcc@8.3.0 failure but spack install cmake@3.13.2%gcc@8.3.1 OK (RHEL compiler)	
        cat >>$HOME/.spack/packages.yaml<<EOF
    cmake:
        compiler: [gcc@8.3.1]
EOF
#    boost:
#        compiler: [gcc@8.3.1]
#    perl:
#        compiler: [gcc@8.3.1]
# Workaround linux-rhel8-cascadelake problem on sulfur for GCC in spack circa 202006
#	cat >>$HOME/.spack/packages.yaml<<EOF
#packages:
#  all:
#    target: [skylake_avx512]
#EOF

# Workaround build error of perl%gcc@installed_by_spack	. perl-data-dumper & perl are e.g. an LLVM build dep
# dnf install perl-data-dumper if missing
	cat >>$HOME/.spack/packages.yaml<<EOF
    perl:
        externals:
        - spec: perl@5.26.3
          prefix: /
          buildable: False
    perl-data-dumper:
        externals:
        - spec: perl-data-dumper@2.167.399
          prefix: /
          buildable: False
EOF
	cat >>$HOME/.spack/packages.yaml<<EOF
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
#DEBUG git clone --depth 1 https://github.com/spack/spack.git 
git clone https://github.com/spack/spack.git 

cd $HOME/apps/spack
# For reproducibility, use a specific version of Spack
# Prefer to use tagged releases https://github.com/spack/spack/releases
#git checkout 0c63c94103acae33c4f3adfe7b90d2ea1165ce46
#Author: Dennis Klein <dennis.klein.github@gmail.com>
#Date:   Fri Jul 24 19:00:55 2020 +0200
#
#    Relax architecture compatibility check (#15972)
#
#    * Relax architecture compatibility check
#    * Add test coverage for the spack.abi module
#cp -p var/spack/repos/builtin/packages/gcc/package.py var/spack/repos/builtin/packages/gcc/package.py_orig
#git checkout 5f9de810ff14455ac6b433398b965d2cf6b4a70f -- var/spack/repos/builtin/packages/gcc/package.py
#ls -l var/spack/repos/builtin/packages/gcc/package.py*

#git checkout v0.16.0 #Works? FAIL 18 Nov.
#git checkout e0ae60edc4ee7dafc5ab99d756d9c63624c281d5 # Works 1 Nov
#git checkout a5faf7d27a3ef1afca28e5d2348130f580be898d # Works 20 Oct
#git checkout 9eb87d1026819875dec8bb7fa00252576fd62006 # Works 23 Sep. 
#git checkout v0.15.4 #works 13 Aug
#git checkout v0.15.0 #works

echo --- Git version and last log entry
git log -1

module() { eval `/usr/bin/modulecmd bash $*`; }

cd bin

export SPACK_ROOT=$HOME/apps/spack
export PATH=$SPACK_ROOT/bin:$PATH
. $SPACK_ROOT/share/spack/setup-env.sh

## IMPORTANT: Install+Use a GCC toolset on Red Hat systems to use
## recent compilers with better architecture support.
## e.g. yum install gcc-toolset-9
#if [ -e /opt/rh/gcc-toolset-9/root/bin/gcc ]; then
#    echo --- Added gcc-toolset-9 to path for RHEL provided GCC9 compilers
#    export PATH=/opt/rh/gcc-toolset-9/root/bin/:$PATH
#fi

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


echo --- START `date`
echo --- diffutils iconv dependency workaround
cp ${SPACK_ROOT}/var/spack/repos/builtin/packages/diffutils/package.py  ${SPACK_ROOT}/var/spack/repos/builtin/packages/diffutils/package.py_orig
sed 's/.*depends_on.*iconv.*//g' ${SPACK_ROOT}/var/spack/repos/builtin/packages/diffutils/package.py_orig >  ${SPACK_ROOT}/var/spack/repos/builtin/packages/diffutils/package.py
#echo --- Debug
#spack install gcc@10.2.0
#spack load gcc@10.2.0
#spack compiler find
#spack install cmake%gcc@10.2.0
#spack load cmake%gcc@10.2.0
#exit
echo --- Convenience
spack install git
echo --- gcc@${gcc_vnew}
spack install gcc@${gcc_vnew}
echo --- load gcc@${gcc_vnew}
spack load gcc@${gcc_vnew}
module list
spack compiler find
echo --- gcc@${gcc_vnew} consumers
spack install cmake@${cmake_vnew}
spack install libxml2@${libxml2_v}%gcc@${gcc_vnew}
#spack install boost@${boost_vnew}%gcc@${gcc_vnew}
spack install boost@${boost_vnew}
spack install python@${python_version}%gcc@${gcc_vnew} # Needed by libflame, good safety measure
spack load python@${python_version}%gcc@${gcc_vnew}
spack install openmpi@${ompi_vnew}%gcc@${gcc_vnew} ^libxml2@${libxml2_v}%gcc@${gcc_vnew}
#spack HDF5 package requires fortran and hl (high level) support to be specifically enabled for use with QE
spack install hdf5@${hdf5_vnew}%gcc@${gcc_vnew} +fortran +hl -mpi #Avoid MPI otherwise nompi build can break
spack install fftw@${fftw_vnew}%gcc@${gcc_vnew} -mpi #Avoid MPI for simplicity
#spack install perl%gcc@${gcc_vnew} # Kludge for broken deps
spack unload python@${python_version}%gcc@${gcc_vnew}
spack unload gcc@${gcc_vnew}
echo --- gcc@${gcc_vold}
spack install gcc@${gcc_vold}
echo --- load gcc@${gcc_vold}
spack load gcc@${gcc_vold}
module list
spack compiler find
spack unload gcc@${gcc_vold}
echo --- gcc@${gcc_vold} consumers
#spack install cmake@${cmake_vold}%gcc@${gcc_vold}
spack install cmake@${cmake_vold} # Some old cmakes do not bootstrap with old gcc
spack install --no-checksum libxml2@${libxml2_v}%gcc@${gcc_vold}
spack install boost@${boost_vold}%gcc@${gcc_vold}
#spack install boost@${boost_vold}
spack install python@${python_version}%gcc@${gcc_vold}
spack load python@${python_version}%gcc@${gcc_vold}
#spack install --no-checksum libxml2@${libxml2_v}%gcc@${gcc_vold}
#spack install openmpi@${ompi_vold}%gcc@${gcc_vold} ^libxml2@${libxml2_v}%gcc@${gcc_vold}
#spack install openmpi@${ompi_vold} ^libxml2@${libxml2_v}%gcc@${gcc_vold}
spack install openmpi@${ompi_vold} ^libxml2@${libxml2_v}%gcc@${gcc_vold}
#spack HDF5 package requires fortran and hl (high level) support to be specifically enabled for use with QE
spack install hdf5@${hdf5_vold}%gcc@${gcc_vold} +fortran +hl -mpi #Avoid MPI otherwise nompi build can break
spack install fftw@${fftw_vold}%gcc@${gcc_vold} -mpi #Avoid MPI for simplicity
spack unload python@${python_version}%gcc@${gcc_vold}
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
spack install llvm@${llvm_vnew} ^python@${python_version}%gcc@${gcc_vnew}
spack load llvm@${llvm_vnew}
spack compiler find
spack unload llvm@${llvm_vnew}
#echo --- llvm@${llvm_vold}
#spack install --no-checksum llvm@${llvm_vold}%gcc@${gcc_vold} ^python@${python_version}%gcc@${gcc_vold}
#spack load llvm@${llvm_vold}%gcc@$gcc_vold
#spack compiler find
#spack unload llvm@${llvm_vold}%gcc@$gcc_vold
echo --- llvm@${llvm_vcuda}
spack install llvm@${llvm_vcuda} ^python@${python_version}%gcc@${gcc_vnew}
spack load llvm@${llvm_vcuda}
spack compiler find
spack unload llvm@${llvm_vcuda}
echo --- BLAS+LAPACK
# Additional development required here due to only partial AMD Rome coverage and
# spack support for optimized BLAS and LAPACK, problems building libflame etc.
#spack install blis%gcc@${gcc_vnew} ^python@${python_version}%gcc@${gcc_vnew}
#if [ "$ourplatform" == "AMD" ]; then
#    spack install amdblis%gcc@${gcc_vnew} ^python@${python_version}%gcc@${gcc_vnew}
#fi

## Spack has no amdlibflame package for LAPACK
## libflame builds fail 20200326 with SHA-1 collisions
##spack install libflame%gcc@${gcc_vnew} ^blis%gcc@${gcc_vnew} # Needs env python to work to build. Python is loaded above.
#spack install netlib-lapack%gcc@${gcc_vnew} # Netlib failback
##spack install openblas%gcc@${gcc_vnew} # Historically crashed in Performance tests with large thread counts on AMD

#spack install openblas%gcc@${gcc_vnew} threads=openmp
#spack install openblas%gcc@${gcc_vold} threads=openmp
spack install openblas threads=openmp

echo --- Modules installed so far
spack find
echo --- Python setup for NEXUS `date`
echo --- New python modules
#spack install py-numpy^blis%gcc@${gcc_vnew} # will pull in libflame (problems 2020-03-27)
spack load gcc@${gcc_vnew}
spack load python@${python_version}%gcc@${gcc_vnew} 
spack install py-numpy%gcc@${gcc_vnew} ^python@${python_version}%gcc@${gcc_vnew} # Will pull in OpenBLAS
spack install py-scipy%gcc@${gcc_vnew} ^python@${python_version}%gcc@${gcc_vnew}
#spack install py-mpi4py%gcc@${gcc_vnew} ^python@${python_version}%gcc@${gcc_vnew} ^openmpi@${ompi_vnew}%gcc@${gcc_vnew} ^libxml2@${libxml2_v}%gcc@${gcc_vnew} 
spack install py-setuptools%gcc@${gcc_vnew} ^python@${python_version}%gcc@${gcc_vnew}
spack install py-mpi4py%gcc@${gcc_vnew} ^openmpi@${ompi_vnew}%gcc@${gcc_vnew} ^py-setuptools%gcc@${gcc_vnew}  ^python@${python_version}%gcc@${gcc_vnew}
spack install py-h5py%gcc@${gcc_vnew}  -mpi ^python@${python_version}%gcc@${gcc_vnew} ^hdf5@${hdf5_vnew}%gcc@${gcc_vnew} +fortran +hl -mpi #Avoid MPI otherwise nompi build can break 
spack install py-pandas%gcc@${gcc_vnew} ^python@${python_version}%gcc@${gcc_vnew}
spack install py-lxml%gcc@${gcc_vnew}  
spack install py-matplotlib%gcc@${gcc_vnew} ^python@${python_version}%gcc@${gcc_vnew}
spack activate py-numpy%gcc@${gcc_vnew} 
spack activate py-scipy%gcc@${gcc_vnew} 
spack activate py-h5py%gcc@${gcc_vnew} 
spack activate py-pandas%gcc@${gcc_vnew} 
spack activate py-lxml%gcc@${gcc_vnew}
spack activate py-matplotlib%gcc@${gcc_vnew}
spack unload python@${python_version}%gcc@${gcc_vnew} 
spack unload gcc@${gcc_vnew}
echo --- Old python modules
# Due to problems installing old versions with old compilers of at least scipy (20201223) we only install numpy.
spack load gcc@${gcc_vold}
spack load python@${python_version}%gcc@${gcc_vold}
spack install py-numpy%gcc@${gcc_vold} ^python@${python_version}%gcc@${gcc_vold} # Will pull in OpenBLAS
#SKIP#BADspack install py-scipy%gcc@${gcc_vold} ^python@${python_version}%gcc@${gcc_vold} # Will pull in py-pybind11 and cmake which won't bootstrap 20201223
#SKIPspack install py-setuptools%gcc@${gcc_vold} ^python@${python_version}%gcc@${gcc_vold}
#SKIP#BAD gives dupe python spack install py-mpi4py%gcc@${gcc_vold} ^python@${python_version}%gcc@${gcc_vold} ^openmpi@${ompi_vold}%gcc@${gcc_vold} ^libxml2@${libxml2_v}%gcc@${gcc_vold}
#SKIPspack install py-mpi4py%gcc@${gcc_vold} ^openmpi@${ompi_vold}%gcc@${gcc_vold} ^py-setuptools%gcc@${gcc_vold}  ^python@${python_version}%gcc@${gcc_vold}
#SKIPspack install py-h5py%gcc@${gcc_vold}  -mpi ^python@${python_version}%gcc@${gcc_vold} ^hdf5@${hdf5_vold}%gcc@${gcc_vold} +fortran +hl -mpi #Avoid MPI otherwise nompi build can break
#SKIPspack install py-pandas%gcc@${gcc_vold} ^python@${python_version}%gcc@${gcc_vold}
#SKIPspack install py-lxml%gcc@${gcc_vold} ^python@${python_version}%gcc@${gcc_vold}
#SKIPspack install py-matplotlib%gcc@${gcc_vold} ^python@${python_version}%gcc@${gcc_vold}
spack activate py-numpy%gcc@${gcc_vold} 
#SKIPspack activate py-scipy%gcc@${gcc_vold} 
#SKIPspack activate py-h5py%gcc@${gcc_vold} 
#SKIPspack activate py-pandas%gcc@${gcc_vold} 
#SKIPspack activate py-lxml%gcc@${gcc_vold} 
#SKIPspack activate py-matplotlib%gcc@${gcc_vold}
spack unload python@${python_version}%gcc@${gcc_vold} 
spack unload gcc@${gcc_vold}
#spack load git
#DEBUG Disable cleanup based on guess there are hidden dependencies
#DEBUGecho --- Remove build dependencies
#DEBUGspack gc --yes-to-all
echo --- PGI setup reminder
echo "To configure the PGI compilers with one of the newly installed C++ libraries:"
echo "spack load gcc@${gcc_vpgi} # For example"
echo "cd /opt/pgi/linux86-64/19.10/bin"
echo "sudo ./makelocalrc -x /opt/pgi/linux86-64/19.10/ -gcc `which gcc` -gpp `which g++` -g77 `which gfortran`"
echo "gcc_vpgi is set to" ${gcc_vpgi}
# sudo ./makelocalrc -x /opt/nvidia/hpc_sdk/Linux_x86_64/20.9/compilers/ -gcc `which gcc` -gpp `which g++` -g77 `which gfortran`

# Bug avoidance 20200902
if [ -e /usr/lib/aomp/bin/clang ]; then
    echo --- AOMP Bug workaround
    echo Change permissions back on AOMP Clang install detected, e.g. sudo chmod og+rx /usr/lib/aomp
fi
echo --- FINISH setup script `date`
