#!/bin/bash

echo --- Script START `date`

localonly=no

if [[ $localonly == "yes" ]]; then
echo --- Local CMake/Make/CTest only. No cdash drop.
fi

# Weekly settings:
export GLOBALTCFG="--timeout 4800 -VV"
export LIMITEDTESTS=""
export LESSLIMITEDTESTS=""
# Nightly settings:
#export GLOBALTCFG="--timeout 360 -VV"
#export LIMITEDTESTS="-R deterministic -LE unstable -E long-"
#export LESSLIMITEDTESTS="-E long- -LE unstable"

# Directory in which to run tests. Should be an absolute path and fastest usable filesystem
test_path=/scratch/${USER}   # RAID FLASH on oxygen

test_dir=${test_path}/QMCPACK_CI_BUILDS_DO_NOT_REMOVE

export QMC_DATA=/data/pk7/QMC_DATA # Route to directory containing performance test files

# CUDA 10 setup
export CUDAVER=10.1
export PATH=/opt/local/bin:/opt/local/sbin:/usr/local/cuda-${CUDAVER}/bin/:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin
export LD_LIBRARY_PATH=/usr/local/cuda-${CUDAVER}/lib64

# Specify GPUs for testing. Obtain device IDs via "nvidia-smi -L"
#export CUDA_VISIBLE_DEVICES=

# PGI2019  setup
export PGI=/opt/pgi
export MANPATH=$MANPATH:$PGI/linux86-64/2019/man
export LM_LICENSE_FILE=$PGI/license.dat
export PATH=$PGI/linux86-64/2019/bin:$PATH

# Intel2019.1 MPI configure setting to avoid MPI crash
# via https://software.intel.com/en-us/forums/intel-clusters-and-hpc-technology/topic/799716
#export FI_PROVIDER=sockets
export I_MPI_FABRICS=shm

module() { eval `/usr/bin/modulecmd bash $*`; }

export SPACK_ROOT=$HOME/apps/spack
. $SPACK_ROOT/share/spack/setup-env.sh

echo --- Spack list
spack find
echo --- Modules list
module list
echo --- End listings

spack load git

module list
if [ -e ${test_path} ]; then

if [ ! -e ${test_dir} ]; then
mkdir ${test_dir}
fi

if [ -e ${test_dir} ]; then
cd ${test_dir}

# Minimize load on GitHub by maintaining a local cloned git used for all builds
if [ ! -e qmcpack ]; then
echo --- Cloning QMCPACK git `date`
git clone https://github.com/QMCPACK/qmcpack.git --depth 1
else
cd qmcpack
echo --- Updating local QMCPACK git `date`
git pull
cd ..
fi

# Sanity check cmake config file present
if [ -e qmcpack/CMakeLists.txt ]; then

export PYTHONPATH=${test_dir}/qmcpack/nexus/lib
echo --- PYTHONPATH=$PYTHONPATH
#
# Quantum Espresso setup/download/build
# Future improvement: use spack version
#

export QE_VERSION=6.4.1
sys=build_intel2020
# QE version 6.x unpacks to qe-; Older versions 5.x uses espresso-
export QE_PREFIX=qe-
export QE_BIN=${test_dir}/${sys}_QE/${QE_PREFIX}${QE_VERSION}/bin
echo --- QE_BIN set to ${QE_BIN}
if [ ! -e ${QE_BIN}/pw.x ]; then
    # Start from clean build if no executable present
    if [ -e ${test_dir}/${sys}_QE ]; then
	rm -r -f ${test_dir}/${sys}_QE
    fi
    mkdir ${test_dir}/${sys}_QE
		
    cd ${test_dir}/${sys}_QE
    cp -p ../qmcpack/external_codes/quantum_espresso/*${QE_VERSION}* .
    ./download_and_patch_qe${QE_VERSION}.sh
    cd ${QE_PREFIX}${QE_VERSION}

(
   source /opt/intel2020/bin/compilervars.sh intel64
    ./configure CC=mpiicc MPIF90=mpiifort F77=mpiifort --with-scalapack=intel --with-hdf5=/home/pk7/apps/hdf5-1.10.1-intel-mpi
    make pwall # No parallel build due to sometimes broken dependencies in QE build system
)
    echo -- New QE executable `ls -l bin/pw.x`
    cd ${test_dir}
else
    echo -- Found existing QE ${QE_VERSION} executable
fi
# Finished with QE


echo --- Starting test builds and tests

#for sys in build_intel2020 build_clangnew_mkl build_gccnew_mkl build_gcccuda build_clangnew_mkl_nompi build_pgi2019_nompi_mkl build_clangnew_mkl_complex build_gcccuda_complex build_gcccuda_full build_gccnew_mkl_aos build_gccnew_mkl_complex build_gccnew_mkl_complex_aos build_intel2020_complex build_intel2020_mixed build_intel2020_complex_mixed build_intel2020_nompi build_intel2019 build_intel2019_complex build_clangold_mkl build_gccold_mkl build_gccold_mkl_complex
for sys in build_intel2020 build_intel2020_complex
do

echo --- START $sys `date`

cd ${test_dir}

if [ -e $sys ]; then
rm -r -f $sys
fi
mkdir $sys
cd $sys


# Set appropriate environment
if [[ $sys == *"gccnew"* ]]; then
ourenv=gccnewbuild
fi
if [[ $sys == *"gccold"* ]]; then
ourenv=gccoldbuild
fi
if [[ $sys == *"gcccuda"* ]]; then
ourenv=gcccudabuild
fi
if [[ $sys == *"clangnew"* ]]; then
ourenv=clangnewbuild
fi
if [[ $sys == *"clangold"* ]]; then
ourenv=clangoldbuild
fi
if [[ $sys == *"clangcuda"* ]]; then
ourenv=clangcudabuild
fi
if [[ $sys == *"intel"* ]]; then
ourenv=gccintelbuild
fi
if [[ $sys == *"pgi2019"* ]]; then
ourenv=gccnewbuild
fi


#
# FUTURE UPGRADE: Use spack environments  
#
# 2019-06 Not possible to install a clean environment with e.g. only a single libxml2
# Load individual modules and treat ambiguities by hand instead 
#
#echo Build is $sys using env $ourenv
#spack env activate $ourenv
#spack load libxml2
#spack load cmake
#spack load openmpi
#spack load hdf5
#spack load boost


#
# Versions should be consistent with setup script
#

# GCC
# Dates at https://gcc.gnu.org/releases.html
gcc_vnew=9.2.0 # 2019-08-12
gcc_vold=7.3.0 # 2018-01-25
#gcc_vold=7.2.0 # 2017-08-14

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

spack load git

# Python version determined by numpy install in setup script
spack load python@3.7.4%gcc@9.2.0 # Has numpy, scipy, h5py, pandas "activated" and available for import

case "$ourenv" in
gccnewbuild) echo $ourenv
	spack load boost@$boost_vnew%gcc@$gcc_vnew
	spack load gcc@$gcc_vnew
	spack load hdf5@$hdf5_vnew%gcc@$gcc_vnew~mpi
	spack load cmake@$cmake_vnew%gcc@$gcc_vnew
	spack load openmpi@$ompi_vnew%gcc@$gcc_vnew
	spack load libxml2@$libxml2_vnew%gcc@$gcc_vnew
	spack load fftw@$fftw_vnew%gcc@$gcc_vnew
;;
gccoldbuild) echo $ourenv
	spack load boost@$boost_vold%gcc@$gcc_vold
	spack load gcc@$gcc_vold
	spack load hdf5@$hdf5_vold%gcc@$gcc_vold~mpi
	spack load cmake@$cmake_vold%gcc@$gcc_vold
	spack load openmpi@$ompi_vold%gcc@$gcc_vold
	spack load libxml2@$libxml2_vold%gcc@$gcc_vold
	spack load fftw@$fftw_vold%gcc@$gcc_vold
;;
gccintelbuild) echo $ourenv
	spack load boost@$boost_vnew%gcc@$gcc_vnew
	spack load gcc@$gcc_vintel # Provides old enough C++ library for Intel compiler
	spack load hdf5@$hdf5_vnew%gcc@$gcc_vnew~mpi
	spack load cmake@$cmake_vnew%gcc@$gcc_vnew
# Use Intel MPI with Intel builds
#	spack load openmpi@$ompi_vnew%gcc@$gcc_vnew 
	spack load libxml2@$libxml2_vnew%gcc@$gcc_vnew
	spack load fftw@$fftw_vnew%gcc@$gcc_vnew
;;
gcccudabuild) echo $ourenv
	spack load boost@$boost_vnew%gcc@$gcc_vnew
	spack load gcc@$gcc_vcuda
	spack load hdf5@$hdf5_vnew%gcc@$gcc_vnew~mpi
	spack load cmake@$cmake_vnew%gcc@$gcc_vnew
	spack load openmpi@$ompi_vnew%gcc@$gcc_vnew
	spack load libxml2@$libxml2_vnew%gcc@$gcc_vnew
	spack load fftw@$fftw_vnew%gcc@$gcc_vnew
;;
clangnewbuild) echo $ourenv
	spack load llvm@$llvm_vnew
	spack load boost@$boost_vnew%gcc@$gcc_vnew
	spack load gcc@$gcc_vnew
	spack load hdf5@$hdf5_vnew%gcc@$gcc_vnew~mpi
	spack load cmake@$cmake_vnew%gcc@$gcc_vnew
	spack load openmpi@$ompi_vnew%gcc@$gcc_vnew
	spack load libxml2@$libxml2_vnew%gcc@$gcc_vnew
	spack load fftw@$fftw_vnew%gcc@$gcc_vnew
;;
clangoldbuild) echo $ourenv
	spack load llvm@$llvm_vold%gcc@$gcc_vold
	spack load boost@$boost_vold%gcc@$gcc_vold
	spack load gcc@$gcc_vold
	spack load hdf5@$hdf5_vold%gcc@$gcc_vold~mpi
	spack load cmake@$cmake_vold%gcc@$gcc_vold
	spack load openmpi@$ompi_vold%gcc@$gcc_vold
	spack load libxml2@$libxml2_vold%gcc@$gcc_vold
	spack load fftw@$fftw_vold%gcc@$gcc_vold
;;
clangcudabuild) echo $ourenv
	spack load llvm@$llvm_vcuda
	spack load boost@$boost_vnew%gcc@$gcc_vnew
	spack load gcc@$gcc_vnew
	spack load hdf5@$hdf5_vnew%gcc@$gcc_vnew~mpi
	spack load cmake@$cmake_vnew%gcc@$gcc_vnew
	spack load openmpi@$ompi_vnew%gcc@$gcc_vnew
	spack load libxml2@$libxml2_vnew%gcc@$gcc_vnew
	spack load fftw@$fftw_vnew%gcc@$gcc_vnew
;;
*) echo "Problems: Unknown build environment"
	exit 1
;;
esac
module list

# Use subshell to allow Intel MKL or Intel compiler setup to contaminate the environment
# No "unload" capability is provided
(

if [[ $sys == *"intel2020"* ]]; then
source /opt/intel2020/bin/compilervars.sh intel64
fi
if [[ $sys == *"intel2019"* ]]; then
source /opt/intel2019/bin/compilervars.sh intel64
fi
if [[ $sys == *"intel2018"* ]]; then
source /opt/intel2018/bin/compilervars.sh intel64
fi

# Construct test name and configure flags
# Compiler and major version, MPI or not
if [[ $sys == *"gcc"* ]]; then
compilerversion=`gcc --version|grep ^gcc|sed 's/^.* //g'|sed 's/\..*//g'`
if [[ $sys == *"nompi"* ]]; then
QMCPACK_TEST_SUBMIT_NAME=GCC${compilerversion}-NoMPI
CTCFG="-DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DQMC_MPI=0"
else
QMCPACK_TEST_SUBMIT_NAME=GCC${compilerversion}
CTCFG="-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DQMC_MPI=1"
export OMPI_CC=gcc
export OMPI_CXX=g++
fi
fi

#Clang/LLVM
if [[ $sys == *"clang"* ]]; then
compilerversion=`clang --version|grep ^clang|sed -e 's/^.* version //g' -e 's/(.*//g'|sed 's/\..*//g'`
if [[ $sys == *"nompi"* ]]; then
QMCPACK_TEST_SUBMIT_NAME=Clang${compilerversion}-NoMPI
CTCFG="-DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DQMC_MPI=0"
else
QMCPACK_TEST_SUBMIT_NAME=Clang${compilerversion}
CTCFG="-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DQMC_MPI=1"
export OMPI_CC=clang
export OMPI_CXX=clang++
fi
fi

# Intel
if [[ $sys == *"intel"* ]]; then
# Assumie the Intel dumpversion string has format AA.B.C.DDD
# AA is historically major year and a good enough unique identifier,
# but Intel 2020 packages currently (12/2019) come with the 19.1 compiler.
# Use 19.1 for this compiler only to avoid breaking test histories.
if [[ $sys == *"intel2020"* ]]; then
compilerversion=`icc -dumpversion|sed -e 's/\([0-9][0-9]\)\.\([0-9]\)\..*/\1.\2/g'` # AA.B
else
compilerversion=`icc -dumpversion|sed -e 's/\..*//g'` # AA year only
fi

if [[ $sys == *"nompi"* ]]; then
QMCPACK_TEST_SUBMIT_NAME=Intel20${compilerversion}-NoMPI
CTCFG="-DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DQMC_MPI=0"
else
QMCPACK_TEST_SUBMIT_NAME=Intel20${compilerversion}
CTCFG="-DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DQMC_MPI=1"
# Add QE to any Intel MPI builds.
# Do not add to non Intel MPI because ctest will use mpirun from "wrong" MPI
CTCFG="$CTCFG -DQE_BIN=${QE_BIN}" 
fi
fi

# PGI
if [[ $sys == *"pgi"* ]]; then
compilerversion=`pgcc -V|grep pgcc|sed 's/^pgcc //g'|sed 's/\..*//g'`
if [[ $sys == *"nompi"* ]]; then
QMCPACK_TEST_SUBMIT_NAME=PGI20${compilerversion}-NoMPI
CTCFG="-DCMAKE_C_COMPILER=pgcc -DCMAKE_CXX_COMPILER=pgc++ -DQMC_MPI=0"
else
CTCFG="-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DQMC_MPI=1"
export OMPI_CC=pgcc
export OMPI_CXX=pgc++
fi
fi

# CUDA
if [[ $sys == *"cuda"* ]]; then
QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-CUDA
CTCFG="$CTCFG -DQMC_CUDA=1"
fi

# MKL
# MKLROOT set in sourced Intel compilervars.sh 
if [[ $sys == *"mkl"* ]]; then
QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-MKL
# MKL setup used by many builds for BLAS, LAPACK etc.
source /opt/intel2020/mkl/bin/mklvars.sh intel64
CTCFG="$CTCFG -DENABLE_MKL=1 -DBLA_VENDOR=Intel10_64lp_seq -DCMAKE_PREFIX_PATH=$MKLROOT/lib"
fi

# Complex
if [[ $sys == *"complex"* ]]; then
QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-Complex
CTCFG="$CTCFG -DQMC_COMPLEX=1"
fi

# Mixed/Full precision
if [[ $sys == *"mixed"* ]]; then
QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-Mixed
CTCFG="$CTCFG -DQMC_MIXED_PRECISION=1"
fi
if [[ $sys == *"full"* ]]; then
QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-Full
CTCFG="$CTCFG -DQMC_MIXED_PRECISION=0"
fi

# SoA/AoS build (label aos only)
if [[ $sys == *"aos"* ]]; then
QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-AoS
CTCFG="$CTCFG -DENABLE_SOA=0"
else
CTCFG="$CTCFG -DENABLE_SOA=1"
fi

# Boilerplate for all tests
CTCFG="$CTCFG -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1"

# Selectively enable AFQMC due to more stringent compiler requirements
if [[ $sys == *"complex"* ]]; then
case "$sys" in
*gccnew*|*clangnew*|*cuda*|*intel*) echo "AFQMC is enabled for this complex build"
CTCFG="$CTCFG -DBUILD_AFQMC=1"
;;
*) echo "AFQMC is disabled for this complex build (either unsupported or unchecked)"
CTCFG="$CTCFG -DBUILD_AFQMC=0"
;;
esac
fi

# Adjust which tests are run to control overall runtime
case "$sys" in
*gccnew_mkl|*gccnew_mkl_complex|*gccnew_mkl_aos|*gccnew_mkl_complex_aos|*intel2020|*intel2020_complex|*clangnew_mkl|*clangnew_mkl_complex|*gcccuda|*gcccuda_complex|*gcccuda_full) echo "Running full ("less limited") test set for $sys"
THETESTS=$LESSLIMITEDTESTS
;;
*) echo "Running limited test set for $sys"
THETESTS=$LIMITEDTESTS
;;
esac

export QMCPACK_TEST_SUBMIT_NAME=${QMCPACK_TEST_SUBMIT_NAME}-Release

echo $QMCPACK_TEST_SUBMIT_NAME
echo $CTCFG
if [[ $localonly == "yes" ]]; then
echo --- START cmake `date` 
cmake ${CTCFG} ${GLOBALTCFG} -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 ../qmcpack/ 
echo --- END cmake `date`
echo --- START make `date` 
make -j 16
echo --- END make `date`
echo --- START ctest `date` 
ctest ${GLOBALTCFG} ${THETESTS}
echo --- END ctest `date`
else
echo --- START ctest `date` 
echo ctest ${CTCFG} ${GLOBALTCFG} -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release ${THETESTS}
ctest ${CTCFG} ${GLOBALTCFG} -DQMC_DATA=${QMC_DATA} -DENABLE_TIMERS=1 -S $PWD/../qmcpack/CMake/ctest_script.cmake,release ${THETESTS}
echo --- END ctest `date`
fi

)

module purge

echo --- END $sys `date`
done

else
echo "ERROR: No CMakeLists. Bad git clone or update"
exit 1
fi

else
echo "ERROR: Unable to make test directory ${test_dir}"
exit 1
fi

else
echo "ERROR: No directory ${test_path}"
exit 1
fi
echo --- Script END `date`
