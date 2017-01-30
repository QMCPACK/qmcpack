#!/bin/bash
# source environment
# MPI wrappers, MKL, and Intel and GCC compiler
source /opt/intel/2017/parallel_studio_xe_2017.1.043/bin/psxevars.sh intel64

# timeout
timeout=1800

# topdir must exist, otherwise it fails
topdir=/data/ci
testdir=${topdir}/scratch/QMCPACK_CI_BUILDS_DO_NOT_REMOVE

if [ -e $topdir ]; then

if [ ! -e $testdir ]; then
mkdir -p $testdir
fi

if [ -e $testdir ]; then
cd $testdir

# Minimize load of GitHub by maintaining a local cloned git used for all builds
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

echo --- Starting test builds and tests

for sys in build_gcc build_gcc_complex build_gcc_mixed build_gcc_complex_mixed build_intel2017 build_intel2017_complex build_intel2017_mixed build_intel2017_complex_mixed
do

echo --- Building for $sys `date`

cd ${testdir}

if [ -e $sys ]; then
rm -rf $sys
fi
mkdir $sys
cd $sys

# export PATH=/opt/local/bin:/opt/local/sbin:/usr/local/cuda/bin/:/usr/lib64/openmpi/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin
# export LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:/usr/local/cuda-7.0/lib64

case $sys in
    "build_gcc")
	export QMCPACK_TEST_SUBMIT_NAME=GCC-MKL-Release
	ctest -DQMC_COMPLEX=0 -DQMC_MIXED_PRECISION=0 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -E 'long' -VV --timeout $timeout
	;;
    "build_gcc_complex")
	export QMCPACK_TEST_SUBMIT_NAME=GCC-MKL-Complex-Release
	ctest -DQMC_COMPLEX=1 -DQMC_MIXED_PRECISION=0 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -E 'long' -VV --timeout $timeout
	;;
    "build_gcc_mixed")
	export QMCPACK_TEST_SUBMIT_NAME=GCC-MKL-Mixed-Release
	ctest -DQMC_COMPLEX=0 -DQMC_MIXED_PRECISION=1 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -E 'long' -VV --timeout $timeout
	;;
    "build_gcc_complex_mixed")
	export QMCPACK_TEST_SUBMIT_NAME=GCC-MKL-Complex-Mixed-Release
	ctest -DQMC_COMPLEX=1 -DQMC_MIXED_PRECISION=1 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -E 'long' -VV --timeout $timeout
	;;
    "build_intel2017")
	source /opt/intel/2017/parallel_studio_xe_2017.1.043/bin/psxevars.sh intel64
	export QMCPACK_TEST_SUBMIT_NAME=Intel2017-Release
	ctest -DQMC_COMPLEX=0 -DQMC_MIXED_PRECISION=0 -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=mpiicpc -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -E 'long' -VV --timeout $timeout
	;;
    "build_intel2017_complex")
	# intel compiler should already be loaded 
	# source /opt/intel/2017/parallel_studio_xe_2017.1.043/bin/psxevars.sh intel64
	export QMCPACK_TEST_SUBMIT_NAME=Intel2017-Complex-Release
	ctest -DQMC_COMPLEX=1 -DQMC_MIXED_PRECISION=0 -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=mpiicpc -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -E 'long' -VV --timeout $timeout
	;;
    "build_intel2017_mixed")
	# intel compiler should already be loaded
	# source /opt/intel/2017/parallel_studio_xe_2017.1.043/bin/psxevars.sh intel64
	export QMCPACK_TEST_SUBMIT_NAME=Intel2017-Mixed-Release
	ctest -DQMC_COMPLEX=0 -DQMC_MIXED_PRECISION=1 -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=mpiicpc -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -E 'long' -VV --timeout $timeout
	;;
    "build_intel2017_complex_mixed")
	# intel compiler should already be loaded
	# source /opt/intel/2017/parallel_studio_xe_2017.1.043/bin/psxevars.sh intel64
	export QMCPACK_TEST_SUBMIT_NAME=Intel2017-Complex-Mixed-Release
	ctest -DQMC_COMPLEX=1 -DQMC_MIXED_PRECISION=1 -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=mpiicpc -S $PWD/../qmcpack/CMake/ctest_script.cmake,release -E 'long' -VV --timeout $timeout
	;;
    *)
	echo "ERROR: Unknown build type $sys"
	;;
esac

done # termination for $sys loop

else
    echo  "ERROR: No CMakeLists. Bad git clone or update."
    exit 1
fi # termination for if block which checks CMakeLists.txt


else
echo "ERROR: Unable to make test directory ${test_dir}"
exit 1
fi # termination for if block which checks testing directory

else
echo "ERROR: No top level directory"
exit 1
fi # termination for if block which checks top level directory



