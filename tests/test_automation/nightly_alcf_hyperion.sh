#!/bin/bash
# source environment
# MPI wrappers, MKL, and Intel and GCC compiler
source /opt/intel/2017/parallel_studio_xe_2017.1.043/bin/psxevars.sh intel64

# timeout
timeout=1800

# topdir must exist, otherwise it fails
topdir=/home/naromero/ci
place=${topdir}/scratch/QMCPACK_CI_BUILDS_DO_NOT_REMOVE

if [ -e $topdir ]; then

if [ ! -e $place ]; then
mkdir -p $place
fi

if [ -e $place ]; then


for sys in build_gcc build_gcc_complex build_gcc_mixed build_gcc_complex_mixed build_intel2017 build_intel2017_complex build_intel2017_mixed build_intel2017_complex_mixed
do

cd $place

if [ -e $sys ]; then
rm -rf $sys
fi
mkdir $sys
cd $sys

echo --- Checkout for $sys `date`
svn checkout https://svn.qmcpack.org/svn/trunk
#svn checkout https://subversion.assembla.com/svn/qmcdev/trunk

if [ -e trunk/CMakeLists.txt ]; then
cd trunk
mkdir $sys
cd $sys
echo --- Building for $sys `date`

# export PATH=/opt/local/bin:/opt/local/sbin:/usr/local/cuda/bin/:/usr/lib64/openmpi/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin
# export LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:/usr/local/cuda-7.0/lib64

case $sys in
    "build_gcc")
	export QMCPACK_TEST_SUBMIT_NAME=GCC-Release
	ctest -DQMC_COMPLEX=0 -DQMC_MIXED_PRECISION=0 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -S $PWD/../CMake/ctest_script.cmake,release -E 'long' -VV --timeout $timeout
	;;
    "build_gcc_complex")
	export QMCPACK_TEST_SUBMIT_NAME=GCC-Complex-Release
	ctest -DQMC_COMPLEX=1 -DQMC_MIXED_PRECISION=0 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -S $PWD/../CMake/ctest_script.cmake,release -E 'long' -VV --timeout $timeout
	;;
    "build_gcc_mixed")
	export QMCPACK_TEST_SUBMIT_NAME=GCC-Mixed-Release
	ctest -DQMC_COMPLEX=0 -DQMC_MIXED_PRECISION=1 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -S $PWD/../CMake/ctest_script.cmake,release -E 'long' -VV --timeout $timeout
	;;
    "build_gcc_complex_mixed")
	export QMCPACK_TEST_SUBMIT_NAME=GCC-Complex-Mixed-Release
	ctest -DQMC_COMPLEX=1 -DQMC_MIXED_PRECISION=1 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -S $PWD/../CMake/ctest_script.cmake,release -E 'long' -VV --timeout $timeout
	;;
    "build_intel2017")
	source /opt/intel/2017/parallel_studio_xe_2017.1.043/bin/psxevars.sh intel64
	export QMCPACK_TEST_SUBMIT_NAME=Intel2017-Release
	ctest -DQMC_COMPLEX=0 -DQMC_MIXED_PRECISION=0 -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=mpiicpc -S $PWD/../CMake/ctest_script.cmake,release -E 'long' -VV --timeout $timeout
	;;
    "build_intel2017_complex")
	# intel compiler should already be loaded 
	# source /opt/intel/2017/parallel_studio_xe_2017.1.043/bin/psxevars.sh intel64
	export QMCPACK_TEST_SUBMIT_NAME=Intel2017-Complex-Release
	ctest -DQMC_COMPLEX=1 -DQMC_MIXED_PRECISION=0 -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=mpiicpc -S $PWD/../CMake/ctest_script.cmake,release -E 'long' -VV --timeout $timeout
	;;
    "build_intel2017_mixed")
	# intel compiler should already be loaded
	# source /opt/intel/2017/parallel_studio_xe_2017.1.043/bin/psxevars.sh intel64
	export QMCPACK_TEST_SUBMIT_NAME=Intel2017-Mixed-Release
	ctest -DQMC_COMPLEX=0 -DQMC_MIXED_PRECISION=1 -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=mpiicpc -S $PWD/../CMake/ctest_script.cmake,release -E 'long' -VV --timeout $timeout
	;;
    "build_intel2017_complex_mixed")
	# intel compiler should already be loaded
	# source /opt/intel/2017/parallel_studio_xe_2017.1.043/bin/psxevars.sh intel64
	export QMCPACK_TEST_SUBMIT_NAME=Intel2017-Complex-Mixed-Release
	ctest -DQMC_COMPLEX=1 -DQMC_MIXED_PRECISION=1 -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=mpiicpc -S $PWD/../CMake/ctest_script.cmake,release -E 'long' -VV --timeout $timeout
	;;
    *)
	echo "ERROR: Unknown build type $sys"
	;;
esac

else
    echo  "ERROR: No CMakeLists. Bad svn checkout."
    exit 1
fi

done

else
echo "ERROR: No directory $place"
exit 1
fi

fi

