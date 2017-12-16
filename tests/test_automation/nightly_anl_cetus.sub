#!/bin/bash
#COBALT -q default
#COBALT -A PSFMat
#COBALT -n 128
#COBALT -t 60
#COBALT -O validation
#COBALT -M xw111luoye@gmail.com

#
# Setup for cetus.alcf.anl.gov
#
# Run the "short" nightlies, requeue if an executable is built
# 
# Location of job script must be set at end ($HOME/.qmcpack_test_jobs/... ) for resubmit
#
# Checkout, build, and run on scratch. Use custom names for checkout
# directory because scratch is shared.
#

#Must be an absolute path
place=/projects/PSFMat/QMCPACK_CI_BUILDS_DO_NOT_REMOVE

#define compiler
#compiler=XL
compiler=Clang++11

if [ ! -e $place ]; then
mkdir $place
fi

if [ -e $place ]; then
cd $place

echo --- Hostname --- $HOSTNAME
echo --- Checkout for $sys `date`

branch=develop
entry=qmcpack-${branch}

if [ ! -e $entry ]; then
echo --- Cloning QMCPACK git `date`
git clone --depth 1 https://github.com/QMCPACK/qmcpack.git $entry
else
echo --- Updating local QMCPACK git `date`
cd $entry
git pull
cd ..
fi

if [ -e $entry/CMakeLists.txt ]; then
cd $entry

git checkout $branch

for sys in build_BGQ_${compiler}_real build_BGQ_${compiler}_complex build_BGQ_${compiler}_MP_real build_BGQ_${compiler}_MP_complex build_BGQ_${compiler}_real_SoA build_BGQ_${compiler}_complex_SoA
do

if [ -e $sys ]; then
rm -r $sys
fi
mkdir $sys
cd $sys

echo --- Building for $sys `date`

# create log file folder if not exist
mydate=`date +%y_%m_%d`
if [ ! -e $place/log/$entry/$mydate ];
then
  mkdir -p $place/log/$entry/$mydate
fi

case $sys in
"build_BGQ_${compiler}_real_SoA")
# Build the real version with ${compiler} compiler on BGQ
    export QMCPACK_TEST_SUBMIT_NAME=BGQ-${compiler}-Real-SoA-Release
    ctest -DCMAKE_TOOLCHAIN_FILE=../config/BGQ_${compiler}_ToolChain.cmake -DQMC_COMPLEX=0 -DENABLE_SOA=1 -S $PWD/../CMake/ctest_script.cmake,release -R "short" -E "dmc|opt|csvmc" -VV &> $place/log/$entry/$mydate/${QMCPACK_TEST_SUBMIT_NAME}.log
    ;;
"build_BGQ_${compiler}_complex_SoA")
# Build the complex version with ${compiler} compiler on BGQ
    export QMCPACK_TEST_SUBMIT_NAME=BGQ-${compiler}-Complex-SoA-Release
    ctest -DCMAKE_TOOLCHAIN_FILE=../config/BGQ_${compiler}_ToolChain.cmake -DQMC_COMPLEX=1 -DENABLE_SOA=1 -S $PWD/../CMake/ctest_script.cmake,release -R "short" -E "dmc|opt|csvmc" -VV &> $place/log/$entry/$mydate/${QMCPACK_TEST_SUBMIT_NAME}.log
    ;;
"build_BGQ_${compiler}_real")
# Build the real version with ${compiler} compiler on BGQ
    export QMCPACK_TEST_SUBMIT_NAME=BGQ-${compiler}-Real-Release
    ctest -DCMAKE_TOOLCHAIN_FILE=../config/BGQ_${compiler}_ToolChain.cmake -DQMC_COMPLEX=0 -S $PWD/../CMake/ctest_script.cmake,release -R short -E "dmc|short-H4-opt" -VV &> $place/log/$entry/$mydate/${QMCPACK_TEST_SUBMIT_NAME}.log
    ;;
"build_BGQ_${compiler}_complex")
# Build the complex version with ${compiler} compiler on BGQ
    export QMCPACK_TEST_SUBMIT_NAME=BGQ-${compiler}-Complex-Release
    ctest -DCMAKE_TOOLCHAIN_FILE=../config/BGQ_${compiler}_ToolChain.cmake -DQMC_COMPLEX=1 -S $PWD/../CMake/ctest_script.cmake,release -R "unit|short|qe-LiH-unpolarized-no-collect-np-16-12-12-nk-16-2-2" -E "dmc|short-H4-opt" -VV &> $place/log/$entry/$mydate/${QMCPACK_TEST_SUBMIT_NAME}.log
    ;;
"build_BGQ_${compiler}_MP_real")
# Build the real version with ${compiler} compiler on BGQ
    export QMCPACK_TEST_SUBMIT_NAME=BGQ-${compiler}-MixedPrec-Real-Release
    ctest -DCMAKE_TOOLCHAIN_FILE=../config/BGQ_${compiler}_ToolChain.cmake -DQMC_COMPLEX=0 -DQMC_MIXED_PRECISION=1 -S $PWD/../CMake/ctest_script.cmake,release -R short -E "dmc|short-H4-opt" -VV &> $place/log/$entry/$mydate/${QMCPACK_TEST_SUBMIT_NAME}.log
    ;;
"build_BGQ_${compiler}_MP_complex")
# Build the complex version with ${compiler} compiler on BGQ
    export QMCPACK_TEST_SUBMIT_NAME=BGQ-${compiler}-MixedPrec-Complex-Release
    ctest -DCMAKE_TOOLCHAIN_FILE=../config/BGQ_${compiler}_ToolChain.cmake -DQMC_COMPLEX=1 -DQMC_MIXED_PRECISION=1 -S $PWD/../CMake/ctest_script.cmake,release -R "unit|short|qe-LiH-unpolarized-no-collect-np-16-12-12-nk-16-2-2" -E "dmc|short-H4-opt" -VV &> $place/log/$entry/$mydate/${QMCPACK_TEST_SUBMIT_NAME}.log
    ;;
*)
    echo "ERROR: Unknown build type $sys"
    ;;
esac

cd ..
echo --- Finished $sys `date`
done

else
echo  "ERROR: No CMakeLists. Bad git clone."
exit 1
fi

else
echo "ERROR: No directory $place"
exit 1
fi
