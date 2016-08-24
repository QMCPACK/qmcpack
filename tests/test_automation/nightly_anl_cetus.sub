#!/bin/bash
#COBALT -q default
#COBALT -A QMCSim
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
place=/projects/QMCSim/QMCPACK_CI_BUILDS_DO_NOT_REMOVE

if [ ! -e $place ]; then
mkdir $place
fi

if [ -e $place ]; then
cd $place

echo --- Hostname --- $HOSTNAME
echo --- Checkout for $sys `date`
svn checkout https://svn.qmcpack.org/svn/trunk
#svn checkout https://subversion.assembla.com/svn/qmcdev/branches/GPU_precision_validation

entry=trunk
#entry=GPU_precision_validation

if [ -e $entry/CMakeLists.txt ]; then
cd $entry

for sys in build_BGQ_XL_real #build_BGQ_XL_complex
do

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
"build_BGQ_XL_real")
# Build the real version with XL compiler on BGQ
    export QMCPACK_TEST_SUBMIT_NAME=BGQ-XL-Real-Release
    ctest -DCMAKE_TOOLCHAIN_FILE=../config/BGQToolChain.cmake -DQMC_COMPLEX=0 -S $PWD/../CMake/ctest_script.cmake,release -R short -VV &> $place/log/$entry/$mydate/${QMCPACK_TEST_SUBMIT_NAME}.log
    ;;
"build_BGQ_XL_complex")
# Build the complex version with XL compiler on BGQ
    export QMCPACK_TEST_SUBMIT_NAME=BGQ-XL-Complex-Release
    ctest -DCMAKE_TOOLCHAIN_FILE=../config/BGQToolChain.cmake -DQMC_COMPLEX=1 -S $PWD/../CMake/ctest_script.cmake,release -R "short|unit|qe-LiH-unpolarized-no-collect" -VV &> $place/log/$entry/$mydate/${QMCPACK_TEST_SUBMIT_NAME}.log
    ;;
*)
    echo "ERROR: Unknown build type $sys"
    ;;
esac

cd ..
echo --- Finished $sys `date`
done

else
echo  "ERROR: No CMakeLists. Bad svn checkout."
exit 1
fi

else
echo "ERROR: No directory $place"
exit 1
fi
