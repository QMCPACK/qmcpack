#!/bin/bash -x

export LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH}

exit_code=0

BUILD_DIR=$(pwd)
echo $BUILD_DIR

QMCNSPACE_FLAG=''
if [[ $1 == 'real' ]]
then
    QMCNSPACE_FLAG="-DQMC_COMPLEX=0"
else
    QMCNSPACE_FLAG="-DQMC_COMPLEX=1"
fi

QMCPRECISION_FLAG=''
if [[ $2 == 'full' ]]
then
    QMCPRECISION_FLAG="-DQMC_MIXED_PRECISION=0"
else
    QMCPRECISION_FLAG="-DQMC_MIXED_PRECISION=1"
fi

echo ""
echo ""
echo "starting build for ${1} ${2} precision"
echo "at $(date)"
echo ""
echo ""

mkdir build_${1}_${2}
cd build_${1}_${2}

cmake ${QMCNSPACE_FLAG} ${QMCPRECISION_FLAG} -DENABLE_SOA=0 -DBUILD_AFQMC=1 -DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpicxx" -DQMC_NO_SLOW_CUSTOM_TESTING_COMMANDS=1 ../.. 2>&1 | tee cmake.out

if [[ $? -ne 0 ]] ; then
  exit 1
fi

make -j 16
if [[ $? -ne 0 ]] ; then
  exit 1
fi

# return the test results to Jenkins
exit ${exit_code}
