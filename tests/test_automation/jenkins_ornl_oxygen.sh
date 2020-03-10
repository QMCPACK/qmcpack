#!/bin/bash --login

#this is will be general move cpu build soon

export LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH}

echo "starting spack using ${SPACK_ENV_FILE}"

# this depends on SPACK_ROOT being set in Jenkinsfile_xxx
# it also supplies QMC_IMMUTABLE_FLAGS which makes it a bit more than the
# environment from the set of loaded spack packages.
. ${SPACK_ENV_FILE}

exit_code=0

BUILD_DIR=$(pwd)
echo $BUILD_DIR

#translate from jenkins matrix options to cmake options
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

rm -rf build_${1}_${2}
mkdir build_${1}_${2}
cd build_${1}_${2}

cmake ${QMCNSPACE_FLAG} ${QMCPRECISION_FLAG} -DENABLE_SOA=1 -DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpicxx" ${QMC_IMMUTABLE_FLAGS} -DQMC_NO_SLOW_CUSTOM_TESTING_COMMANDS=1 ../.. 2>&1 | tee cmake.out

if [[ $? -ne 0 ]] ; then
  exit 1
fi

make -j ${JNK_THREADS}
if [[ $? -ne 0 ]] ; then
  exit 1
fi

# return the test results to Jenkins
exit ${exit_code}
