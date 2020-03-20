#!/bin/bash -x

exit_code=0

echo "starting spack using ${SPACK_ENV_FILE}"

# this depends on SPACK_ROOT being set in Jenkinsfile_xxx
# it also supplies QMC_IMMUTABLE_FLAGS which makes it a bit more than the
# environment from the set of loaded spack packages.
. ${SPACK_ENV_FILE}

cd build_${1}_${2}
BUILD_DIR=$(pwd)
echo "BUILD_DIR: ${BUILD_DIR}"

export LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH}

module list

# this keeps tee from eating the exit status
set -o pipefail

# exclude failing hamiltonian test
EXCLUDE_FLAG=''
if [[ $2 == 'mixed' ]]; then
    EXCLUDE_FLAG='-E unit_test_hamiltonian'
fi


ctest -j${JNK_THREADS} -L unit ${EXCLUDE_FLAG} --output-on-failure --timeout 120 2>&1 | tee ${1}_${2}_ctest.out
ret=$?
if [[ ${ret} -ne 0 ]] ; then
  exit_code=${ret}
fi

exit ${exit_code}
