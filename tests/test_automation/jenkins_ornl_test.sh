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

ctest -j${JNK_THREADS} -L unit --output-on-failure --timeout 120 2>&1 | tee ${1}_${2}_ctest.out
ret=$?
if [[ ${ret} -ne 0 ]] ; then
  exit_code=${ret}
fi

exit ${exit_code}
