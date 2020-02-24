#!/bin/bash -x

exit_code = 0

cd build_${ARGV}[1]_${ARGV}[2]
BUILD_DIR=$(pwd)
echo $BUILD_DIR


ctest -j8 -L unit --output-on-failure --timeout 120
ret=$?
if [[ ${ret} -ne 0 ]] ; then
  exit_code=${ret}
fi

exit ${exit_code}
