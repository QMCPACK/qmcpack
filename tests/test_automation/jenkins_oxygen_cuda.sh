#!/bin/bash -x

export PATH=/usr/local/cuda/bin:${PATH}
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:/usr/local/lib:${LD_LIBRARY_PATH}

exit_code=0

BUILD_DIR=$(pwd)
echo $BUILD_DIR

echo ""
echo ""
echo "cache source in tmpfs"
echo "at $(date)"
echo ""
echo ""

rm -rf /dev/shm/${BUILD_TAG}-src
cp -R ${BUILD_DIR} /dev/shm/${BUILD_TAG}-src


echo ""
echo ""
echo "starting new test for real full precision"
echo "at $(date)"
echo ""
echo ""

mkdir -p /dev/shm/${BUILD_TAG}-build
cd /dev/shm/${BUILD_TAG}-build

time cmake -DQMC_COMPLEX=0 -DQMC_MIXED_PRECISION=0 -DENABLE_SOA=0 -DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpicxx" -DENABLE_CUDA=1 -DQMC_NO_SLOW_CUSTOM_TESTING_COMMANDS=1 /dev/shm/${BUILD_TAG}-src 2>&1 | tee cmake.out
if [[ $? -ne 0 ]] ; then
  rm -rf /dev/shm/${BUILD_TAG}-build
  rm -rf /dev/shm/${BUILD_TAG}-src
  exit 1
fi

make -j 8
if [[ $? -ne 0 ]] ; then
  rm -rf /dev/shm/${BUILD_TAG}-build
  rm -rf /dev/shm/${BUILD_TAG}-src
  exit 1
fi

ctest -L unit --output-on-failure --timeout 120
ret=$?
if [[ ${ret} -ne 0 ]] ; then
  exit_code=${ret}
fi

echo ""
echo ""
echo "starting new test for real mixed precision"
echo "at $(date)"
echo ""
echo ""

cd ../
rm -rf ./${BUILD_TAG}-build
mkdir -p ${BUILD_TAG}-build
cd ${BUILD_TAG}-build

time cmake -DQMC_COMPLEX=0 -DQMC_MIXED_PRECISION=1 -DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpicxx" -DENABLE_CUDA=1 -DQMC_NO_SLOW_CUSTOM_TESTING_COMMANDS=1/dev/shm/${BUILD_TAG}-src 2>&1 | tee cmake.out
if [[ $? -ne 0 ]] ; then
  rm -rf ./${BUILD_TAG}-build
  rm -rf ./${BUILD_TAG}-src
  exit 1
fi

make -j 8
if [[ $? -ne 0 ]] ; then
  rm -rf /dev/shm/${BUILD_TAG}-build
  rm -rf /dev/shm/${BUILD_TAG}-src
  exit 1
fi

ctest -L unit --output-on-failure --timeout 120
ret=$?
if [[ ${ret} -ne 0 ]] ; then
  exit_code=${ret}
fi

echo ""
echo ""
echo "starting new test for complex full precision"
echo "at $(date)"
echo ""
echo ""

cd ../
rm -rf ./${BUILD_TAG}-build
mkdir -p ${BUILD_TAG}-build
cd ${BUILD_TAG}-build

time cmake -DQMC_COMPLEX=1 -DQMC_MIXED_PRECISION=0 -DENABLE_SOA=0 -DBUILD_AFQMC=1 -DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpicxx" -DENABLE_CUDA=1 -DQMC_NO_SLOW_CUSTOM_TESTING_COMMANDS=1 /dev/shm/${BUILD_TAG}-src 2>&1 | tee cmake.out
if [[ $? -ne 0 ]] ; then
  rm -rf /dev/shm/${BUILD_TAG}-build
  rm -rf /dev/shm/${BUILD_TAG}-src
  exit 1
fi

make -j 8
if [[ $? -ne 0 ]] ; then
  rm -rf /dev/shm/${BUILD_TAG}-build
  rm -rf /dev/shm/${BUILD_TAG}-src
  exit 1
fi

ctest -L unit --output-on-failure --timeout 120
ret=$?
if [[ ${ret} -ne 0 ]] ; then
  exit_code=${ret}
fi

echo ""
echo ""
echo "starting new test for complex mixed precision"
echo "at $(date)"
echo ""
echo ""

cd ../
rm -rf ./${BUILD_TAG}-build
mkdir -p ${BUILD_TAG}-build
cd ${BUILD_TAG}-build

time cmake -DQMC_COMPLEX=1 -DQMC_MIXED_PRECISION=1 -DBUILD_AFQMC=1 -DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpicxx" -DENABLE_CUDA=1 -DQMC_NO_SLOW_CUSTOM_TESTING_COMMANDS=1 /dev/shm/${BUILD_TAG}-src 2>&1 | tee cmake.out
if [[ $? -ne 0 ]] ; then
  rm -rf /dev/shm/${BUILD_TAG}-build
  rm -rf /dev/shm/${BUILD_TAG}-src
  exit 1
fi

make -j 8
if [[ $? -ne 0 ]] ; then
  rm -rf /dev/shm/${BUILD_TAG}-build
  rm -rf /dev/shm/${BUILD_TAG}-src
  exit 1
fi

ctest -L unit --output-on-failure --timeout 120
ret=$?
if [[ ${ret} -ne 0 ]] ; then
  exit_code=${ret}
fi

# final cleanup of tmpfs
cd ../
rm -rf ./${BUILD_TAG}-build
rm -rf ./${BUILD_TAG}-src

# return the test results to Jenkins
exit ${exit_code}
