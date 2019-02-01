#!/bin/bash -x

BUILD_DIR=$(pwd)
echo $BUILD_DIR

cat > $BUILD_TAG.pbs << EOF
#PBS -A MAT151ci
#PBS -N $BUILD_TAG
#PBS -j oe
#PBS -l walltime=1:00:00,nodes=1
#PBS -d $BUILD_DIR
#PBS -l partition=gpu

source /sw/rhea/environment-modules/3.2.10/rhel6.7_gnu4.4.7/init/bash

module unload PE-intel
module load PE-gnu/5.3.0-1.10.2
module load fftw
export FFTW_HOME=\$FFTW3_DIR
module load hdf5
module load git
module load cudatoolkit/8.0.44
module load cmake/3.6.1
module load boost/1.61.0

env
module list

echo ""
echo ""
echo "cache source in tmpfs"
echo "at $(date)"
echo ""
echo ""

rm -rf /tmp/${BUILD_TAG}-src
cp -R ${BUILD_DIR} /tmp/${BUILD_TAG}-src


echo ""
echo ""
echo "starting new test for real full precision"
echo "at $(date)"
echo ""
echo ""

mkdir -p /tmp/${BUILD_TAG}-build
cd /tmp/${BUILD_TAG}-build

time cmake -DQMC_COMPLEX=0 -DQMC_MIXED_PRECISION=0 -DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpicxx" -DCMAKE_CXX_FLAGS="-mno-bmi2 -mno-avx2" -DCMAKE_C_FLAGS="-mno-bmi2 -mno-avx2" -DBLAS_blas_LIBRARY="/usr/lib64/libblas.so.3" -DLAPACK_lapack_LIBRARY="/usr/lib64/atlas/liblapack.so.3" -DHDF5_INCLUDE_DIR="/sw/rhea/hdf5/1.8.11/rhel6.6_intel14.0.4/include" -DQMC_CUDA=1 /tmp/${BUILD_TAG}-src 2>&1 | tee cmake.out

# hacky way to check on cmake. works for now
if ! ( grep -- '-- The C compiler identification is GNU 5.3.0' cmake.out && \
       grep -- '-- The CXX compiler identification is GNU 5.3.0' cmake.out ) ;
then
  echo "compiler version mismatch. exiting."
  exit 1
fi

# because Andreas tells me (and I observe) that GPU builds are unstable with Cmake
time make -j 24
time make -j 24

time ctest -L unit --output-on-failure



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

time cmake -DQMC_COMPLEX=0 -DQMC_MIXED_PRECISION=1 -DENABLE_SOA=1 -DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpicxx" -DCMAKE_CXX_FLAGS="-mno-bmi2 -mno-avx2" -DCMAKE_C_FLAGS="-mno-bmi2 -mno-avx2" -DBLAS_blas_LIBRARY="/usr/lib64/libblas.so.3" -DLAPACK_lapack_LIBRARY="/usr/lib64/atlas/liblapack.so.3" -DHDF5_INCLUDE_DIR="/sw/rhea/hdf5/1.8.11/rhel6.6_intel14.0.4/include" -DQMC_CUDA=1 /tmp/${BUILD_TAG}-src 2>&1 | tee cmake.out

time make -j 24
time make -j 24

time ctest -L unit --output-on-failure



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

time cmake -DQMC_COMPLEX=1 -DQMC_MIXED_PRECISION=0 -DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpicxx" -DCMAKE_CXX_FLAGS="-mno-bmi2 -mno-avx2" -DCMAKE_C_FLAGS="-mno-bmi2 -mno-avx2" -DBLAS_blas_LIBRARY="/usr/lib64/libblas.so.3" -DLAPACK_lapack_LIBRARY="/usr/lib64/atlas/liblapack.so.3" -DHDF5_INCLUDE_DIR="/sw/rhea/hdf5/1.8.11/rhel6.6_intel14.0.4/include" -DQMC_CUDA=1 /tmp/${BUILD_TAG}-src 2>&1 | tee cmake.out

time make -j 24
time make -j 24

time ctest -L unit --output-on-failure


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

time cmake -DQMC_COMPLEX=1 -DQMC_MIXED_PRECISION=1 -DENABLE_SOA=1 -DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpicxx" -DCMAKE_CXX_FLAGS="-mno-bmi2 -mno-avx2" -DCMAKE_C_FLAGS="-mno-bmi2 -mno-avx2" -DBLAS_blas_LIBRARY="/usr/lib64/libblas.so.3" -DLAPACK_lapack_LIBRARY="/usr/lib64/atlas/liblapack.so.3" -DHDF5_INCLUDE_DIR="/sw/rhea/hdf5/1.8.11/rhel6.6_intel14.0.4/include" -DQMC_CUDA=1 /tmp/${BUILD_TAG}-src 2>&1 | tee cmake.out

time make -j 24
time make -j 24

time ctest -L unit --output-on-failure

# final cleanup of tmpfs
cd ../
rm -rf ./${BUILD_TAG}-build
rm -rf ./${BUILD_TAG}-src

EOF

/home/mat151ci_auser/blocking_qsub $BUILD_DIR $BUILD_TAG.pbs

cp $BUILD_DIR/$BUILD_TAG.o* ../

## this end of job logic could probably be more elegant
## hacks to get us going

cp $BUILD_DIR/$BUILD_TAG.o* ../

# check for correct test output from all builds
if [ $(grep '100% tests passed, 0 tests failed out of [0-9]*' ../$BUILD_TAG.o* | wc -l) -ne 4 ]
then
   echo; echo
   echo One or more build variants failed. Check the build log for details.
   echo; echo
fi

# set the return code for the script
[ $(grep '100% tests passed, 0 tests failed out of [0-9]*' ../$BUILD_TAG.o* | wc -l) -eq 4 ]
