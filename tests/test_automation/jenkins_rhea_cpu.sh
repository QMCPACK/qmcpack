#!/bin/bash -x

BUILD_DIR=$(pwd)
echo $BUILD_DIR

cat > $BUILD_TAG.pbs << EOF
#PBS -A MAT151
#PBS -N $BUILD_TAG
#PBS -j oe
#PBS -l walltime=1:00:00,nodes=1
#PBS -d $BUILD_DIR
#PBS -l partition=rhea

cd $BUILD_DIR

source /sw/rhea/environment-modules/3.2.10/rhel6.7_gnu4.4.7/init/bash

module swap PE-intel PE-gnu/5.3.0-1.10.2
module load fftw
export FFTW_HOME=\$FFTW3_DIR
module load hdf5
module load git

env

module list

mkdir -p build

cd build 

cmake -DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpicxx" -DCMAKE_CXX_FLAGS="-std=c++11" -DBLAS_blas_LIBRARY="/usr/lib64/libblas.so.3" -DLAPACK_lapack_LIBRARY="/usr/lib64/atlas/liblapack.so.3" -DHDF5_INCLUDE_DIR="/sw/rhea/hdf5/1.8.11/rhel6.6_gnu4.8.2/include" ..

make -j 24

ctest -L unit
#ctest -R short-LiH_dimer_ae-vmc_hf_noj-16-1

EOF

cp $BUILD_TAG.pbs $BUILD_DIR

cd $BUILD_DIR

source scl_source enable rh-python35
which python 

$BUILD_DIR/../../../scripts/blocking_qsub.py $BUILD_DIR $BUILD_TAG.pbs

## this end of job logic could probably be more elegant
## hacks to get us going

cp $BUILD_DIR/$BUILD_TAG.o* ../

# explicitly check for correct test output

grep '100% tests passed, 0 tests failed out of [0-9]*' ../$BUILD_TAG.o*
