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

module unload PE-intel
module load PE-gnu/5.3.0-1.10.2
module load fftw
export FFTW_HOME=\$FFTW3_DIR
module load hdf5
module load git

env
module list


echo ""
echo ""
echo "starting new test for real full precision"
echo ""
echo ""

sleep 10
EOF

/home/bgl/blocking_qsub $BUILD_DIR $BUILD_TAG.pbs

cp $BUILD_DIR/$BUILD_TAG.o* ../

## this end of job logic could probably be more elegant
## hacks to get us going

cp $BUILD_DIR/$BUILD_TAG.o* ../

# explicitly check for correct test output from all builds
[ $(grep '100% tests passed, 0 tests failed out of [0-9]*' ../$BUILD_TAG.o* | wc -l) -eq 4 ]
