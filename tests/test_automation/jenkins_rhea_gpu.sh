#!/bin/bash -x

BUILD_DIR=$(pwd)
echo $BUILD_DIR

cat > $BUILD_TAG.pbs << EOF

EOF

cp $BUILD_TAG.pbs $BUILD_DIR

cd $BUILD_DIR

source scl_source enable rh-python35
which python

$BUILD_DIR/../../../scripts/blocking_qsub.py $BUILD_DIR $BUILD_TAG.pbs

## this end of job logic could probably be more elegant
## hacks to get us going

cp $BUILD_DIR/$BUILD_TAG.o* ../

# explicitly check for correct test output from all builds
[ $(grep '100% tests passed, 0 tests failed out of [0-9]*' ../$BUILD_TAG.o* | wc -l) -eq 4 ]
