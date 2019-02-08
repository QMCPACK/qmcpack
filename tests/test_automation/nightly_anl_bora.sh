#!/bin/bash
#
# Setup for bora.alcf.anl.gov
#
# Run the "short" nightlies
# 

export TEST_SITE_NAME=bora.alcf.anl.gov
export N_PROCS_BUILD=24
export N_PROCS=32
export CC=mpicc
export CXX=mpicxx
export BOOST_ROOT=/sandbox/opt/boost_1_61_0

QE_BIN=/sandbox/opt/qe-stable/qe-6.3/bin
QMC_DATA=/sandbox/opt/h5data

#Must be an absolute path
place=/sandbox/QMCPACK_CI_BUILDS_DO_NOT_REMOVE

#define and load compiler
compiler=Intel2019

if [ ! -e $place ]; then
mkdir $place
fi

if [ -e $place ]; then
cd $place

echo --- Hostname --- $HOSTNAME
echo --- Checkout for $sys `date`

branch=develop
entry=qmcpack-${branch}

if [ ! -e $entry ]; then
echo --- Cloning QMCPACK git `date`
git clone --depth 1 https://github.com/QMCPACK/qmcpack.git $entry
else
echo --- Updating local QMCPACK git `date`
cd $entry
git pull
cd ..
fi

if [ -e $entry/CMakeLists.txt ]; then
cd $entry

git checkout $branch

for sys in Real-SoA Real-Mixed-SoA Complex-SoA Complex-Mixed-SoA Real Real-Mixed Complex Complex-Mixed
do

folder=build_$compiler_$sys

if [ -e $folder ]; then
rm -r $folder
fi
mkdir $folder
cd $folder

echo --- Building for $sys `date`

# create log file folder if not exist
mydate=`date +%y_%m_%d`
if [ ! -e $place/log/$entry/$mydate ];
then
  mkdir -p $place/log/$entry/$mydate
fi

CTEST_FLAGS="-D QE_BIN=$QE_BIN -D QMC_DATA=$QMC_DATA -D ENABLE_TIMERS=1 -D C_FLAGS=-xCOMMON-AVX512 -D CXX_FLAGS=-xCOMMON-AVX512"

if [[ $sys == *"Complex"* ]]; then
  CTEST_FLAGS="$CTEST_FLAGS -D QMC_COMPLEX=1"
fi

if [[ $sys == *"-SoA"* ]]; then
  CTEST_FLAGS="$CTEST_FLAGS -D ENABLE_SOA=1"
fi

if [[ $sys == *"-Mixed"* ]]; then
  CTEST_FLAGS="$CTEST_FLAGS -D QMC_MIXED_PRECISION=1"
fi

export QMCPACK_TEST_SUBMIT_NAME=${compiler}-${sys}-Release

ctest $CTEST_FLAGS -S $PWD/../CMake/ctest_script.cmake,release -VV -E 'long' --timeout 800 &> $place/log/$entry/$mydate/${QMCPACK_TEST_SUBMIT_NAME}.log

cd ..
echo --- Finished $sys `date`
done

else
echo  "ERROR: No CMakeLists. Bad git clone."
exit 1
fi

else
echo "ERROR: No directory $place"
exit 1
fi
