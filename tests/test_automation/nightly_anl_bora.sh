#!/bin/bash
#
# Setup for bora.alcf.anl.gov
#
# Run the "short" nightlies
# 

export N_PROCS=32
export CC=mpicc
export CXX=mpicxx
export BOOST_ROOT=/sandbox/opt/qmcdev/trunk/external_codes/boost_1_55_0

QE_BIN=/sandbox/opt/qe-6.2.1/bin
QMC_DATA=/sandbox/yeluo/benchmark
site_name=bora.alcf.anl.gov

#Must be an absolute path
place=/sandbox/QMCPACK_CI_BUILDS_DO_NOT_REMOVE

#define and load compiler
compiler=Intel2018

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

CTEST_FLAGS="-D QE_BIN=$QE_BIN -D QMC_DATA=$QMC_DATA -D CTEST_SITE=$site_name"

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

ctest $CTEST_FLAGS -S $PWD/../CMake/ctest_script.cmake,release -VV -E 'long' --timeout 600 &> $place/log/$entry/$mydate/${QMCPACK_TEST_SUBMIT_NAME}.log

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
