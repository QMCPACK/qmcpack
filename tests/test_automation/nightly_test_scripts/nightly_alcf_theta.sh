#!/bin/bash
#COBALT -q debug-cache-quad
#COBALT -A CSC249ADSE09
#COBALT -n 1
#COBALT -t 60
#COBALT -O nightly
#COBALT --attrs mcdram=cache:numa=quad

#
# Setup for theta.alcf.anl.gov
# 

module unload cray-libsci
module load cray-hdf5-parallel
module load gcc
module load cmake/3.14.5
module load miniconda-3/latest

export TEST_SITE_NAME=theta.alcf.anl.gov
export N_PROCS=16
export N_PROCS_BUILD=16
export N_CONCURRENT_TESTS=1

#QE_BIN=/scratch/opt/qe-stable/qe-6.4.1/bin
QMC_DATA=/projects/catalyst/yeluo/benchmark/h5data

#Must be an absolute path
place=/projects/CSC249ADSE09/QMCPACK_CI_BUILDS_DO_NOT_REMOVE

#define and load compiler
compiler=Intel19.0.5
CC=cc
CXX=CC

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

for sys in Real Complex # Real-Mixed Complex-Mixed
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

CTEST_FLAGS="-D CMAKE_C_COMPILER=$CC -D CMAKE_CXX_COMPILER=$CXX -D QMC_DATA=$QMC_DATA -D ENABLE_TIMERS=1
             -DQMC_OPTIONS='-DMPIEXEC_EXECUTABLE=/bin/sh;-DMPIEXEC_NUMPROC_FLAG=$place/$entry/tests/scripts/aprunhelper.sh'"

if [[ $sys == *"Complex"* ]]; then
  CTEST_FLAGS="$CTEST_FLAGS -D QMC_COMPLEX=1"
fi

if [[ $sys == *"-Mixed"* ]]; then
  CTEST_FLAGS="$CTEST_FLAGS -D QMC_MIXED_PRECISION=1"
fi

export QMCPACK_TEST_SUBMIT_NAME=${compiler}-${sys}-Release

ctest $CTEST_FLAGS -S $PWD/../CMake/ctest_script.cmake,release \
      -VV --stop-time $(date --date=now+28mins +%H:%M:%S) \
      -R 'deterministic|performance-NiO-cpu-a128-e1536' \
      --timeout 400 &> $place/log/$entry/$mydate/${QMCPACK_TEST_SUBMIT_NAME}.log

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
