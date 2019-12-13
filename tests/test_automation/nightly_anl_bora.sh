#!/bin/bash
#
# Setup for bora.alcf.anl.gov
#
# Run the "short" nightlies
# 

# load necessary modules
source /etc/profile.d/z00_lmod.sh
if [ -d /scratch/packages/modulefiles ]; then
  module use /scratch/packages/modulefiles
fi

module load intel-mkl intel/18.4
module load openmpi/4.0.2-intel
module load cuda/10.1

export TEST_SITE_NAME=bora.alcf.anl.gov
export N_PROCS=16
export N_PROCS_BUILD=16
export N_CONCURRENT_TESTS=16
export CC=mpicc
export CXX=mpicxx

# run on socket 1
NUMA_ID=1

QE_BIN=/scratch/opt/qe-stable/qe-6.4.1/bin
QMC_DATA=/scratch/opt/h5data

#Must be an absolute path
place=/scratch/QMCPACK_CI_BUILDS_DO_NOT_REMOVE

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

for sys in Real-SoA Real-Mixed-SoA Complex-SoA Complex-Mixed-SoA Real Real-Mixed Complex Complex-Mixed \
           Real-Mixed-SoA-CUDA2 Complex-Mixed-SoA-CUDA2
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

CTEST_FLAGS="-D QMC_DATA=$QMC_DATA -D ENABLE_TIMERS=1 -D C_FLAGS=-xCOMMON-AVX512 -D CXX_FLAGS=-xCOMMON-AVX512"

if [[ $sys == *"-CUDA2"* ]]; then
  CTEST_FLAGS="$CTEST_FLAGS -D ENABLE_CUDA=1 -D CUDA_ARCH=sm_61 -L 'deterministic|performance' -LE unstable"
else
  CTEST_FLAGS="$CTEST_FLAGS -D QE_BIN=$QE_BIN"
fi

if [[ $sys == *"Complex"* ]]; then
  CTEST_FLAGS="$CTEST_FLAGS -D QMC_COMPLEX=1"
fi

if [[ $sys == *"-SoA"* ]]; then
  CTEST_FLAGS="$CTEST_FLAGS -D ENABLE_SOA=1"
else
  CTEST_FLAGS="$CTEST_FLAGS -D ENABLE_SOA=0"
fi

if [[ $sys == *"-Mixed"* ]]; then
  CTEST_FLAGS="$CTEST_FLAGS -D QMC_MIXED_PRECISION=1"
fi

export QMCPACK_TEST_SUBMIT_NAME=${compiler}-${sys}-Release

numactl -N $NUMA_ID \
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
