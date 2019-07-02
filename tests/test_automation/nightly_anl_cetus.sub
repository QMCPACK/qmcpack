#!/bin/bash
#COBALT -q default
#COBALT -A PSFMat
#COBALT -n 128
#COBALT -t 60
#COBALT -O validation

#
# Setup for cetus.alcf.anl.gov
#
# Run the "short" nightlies, requeue if an executable is built
# 

export TEST_SITE_NAME=cetus.alcf.anl.gov
export N_PROCS_BUILD=24
export N_CONCURRENT_TESTS=1

#Must be an absolute path
place=/projects/PSFMat/QMCPACK_CI_BUILDS_DO_NOT_REMOVE

#define and load compiler
#compiler=XL
compiler=Clang

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

for sys in Real-SoA Real-Mixed-SoA Complex-SoA Complex-Mixed-SoA
do

folder=build_BGQ_$sys

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

CTEST_FLAGS="-DCMAKE_TOOLCHAIN_FILE=../config/BGQ_${compiler}_ToolChain.cmake"

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

ctest $CTEST_FLAGS -S $PWD/../CMake/ctest_script.cmake,release --stop-time `date --date=now+55mins +%H:%M:%S` -R "deterministic" -VV &> $place/log/$entry/$mydate/${QMCPACK_TEST_SUBMIT_NAME}.log

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
