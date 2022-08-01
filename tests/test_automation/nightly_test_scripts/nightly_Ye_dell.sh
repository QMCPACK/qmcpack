#!/bin/bash
#
# Setup for dell-laptop with gen9
#
# Run the "short" nightlies
#
 
export TEST_SITE_NAME=Ye-Dell-Laptop
export N_PROCS_BUILD=12
export N_PROCS=12

#Must be an absolute path
place=/home/QMCPACK_NIGHTLY_BUILDS

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

for sys in OneAPI-Offload-Real OneAPI-Offload-Real-Mixed OneAPI-Real OneAPI-Real-Mixed
do

folder=build_$sys

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

CTEST_FLAGS=""

# compiler dependent options
if [[ $sys == *"OneAPI"* ]]; then
  #define and load compiler
  export LIBOMPTARGET_PLUGIN=OPENCL
  source /opt/intel/oneapi/setvars.sh 
  export CC=icx
  export CXX=icpx

  CTEST_FLAGS="$CTEST_FLAGS -DQMC_MPI=0"
  if [[ $sys == *"Offload"* ]]; then
    CTEST_FLAGS="$CTEST_FLAGS -DQMC_OPTIONS='-DENABLE_OFFLOAD=ON;-DOFFLOAD_TARGET=spir64;-DQMC_PERFORMANCE_NIO_MAX_ATOMS=64'"
  else
    CTEST_FLAGS="$CTEST_FLAGS -DQMC_OPTIONS='-DQMC_PERFORMANCE_NIO_MAX_ATOMS=64'"
  fi

  CTEST_LABELS="-L 'deterministic' -LE unstable"
  export N_CONCURRENT_TESTS=1
elif [[ $sys == *"GCC"* ]]; then
  export CC=mpicc
  export CXX=mpicxx

  CTEST_FLAGS="$CTEST_FLAGS -DQMC_OPTIONS='-DQMC_PERFORMANCE_NIO_MAX_ATOMS=64'"
  CTEST_LABELS="-L 'deterministic' -LE unstable"
  export N_CONCURRENT_TESTS=16
fi

if [[ $sys == *"Complex"* ]]; then
  CTEST_FLAGS="$CTEST_FLAGS -D QMC_COMPLEX=1"
fi

if [[ $sys == *"-Mixed"* ]]; then
  CTEST_FLAGS="$CTEST_FLAGS -D QMC_MIXED_PRECISION=1"
fi

export QMCPACK_TEST_SUBMIT_NAME=${sys}-Release

ctest -DCMAKE_C_FLAGS="$CMAKE_C_FLAGS" -DCMAKE_CXX_FLAGS="$CMAKE_CXX_FLAGS" \
      $CTEST_FLAGS $CTEST_LABELS -S $PWD/../CMake/ctest_script.cmake,release -VV --timeout 600 &> $place/log/$entry/$mydate/${QMCPACK_TEST_SUBMIT_NAME}.log

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
