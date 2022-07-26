#!/bin/bash
#
# Setup for Ye-Ryzen-Box
#
# Run the "short" nightlies
#
 
export TEST_SITE_NAME=Ye-Ryzen-Box
export N_PROCS_BUILD=16
export N_PROCS=16

export PATH=/opt/cmake/current/bin:$PATH
source /usr/share/lmod/lmod/init/bash
export MODULEPATH=/home/packages/modules:$MODULEPATH

export BOOST_ROOT=/home/packages/math-libraries/boost/current
export FFTW_HOME=/home/packages/math-libraries/fftw/current

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

for sys in GCC9-Real GCC9-Real-Mixed GCC9-Complex GCC9-Complex-Mixed \
           AOMP-Offload-Real AOMP-Offload-Complex AOMP-Offload-Real-Mixed AOMP-Offload-Complex-Mixed
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

CTEST_FLAGS="-D CMAKE_PREFIX_PATH=/home/packages/math-libraries/OpenBLAS/current -D QMC_DATA=/home/yeluo/opt/QMCDATA"

# compiler dependent options
if [[ $sys == *"AOMP"* ]]; then
  #define and load compiler
  module load aomp/aomp-repo
  export CC=mpicc
  export CXX=mpicxx

  if [[ $sys == *"Offload"* ]]; then
    CTEST_FLAGS="$CTEST_FLAGS -DQMC_OPTIONS='-DENABLE_OFFLOAD=ON;-DOFFLOAD_TARGET=amdgcn-amd-amdhsa;-DOFFLOAD_ARCH=gfx906;-DQMC_PERFORMANCE_NIO_MAX_ATOMS=64'"
  else
    CTEST_FLAGS="$CTEST_FLAGS -DQMC_OPTIONS='-DQMC_PERFORMANCE_NIO_MAX_ATOMS=64'"
  fi

  CTEST_LABELS="-L 'deterministic|performance' -LE unstable"
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
