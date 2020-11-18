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

module load cmake
module load intel-mkl
module load cuda/10.1

export TEST_SITE_NAME=bora.alcf.anl.gov
export N_PROCS=16
export N_PROCS_BUILD=16

# run on socket 1
NUMA_ID=0

QE_BIN=/scratch/opt/qe-stable/qe-6.4.1/bin
QMC_DATA=/scratch/opt/h5data

#Must be an absolute path
place=/scratch/QMCPACK_CI_BUILDS_DO_NOT_REMOVE

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

for sys in ClangDev-Offload-CUDA-Real ClangDev-Offload-CUDA-Complex ClangDev-Offload-CUDA-Real-Mixed ClangDev-Offload-CUDA-Complex-Mixed \
           Intel18-Real Intel18-Real-Mixed Intel18-Complex Intel18-Complex-Mixed \
           Intel18-Real-Mixed-CUDA2 Intel18-Complex-Mixed-CUDA2 Intel18-Real-Mixed-legacy-CUDA Intel18-Complex-Mixed-legacy-CUDA
do

echo --- Building for $sys `date`

# create log file folder if not exist
mydate=`date +%y_%m_%d`
if [ ! -e $place/log/$entry/$mydate ];
then
  mkdir -p $place/log/$entry/$mydate
fi

# options common to all cases
CTEST_FLAGS="-DQMC_DATA=$QMC_DATA"

# compiler dependent options
if [[ $sys == *"ClangDev"* ]]; then
  #define and load compiler
  module load llvm/master-latest
  module load openmpi/4.0.2-llvm
  export CC=mpicc
  export CXX=mpicxx

  CTEST_FLAGS="$CTEST_FLAGS -DCMAKE_C_FLAGS=-march=skylake -DCMAKE_CXX_FLAGS=-march=skylake"
  if [[ $sys == *"Offload-CUDA"* ]]; then
    CTEST_FLAGS="$CTEST_FLAGS -DQMC_OPTIONS='-DENABLE_OFFLOAD=ON;-DUSE_OBJECT_TARGET=ON;-DCUDA_HOST_COMPILER=`which gcc`;-DCUDA_PROPAGATE_HOST_FLAGS=OFF;-DQMC_NIO_MAX_SIZE=16'"
    CTEST_FLAGS="$CTEST_FLAGS -DENABLE_CUDA=1 -DCUDA_ARCH=sm_61"
    CTEST_LABELS="-L deterministic -LE unstable"
    export N_CONCURRENT_TESTS=4
  else
    CTEST_FLAGS="$CTEST_FLAGS -DQE_BIN=$QE_BIN"
    CTEST_LABELS="-L 'deterministic|performance' -LE unstable"
    export N_CONCURRENT_TESTS=16
  fi
elif [[ $sys == *"Intel"* ]]; then
  #define and load compiler
  module load intel/18.4
  module load openmpi/4.0.2-intel
  export CC=mpicc
  export CXX=mpicxx

  CTEST_FLAGS="$CTEST_FLAGS -DCMAKE_C_FLAGS=-xCOMMON-AVX512 -DCMAKE_CXX_FLAGS=-xCOMMON-AVX512"
  if [[ $sys == *"-CUDA2"* ]]; then
    CTEST_FLAGS="$CTEST_FLAGS -DENABLE_CUDA=1 -DCUDA_ARCH=sm_61 -DQMC_OPTIONS='-DQMC_NIO_MAX_SIZE=16'"
    CTEST_LABELS="-L 'deterministic|performance' -LE unstable"
    export N_CONCURRENT_TESTS=4
  elif [[ $sys == *"-legacy-CUDA"* ]]; then
    CTEST_FLAGS="$CTEST_FLAGS -DQMC_CUDA=1 -DCUDA_ARCH=sm_61 -DQMC_OPTIONS='-DQMC_NIO_MAX_SIZE=16'"
    CTEST_LABELS="-L 'deterministic|performance' -LE unstable"
    export N_CONCURRENT_TESTS=1
  else
    CTEST_FLAGS="$CTEST_FLAGS -DQE_BIN=$QE_BIN"
    export N_CONCURRENT_TESTS=16
  fi
fi

# compiler independent options
if [[ $sys == *"-Complex"* ]]; then
  CTEST_FLAGS="$CTEST_FLAGS -DQMC_COMPLEX=1"
fi

if [[ $sys == *"-Mixed"* ]]; then
  CTEST_FLAGS="$CTEST_FLAGS -DQMC_MIXED_PRECISION=1"
fi

export QMCPACK_TEST_SUBMIT_NAME=${sys}-Release

folder=build_${sys}
if [ -e $folder ]; then
rm -r $folder
fi
mkdir $folder
cd $folder

numactl -N $NUMA_ID \
ctest $CTEST_FLAGS $CTEST_LABELS -S $PWD/../CMake/ctest_script.cmake,release -VV -E 'long' --timeout 800 &> $place/log/$entry/$mydate/${QMCPACK_TEST_SUBMIT_NAME}.log

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
