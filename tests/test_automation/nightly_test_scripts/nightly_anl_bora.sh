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

module use /nfs/gce/projects/QMCPACK_dev/spack/share/spack/lmod/linux-ubuntu18.04-x86_64/Core
module load cuda/11.2.2-lo3x6k7

export TEST_SITE_NAME=bora.alcf.anl.gov
export N_PROCS=16
export N_PROCS_BUILD=16
export N_CONCURRENT_TESTS=16

# run on socket 1
NUMA_ID=0

QE_BIN=/scratch/opt/qe-stable/qe-6.4.1/bin
QMC_DATA=/scratch/opt/h5data

#Must be an absolute path
place=/scratch2/QMCPACK_CI_BUILDS_DO_NOT_REMOVE

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

for sys in ROCm-CUDA2HIP-Real ROCm-CUDA2HIP-Complex ROCm-CUDA2HIP-Real-Mixed ROCm-CUDA2HIP-Complex-Mixed \
           Intel19-Real Intel19-Real-Mixed Intel19-Complex Intel19-Complex-Mixed \
           ClangDev-Offload-Real ClangDev-Offload-Complex ClangDev-Offload-Real-Mixed ClangDev-Offload-Complex-Mixed \
           ClangDev-Offload-CUDA-Real ClangDev-Offload-CUDA-Complex ClangDev-Offload-CUDA-Real-Mixed ClangDev-Offload-CUDA-Complex-Mixed \
           Intel19-CUDA2-Real-Mixed Intel19-CUDA2-Complex-Mixed Intel19-legacy-CUDA-Real-Mixed Intel19-legacy-CUDA-Complex-Mixed
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

  if [[ $sys == *"Offload-CUDA"* ]]; then
    CTEST_FLAGS="$CTEST_FLAGS -DQMC_OPTIONS='-DENABLE_OFFLOAD=ON;-DUSE_OBJECT_TARGET=ON;-DQMC_PERFORMANCE_NIO_MAX_ATOMS=32'"
    CTEST_FLAGS="$CTEST_FLAGS -DENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=61 -DBUILD_AFQMC=ON"
  elif [[ $sys == *"Offload"* ]]; then
    CTEST_FLAGS="$CTEST_FLAGS -DQMC_OPTIONS='-DENABLE_OFFLOAD=ON;-DOFFLOAD_ARCH=sm_61;-DUSE_OBJECT_TARGET=ON;-DQMC_PERFORMANCE_NIO_MAX_ATOMS=32'"
    CTEST_FLAGS="$CTEST_FLAGS -DBUILD_AFQMC=ON"
  else
    CTEST_FLAGS="$CTEST_FLAGS -DQE_BIN=$QE_BIN"
  fi

  export LIBOMP_USE_HIDDEN_HELPER_TASK=OFF
  if [[ $sys == *"Real-Mixed"* ]]; then
    #CTEST_LABELS="-E 'long|short-LiH|short-bccH_1x1x1|short-H2O_dimer|short-diamondC_1x1x1|short-diamondC_2x1x1_pp-dmc|short-li2_sto|short-NiO_a4_e48_pp'"
    CTEST_LABELS="-E 'restart|save|long|developer|He_param|example|cpu_limit|short-LiH|short-bccH_1x1x1|short-H2O_dimer|short-diamondC_1x1x1|short-diamondC_2x1x1_pp-dmc|short-li2_sto|short-NiO_a4_e48_pp|short-C'"
  else
    CTEST_LABELS="-L 'deterministic|QMCPACK' -LE unstable -E long"
  fi
elif [[ $sys == *"ROCm"* ]]; then
  #define and load compiler
  module load rocm/afar001-432 aomp/afar001-432
  module load openmpi/4.0.2-llvm

  if [[ $sys == *"Offload-CUDA2HIP"* ]]; then
    CTEST_FLAGS="$CTEST_FLAGS -DQMC_OPTIONS='-DCMAKE_HIP_ARCHITECTURES=gfx1030;-DENABLE_OFFLOAD=ON;-DQMC_PERFORMANCE_NIO_MAX_ATOMS=32'"
    CTEST_FLAGS="$CTEST_FLAGS -DENABLE_CUDA=ON -DQMC_CUDA2HIP=ON"
  elif [[ $sys == *"CUDA2HIP"* ]]; then
    CTEST_FLAGS="$CTEST_FLAGS -DQMC_OPTIONS='-DCMAKE_HIP_ARCHITECTURES=gfx1030;-DQMC_PERFORMANCE_NIO_MAX_ATOMS=32'"
    CTEST_FLAGS="$CTEST_FLAGS -DENABLE_CUDA=ON -DQMC_CUDA2HIP=ON"
  elif [[ $sys == *"Offload"* ]]; then
    CTEST_FLAGS="$CTEST_FLAGS -DQMC_OPTIONS='-DENABLE_OFFLOAD=ON;-DOFFLOAD_TARGET=amdgcn-amd-amdhsa;-DOFFLOAD_ARCH=gfx1030;-DQMC_PERFORMANCE_NIO_MAX_ATOMS=32'"
  else
    CTEST_FLAGS="$CTEST_FLAGS -DQE_BIN=$QE_BIN"
  fi

  if [[ $sys == *"Real-Mixed"* ]]; then
    #CTEST_LABELS="-E 'long|short-LiH|short-bccH_1x1x1|short-H2O_dimer|short-diamondC_1x1x1|short-diamondC_2x1x1_pp-dmc|short-li2_sto|short-NiO_a4_e48_pp'"
    CTEST_LABELS="-E 'restart|save|long|developer|He_param|example|cpu_limit|short-LiH|short-bccH_1x1x1|short-H2O_dimer|short-diamondC_1x1x1|short-diamondC_2x1x1_pp-dmc|short-li2_sto|short-NiO_a4_e48_pp|short-C'"
  else
    CTEST_LABELS="-L 'deterministic|QMCPACK' -LE unstable -E long"
  fi
elif [[ $sys == *"Intel"* ]]; then
  #define and load compiler
  module load intel/20.2
  module load openmpi/4.0.2-intel

  CTEST_FLAGS="$CTEST_FLAGS -DCMAKE_C_FLAGS=-xCOMMON-AVX512 -DCMAKE_CXX_FLAGS=-xCOMMON-AVX512"
  if [[ $sys == *"-CUDA2"* ]]; then
    CTEST_FLAGS="$CTEST_FLAGS -DENABLE_CUDA=1 -DCMAKE_CUDA_ARCHITECTURES=61 -DQMC_OPTIONS='-DQMC_PERFORMANCE_NIO_MAX_ATOMS=64'"
    CTEST_LABELS="-L 'deterministic|performance' -LE unstable -E long"
    export N_CONCURRENT_TESTS=4
  elif [[ $sys == *"-legacy-CUDA"* ]]; then
    CTEST_FLAGS="$CTEST_FLAGS -DQMC_CUDA=1 -DCMAKE_CUDA_ARCHITECTURES=61 -DQMC_OPTIONS='-DQMC_PERFORMANCE_NIO_MAX_ATOMS=64'"
    CTEST_LABELS="-L 'deterministic|performance' -LE unstable -E long"
    export N_CONCURRENT_TESTS=1
  else
    CTEST_FLAGS="$CTEST_FLAGS -DQE_BIN=$QE_BIN"
    CTEST_LABELS="-E long"
    export N_CONCURRENT_TESTS=16
  fi
fi

CTEST_FLAGS="-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx $CTEST_FLAGS"

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

logfile=$place/log/$entry/$mydate/${QMCPACK_TEST_SUBMIT_NAME}.log
echo "$CTEST_FLAGS $CTEST_LABELS" | tee $logfile
echo PATH=$PATH >> $logfile
echo OMPI_CXX=${OMPI_CXX} >> $logfile
mpicxx -v >> $logfile
echo >> $logfile

numactl -N $NUMA_ID \
ctest $CTEST_FLAGS $CTEST_LABELS -S $PWD/../CMake/ctest_script.cmake,release -VV --timeout 800 &>> $logfile

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
