#!/bin/bash
#SBATCH -A mat189
#SBATCH -J nightly_spock
#SBATCH -o nightly_spock.%j
#SBATCH -e nightly_spock.%j
#SBATCH -t 00:40:00
#SBATCH -p ecp
#SBATCH -N 1

base_dir=/gpfs/alpine/proj-shared/mat189/wgodoy/nightly_olcf_spock
qmc_data_dir=/gpfs/alpine/mat189/proj-shared/qmc_data/Benchmark

cd ${base_dir}

# setup proxy
export https_proxy=http://proxy.ccs.ornl.gov:3128/

# create log file directory if it doesn't exist
log_dir=${base_dir}/log/$(date +%y_%m_%d)
mkdir -p ${log_dir}

# required to load openblas
module load PrgEnv-gnu
# module load llvm # has a conflict with serial cray-hdf5
module load gcc
module load cmake
module load git

# QMCPACK dependencies
module load libxml2/2.9.12
module load fftw
module load cray-hdf5 # require serial hdf5
module load boost
module load openblas/0.3.17-omp
module load rocm

module list

export TEST_SITE_NAME=spock.olcf.ornl.gov
export N_PROCS_BUILD=64
export N_CONCURRENT_TESTS=1


# implement git clone logic
mkdir -p qmcpack
rm -fr qmcpack
git clone --branch develop --depth 1 https://github.com/QMCPACK/qmcpack.git
cd qmcpack/build

# Start real build test
now=$(date +"%T")
echo "Start GCC10-NoMPI-CUDA2HIP-Release-Real test ${now}"
export QMCPACK_TEST_SUBMIT_NAME=GCC10-NoMPI-CUDA2HIP-Real-Release

CTEST_FLAGS="-DCMAKE_C_COMPILER=gcc \
      -DCMAKE_CXX_COMPILER=g++ \
      -DQMC_MPI=0 \
      -DENABLE_CUDA=ON \
      -DQMC_CUDA2HIP=ON \
      -DQMC_COMPLEX=0 \
      -DQMC_OPTIONS='-DQMC_DATA=${qmc_data_dir};-DQMC_PERFORMANCE_NIO_MAX_ATOMS=32'"

ctest ${CTEST_FLAGS} \
      -S $(pwd)/../CMake/ctest_script.cmake,release \
      --stop-time $(date --date=now+20mins +%H:%M:%S) \
      -VV -R 'deterministic|performance-NiO' --timeout 600 &> \
      ${log_dir}/${QMCPACK_TEST_SUBMIT_NAME}.log

unset QMCPACK_TEST_SUBMIT_NAME

# Start complex build test
now=$(date +"%T")
echo "Start GCC10-NoMPI-CUDA2HIP-Release-Complex test ${now}"
export QMCPACK_TEST_SUBMIT_NAME=GCC10-NoMPI-CUDA2HIP-Complex-Release

cd ${base_dir}/qmcpack/build
rm -fr *

CTEST_FLAGS="-DCMAKE_C_COMPILER=gcc \
      -DCMAKE_CXX_COMPILER=g++ \
      -DQMC_MPI=0 \
      -DENABLE_CUDA=ON \
      -DQMC_CUDA2HIP=ON \
      -DQMC_COMPLEX=1 \
      -DQMC_OPTIONS='-DQMC_DATA=${qmc_data_dir};-DQMC_PERFORMANCE_NIO_MAX_ATOMS=32'"

ctest ${CTEST_FLAGS} \
      -S $(pwd)/../CMake/ctest_script.cmake,release \
      --stop-time $(date --date=now+20mins +%H:%M:%S) \
      -VV -R 'deterministic|performance-NiO' --timeout 600 &> \
      ${log_dir}/${QMCPACK_TEST_SUBMIT_NAME}.log

unset QMCPACK_TEST_SUBMIT_NAME

# make it recursive to run nightly
# create next log file folder
log_dir=${base_dir}/log/$(date --date="tomorrow" +%y_%m_%d)
mkdir -p ${log_dir}

# go there to capture the job output file
cd ${log_dir}

echo "Submit sbatch job for 22:00 nightly"
sbatch -b 22:00 ${base_dir}/qmcpack/tests/test_automation/nightly_test_scripts/nightly_olcf_spock.sh 

