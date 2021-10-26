#!/bin/bash
#SBATCH -A mat189
#SBATCH -J nightly_spock
#SBATCH -o nightly_spock.%j
#SBATCH -e nightly_spock.%j
#SBATCH -t 00:15:00
#SBATCH -p ecp
#SBATCH -N 1

base_dir=/gpfs/alpine/proj-shared/mat189/wgodoy/nightly_olcf_spock

cd ${base_dir}

# setup proxy
export https_proxy=http://proxy.ccs.ornl.gov:3128/

# create log file directory if it doesn't exist
log_dir=${base_dir}/log/$(date +%y_%m_%d)
mkdir -p ${log_dir}

# module load llvm # has a conflict with serial cray-hdf5/1.12.0.6
module load gcc
module load cmake
module load git

# QMCPACK dependencies
module load libxml2/2.9.12
module load fftw
module load cray-hdf5/1.12.0.6 # require serial hdf5
module load boost
module load openblas
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
echo "Start gcc10-nompi-cuda2hip-release-real test"
export QMCPACK_TEST_SUBMIT_NAME=gcc10-nompi-cuda2hip-release-real

CTEST_FLAGS="-DCMAKE_C_COMPILER=gcc \
      -DCMAKE_CXX_COMPILER=g++ \
      -DQMC_MPI=0 \
      -DENABLE_CUDA=ON \
      -DQMC_CUDA2HIP=ON \
      -DQMC_COMPLEX=0"

ctest ${CTEST_FLAGS} \
      -S $(pwd)/../CMake/ctest_script.cmake,release \
      --stop-time $(date --date=now+15mins +%H:%M:%S) \
      -VV -L 'deterministic' --timeout 600 &> \
      ${log_dir}/${QMCPACK_TEST_SUBMIT_NAME}.log

unset QMCPACK_TEST_SUBMIT_NAME

# Start complex build test
echo "Start gcc10-nompi-cuda2hip-release-complex test"
export QMCPACK_TEST_SUBMIT_NAME=gcc10-nompi-cuda2hip-release-complex

cd ${base_dir}/qmcpack/build
rm -fr *

CTEST_FLAGS="-DCMAKE_C_COMPILER=gcc \
      -DCMAKE_CXX_COMPILER=g++ \
      -DQMC_MPI=0 \
      -DENABLE_CUDA=ON \
      -DQMC_CUDA2HIP=ON \
      -DQMC_COMPLEX=1"

ctest ${CTEST_FLAGS} \
      -S $(pwd)/../CMake/ctest_script.cmake,release \
      --stop-time $(date --date=now+15mins +%H:%M:%S) \
      -VV -L 'deterministic' --timeout 600 &> \
      ${log_dir}/${QMCPACK_TEST_SUBMIT_NAME}.log

unset QMCPACK_TEST_SUBMIT_NAME

# make it recursive to run nightly
# create next log file folder
log_dir=${base_dir}/log/$(date --date="tomorrow" +%y_%m_%d)
mkdir -p ${log_dir}

# go there to capture the job output file
cd ${log_dir}

echo "Submit sbatch job for 22:00"
sbatch -b 22:00 ${base_dir}/qmcpack/tests/test_automation/nightly_test_scripts/nightly_olcf_spock.sh 

