#!/bin/bash

set -x
HOST_NAME=$(hostname -s)

case "$1" in 

  # Configure qmcpack using cmake out-of-source builds 
  configure)
  
    echo "Use recent CMake v3.26.3"
    export PATH=$HOME/opt/cmake/3.26.3/bin:$PATH
    # Make current environment variables available to subsequent steps, ctest
    echo "PATH=$PATH" >> $GITHUB_ENV

    QMC_DATA_DIR=/scratch/ci/QMC_DATA_FULL

    if [ -d ${GITHUB_WORKSPACE}/../qmcpack-build ]
    then
      echo "Found existing out-of-source build directory ${GITHUB_WORKSPACE}/../qmcpack-build, removing"
      rm -fr ${GITHUB_WORKSPACE}/../qmcpack-build
    fi

    echo "Creating new out-of-source build directory ${GITHUB_WORKSPACE}/../qmcpack-build"
    cd ${GITHUB_WORKSPACE}/.. && mkdir qmcpack-build && cd qmcpack-build
    
    # Build variants
    # Real or Complex configuration
    case "${GH_JOBNAME}" in
      *"Real"*)
        echo 'Configure for real build -DQMC_COMPLEX=0'
        IS_COMPLEX=0
      ;;
      *"Complex"*)
        echo 'Configure for complex build -DQMC_COMPLEX=1'
        IS_COMPLEX=1
      ;; 
    esac

    # Mixed or Non-Mixed (default, full) precision, used with GPU code
    case "${GH_JOBNAME}" in
      *"Mixed"*)
        echo 'Configure for mixed precision build -DQMC_MIXED_PRECISION=1'
        IS_MIXED_PRECISION=1
      ;; 
      *)
        IS_MIXED_PRECISION=0
      ;;
    esac

    case "${GH_JOBNAME}" in
      *"GCC11-NoMPI-MKL-"*)
        echo 'Configure for building with GCC and Intel MKL'

        source /opt/intel/oneapi/setvars.sh

        cmake -GNinja \
              -DBLA_VENDOR=Intel10_64lp \
              -DQMC_MPI=0 \
              -DQMC_COMPLEX=$IS_COMPLEX \
              -DQMC_MIXED_PRECISION=$IS_MIXED_PRECISION \
              -DCMAKE_BUILD_TYPE=RelWithDebInfo \
              -DQMC_DATA=$QMC_DATA_DIR \
              ${GITHUB_WORKSPACE}

      ;;
    esac
    ;;  

  build)
    cd ${GITHUB_WORKSPACE}/../qmcpack-build
    ninja
    ;;
   
  test)
    source /opt/intel/oneapi/setvars.sh
    echo "Running deterministic tests"
    cd ${GITHUB_WORKSPACE}/../qmcpack-build
    ctest --output-on-failure -L deterministic -j 32
    ;;
    
esac
