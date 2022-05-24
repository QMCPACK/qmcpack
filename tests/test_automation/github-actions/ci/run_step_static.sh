#!/bin/bash

set -x
HOST_NAME=$(hostname -s)

case "$1" in 

  # Configure qmcpack using cmake out-of-source builds 
  configure)
    
    if [ -d ${GITHUB_WORKSPACE}/../qmcpack-build ]
    then
      echo "Found existing out-of-source build directory ${GITHUB_WORKSPACE}/../qmcpack-build, removing"
      rm -fr ${GITHUB_WORKSPACE}/../qmcpack-build
    fi
    
    echo "Creating new out-of-source build directory ${GITHUB_WORKSPACE}/../qmcpack-build"
    cd ${GITHUB_WORKSPACE}/..
    mkdir qmcpack-build
    cd qmcpack-build
    
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
      *"ClangTidy14-NoMPI-"*)
        echo 'Configure for debug mode'
        cmake -GNinja \
              -DCMAKE_C_COMPILER=clang \
              -DCMAKE_CXX_COMPILER=clang++ \
              -DQMC_MIXED_PRECISION=$IS_MIXED_PRECISION \
              -DCMAKE_BUILD_TYPE=Debug \
              -DCMAKE_CXX_CLANG_TIDY="clang-tidy" \
              -DQMC_COMPLEX=$IS_COMPLEX \
              ${GITHUB_WORKSPACE}
      ;;
    esac
    ;;

  # Build using ninja (~ 2 hours 15 minutes on GitHub-hosted runner)
  build)
    cd ${GITHUB_WORKSPACE}/../qmcpack-build
    ninja
    ;;
  *)
    echo " Invalid step" "$1"
    exit -1
    ;;
esac
