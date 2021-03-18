#!/bin/bash

case "$1" in 

  # Configure qmcpack using cmake out-of-source builds 
  configure)
   
    cd ${GITHUB_WORKSPACE}/..
    mkdir qmcpack-build
    cd qmcpack-build
    cmake -GNinja -DCMAKE_INSTALL_PREFIX=${GITHUB_WORKSPACE}/../qmcpack-install ${GITHUB_WORKSPACE}
    ;;

  # Build using ninja
  build)
    cd ${GITHUB_WORKSPACE}/../qmcpack-build
    ninja
    ;;

  # Run unit and deterministic tests
  test)
    cd ${GITHUB_WORKSPACE}/../qmcpack-build
    ctest -L unit
    ctest -L deterministic
    ;;

  # Install the library
  install)
    cd ${GITHUB_WORKSPACE}/../qmcpack-build
    ninja install
    ;;

  *)
    echo " Invalid step" "$1"
    exit -1
    ;;
esac
