#!/bin/bash

# This script is used to configure the ORNL nitrogen cluster for CI testing
# of the github-actions/ci workflow.
# Run manually to reproduce installation of dependencies in $PREFIX
# $SOFT_TMP keeps sources, installers, etc. Can be removed.

# create directories, change PREFIX and SOFT_TMP to a location of your choice
PREFIX=$HOME/opt
SOFT_TMP=$HOME/software
mkdir -p $PREFIX
mkdir -p $SOFT_TMP

# A recent CMake ( > 3.21) is required to build the latest version of QMCPACK
# CMake installer
cd $PREFIX && mkdir -p cmake
cd $SOFT_TMP
CMAKE_VERSION=3.24.3
wget https://github.com/Kitware/CMake/releases/download/v$CMAKE_VERSION/cmake-$CMAKE_VERSION-linux-x86_64.sh
sh cmake-$CMAKE_VERSION-linux-x86_64.sh --skip-license --prefix=$PREFIX/cmake/$CMAKE_VERSION
export PATH=$PREFIX/cmake/$CMAKE_VERSION/bin:$PATH
