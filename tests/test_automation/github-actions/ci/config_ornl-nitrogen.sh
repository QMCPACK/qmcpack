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

# Recent CMake and LLVM are required to build the latest version of QMCPACK
# CMake installer
cd $PREFIX && mkdir -p cmake
cd $SOFT_TMP
CMAKE_VERSION=3.24.3
wget https://github.com/Kitware/CMake/releases/download/v$CMAKE_VERSION/cmake-$CMAKE_VERSION-linux-x86_64.sh
sh cmake-$CMAKE_VERSION-linux-x86_64.sh --skip-license --prefix=$PREFIX/cmake/$CMAKE_VERSION
export PATH=$PREFIX/cmake/$CMAKE_VERSION/bin:$PATH

# LLVM 
cd $SOFT_TMP
LLVM_VERSION=15.0.7
LLVM_TAG=llvmorg-$LLVM_VERSION
git clone --depth=1 --branch $LLVM_TAG https://github.com/llvm/llvm-project.git
cd llvm-project

cmake -S llvm -B build -G Ninja \
      -DLLVM_ENABLE_PROJECTS="clang;clang-tools-extra;lld" \
      -DLLVM_ENABLE_RUNTIMES="libcxx;libcxxabi;openmp" \
      -DCLANG_OPENMP_NVPTX_DEFAULT_ARCH=sm_70 \
      -DLIBOMPTARGET_NVPTX_COMPUTE_CAPABILITIES=37,60,70 \
      -DCMAKE_INSTALL_PREFIX=$HOME/opt/llvm/$LLVM_VERSION \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_C_COMPILER=/usr/bin/gcc \
      -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
      -DGCC_INSTALL_PREFIX=/usr 

cd build && ninja && ninja install
