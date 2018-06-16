#!/bin/bash

# Xray instrumentation is a feature of relatively modern llvm
# It is only supported on linux as of clang 4.0.1
# While nvcc seems incomplatible with llvm > 5.0.1 even with
# a rewritten llvm version string
# Xray leaves entrypoints in the instrumented code that at runtime
# can be patched to call profiling routines in the compiler_rt library.
# It appears to be much less of a runtime drag that many profilers.

# Requirements:
# It should be sufficient to install spack packages loaded
# below with the same specs.
# Additionally if you build llvm with your system gcc you need libatomic installed.
# this should be done with the same package manager that installed you system gcc

# DMKL_ROOT should be defined to your mkl path
# DCUDA_TOOLKIT_ROOTDIR should point at the root of your CUDA install
# DCUDA_NVCC_FLAGS should have -arch appropriate to your CUDA version and hardware
#
# nvcc and libc++ don't get along hence -stdlib=libstdc++
# It is both a compiling and linking flag

spack load mpich%clang@4.0.1
spack load hdf5%clang@4.0.1
spack load ninja
spack load llvm@4.0.1 ^ncurses+termlib
spack load zlib%clang@4.0.1
 
CXXFLAGS="-stdlib=libstdc++" LDFLAGS="-stdlib=libstdc++" cmake -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH} -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_C_COMPILER=mpicc -DENABLE_MKL=1 -DMKL_ROOT="/opt/intel2018/mkl" -GNinja -DQMC_MPI=1 -DQMC_CUDA=1 -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-9.1 -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DCUDA_NVCC_FLAGS="-arch=sm_60;-Drestrict=__restrict__;-DNO_CUDA_MAIN;-O3" -DHDF5_ROOT=/data/epd/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/hdf5-1.10.2-qqmot24bg6uetn3xhpxjlwafvxr4p5pp/ -DXRAY_PROFILE=1 -DXRAY_INSTRUCTION_THRESHOLD=50 -DXRAY_GPU_MOST=1 ..

ninja

#you will see numerous warnings about loop vectorizations.
