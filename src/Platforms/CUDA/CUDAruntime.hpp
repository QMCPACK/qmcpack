//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/**@file CUDAruntime.hpp
 *@brief handle CUDA/HIP runtime selection.
 */
#ifndef QMCPLUSPLUS_CUDA_RUNTIME_H
#define QMCPLUSPLUS_CUDA_RUNTIME_H

#include <cstddef>
#include <string> // Positioned here to avoid conflict between CUDA and GCC >= 12 header files. https://github.com/QMCPACK/qmcpack/pull/4814
#include "config.h"
#ifndef QMC_CUDA2HIP
#include <cuda_runtime.h>
#else
#include <hip/hip_runtime.h>
#include "Platforms/ROCm/cuda2hip.h"
#endif

#define cudaErrorCheck(ans, cause)                \
  {                                               \
    cudaAssert((ans), cause, __FILE__, __LINE__); \
  }

// If the cause is largely redundant with the __FILE__ and __LINE__ information
#define cudaCheck(ans)                         \
  {                                            \
    cudaAssert((ans), "", __FILE__, __LINE__); \
  }

/// prints CUDA error messages. Always use cudaErrorCheck macro.
void cudaAssert(cudaError_t code, const std::string& cause, const char* filename, int line, bool abort = true);

size_t getCUDAdeviceFreeMem();

#endif
