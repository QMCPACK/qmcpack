//////////////////////////////////////////////////////////////////////////////////////
//// This file is distributed under the University of Illinois/NCSA Open Source License.
//// See LICENSE file in top directory for details.
////
//// Copyright (c) 2019 QMCPACK developers.
////
//// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
////
//// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
////////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_CUDA_ERROR_H
#define QMCPLUSPLUS_CUDA_ERROR_H

#ifdef QMC_CUDA
#ifndef QMC_CUDA2HIP
#include <cuda_runtime_api.h>
#else
#include <hip/hip_runtime.h>
#include "Platforms/ROCm/cuda2hip.h"
#endif
#endif

#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>

#define cudaErrorCheck(ans, cause)                \
  {                                               \
    cudaAssert((ans), cause, __FILE__, __LINE__); \
  }

// If the cause is largely redundant with the __FILE__ and __LINE__ information
#define cudaCheck(ans)                         \
  {                                            \
    cudaAssert((ans), "", __func__, __FILE__, __LINE__); \
  }

/// prints CUDA error messages. Always use cudaErrorCheck macro.
inline void cudaAssert(cudaError_t code, const std::string& cause, const char* function, const char* filename, int line, bool abort = true)
{
  if (code != cudaSuccess)
  {
    std::ostringstream err;
    err << "cudaAssert: " << cudaGetErrorName(code) << " " << cudaGetErrorString(code)
        << ", function " << function
        << ", file " << filename
        << ", line " << line << std::endl
        << cause << std::endl;
    std::cerr << err.str();
    if (abort)
      throw std::runtime_error(cause);
  }
}

#endif
