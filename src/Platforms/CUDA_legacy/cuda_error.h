//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
// Copyright(C) 2022 Advanced Micro Devices, Inc. All rights reserved.
//
// File developed by: Jakub Kurzak, jakurzak@amd.com, Advanced Micro Devices, Inc.
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jakub Kurzak, jakurzak@amd.com, Advanced Micro Devices, Inc.
//////////////////////////////////////////////////////////////////////////////////////


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

#define cudaCheck(call) \
{ \
    cudaError_t code = call; \
    if (code != cudaSuccess) \
      cudaThrow(#call, code, __func__, __FILE__, __LINE__); \
}

inline void cudaThrow(const char* call, cudaError_t code,
                      const char* func, const char* file, int line)
{
  char const* name = cudaGetErrorName(code);
  char const* string = cudaGetErrorString(code);
  std::cerr << call << " returned " << name << " (" << string << ")." << std::endl;
  std::cerr << "func: " << func << std::endl;
  std::cerr << "file: " << file << std::endl;
  std::cerr << "line: " << line << std::endl;
  throw std::runtime_error("");
}

#define cudaCheckMalloc(call, ...) \
{ \
    cudaError_t code = call; \
    if (code != cudaSuccess) \
      cudaThrowMalloc(#call, code, ##__VA_ARGS__, __func__, __FILE__, __LINE__); \
}

inline void cudaThrowMalloc(const char* call, cudaError_t code,
                            const char* func, const char* file, int line)
{
  char const* name = cudaGetErrorName(code);
  char const* string = cudaGetErrorString(code);
  std::cerr << call << " returned " << name << " (" << string << ")." << std::endl;
  std::cerr << "func: " << func << std::endl;
  std::cerr << "file: " << file << std::endl;
  std::cerr << "line: " << line << std::endl;
  throw std::runtime_error("");
}

inline void cudaThrowMalloc(const char* call, cudaError_t code,
                            std::size_t size, const char* purpose,
                            const char* func, const char* file, int line)
{
  char const* name = cudaGetErrorName(code);
  char const* string = cudaGetErrorString(code);
  std::cerr << call << " returned " << name << " (" << string << ")." << std::endl;
  std::cerr << "func: " << func << std::endl;
  std::cerr << "file: " << file << std::endl;
  std::cerr << "line: " << line << std::endl;
  std::cerr << "size: " << size << std::endl;
  std::cerr << "purpose: " << purpose << std::endl;
  throw std::runtime_error("");
}

#endif
