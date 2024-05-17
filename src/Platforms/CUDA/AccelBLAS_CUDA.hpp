//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
// File refactored from: MatrixDelayedUpdateCUDA.h
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_CUDA_ACCELBLAS_CUDA_H
#define QMCPLUSPLUS_CUDA_ACCELBLAS_CUDA_H

#include "AccelBLAS.hpp"
#include "CUDA/CUDAruntime.hpp"
#include "CUDA/QueueCUDA.hpp"
#include "CUDA/cuBLAS.hpp"

namespace qmcplusplus
{
namespace compute
{
template<>
class BLASHandle<PlatformKind::CUDA>
{
  public:
  // cuda stream, not owned, reference-only
  cudaStream_t h_stream;
  // cublas handle
  cublasHandle_t h_cublas;

  BLASHandle(Queue<PlatformKind::CUDA>& queue): h_stream(queue.getNative())
  {
    cublasErrorCheck(cublasCreate(&h_cublas), "cublasCreate failed!");
    cublasErrorCheck(cublasSetStream(h_cublas, h_stream), "cublasSetStream failed!");
  }

  ~BLASHandle()
  {
    cublasErrorCheck(cublasDestroy(h_cublas), "cublasDestroy failed!");
  }
};
}
}
#endif
