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

#ifndef QMCPLUSPLUS_CUDA_LINEAR_ALGEBRA_HANDLES_H
#define QMCPLUSPLUS_CUDA_LINEAR_ALGEBRA_HANDLES_H

#include "CUDA/CUDAruntime.hpp"
#include "CUDA/cuBLAS.hpp"

namespace qmcplusplus
{
struct CUDALinearAlgebraHandles
{
  // cuda stream, not owned, reference-only
  cudaStream_t h_stream;
  // cublas handle
  cublasHandle_t h_cublas;

  CUDALinearAlgebraHandles(cudaStream_t stream): h_stream(stream)
  {
    cublasErrorCheck(cublasCreate(&h_cublas), "cublasCreate failed!");
    cublasErrorCheck(cublasSetStream(h_cublas, stream), "cublasSetStream failed!");
  }

  ~CUDALinearAlgebraHandles()
  {
    cublasErrorCheck(cublasDestroy(h_cublas), "cublasDestroy failed!");
  }
};
}
#endif
