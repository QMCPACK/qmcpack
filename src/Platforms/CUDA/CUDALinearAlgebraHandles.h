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

#include "ResourceCollection.h"
#include "CUDA/CUDAruntime.hpp"
#include "CUDA/cuBLAS.hpp"

namespace qmcplusplus
{
struct CUDALinearAlgebraHandles : public Resource
{
  // CUDA specific variables
  cudaStream_t hstream;
  cublasHandle_t h_cublas;

  CUDALinearAlgebraHandles() : Resource("CUDALinearAlgebraHandles")
  {
    cudaErrorCheck(cudaStreamCreate(&hstream), "cudaStreamCreate failed!");
    cublasErrorCheck(cublasCreate(&h_cublas), "cublasCreate failed!");
    cublasErrorCheck(cublasSetStream(h_cublas, hstream), "cublasSetStream failed!");
  }

  CUDALinearAlgebraHandles(const CUDALinearAlgebraHandles&) : CUDALinearAlgebraHandles() {}

  ~CUDALinearAlgebraHandles()
  {
    cublasErrorCheck(cublasDestroy(h_cublas), "cublasDestroy failed!");
    cudaErrorCheck(cudaStreamDestroy(hstream), "cudaStreamDestroy failed!");
  }

  std::unique_ptr<Resource> makeClone() const override { return std::make_unique<CUDALinearAlgebraHandles>(*this); }
};
}
#endif
