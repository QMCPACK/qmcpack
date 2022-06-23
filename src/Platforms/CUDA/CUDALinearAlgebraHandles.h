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
class CUDALinearAlgebraHandles : public Resource
{
  struct DeleteCUDAStream
  {
    void operator()(cudaStream_t* stream_ptr)
    {
      if (stream_ptr)
      {
        cudaErrorCheck(cudaStreamDestroy(*stream_ptr), "cudaStreamDestroy failed!");
        delete stream_ptr;
        std::cout << "CUDAStream destroyed\n";
      }
    }
  };
  struct DeleteCuBLASHandle
  {
    void operator()(cublasHandle_t* cublas_handle_ptr)
    {
      if (cublas_handle_ptr)
      {
        cublasErrorCheck(cublasDestroy(*cublas_handle_ptr), "cublasDestroy failed!");
        delete cublas_handle_ptr;
      }
    }
  };

  // CUDA specific variables
  std::shared_ptr<cudaStream_t> hstream_ptr{nullptr, DeleteCUDAStream()};
  std::shared_ptr<cublasHandle_t> h_cublas_ptr{nullptr, DeleteCuBLASHandle()};

public:
  CUDALinearAlgebraHandles() : Resource("CUDALinearAlgebraHandles")
  {
    hstream_ptr  = std::shared_ptr<cudaStream_t>(new cudaStream_t{}, DeleteCUDAStream());
    h_cublas_ptr = std::shared_ptr<cublasHandle_t>(new cublasHandle_t{}, DeleteCuBLASHandle());
    cudaErrorCheck(cudaStreamCreate(hstream_ptr.get()), "cudaStreamCreate failed!");
    cublasErrorCheck(cublasCreate(h_cublas_ptr.get()), "cublasCreate failed!");
    cublasErrorCheck(cublasSetStream(*h_cublas_ptr, *hstream_ptr), "cublasSetStream failed!");
  }

  cublasHandle_t getCuBLAS() const { return *h_cublas_ptr; }
  cudaStream_t getStream() const { return *hstream_ptr; }

  Resource* makeClone() const override { return new CUDALinearAlgebraHandles(); }
};
} // namespace qmcplusplus
#endif
