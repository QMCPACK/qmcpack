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

namespace BLAS
{
template<>
inline void gemm<PlatformKind::CUDA, float>(BLASHandle<PlatformKind::CUDA>& handle,
                           const char transa,
                           const char transb,
                           int m,
                           int n,
                           int k,
                           const float* alpha,
                           const float* A,
                           int lda,
                           const float* B,
                           int ldb,
                           const float* beta,
                           float* C,
                           int ldc)
{
  cublasErrorCheck(cublasSgemm(handle.h_cublas, cuBLAS::convertOperation(transa), cuBLAS::convertOperation(transb), m, n, k, alpha, A, lda, B, ldb, beta, C, ldc),
"cublasSgemm failed!");
}

template<>
inline void gemm<PlatformKind::CUDA, double>(BLASHandle<PlatformKind::CUDA>& handle,
                           const char transa,
                           const char transb,
                           int m,
                           int n,
                           int k,
                           const double* alpha,
                           const double* A,
                           int lda,
                           const double* B,
                           int ldb,
                           const double* beta,
                           double* C,
                           int ldc)
{
  cublasErrorCheck(cublasDgemm(handle.h_cublas, cuBLAS::convertOperation(transa), cuBLAS::convertOperation(transb), m, n, k, alpha, A, lda, B, ldb, beta, C, ldc),
"cublasDgemm failed!");
}

template<>
inline void gemm<PlatformKind::CUDA, std::complex<float>>(BLASHandle<PlatformKind::CUDA>& handle,
                           const char transa,
                           const char transb,
                           int m,
                           int n,
                           int k,
                           const std::complex<float>* alpha,
                           const std::complex<float>* A,
                           int lda,
                           const std::complex<float>* B,
                           int ldb,
                           const std::complex<float>* beta,
                           std::complex<float>* C,
                           int ldc)
{
  cublasErrorCheck(cublasCgemm(handle.h_cublas, cuBLAS::convertOperation(transa), cuBLAS::convertOperation(transb), m, n, k, castCUDAType(alpha), castCUDAType(A), lda, castCUDAType(B), ldb, castCUDAType(beta), castCUDAType(C), ldc),
"cublasCgemm failed!");
}

template<>
inline void gemm<PlatformKind::CUDA, std::complex<double>>(BLASHandle<PlatformKind::CUDA>& handle,
                           const char transa,
                           const char transb,
                           int m,
                           int n,
                           int k,
                           const std::complex<double>* alpha,
                           const std::complex<double>* A,
                           int lda,
                           const std::complex<double>* B,
                           int ldb,
                           const std::complex<double>* beta,
                           std::complex<double>* C,
                           int ldc)
{
  cublasErrorCheck(cublasZgemm(handle.h_cublas, cuBLAS::convertOperation(transa), cuBLAS::convertOperation(transb), m, n, k, castCUDAType(alpha), castCUDAType(A), lda, castCUDAType(B), ldb, castCUDAType(beta), castCUDAType(C), ldc),
"cublasZgemm failed!");
}
}
}
}
#endif
