//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SYCL_ACCELBLAS_SYCL_H
#define QMCPLUSPLUS_SYCL_ACCELBLAS_SYCL_H

#include "AccelBLASHandle.hpp"
#include "SYCL/QueueSYCL.hpp"
#include "SYCL/syclBLAS.hpp"

namespace qmcplusplus
{
namespace compute
{
template<>
class BLASHandle<PlatformKind::SYCL>
{
public:
  BLASHandle(Queue<PlatformKind::SYCL>& queue) : queue_(queue.getNative()) {}
  // sycl queue, not owned, reference-only
  sycl::queue& queue_;
};

namespace BLAS
{
template<typename T>
inline void gemm(BLASHandle<PlatformKind::SYCL>& handle,
                 const char transa,
                 const char transb,
                 int m,
                 int n,
                 int k,
                 const T& alpha,
                 const T* A,
                 int lda,
                 const T* B,
                 int ldb,
                 const T& beta,
                 T* C,
                 int ldc)
{
  oneapi::mkl::blas::gemm(handle.queue_, syclBLAS::convertTransEnum(transa), syclBLAS::convertTransEnum(transb), m, n, k, alpha, A, lda, B, ldb,
                                 beta, C, ldc);
}

template<typename T>
inline void gemv_batched(BLASHandle<PlatformKind::SYCL>& handle,
                         const char trans,
                         const int m,
                         const int n,
                         const T* alpha,
                         const T* const A[],
                         const int lda,
                         const T* const x[],
                         const int incx,
                         const T* beta,
                         T* const y[],
                         const int incy,
                         const int batch_count)
{
}

template<typename T>
inline void ger_batched(BLASHandle<PlatformKind::SYCL>& handle,
                        const int m,
                        const int n,
                        const T* alpha,
                        const T* const x[],
                        const int incx,
                        const T* const y[],
                        const int incy,
                        T* const A[],
                        const int lda,
                        const int batch_count)
{
}

template<typename T>
inline void copy_batched(BLASHandle<PlatformKind::SYCL>& handle,
                         const int n,
                         const T* const in[],
                         const int incx,
                         T* const out[],
                         const int incy,
                         const int batch_count)
{
}

inline void gemm_batched(BLASHandle<PlatformKind::SYCL>& handle,
                         const char transa,
                         const char transb,
                         int m,
                         int n,
                         int k,
                         const float& alpha,
                         const float* const A[],
                         int lda,
                         const float* const B[],
                         int ldb,
                         const float& beta,
                         float* const C[],
                         int ldc,
                         int batchCount)
{
}

inline void gemm_batched(BLASHandle<PlatformKind::SYCL>& handle,
                         const char transa,
                         const char transb,
                         int m,
                         int n,
                         int k,
                         const std::complex<float>& alpha,
                         const std::complex<float>* const A[],
                         int lda,
                         const std::complex<float>* const B[],
                         int ldb,
                         const std::complex<float>& beta,
                         std::complex<float>* const C[],
                         int ldc,
                         int batchCount)
{
}

inline void gemm_batched(BLASHandle<PlatformKind::SYCL>& handle,
                         const char transa,
                         const char transb,
                         int m,
                         int n,
                         int k,
                         const double& alpha,
                         const double* const A[],
                         int lda,
                         const double* const B[],
                         int ldb,
                         const double& beta,
                         double* const C[],
                         int ldc,
                         int batchCount)
{
}

inline void gemm_batched(BLASHandle<PlatformKind::SYCL>& handle,
                         const char transa,
                         const char transb,
                         int m,
                         int n,
                         int k,
                         const std::complex<double>& alpha,
                         const std::complex<double>* const A[],
                         int lda,
                         const std::complex<double>* const B[],
                         int ldb,
                         const std::complex<double>& beta,
                         std::complex<double>* const C[],
                         int ldc,
                         int batchCount)
{
}

} // namespace BLAS
} // namespace compute
} // namespace qmcplusplus
#undef castNativeType
#endif
