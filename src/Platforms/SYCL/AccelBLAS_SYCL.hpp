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
                         const size_t batch_count)
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
                        const size_t batch_count)
{
}

template<typename T>
inline void copy_batched(BLASHandle<PlatformKind::SYCL>& handle,
                         const int n,
                         const T* const in[],
                         const int incx,
                         T* const out[],
                         const int incy,
                         const size_t batch_count)
{
}

template<typename T>
inline void gemm_batched(BLASHandle<PlatformKind::SYCL>& handle,
                         const char transa,
                         const char transb,
                         syclBLAS::syclBLAS_int m,
                         syclBLAS::syclBLAS_int n,
                         syclBLAS::syclBLAS_int k,
                         const T& alpha,
                         const T* const A[],
                         syclBLAS::syclBLAS_int lda,
                         const T* const B[],
                         syclBLAS::syclBLAS_int ldb,
                         const T& beta,
                         T* const C[],
                         syclBLAS::syclBLAS_int ldc,
                         const size_t batch_count)
{
  auto trans_a = syclBLAS::convertTransEnum(transa);
  auto trans_b = syclBLAS::convertTransEnum(transa);
oneapi::mkl::blas::gemm_batch(handle.queue_, sycl::span{&trans_a, 1}, sycl::span{&trans_b, 1}, sycl::span{&m, 1}, sycl::span{&n, 1}, sycl::span{&k, 1},
           sycl::span{const_cast<T*>(&alpha), 1}, sycl::span{const_cast<const T**>(A), batch_count}, sycl::span{&lda, 1}, sycl::span{const_cast<const T**>(B), batch_count},
	sycl::span{&ldb, 1}, sycl::span{const_cast<T*>(&beta), 1}, sycl::span{const_cast<T**>(C), batch_count}, sycl::span{&ldc, 1}, 1, sycl::span{const_cast<size_t*>(&batch_count), 1});
  //syclBLAS::syclBLAS_int bc = batch_count;
  //oneapi::mkl::blas::gemm_batch(handle.queue_, &trans_a, &trans_b, &m, &n, &k, const_cast<const T*>(&alpha), const_cast<const T**>(A), &lda, const_cast<const T**>(B), &ldb, const_cast<const T*>(&beta), const_cast<T**>(C), &ldc , 1, &bc);
}

} // namespace BLAS
} // namespace compute
} // namespace qmcplusplus
#undef castNativeType
#endif
