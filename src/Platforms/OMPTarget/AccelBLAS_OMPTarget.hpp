//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_ACCELBLAS_OMPTARGET_H
#define QMCPLUSPLUS_ACCELBLAS_OMPTARGET_H

#include "AccelBLASHandle.hpp"
#include "QueueOMPTarget.hpp"
#include "ompBLAS.hpp"

namespace qmcplusplus
{
namespace compute
{
template<>
class BLASHandle<PlatformKind::OMPTARGET>
{
public:
  ompBLAS::ompBLAS_handle h_ompblas;

  BLASHandle(Queue<PlatformKind::OMPTARGET>& queue) : h_ompblas(0) {}
};

namespace BLAS
{

template<typename T>
inline void gemm(BLASHandle<PlatformKind::OMPTARGET>& handle,
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
  if (ompBLAS::gemm(handle.h_ompblas, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) != 0)
    throw std::runtime_error("ompBLAS::gemm failed!");
}

template<typename T>
inline void gemm_batched(BLASHandle<PlatformKind::OMPTARGET>& handle,
                         const char transa,
                         const char transb,
                         int m,
                         int n,
                         int k,
                         const T& alpha,
                         const T* const A[],
                         int lda,
                         const T* const B[],
                         int ldb,
                         const T& beta,
                         T* const C[],
                         int ldc,
                         int batchCount)
{
  if (ompBLAS::gemm_batched(handle.h_ompblas, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc,
                            batchCount) != 0)
    throw std::runtime_error("ompBLAS::gemm_batched failed!");
}


template<typename T>
inline void gemv(BLASHandle<PlatformKind::OMPTARGET>& handle,
                         const char trans,
                         const int m,
                         const int n,
                         const T& alpha,
                         const T* const A,
                         const int lda,
                         const T* const x,
                         const int incx,
                         const T& beta,
                         T* const y,
                         const int incy)
{
  if (ompBLAS::gemv(handle.h_ompblas, trans, m, n, alpha, A, lda, x, incx, beta, y, incy) != 0)
    throw std::runtime_error("ompBLAS::gemv_batched failed!");
}

template<typename T>
inline void gemv_batched(BLASHandle<PlatformKind::OMPTARGET>& handle,
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
  if (ompBLAS::gemv_batched(handle.h_ompblas, trans, m, n, alpha, A, lda, x, incx, beta, y, incy, batch_count) != 0)
    throw std::runtime_error("ompBLAS::gemv_batched failed!");
}

template<typename T>
inline void ger(BLASHandle<PlatformKind::OMPTARGET>& handle,
                        const int m,
                        const int n,
                        const T& alpha,
                        const T* const x,
                        const int incx,
                        const T* const y,
                        const int incy,
                        T* const A,
                        const int lda)
{
  if (ompBLAS::ger(handle.h_ompblas, m, n, alpha, x, incx, y, incy, A, lda) != 0)
    throw std::runtime_error("ompBLAS::ger_batched failed!");
}

template<typename T>
inline void ger_batched(BLASHandle<PlatformKind::OMPTARGET>& handle,
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
  if (ompBLAS::ger_batched(handle.h_ompblas, m, n, alpha, x, incx, y, incy, A, lda, batch_count) != 0)
    throw std::runtime_error("ompBLAS::ger_batched failed!");
}

template<typename T>
inline void copy_batched(BLASHandle<PlatformKind::OMPTARGET>& handle,
                         const int n,
                         const T* const x[],
                         const int incx,
                         T* const y[],
                         const int incy,
                         const int batch_count)
{
  if (ompBLAS::copy_batched(handle.h_ompblas, n, x, incx, y, incy, batch_count) != 0)
    throw std::runtime_error("ompBLAS::copy_batched failed!");
}

} // namespace BLAS
} // namespace compute
} // namespace qmcplusplus
#undef castNativeType
#endif
